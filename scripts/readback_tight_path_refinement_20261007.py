"""Independent readback of a terminal smaller-step tight-path recovery window.

Only saved arrays, graph certificates and arithmetic ledgers are read. This
module does not import a process controller, mutate another reader's globals,
call NNPS, evaluate material/contact, or solve a physical step. A failed coarse
run remains failed even when its fully qualified prefix supplies references.
"""
import argparse
import math
from pathlib import Path
import sys
import time
import traceback
import json
import numpy as np
import readback_visible_support_process as old

ROOT = old.ROOT
CLASS = old.CLASS
COARSE_DT, FINE_DT, MAXIMUM_PAIRS = .000015625, .0000078125, 16
PASS_COARSE = 'QUALIFIED_CONTINUATION_COARSE_NOT_ACCEPTED'
PASS_FINE = 'PASS_VISIBLE_CONTINUATION_WINDOW_NOT_DEPOSITION_VALIDATION'
PASS_PREFIX = 'PASS_QUALIFIED_FINE_PREFIX_OF_FAILED_COARSE'
PASS_READBACK = 'PASS_TIGHT_PATH_REFINEMENT_READBACK_NOT_PHYSICAL_VALIDATION'
FIRST = ROOT/'outputs/tight_path_process_v2_20261004/window_002/fine/step_008/accepted.npz'
FIRST_SHA = '857d1fa88b2b74939d9c308306506918611c823847fde13f04684c24324ae408'
FIRST_AUDIT_SHA = '69984500c9337f4f2b4eb4e089755a38290c9b7fa48217e8dafd7f197b932212'
COMPONENTS = ROOT/'outputs/visible_top_components_20261004/v1/summary.json'
COMPONENTS_SHA = '7fbb5fdca371335ff528f4132fb840948dc3239de8eb50e420b3533014505cab'
CONTRACT = ROOT/'docs/TIGHT_PATH_REFINEMENT_CONTRACT_20261007.md'
PATH_COMPONENTS = ROOT/'outputs/visible_path_tight_20261004/qualification_v1/summary.json'
PATH_COMPONENTS_SHA = '3d461593a2fab30900a83bcde5b2835b353bec3d66783393b3500ce493b750bc'
PATH_QUALIFICATION = 'FIXED_RETAINED_GRAPH_RODRIGUES_CANCELLATION_ENCLOSURE_V1'
QUALIFICATION = 'VISIBLE_SUPPORT_SHOULDER_AND_TOP_EVENT_SCREEN_V1'
ALLOWED_FEATURES = ['SHOULDER', 'TOP']
require, clean, digest = old.require, old.clean, old.digest
positions, pose = old.positions, old.pose
Reader = old.Reader


def project_path(path):
    path = Path(path)
    path = (ROOT/path).resolve() if not path.is_absolute() else path.resolve()
    require(path.is_relative_to(ROOT), 'path escaped project')
    return path


def finite_properties(pa):
    for key, prop in pa.properties.items():
        require(np.isfinite(prop.get_npy_array()).all(), 'nonfinite particle field: '+key)


def validate_state(pa, state, context, number, dt, restart):
    initial, parent, p = context['initial'], context['parent'], context['p']
    require(state['classification'] == CLASS and state['restart_allowed'] is restart,
            'model/restart flag differs')
    require(state['step'] == parent['step']+number and
            abs(state['time_s']-(parent['time_s']+number*dt)) < 1e-13,
            'dynamic start time/global step differs')
    require(clean(state['spatial_model_transition']) == clean(parent['spatial_model_transition']),
            'continuation changed or repeated the spatial model transition')
    require(state['parameters'] == p and state['parameter_provenance'] == parent['parameter_provenance'],
            'parameters or provenance changed')
    require(state['source_sha256'] == parent['source_sha256'], 'underlying model source changed')
    require(set(pa.properties) == set(initial.properties), 'particle field schema changed')
    for key in old.INVARIANTS:
        require(np.array_equal(getattr(pa, key), getattr(initial, key)), 'identity/mass changed: '+key)
    require(np.array_equal(state['free'], parent['free']), 'prescribed/free axes changed')
    finite_properties(pa)
    require(np.all(pa.rho > 0) and np.all(pa.ep >= 0) and
            np.all(pa.temperature >= p['room']) and np.all(pa.temperature < p['melt']),
            'stored material/temperature domain failed')
    old.validate_interval_energy_state(state['interval_energy'])
    old.validate_bore_crossing_state(state['bore_observation'])
    records = state['interval_energy']['records']
    require(len(records) == state['step'] and records[:parent['step']] == context['prefix'],
            'complete inherited energy records were lost or reset')
    require(state['crossing_ids'] == state['bore_observation']['crossing_ids'] and
            state['crossing_mass_kg'] == state['bore_observation']['crossing_mass_kg'],
            'observer IDs/mass differ')


def correction_audit(saved):
    """Replay the saved graph's moment identity, not a force/material map."""
    x, i, j, g = saved['position1'], saved['i'], saved['j'], saved['g']
    v, correction = saved['volume'], saved['correction']
    require(correction.shape == (len(x), 3, 3) and np.isfinite(correction).all(),
            'invalid saved correction matrices')
    moment = np.zeros_like(correction)
    pair = (x[j]-x[i])[:, :, None]*g[:, None, :]
    for target, neighbor in ((i, j), (j, i)):
        np.add.at(moment, target, v[neighbor, None, None]*pair)
    ev = np.linalg.eigvalsh(moment)
    require(np.all(ev[:, 0] > 0) and np.all(ev[:, -1] <= 1e6*ev[:, 0]),
            'saved visible graph lost admissible moment support')
    errors = np.max(np.abs(moment@correction-np.eye(3)), axis=(1, 2))
    limits = 128*np.finfo(float).eps*np.maximum(1.,
        np.linalg.norm(moment, ord=np.inf, axis=(1, 2))*
        np.linalg.norm(correction, ord=np.inf, axis=(1, 2)))
    error = float(errors.max())
    require(np.all(errors <= limits), 'saved correction is inconsistent with the actual kept graph')
    return dict(moment_inverse_identity_max_absolute_error=error,
                moment_inverse_scaled_roundoff_max=float(np.max(errors/limits)),
                moment_condition_number_max=float(np.max(ev[:, -1]/ev[:, 0])))


def tight_path_arithmetic_audit(report, saved):
    require(report['classification'] == PATH_QUALIFICATION and report['nnps_evaluations'] == 0 and
            report['edge_deletions'] == 0, 'path used an unregistered qualification or changed the graph')
    bound = report['speed_bound']
    require(bound['minimum_tight_angle_rad'] == 1e-6 and
            bound['maximum_tight_angle_rad'] == np.pi/2 and
            bound['roundoff_factor'] == 8192*np.finfo(float).eps and
            report['roundoff_margin_m'] >= bound['original_displacement_roundoff_margin_m'],
            'registered path-bound domain or roundoff allowance differs')
    if bound['method'] == 'INHERITED_TRIANGLE_BOUND_FALLBACK':
        require(not 1e-6 <= report['rotation_angle_rad'] <= np.pi/2 and
                bound['fallback_reason'] == 'ROTATION_ANGLE_OUTSIDE_REGISTERED_TIGHT_DOMAIN' and
                bound['tightened_particle_count'] == 0,
                'unregistered tight-bound fallback')
        return
    require(bound['method'] == 'MIN_INHERITED_AND_MIDPOINT_DERIVATIVE_PLUS_HALF_CURVATURE' and
            bound['fallback_reason'] is None and 1e-6 <= report['rotation_angle_rad'] <= np.pi/2,
            'tight-bound method/domain mismatch')
    arrays = {key:np.asarray(bound[key]) for key in ('per_particle_old_speed_m',
        'per_particle_tight_speed_m', 'per_particle_selected_speed_m',
        'per_particle_midpoint_derivative_m', 'per_particle_curvature_upper_m',
        'per_particle_absolute_derivative_roundoff_m')}
    require(all(a.shape == (len(saved['particle_id']),) and np.isfinite(a).all() and np.all(a >= 0.)
                for a in arrays.values()), 'invalid saved per-particle bound arrays')
    oldspeed, tight, selected = (arrays[k] for k in ('per_particle_old_speed_m',
        'per_particle_tight_speed_m', 'per_particle_selected_speed_m'))
    reconstructed = np.nextafter((arrays['per_particle_midpoint_derivative_m']+
        .5*arrays['per_particle_curvature_upper_m']+arrays['per_particle_absolute_derivative_roundoff_m'])*
        (1.+8192*np.finfo(float).eps), np.inf)
    require(np.array_equal(tight, reconstructed) and np.array_equal(selected, np.minimum(oldspeed,tight)) and
            bound['tightened_particle_count'] == int(np.count_nonzero(tight<oldspeed)),
            'saved derivative/curvature/enlargement arithmetic differs')
    margin = float(np.nextafter(bound['original_displacement_roundoff_margin_m']+
        2*bound['uniform_position_evaluation_error_m'], np.inf))
    require(report['roundoff_margin_m'] == margin and
            bound['extra_displacement_roundoff_margin_m'] == margin-bound['original_displacement_roundoff_margin_m'],
            'path allowance lost its explicit uniform position error')
    speeds = np.maximum(selected[saved['i']],selected[saved['j']])
    require(report['maximum_edge_body_speed_per_fraction_m'] == float(speeds.max(initial=0.)),
            'edge-speed bound does not match the supplied retained graph')
    for cert in report['certificates']:
        lo,hi = cert['start_fraction'],cert['end_fraction']
        edges = np.asarray(cert['edge_indices'],dtype=int)
        require(cert['reference_fraction'] == .5*(lo+hi) and
                cert['displacement_plus_roundoff_bound_m'] >=
                .5*(hi-lo)*float(speeds[edges].max())+margin,
                'certificate bound does not enclose its own retained edges')


def candidate_audit(reader, stage, directory, row, previous, previous_state, context, frozen):
    path = reader.bind(directory/row['checkpoint'], row['sha256'])
    pa, state = reader.checkpoint(path)
    require(path.parent == directory/f"step_{row['step']:03d}" and
            path.name == 'NOT_RESTARTABLE_candidate.npz', 'candidate path is not the registered uncommitted file')
    required_artifacts = ('NOT_RESTARTABLE_candidate.npz', 'FULL_UNCOMMITTED_candidate.npz',
        'FULL_UNCOMMITTED_candidate.json', 'candidate_evidence.npz', 'energy.json',
        'visible_graph.npz', 'visible_graph.json', 'visible_path_audit.json',
        'continuous_shell_audit.json', 'geometry_audit.json', 'scalar_batch_equivalence.json',
        'sampled_nonbore_events.json', 'contact_event_qualification.json')
    for name in required_artifacts:
        require((path.parent/name).relative_to(directory).as_posix() in row['artifacts'],
                'candidate lacks required artifact hash: '+name)
    for relative, expected in row['artifacts'].items():
        actual = (directory/relative).resolve()
        require(actual.is_relative_to(path.parent.resolve()), 'candidate artifact escaped its step')
        reader.bind(actual, expected)
    dt = FINE_DT if stage == 'fine' else COARSE_DT
    require(row['dt_s'] == dt and row['global_step'] == context['parent']['step']+row['step'] and
            row['time_s'] == state['time_s'], 'candidate history time/step differs')
    require(abs(state['time_s']-previous_state['time_s']-dt) < 1e-15,
            'candidate actual immediate-parent time differs')
    validate_state(pa, state, context, row['step'], dt, False)
    require(state['implementation_source_sha256'] == frozen, 'candidate implementation source differs')
    require(state['interval_energy']['records'][:-1] == previous_state['interval_energy']['records'],
            'candidate interval direct parent differs')
    sd = path.parent
    graph_meta = reader.read(sd/'visible_graph.json')
    require(graph_meta['parent_time_s'] == previous_state['time_s'] and
            graph_meta['endpoint_time_s'] == state['time_s'], 'graph path time differs')
    saved, graph = old.graph_audit(reader, sd, pa, state, previous, context['p'], row)
    graph.update(correction_audit(saved))
    path_report = reader.read(sd/'visible_path_audit.json')
    tight_path_arithmetic_audit(path_report, saved)
    require(path_report['analytic_batch_queries'] == len(path_report['query_log']) and
            path_report['analytic_batch_queries'] == path_report['strict_geometry_batch_queries']+
                path_report['enclosure_batch_queries'] and
            path_report['analytic_batch_queries'] <= path_report['maximum_queries'] <= 256 and
            path_report['deepest_subdivision'] <= path_report['maximum_depth'] <= 32 and
            path_report['edge_deletions'] == 0 and path_report['physical_solver_evaluations'] == 0,
            'visible path certificate scope/budget differs')
    require(path_report['certified_edge_interval_count'] ==
            sum(len(c['edge_indices']) for c in path_report['certificates']),
            'certificate interval count differs')
    journal = reader.npz(sd/'FULL_UNCOMMITTED_candidate.npz')
    j = reader.read(sd/'FULL_UNCOMMITTED_candidate.json')
    energy = reader.read(sd/'energy.json')
    require(journal['accepted'].item() is False and journal['restart_allowed'].item() is False and
            j['accepted'] is False and j['restart_allowed'] is False, 'journal became accepted')
    values = dict(particle_id=pa.particle_id, position_m=positions(pa),
        velocity_m_s=np.column_stack((pa.u, pa.v, pa.w)), density_kg_m3=pa.rho,
        stress_Pa=pa.S.reshape(-1, 3, 3), temperature_K=pa.temperature, plastic_strain=pa.ep,
        contact_history=state['contact_history'])
    for key, value in values.items():
        require(np.array_equal(journal[key], value), 'full journal field differs: '+key)
    for key, value in journal.items():
        if value.dtype.kind in 'fc': require(np.isfinite(value).all(), 'nonfinite journal: '+key)
    require(np.allclose(journal['displacement_m'], positions(pa)-positions(previous), rtol=0, atol=1e-18),
            'journal displacement differs')
    require(np.array_equal(journal['force_N'], journal['internal_force_N']+journal['contact_force_N']),
            'complete force journal components do not sum')
    require(np.array_equal(journal['contact_force_N'], np.asarray(j['contact']['force_N'])) and
            np.array_equal(journal['contact_history'], np.asarray(j['contact']['history'])),
            'contact arrays differ between complete journals')
    require(clean(energy['mechanical']) == clean(state['energy']) == clean(j['energy']) and
            clean(energy['thermal']) == clean(j['thermal']) and
            clean(energy['interval']) == clean(state['interval_energy']), 'energy/thermal journal differs')
    actual = float(np.max(np.abs(journal['force_N'][np.asarray(state['free'], bool)])))
    require(actual <= j['metrics']['force_limit_N'] and j['metrics']['scaled_residual'] <= 1,
            'saved complete residual did not converge')
    e, th, p = energy['mechanical'], energy['thermal'], context['p']
    last = state['interval_energy']['records'][-1]
    require(last == {k:e[k] for k in last}, 'interval last record differs from raw energy')
    thermal = math.fsum(pa.m*p['cp']*(pa.temperature-previous.temperature))
    plastic = math.fsum(pa.plastic_work-previous.plastic_work)
    limit = 1e-12+1e-10*abs(th['total_heat_input_J'])
    require(abs(thermal-th['measured_thermal_increment_J']) <= limit and
            abs(thermal-th['total_heat_input_J']) <= limit, 'saved heat capacity balance failed')
    require(abs(plastic-th['plastic_work_J']) <= limit and
            abs(th['plastic_heat_J']-p['plastic_heat_fraction']*plastic) <= limit,
            'actual plastic work/heat differ')
    require(abs(th['internal_conduction_J']) <= limit, 'internal pair heat is not conservative')
    observation = old.observe_square_bore_crossings(previous_state['bore_observation'],
        positions(previous), positions(pa), pa.particle_id, pa.component_id, pa.m,
        previous_state['time_s'], state['time_s'], lambda t:pose(p, t))
    require(clean(observation['next_state']) == clean(state['bore_observation']),
            'strict crossing observer does not replay')
    shell = reader.read(sd/'continuous_shell_audit.json')
    geometry = reader.read(sd/'geometry_audit.json')
    require(shell['status'] == 'CERTIFIED_NO_HIDDEN_SHELL_CONTACT' and not shell['hidden_contacts'],
            'continuous particle shell certificate failed')
    require(not geometry['transient_contacts_missed_by_endpoints'] and
            set(geometry['new_non_bore_features']) <= set(ALLOWED_FEATURES), 'unqualified geometry event')
    require(reader.read(sd/'scalar_batch_equivalence.json')['status'] == 'PASS_SAME_FILTERED_COMPLETE_MAP',
            'same-filtered scalar/batch qualification absent')
    require(reader.read(sd/'sampled_nonbore_events.json') == row['events'] and
            set(e['feature'] for e in row['events']) <= set(ALLOWED_FEATURES), 'event history differs')
    qualification = reader.read(sd/'contact_event_qualification.json')
    require(qualification == dict(qualification=QUALIFICATION,
        allowed_nonbore_features=ALLOWED_FEATURES,
        actual_nonbore_features=sorted(set(e['feature'] for e in row['events'])),
        status='CANDIDATE_GATES_PASS_REQUIRES_SAME_PARENT_PAIR_SCREEN',
        physics_changed=False, candidate_sha256=row['sha256'], time_s=state['time_s']),
        'candidate TOP domain certificate differs from actual registered candidate')
    # Bind sampled endpoint events to the saved complete contact response;
    # this arithmetic check never recalculates geometry or contact forces.
    features = j['contact']['features']
    endpoint_load = np.asarray(j['contact']['endpoint_normal_load_N'])
    require(endpoint_load.shape == (len(pa.particle_id), len(features)),
            'contact endpoint load shape differs')
    observed_endpoint = {(int(e['particle_id']), e['feature']) for e in row['events'] if e['fraction'] == 1.}
    actual_endpoint = {(int(pa.particle_id[k]), features[f])
        for k, f in np.argwhere(endpoint_load > 0) if features[f] in ALLOWED_FEATURES}
    require(observed_endpoint == actual_endpoint, 'sampled endpoint contact events do not match raw contact load')
    top_index = features.index('TOP')
    top_contact = {}
    for key in ('normal_potential0_J', 'normal_potential1_J', 'tangential_potential0_J',
            'tangential_potential1_J', 'friction_dissipation_J', 'algorithmic_dissipation_J',
            'release_energy_J', 'normal_particle_work_J', 'normal_wall_work_J',
            'particle_work_J', 'wall_work_J'):
        values = np.asarray(j['contact'][key])
        require(values.shape == endpoint_load.shape and np.isfinite(values).all(),
                'invalid saved contact feature field: '+key)
        top_contact[key] = math.fsum(values[:, top_index])
    return dict(pa=pa, state=state, graph=saved, row=row,
        summary=dict(**old.state_stats(pa, state), **graph,
            checkpoint=path.relative_to(ROOT).as_posix(), checkpoint_sha256=row['sha256'],
            top_contact=top_contact,
            accepted=False, wall_seconds=row['wall_seconds'], force_residual_N=actual,
            global_raw_residual_J=e['global_raw_residual_J'], absolute_work_scale_J=e['absolute_work_scale_J'],
            external_work_J=e['prescribed_boundary_work_J']+e['tool_work_J'], plastic_work_J=plastic,
            thermal_capacity_increment_J=thermal, heat_balance_residual_J=thermal-th['total_heat_input_J']))


def replay_screens(reader, directory, cases, context):
    seeds = set(context['seeds'])
    fine, coarse, reports = cases.get('fine', []), cases.get('coarse', []), []
    for path in sorted((directory/'fine').glob('pair_*_temporal_screen.json')):
        screen = reader.read(path)
        number = screen['pair']
        require(number == len(reports)+1 and 2*number <= len(fine) and number <= len(coarse),
                'screen is not bound to complete candidate pairs')
        a, b = fine[2*number-2:2*number]
        c = coarse[number-1]
        require(screen['parent_sha256'] == context['parent_sha256'] and
                screen['coarse_checkpoint_sha256'] == c['row']['sha256'] and
                screen['fine_candidate_sha256'] == [a['row']['sha256'], b['row']['sha256']] and
                abs(screen['time_s']-b['state']['time_s']) < 1e-15 and
                abs(screen['time_s']-c['state']['time_s']) < 1e-15,
                'same-origin equal-time comparison binding differs')
        for item in (a, b, c):
            seeds.update(e['particle_id'] for e in item['row']['events'])
            seeds.update(item['row']['blocked_particle_ids'])
            seeds.update(item['state']['bore_observation']['crossing_ids'])
        require(screen['cumulative_event_and_blocked_ids'] == sorted(seeds),
                'cumulative event/crossing/blocked seeds differ')
        selected = set(seeds)
        for item in (a, b, c):
            g = item['graph']; ids = g['particle_id']
            for left, right in (('i', 'j'), ('raw_i', 'raw_j')):
                x, y = ids[g[left]], ids[g[right]]
                mask = np.isin(x, list(seeds)) | np.isin(y, list(seeds))
                selected.update(np.r_[x[mask], y[mask]].tolist())
        glob, local = old.compare_arrays(b['pa'], c['pa'], sorted(selected), context['p'])
        require(clean(glob) == screen['global_screen'] and clean(local) == screen['local_screen'],
                'saved global/local array comparison does not replay')
        match = b['state']['bore_observation']['crossing_ids'] == c['state']['bore_observation']['crossing_ids']
        require(screen['crossing_ids_match'] == match and
                screen['fine_crossing_ids'] == b['state']['bore_observation']['crossing_ids'] and
                screen['coarse_crossing_ids'] == c['state']['bore_observation']['crossing_ids'],
                'strict crossing comparison differs')
        same = b['row']['graph_fingerprint'] == c['row']['graph_fingerprint']
        require(screen['endpoint_graphs_match'] == same, 'endpoint graph comparison differs')
        reports.append(dict(pair=number, time_s=screen['time_s'], passed=glob['pass'] and local['passed'] and match,
            sha256=digest(path), global_screen=glob, local_particle_count=len(selected),
            local_max_stress_relative=max(x['stress_relative'] for x in local['rows']),
            local_max_plastic_relative=max(x['plastic_relative'] for x in local['rows']),
            crossing_ids_match=match, endpoint_graphs_match=same,
            cumulative_event_and_blocked_ids=sorted(seeds)))
    return reports


def path_component_sources(reader):
    reader.bind(PATH_COMPONENTS, PATH_COMPONENTS_SHA)
    receipt = reader.read(PATH_COMPONENTS)
    require(receipt['status'] == 'PASS_VISIBLE_PATH_TIGHT_COMPONENTS_AND_SAVED_REPLAY_NOT_PHYSICAL_VALIDATION' and
            receipt['qualification_passed'] is True and receipt['physical_solver_evaluations'] == 0 and
            receipt['nnps_evaluations'] == 0 and receipt['accepted_time_added_s'] == 0. and
            receipt['original_failure_preserved'] is True and receipt['original_candidate_restart_allowed'] is False and
            receipt['before_after_source_input_hashes_match'] is True,
            'tight-path prerequisite is not the independently qualified diagnostic receipt')
    protocol_path = PATH_COMPONENTS.parent/'protocol.json'
    reader.bind(protocol_path,receipt['protocol_sha256']); protocol=reader.read(protocol_path)
    require(protocol['registered_before_workers'] is True and protocol['maximum_queries'] == 256 and
            protocol['maximum_depth'] == 32 and protocol['expected_replay_edge_count'] == 7789 and
            protocol['no_force_nnps_or_real_process_evaluations'] is True and
            protocol['no_retroactive_acceptance_or_restart'] is True and
            protocol['contract'] == receipt['contract'] and
            protocol['source_and_input_sha256'] == receipt['source_and_input_sha256'],
            'tight-path qualification protocol differs')
    components=reader.read(PATH_COMPONENTS.parent/'components.json')
    require(components == receipt['components'] and components['status'] == 'PASS_COMPONENTS' and
            components['tests_run'] == 26 and all(components[k] == 0 for k in
                ('failures','errors','skipped','physical_solver_evaluations','nnps_evaluations')),
            'tight-path independent component result differs')
    evidence=dict(receipt['source_and_input_sha256'])
    evidence[receipt['replay_file']]=receipt['replay_sha256']
    require([v['worker'] for v in receipt['workers']] == ['components','replay'], 'prerequisite worker order differs')
    for worker in receipt['workers']:
        require(worker['exit_code'] == 0 and worker['external_limit_seconds'] == 120 and
                worker['external_wall_seconds'] <= 120 and
                reader.read(PATH_COMPONENTS.parent/(worker['worker']+'_supervision.json')) == worker,
                'tight-path prerequisite execution failed or differs')
        evidence[worker['log']]=worker['log_sha256']
    reader.sources(evidence)
    replay=reader.read(project_path(receipt['replay_file']))
    reader.bind(project_path(replay['original_result']),replay['original_result_sha256'])
    original=reader.read(project_path(replay['original_result']))
    require(receipt['replay_status'] == replay['status'] == 'CERTIFIED_VISIBLE_SUPPORT_PATH' and
            replay['classification'] == PATH_QUALIFICATION and replay['diagnostic_only'] is True and
            replay['original_candidate_restart_allowed'] is False and replay['original_failure_preserved'] is True and
            replay['accepted_time_added_s'] == 0. and replay['edge_count'] == 7789 and
            original['status'] == 'UNRESOLVED_VISIBLE_SUPPORT_PATH' and
            original['failure_kind'] == 'MAXIMUM_ANALYTIC_BATCH_QUERIES' and
            original['input_sha256'] == replay['input_sha256'] and
            original['graph_sha256'] == replay['graph_sha256'] and original['candidate_sha256'] == replay['candidate_sha256'],
            'saved geometry replay changed or retrospectively promoted the original failure')
    candidate_path=project_path(replay['original_result']).parent/'NOT_RESTARTABLE_candidate.npz'
    reader.bind(candidate_path,replay['candidate_sha256'])
    require(reader.checkpoint(candidate_path)[1]['restart_allowed'] is False,
            'original failed candidate became restartable')
    return evidence


def parent_context(reader, directory, meta):
    number = int(directory.name.removeprefix('window_'))
    require(directory.name == f'window_{number:03d}' and 1 <= number <= 4,
            'unregistered window name/index')
    parent_path = project_path(meta['parent_checkpoint'])
    reader.bind(parent_path, meta['parent_sha256'])
    initial, parent = reader.checkpoint(parent_path)
    require(parent['classification'] == CLASS and parent['restart_allowed'] is True and
            parent_path.name == 'accepted.npz', 'true restartable visible parent required')
    old.validate_interval_energy_state(parent['interval_energy'])
    old.validate_bore_crossing_state(parent['bore_observation'])
    require(len(parent['interval_energy']['records']) == parent['step'] and
            parent['step'] == meta['start_global_step'] and parent['time_s'] == meta['start_time_s'],
            'parent time/global step/full energy prefix differs')
    transition = parent['spatial_model_transition']
    require(transition['to_classification'] == CLASS and transition['diagnostic_only'] is False and
            transition['fields_reset'] is False and transition['interval_energy_reset'] is False and
            transition['observer_reset'] is False and
            transition['qualification'] == 'PASS_PARENT_MAP_TANGENT_AND_COMPONENT_GATES',
            'parent lacks qualified preserved spatial-model transition')
    finite_properties(initial)
    for key in ('source_sha256', 'implementation_source_sha256', 'controller_source_sha256'):
        reader.sources(parent[key])
    audit_path = project_path(meta['parent_readback_summary'])
    reader.bind(audit_path, meta['parent_readback_sha256'])
    prior = reader.read(audit_path)
    require(prior['status'] in (PASS_READBACK, 'PASS_TIGHT_PATH_PROCESS_READBACK_NOT_PHYSICAL_VALIDATION') and
            prior['accepted_history'][-1]['sha256'] == meta['parent_sha256'] and
            prior['accepted_global_steps'] == parent['step'] and prior['accepted_time_s'] == parent['time_s'],
            'independent parent audit does not end at the true parent')
    reader.sources(prior['evidence_sha256'])
    parent_folder = parent_path.parent.parent
    parent_manifest = reader.read(parent_folder/'manifest.json')
    parent_history = reader.read(parent_folder/'accepted_history.json')
    require(parent_history and parent_history[-1]['sha256'] == meta['parent_sha256'] and
            parent_folder/parent_history[-1]['checkpoint'] == parent_path and
            parent_history[-1]['global_step'] == parent['step'], 'parent is not its accepted branch endpoint')
    if number == 1:
        require(parent_path == FIRST and digest(parent_path) == FIRST_SHA and
                digest(audit_path) == FIRST_AUDIT_SHA and parent['step'] == 695 and
                abs(parent['time_s']-.0511875) < 1e-14, 'first continuation parent differs')
        require(parent_manifest['status'] == 'FAILED_WITH_ACCEPTED_PREFIX' and len(parent_history) == 8,
                'first parent is not the preserved audited 8-step prefix of the failed fine branch')
        require(prior['status'] == 'PASS_TIGHT_PATH_PROCESS_READBACK_NOT_PHYSICAL_VALIDATION' and
                prior['window_completed'] is False and prior['added_accepted_steps'] == 8 and
                prior['top_event_process_screen_passed'] is True,
                'refinement parent lost the audited incomplete window or inherited qualified TOP event')
    else:
        require(parent_folder == directory.parent/f'window_{number-1:03d}'/'fine' and
                parent_manifest['status'] == PASS_FINE and len(parent_history) == 2*MAXIMUM_PAIRS,
                'later window did not inherit the previous complete accepted fine branch')
        require(prior['status'] == PASS_READBACK and prior['top_event_process_screen_passed'] is True,
                'later window lacks the independently qualified inherited TOP process event')
    last_row = parent_history[-1]
    parent_screen = parent_folder/f"pair_{(last_row['step']+1)//2:03d}_temporal_screen.json"
    reader.bind(parent_screen, parent['window_temporal_screen_sha256'])
    screen = reader.read(parent_screen)
    require(screen['global_screen']['pass'] and screen['local_screen']['passed'] and
            screen['crossing_ids_match'], 'inherited locality screen is unqualified')
    seeds = set(screen['cumulative_event_and_blocked_ids']) | set(parent['bore_observation']['crossing_ids']) | {270, 339, 341}
    require(meta['inherited_seed_ids'] == sorted(seeds), 'window forgot prior event/blocked/crossing seeds')
    expected_source = dict(parent['implementation_source_sha256'])
    reader.bind(COMPONENTS, COMPONENTS_SHA)
    components = reader.read(COMPONENTS)
    require(components['status'] == 'PASS_TOP_CONTACT_COMPONENTS_NOT_PROCESS_VALIDATION' and
            components['exit_code'] == 0 and components['tests_run'] == 15 and
            components['physical_solver_evaluations'] == 0 and not components['changed_sources'],
            'TOP component prerequisite receipt is not the qualified 15-test evidence')
    for key in ('source_sha256', 'output_sha256'):
        reader.sources(components[key])
        expected_source.update(components[key])
    expected_source.update(path_component_sources(reader))
    for path in (ROOT/'scripts/run_tight_path_refinement_20261007.py', ROOT/'scripts/visible_tight_path_candidate.py',
                 CONTRACT, COMPONENTS, PATH_COMPONENTS, ROOT/'scripts/visible_support_path_tight.py',
                 parent_path, audit_path, ROOT/'scripts/run_visible_support_process.py'):
        reader.bind(path)
        expected_source[path.relative_to(ROOT).as_posix()] = digest(path)
    require(meta['source_sha256'] == expected_source, 'continuation source collection changed beyond registered controller links')
    reader.sources(expected_source)
    require(meta['branch_origin_sha256'] == parent['branch_origin_sha256'] and
            meta['window_origin_sha256'] == meta['parent_sha256'], 'branch/window origin conflated')
    return dict(initial=initial, parent=parent, p=parent['parameters'],
        prefix=parent['interval_energy']['records'], parent_path=parent_path,
        parent_sha256=meta['parent_sha256'], audit_path=audit_path,
        audit_sha256=meta['parent_readback_sha256'], seeds=seeds, frozen=expected_source, number=number,
        inherited_top_event_passed=bool(prior['top_event_process_screen_passed']))


def stage_audit(reader, stage, directory, context):
    folder = directory/stage
    if not folder.exists():
        return [], dict(status='NOT_RUN', accepted_steps=0, candidate_count=0, history=[])
    meta = reader.read(folder/'manifest.json')
    supervision = reader.read(project_path(meta['supervision_receipt']))
    require(supervision['limit_s'] == 600 and supervision['external_wall_s'] <= 605,
            'external supervision budget differs')
    timed_out = supervision['exit_code'] == 124
    require(meta['status'] != 'RUNNING' or timed_out, 'stage is still running')
    if meta['status'] == 'RUNNING':
        reader.read(folder/'HARD_TIMEOUT.json')
    successes = (PASS_COARSE,) if stage == 'coarse' else (PASS_FINE, PASS_PREFIX)
    successful = meta['status'] in successes
    require((supervision['exit_code'] == 0) == successful, 'terminal manifest/exit code differs')
    require(meta['classification'] == CLASS and meta['stage'] == stage and
            meta['window'] == context['number'], 'stage/model/window differs')
    parent = context['parent']
    require(meta['parent_sha256'] == context['parent_sha256'] and
            project_path(meta['parent_checkpoint']) == context['parent_path'] and
            meta['start_time_s'] == parent['time_s'] and meta['start_global_step'] == parent['step'] and
            meta['window_origin_sha256'] == context['parent_sha256'] and
            meta['branch_origin_sha256'] == parent['branch_origin_sha256'], 'stage shared origin differs')
    require(meta['parent_readback_sha256'] == context['audit_sha256'] and
            project_path(meta['parent_readback_summary']) == context['audit_path'], 'stage parent audit differs')
    require(meta['source_sha256'] == context['frozen'] and meta['provenance'] == parent['parameter_provenance'] and
            meta['inherited_seed_ids'] == sorted(context['seeds']), 'stage source/inputs/inherited seeds differ')
    require(meta['contact_event_qualification'] == QUALIFICATION and
            meta['allowed_nonbore_features'] == ALLOWED_FEATURES and
            project_path(meta['contact_event_contract']) == CONTRACT and
            meta['contact_event_contract_sha256'] == digest(CONTRACT) and
            project_path(meta['component_receipt']) == COMPONENTS and
            meta['component_receipt_sha256'] == COMPONENTS_SHA,
            'stage changed its registered TOP event/component qualification')
    require(meta['path_certificate_qualification'] == PATH_QUALIFICATION and
            project_path(meta['path_component_receipt']) == PATH_COMPONENTS and
            meta['path_component_receipt_sha256'] == PATH_COMPONENTS_SHA,
            'stage changed its registered path qualification')
    reader.sources(meta['source_sha256'])
    count = meta['reference_pair_count']
    require(meta['registered_pair_count'] == MAXIMUM_PAIRS and 1 <= count <= MAXIMUM_PAIRS and
            abs(meta['registered_end_s']-(parent['time_s']+MAXIMUM_PAIRS*COARSE_DT)) < 1e-15 and
            abs(meta['reference_end_s']-(parent['time_s']+count*COARSE_DT)) < 1e-15,
            'registered/reference window lengths differ')
    require(meta['forecast_wall_s'] <= 450 and
            abs(meta['forecast_wall_s']-(30+1.5*count*(2 if stage == 'fine' else 1)*
                meta['complete_stage_cost_per_step_s'])) < 1e-10, 'registered complete-stage forecast differs')
    if not timed_out:
        require(meta['wall_seconds'] <= 540.5, 'stage exceeded internal time cap')
    rows = reader.read(folder/'candidate_history.json')
    history = reader.read(folder/'accepted_history.json')
    require(project_path(meta['candidate_history_path']) == folder/'candidate_history.json',
            'candidate history path differs')
    transfer, ts = reader.checkpoint(folder/'CONTINUATION_PARENT_NOT_RESTARTABLE.npz')
    old.properties_equal(transfer, context['initial'])
    require(clean(ts) == clean(dict(parent, restart_allowed=False)),
            'continuation transfer is not an exact parent copy with restart disabled')
    previous, previous_state, direct = context['initial'], parent, context['parent_sha256']
    records = []
    for index, row in enumerate(rows, 1):
        require(row['step'] == index and row['parent_sha256'] == direct,
                'candidate chain does not follow its actual immediate parent')
        item = candidate_audit(reader, stage, folder, row, previous, previous_state, context, context['frozen'])
        records.append(item)
        previous, previous_state, direct = item['pa'], item['state'], row['sha256']
        if stage == 'fine' and index % 2 == 0:
            published = folder/f'step_{index:03d}'/'accepted.npz'
            if published.exists(): direct = digest(reader.bind(published))
    require(len(rows) <= count*(2 if stage == 'fine' else 1), 'candidate rows exceed registered reference scope')
    if successful:
        require(meta['candidate_steps'] == len(rows) == count*(2 if stage == 'fine' else 1),
                'successful stage candidate count is incomplete')
        require(abs(meta['candidate_time_s']-meta['reference_end_s']) < 1e-14,
                'successful stage stopped before its reference endpoint')
    else:
        require(0 <= meta['candidate_steps'] <= len(rows), 'manifest invents complete candidates')
    if rows and meta['candidate_steps'] == len(rows):
        require(meta['candidate_time_s'] == rows[-1]['time_s'], 'manifest candidate endpoint differs')
    if successful:
        top_ids = sorted({e['particle_id'] for row in rows for e in row['events'] if e['feature'] == 'TOP'})
        require(meta['top_event_observed_in_candidates'] == bool(top_ids) and
                meta['top_particle_ids'] == top_ids, 'manifest invents or omits observed TOP candidates')
    published = sorted(folder.glob('step_*/accepted.npz'))
    if stage == 'coarse':
        require(not published and not history and meta['accepted_steps'] == 0 and
                meta['accepted_time_s'] == parent['time_s'] and count == MAXIMUM_PAIRS,
                'coarse became restartable or changed its registered scope')
    known = {Path(row['checkpoint']).parent.name for row in rows}
    extras = []
    for sd in sorted(folder.glob('step_*')):
        if sd.name not in known:
            files = [reader.bind(p).relative_to(ROOT).as_posix() for p in sorted(sd.rglob('*')) if p.is_file()]
            require(not (sd/'accepted.npz').exists(), 'unqualified partial candidate was published')
            extras.append(dict(directory=sd.relative_to(ROOT).as_posix(),
                qualification='UNREGISTERED_FAILED_OR_PARTIAL_CANDIDATE_NOT_ACCEPTED', evidence_files=files))
    return records, dict(status=meta['status'], timed_out=timed_out, exit_code=supervision['exit_code'],
        wall_seconds=meta.get('wall_seconds'), external_wall_seconds=supervision['external_wall_s'],
        forecast_wall_s=meta['forecast_wall_s'], candidate_count=len(rows),
        declared_candidate_steps=meta['candidate_steps'], accepted_steps=meta['accepted_steps'],
        accepted_time_s=meta['accepted_time_s'], reference_pair_count=count,
        coarse_original_status=meta['coarse_original_status'], candidates=[v['summary'] for v in records],
        unregistered_partial_evidence=extras, manifest_sha256=digest(folder/'manifest.json'), history=history)


def audit_accepted(reader, directory, cases, stages, screens, context):
    accepted, direct = [], context['parent_sha256']
    folder = directory/'fine'
    history = stages['fine']['history']
    publishing = {'restart_allowed', 'implementation', 'source_sha256', 'implementation_source_sha256',
        'controller_source_sha256', 'continuation_parent_sha256', 'branch_origin_sha256',
        'window_origin_sha256', 'window_temporal_screen_sha256', 'window_coarse_reference_sha256',
        'event_candidate_accepted', 'contact_event_qualification', 'path_certificate_qualification'}
    for index, path in enumerate(sorted(folder.glob('step_*/accepted.npz')), 1):
        require(path.parent.name == f'step_{index:03d}' and index <= len(cases['fine']),
                'accepted checkpoints are not a contiguous qualified prefix')
        candidate = cases['fine'][index-1]
        pa, state = reader.checkpoint(path)
        old.properties_equal(pa, candidate['pa'])
        validate_state(pa, state, context, index, FINE_DT, True)
        require(state['path_certificate_qualification'] == dict(qualification=PATH_QUALIFICATION,
            contract_sha256=digest(CONTRACT), component_receipt_sha256=PATH_COMPONENTS_SHA,physics_changed=False),
            'published path qualification differs from independently registered certificate')
        require(state['contact_event_qualification'] == dict(qualification=QUALIFICATION,
            allowed_nonbore_features=ALLOWED_FEATURES, contract_sha256=digest(CONTRACT),
            component_receipt_sha256=COMPONENTS_SHA, physics_changed=False),
            'accepted qualification metadata differs from the registered unchanged physical map')
        require(clean({k:v for k,v in state.items() if k not in publishing}) ==
                clean({k:v for k,v in candidate['state'].items() if k not in publishing}),
                'publication changed non-controller state fields')
        pair = (index+1)//2
        require(pair <= len(screens) and screens[pair-1]['passed'], 'accepted half lacks a passing complete pair screen')
        require(state['continuation_parent_sha256'] == direct and
                state['branch_origin_sha256'] == context['parent']['branch_origin_sha256'] and
                state['window_origin_sha256'] == context['parent_sha256'] and
                state['implementation_source_sha256'] == context['frozen'] and
                state['controller_source_sha256'] == context['frozen'] and
                state['window_temporal_screen_sha256'] == screens[pair-1]['sha256'] and
                state['window_coarse_reference_sha256'] == cases['coarse'][pair-1]['row']['sha256'],
                'published direct-parent/source/branch/window/screen binding differs')
        old.properties_equal(pa, reader.checkpoint(path.parent/'UNPUBLISHED_candidate.npz')[0])
        staged = reader.checkpoint(path.parent/'UNPUBLISHED_candidate.npz')[1]
        require(clean(staged) == clean(dict(state, restart_allowed=False)), 'staged publication metadata differs')
        require(reader.read(path.parent/'transport.json') == clean(state['bore_observation']) and
                reader.read(path.parent/'tool_pose.json') == clean(pose(context['p'], state['time_s'])),
                'published observer/tool pose output differs')
        reader.bind(path.parent/'fields.vtk')
        previous_hash = direct
        direct = digest(path)
        entry = dict(step=index, global_step=state['step'], time_s=state['time_s'], sha256=direct,
            parent_sha256=previous_hash, checkpoint=path.relative_to(ROOT).as_posix())
        accepted.append(entry)
        if index <= len(history):
            row = history[index-1]
            require(row['step'] == index and row['global_step'] == state['step'] and row['dt_s'] == FINE_DT and
                    row['time_s'] == state['time_s'] and row['sha256'] == direct and
                    folder/row['checkpoint'] == path and row['parent_sha256'] == previous_hash and
                    row['temporal_screen_sha256'] == state['window_temporal_screen_sha256'] and
                    row['coarse_reference_sha256'] == state['window_coarse_reference_sha256'],
                    'accepted history row differs from actual checkpoint')
    require(len(history) <= len(accepted), 'history lists missing accepted checkpoint')
    fine = stages['fine']
    if fine['status'] != 'NOT_RUN':
        require(0 <= fine['accepted_steps'] <= len(history), 'manifest invents accepted states')
        if fine['accepted_steps']:
            require(fine['accepted_time_s'] == history[fine['accepted_steps']-1]['time_s'],
                    'manifest accepted endpoint differs')
        else:
            require(fine['accepted_time_s'] == context['parent']['time_s'], 'empty fine branch claims time advance')
    if fine['status'] in (PASS_FINE, PASS_PREFIX):
        count = fine['reference_pair_count']
        require(len(accepted) == len(history) == fine['accepted_steps'] == 2*count and
                len(screens) == count and all(v['passed'] for v in screens), 'successful fine scope is incomplete')
        require((fine['status'] == PASS_FINE) == (count == MAXIMUM_PAIRS), 'full-window/prefix status conflated')
    return accepted


def export_accepted(reader, output, accepted, context):
    initial, parent, p = context['initial'], context['parent'], context['p']
    frames = [(initial, parent)]+[reader.checkpoint(ROOT/row['checkpoint']) for row in accepted]
    endpa, endstate = frames[-1]
    times = np.array([state['time_s'] for _, state in frames])
    poses = [pose(p, t) for t in times]
    fields = {'particle_field__'+name:np.array([pa.properties[name].get_npy_array() for pa, _ in frames])
              for name in initial.properties}
    trajectory = output/'trajectory.npz'
    np.savez_compressed(trajectory, time_s=times,
        global_step=np.array([state['step'] for _, state in frames], dtype=np.int64),
        particle_id=initial.particle_id, position_m=np.array([positions(pa) for pa, _ in frames]),
        velocity_m_s=np.array([np.column_stack((pa.u, pa.v, pa.w)) for pa, _ in frames]),
        density_kg_m3=np.array([pa.rho for pa, _ in frames]),
        temperature_K=np.array([pa.temperature for pa, _ in frames]),
        stress_Pa=np.array([pa.S.reshape(-1, 3, 3) for pa, _ in frames]),
        plastic_strain=np.array([pa.ep for pa, _ in frames]),
        plastic_work_J=np.array([pa.plastic_work for pa, _ in frames]),
        contact_history=np.array([state['contact_history'] for _, state in frames]),
        tool_origin_m=np.array([v['origin'] for v in poses]), tool_rotation=np.array([v['rotation'] for v in poses]),
        mass_kg=initial.m, component_id=initial.component_id, material_id=initial.material_id,
        tracer_type=initial.tracer_type, wholebar_tag=initial.wholebar_tag, **fields)
    ledger = output/'complete_interval_energy.json'
    ledger.write_text(json.dumps(clean(endstate['interval_energy']), indent=2, allow_nan=False)+'\n', encoding='utf-8')
    before = (positions(initial)-poses[0]['origin'])@poses[0]['rotation']
    after = (positions(endpa)-poses[-1]['origin'])@poses[-1]['rotation']
    response = []
    for pid in endstate['bore_observation']['crossing_ids']:
        k = int(np.flatnonzero(initial.particle_id == pid)[0])
        r0, r1 = float(np.linalg.norm(before[k, :2])), float(np.linalg.norm(after[k, :2]))
        response.append(dict(particle_id=pid, mass_kg=float(endpa.m[k]),
            already_crossed_at_window_parent=pid in parent['bore_observation']['crossing_ids'],
            plastic_strain_start=float(initial.ep[k]), plastic_strain_end=float(endpa.ep[k]),
            plastic_strain_increment=float(endpa.ep[k]-initial.ep[k]),
            plastic_work_increment_J=float(endpa.plastic_work[k]-initial.plastic_work[k]),
            temperature_start_K=float(initial.temperature[k]), temperature_end_K=float(endpa.temperature[k]),
            temperature_increment_K=float(endpa.temperature[k]-initial.temperature[k]),
            tool_radius_start_m=r0, tool_radius_end_m=r1, tool_radius_increment_m=r1-r0,
            tool_z_start_m=float(before[k, 2]), tool_z_end_m=float(after[k, 2]),
            tool_z_increment_m=float(after[k, 2]-before[k, 2]),
            strict_square_bore_projection_start=bool(np.all(np.abs(before[k, :2]) < p['bore_width']/2)),
            strict_square_bore_projection_end=bool(np.all(np.abs(after[k, :2]) < p['bore_width']/2)),
            square_bore_projection_margin_end_m=float(p['bore_width']/2-np.max(np.abs(after[k, :2])))))
    return endpa, endstate, response, dict(path=trajectory.relative_to(ROOT).as_posix(),
        sha256=digest(trajectory), frame_count=len(frames), saved_particle_property_count=len(initial.properties),
        scope='ACTUAL_WINDOW_TRUE_PARENT_PLUS_NEW_FINE_ACCEPTED_ONLY_NO_INTERPOLATION_NO_COARSE',
        inherited_readback_summary=context['audit_path'].relative_to(ROOT).as_posix(),
        inherited_readback_sha256=context['audit_sha256'],
        complete_interval_energy_path=ledger.relative_to(ROOT).as_posix(), complete_interval_energy_sha256=digest(ledger))


def top_event_observation(cases, accepted, screens, context):
    """Observe saved per-feature work and sampled contacts; no new tolerance."""
    branches = {}
    for name, records in (('coarse_qualified_reference', cases['coarse']),
                           ('fine_accepted_only', cases['fine'][:len(accepted)])):
        first = next((item for item in records if any(e['feature'] == 'TOP' for e in item['row']['events'])), None)
        first_interval = None
        if first is not None:
            events = [e for e in first['row']['events'] if e['feature'] == 'TOP']
            end = first['state']['time_s']; dt = first['row']['dt_s']
            first_interval = dict(step=first['row']['step'], start_time_s=end-dt, end_time_s=end,
                particle_ids=sorted({e['particle_id'] for e in events}),
                first_sample_fraction=min(e['fraction'] for e in events),
                sampled_event_time_is_not_exact_contact_onset=True)
        ids = sorted({e['particle_id'] for item in records for e in item['row']['events'] if e['feature'] == 'TOP'})
        endpoint_ids = sorted({e['particle_id'] for e in records[-1]['row']['events']
                              if e['feature'] == 'TOP' and e['fraction'] == 1.}) if records else []
        energy = {}
        for key in ('friction_dissipation_J', 'algorithmic_dissipation_J', 'release_energy_J',
                    'normal_particle_work_J', 'normal_wall_work_J', 'particle_work_J', 'wall_work_J'):
            energy[key] = math.fsum(item['summary']['top_contact'][key] for item in records)
        for kind in ('normal', 'tangential'):
            start = records[0]['summary']['top_contact'][kind+'_potential0_J'] if records else None
            end = records[-1]['summary']['top_contact'][kind+'_potential1_J'] if records else None
            energy[kind+'_potential'] = dict(start_J=start, end_J=end,
                change_J=end-start if records else None)
        branches[name] = dict(candidate_count=len(records), first_sampled_TOP_increment=first_interval,
            observed_particle_ids=ids, endpoint_particle_ids=endpoint_ids,
            endpoint_time_s=records[-1]['state']['time_s'] if records else context['parent']['time_s'],
            whole_branch_TOP_energy=energy)
    coarse, fine = branches['coarse_qualified_reference'], branches['fine_accepted_only']
    observed_both = bool(coarse['observed_particle_ids'] and fine['observed_particle_ids'])
    passed = context['inherited_top_event_passed'] or observed_both
    return dict(qualification=QUALIFICATION, process_screen_passed=passed,
        inherited_process_screen_passed=context['inherited_top_event_passed'],
        actual_TOP_observed_in_coarse_and_accepted_fine=observed_both,
        branches=branches,
        equal_endpoint_time=coarse['endpoint_time_s'] == fine['endpoint_time_s'],
        endpoint_TOP_ids_equal=coarse['endpoint_particle_ids'] == fine['endpoint_particle_ids'],
        accepted_pair_count=len(accepted)//2,
        all_accepted_pairs_passed=all(s['passed'] for s in screens[:len(accepted)//2]),
        contact_history_temporal_convergence='NOT_INDEPENDENTLY_CERTIFIED',
        scope='SAVED_CONTACT_OBSERVATION_AND_REGISTERED_PROCESS_NUMERICAL_SCREEN_NOT_PHYSICAL_VALIDATION')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--window', required=True, help='existing terminal window directory')
    parser.add_argument('--output', required=True, help='new readback output directory')
    args = parser.parse_args()
    directory, output = project_path(args.window), project_path(args.output)
    require(directory.parent == ROOT/'outputs/tight_path_refinement_20261007', 'unexpected tight-path process series')
    require(directory.is_dir() and not output.exists(), 'window must exist and output must be new')
    require(not output.is_relative_to(directory/'coarse') and not output.is_relative_to(directory/'fine'),
            'readback outputs cannot enter process evidence directories')
    output.mkdir(parents=True, exist_ok=False)
    begin = time.perf_counter()
    reader = Reader()
    report = dict(status='READBACK_RUNNING', classification='SAVED_TIGHT_PATH_REFINEMENT_EVIDENCE_AUDIT',
        command=[sys.executable, *sys.argv], physical_solver_evaluations=0, nnps_evaluations=0,
        geometry_certificate_evaluations=0, production_authorized=False, deposition_validated=False,
        window_directory=directory.relative_to(ROOT).as_posix())
    try:
        reader.bind(Path(__file__))
        reader.bind(Path(old.__file__))
        coarse_meta = reader.read(directory/'coarse/manifest.json')
        context = parent_context(reader, directory, coarse_meta)
        cases, stages = {}, {}
        for stage in ('coarse', 'fine'):
            cases[stage], stages[stage] = stage_audit(reader, stage, directory, context)
        coarse, fine = stages['coarse'], stages['fine']
        if fine['status'] != 'NOT_RUN':
            require(coarse['status'] in (PASS_COARSE, 'FAILED_NO_NEW_ACCEPTED_STATE') and
                    0 < len(cases['coarse']) <= MAXIMUM_PAIRS and
                    fine['reference_pair_count'] == len(cases['coarse']) and
                    fine['coarse_original_status'] == coarse['status'],
                    'fine does not use the declared fully audited coarse prefix')
            if fine['status'] == PASS_PREFIX:
                require(coarse['status'] == 'FAILED_NO_NEW_ACCEPTED_STATE' and
                        len(cases['coarse']) < MAXIMUM_PAIRS, 'partial-prefix success lost original coarse failure')
        screens = replay_screens(reader, directory, cases, context)
        accepted = audit_accepted(reader, directory, cases, stages, screens, context)
        endpa, endstate, response, trajectory = export_accepted(reader, output, accepted, context)
        energy = endstate['interval_energy']
        newrecords = energy['records'][context['parent']['step']:]
        new = cases['fine'][:len(accepted)]
        full_window = fine['status'] == PASS_FINE and len(accepted) == 2*MAXIMUM_PAIRS
        top_observation = top_event_observation(cases, accepted, screens, context)
        if full_window:
            require(top_observation['process_screen_passed'],
                    'successful TOP window did not observe actual TOP in coarse and accepted fine or inherit its qualification')
        report.update(status=PASS_READBACK, stages=stages, pair_screens=screens,
            path_certificate_qualification=PATH_QUALIFICATION,
            path_component_receipt_sha256=PATH_COMPONENTS_SHA,
            original_failed_candidate_not_promoted=True,
            top_event_process_screen_passed=top_observation['process_screen_passed'],
            top_event_observation=top_observation,
            window_completed=full_window,
            window_outcome=('FULL_REGISTERED_WINDOW_ACCEPTED_NOT_DEPOSITION' if full_window else
                'PARTIAL_ACCEPTED_PREFIX_OR_NO_ADVANCE_ORIGINAL_FAILURES_PRESERVED'),
            parent_checkpoint=context['parent_path'].relative_to(ROOT).as_posix(),
            parent_sha256=context['parent_sha256'], parent_readback_summary=context['audit_path'].relative_to(ROOT).as_posix(),
            parent_readback_sha256=context['audit_sha256'], branch_origin_sha256=context['parent']['branch_origin_sha256'],
            window_origin_sha256=context['parent_sha256'],
            initial_accepted_global_steps=context['parent']['step'], added_accepted_steps=len(accepted),
            accepted_global_steps=endstate['step'], accepted_time_s=endstate['time_s'],
            accepted_time_advanced_s=endstate['time_s']-context['parent']['time_s'],
            actual_accepted_frame_count_including_independent_initial=endstate['step']+1,
            inherited_interval_never_reset=True, repeated_model_transition=False, accepted_history=accepted,
            unrecorded_published_halves_reconciled=len(accepted)-len(fine['history']),
            final_accepted_fields=old.state_stats(endpa, endstate), crossing_events=endstate['bore_observation']['events'],
            crossed_particle_response=response, trajectory=trajectory,
            complete_prefix_energy={k:v for k,v in energy.items() if k != 'records'},
            new_accepted_energy=dict(signed_raw_residual_J=math.fsum(r['global_raw_residual_J'] for r in newrecords),
                absolute_step_raw_residual_J=math.fsum(abs(r['global_raw_residual_J']) for r in newrecords),
                absolute_work_scale_J=math.fsum(r['absolute_work_scale_J'] for r in newrecords),
                external_work_J=math.fsum(v['summary']['external_work_J'] for v in new),
                plastic_work_J=math.fsum(v['summary']['plastic_work_J'] for v in new),
                measured_thermal_increment_J=math.fsum(v['summary']['thermal_capacity_increment_J'] for v in new)),
            limitations=['PASS certifies saved evidence consistency, not a successful full process window or deposition.',
                'Local 15.625/7.8125us screening preserves prior failures; it does not certify convergence of the inherited 0--51.1875ms trajectory.',
                'Every inherited energy prefix is validated; prior spatial/field history is bound through the parent readback receipt.',
                'Actual stored graph/correction, path inputs and whole [0,1] edge coverage are checked without NNPS or geometry reevaluation.',
                'All saved particle fields survive candidate publication unchanged; journals, heat and strict-cross arithmetic are replayed without a fresh force/material/thermal solve.',
                'Failed coarse status and partial candidate artifacts remain intact; only qualified rows can supply equal-time references.',
                'The NPZ contains the true window parent and actual newly published fine frames, never coarse states or interpolation.',
                'Crossed-particle plastic response and bore projection are observations, not deposition or experimental validation.'])
        reader.finish()
        report['evidence_sha256'] = reader.evidence
        report['source_entry_checks_including_repeats'] = reader.source_entries
        code = 0
    except Exception:
        report.update(status='FAIL_TIGHT_PATH_REFINEMENT_READBACK', error=traceback.format_exc(),
                      evidence_sha256=reader.evidence)
        code = 2
    report.update(exit_code=code, wall_seconds=time.perf_counter()-begin)
    (output/'summary.json').write_text(json.dumps(clean(report), indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k:report.get(k) for k in ('status', 'exit_code', 'window_completed',
        'added_accepted_steps', 'accepted_global_steps', 'accepted_time_s', 'wall_seconds', 'error')}, indent=2))
    return code


if __name__ == '__main__':
    raise SystemExit(main())
