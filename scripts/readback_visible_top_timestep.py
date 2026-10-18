"""Independent saved-evidence audit of nonrestartable TOP time-step trials.

Replays arithmetic, references, certificates and local comparisons only. A
consistent failed numerical screen remains failed; no candidate is published.
"""
import argparse
import json
import math
from pathlib import Path
import sys
import time
import traceback
import numpy as np
import readback_visible_top_event as reference

old = reference.old
ROOT, CLASS = reference.ROOT, reference.CLASS
BASE = ROOT/'outputs/visible_top_timestep_20261004'
REFERENCE = ROOT/'outputs/visible_top_event_20261004/window_001'
CONTRACT = ROOT/'docs/VISIBLE_TOP_TIMESTEP_SCREEN_CONTRACT_20261004.md'
CONTROLLER = ROOT/'scripts/screen_visible_top_timestep.py'
PASS_PROCESS = 'PASS_LARGER_TIMESTEP_REFERENCE_SCREEN_NOT_RESTARTABLE'
FAIL_PROCESS = 'FAILED_NONRESTARTABLE_TIMESTEP_SCREEN'
require, clean, digest = old.require, old.clean, old.digest
positions, pose, Reader = old.positions, old.pose, old.Reader
project_path, validate_state = reference.project_path, reference.validate_state
correction_audit = reference.correction_audit
QUALIFICATION, ALLOWED_FEATURES = reference.QUALIFICATION, reference.ALLOWED_FEATURES


def candidate_audit(reader, dt, directory, row, previous, previous_state, context, frozen):
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


def trial_context(reader, meta, multiplier):
    refaudit = REFERENCE/'readback/summary.json'
    reader.bind(refaudit, meta['reference_readback_sha256'])
    audit = reader.read(refaudit)
    require(project_path(meta['reference_readback']) == refaudit and
            audit['status'] == reference.PASS_READBACK and audit['window_completed'] and
            audit['top_event_process_screen_passed'] and audit['added_accepted_steps'] == 32,
            'actual 62.5us reference was not independently qualified')
    reader.sources(audit['evidence_sha256'])
    context = reference.parent_context(reader, REFERENCE, reader.read(REFERENCE/'coarse/manifest.json'))
    require(project_path(meta['parent_checkpoint']) == context['parent_path'] and
            meta['parent_sha256'] == context['parent_sha256'] == audit['parent_sha256'] and
            meta['start_time_s'] == context['parent']['time_s'] and
            meta['start_global_step'] == context['parent']['step'], 'diagnostic true parent differs')
    folder = REFERENCE/'fine'
    manifest = reader.read(folder/'manifest.json')
    history = reader.read(folder/'accepted_history.json')
    rows = reader.read(folder/'candidate_history.json')
    require(manifest['status'] == reference.PASS_FINE and len(rows) == len(history) == 32,
            '62.5us reference is incomplete')
    refs, states, direct = {}, [], context['parent_sha256']
    for number, (row, entry) in enumerate(zip(rows, history), 1):
        path = folder/entry['checkpoint']
        reader.bind(path, entry['sha256'])
        pa, state = reader.checkpoint(path)
        validate_state(pa, state, context, number, reference.FINE_DT, True)
        require(entry['step'] == row['step'] == number and entry['parent_sha256'] == direct and
                state['continuation_parent_sha256'] == direct and
                entry['dt_s'] == reference.FINE_DT and entry['time_s'] == state['time_s'] and
                entry['sha256'] == audit['accepted_history'][number-1]['sha256'],
                'reference accepted chain differs from independent audit')
        direct = entry['sha256']; refs[path.relative_to(ROOT).as_posix()] = direct
        for name, value in row['artifacts'].items():
            path = (folder/name).resolve()
            require(path.is_relative_to(folder.resolve()), 'reference artifact escaped registered branch')
            reader.bind(path, value); refs[path.relative_to(ROOT).as_posix()] = value
        states.append((pa, state))
    for path in (refaudit, folder/'manifest.json', folder/'accepted_history.json', folder/'candidate_history.json'):
        reader.bind(path); refs[path.relative_to(ROOT).as_posix()] = digest(path)
    require(meta['reference_sha256'] == refs, 'complete frozen reference set differs')
    frozen = dict(context['frozen'])
    for path in (CONTROLLER, CONTRACT, refaudit):
        reader.bind(path); frozen[path.relative_to(ROOT).as_posix()] = digest(path)
    require(meta['source_sha256'] == frozen, 'diagnostic source collection changed')
    reader.sources(frozen)
    count, dt = 32//multiplier, reference.FINE_DT*multiplier
    cost = max(float(manifest['wall_seconds'])/32, max(row['wall_seconds'] for row in rows))*multiplier
    require(meta['candidate_count'] == count and meta['dt_s'] == dt and
            abs(meta['forecast_wall_s']-(30+1.5*count*cost)) < 1e-10 and
            meta['forecast_wall_s'] <= 450, 'registered step count, dt or complete-stage forecast differs')
    context.update(frozen=frozen, reference_rows=rows, reference_history=history,
        reference_states=states, count=count, dt=dt, reference_audit_sha256=digest(refaudit))
    return context


def replay_screen(reader, directory, item, context, number, multiplier, seeds):
    row, state = item['row'], item['state']
    screen = reader.read(directory/f'step_{number:03d}/temporal_screen.json')
    first, last = (number-1)*multiplier, number*multiplier
    refpa, refstate = context['reference_states'][last-1]
    entry = context['reference_history'][last-1]
    refpath = REFERENCE/'fine'/entry['checkpoint']
    require(screen['step'] == number and screen['time_s'] == state['time_s'] and
            abs(state['time_s']-refstate['time_s']) < 1e-14 and
            screen['candidate_sha256'] == row['sha256'] and
            project_path(screen['reference_checkpoint']) == refpath and
            screen['reference_sha256'] == entry['sha256'] and
            screen['reference_fine_steps'] == list(range(first+1, last+1)),
            'shared-endpoint screen is not bound to its complete reference interval')
    seeds.update(e['particle_id'] for e in row['events'])
    seeds.update(row['blocked_particle_ids']); seeds.update(state['bore_observation']['crossing_ids'])
    graphs = [item['graph']]
    for refrow in context['reference_rows'][first:last]:
        seeds.update(e['particle_id'] for e in refrow['events']); seeds.update(refrow['blocked_particle_ids'])
        graphs.append(reader.npz(REFERENCE/'fine'/f"step_{refrow['step']:03d}"/'visible_graph.npz'))
    seeds.update(refstate['bore_observation']['crossing_ids'])
    require(screen['cumulative_event_and_blocked_ids'] == sorted(seeds),
            'event, blocked-edge or crossed-particle cumulative seeds were dropped')
    selected = set(seeds)
    for graph in graphs:
        ids = graph['particle_id']
        for left, right in (('i', 'j'), ('raw_i', 'raw_j')):
            x, y = ids[graph[left]], ids[graph[right]]
            mask = np.isin(x, list(seeds)) | np.isin(y, list(seeds))
            selected.update(np.r_[x[mask], y[mask]].tolist())
    # Reference is deliberately first: existing fine-state stress/ep norms and
    # floors define the registered denominators, including for failed screens.
    glob, local = old.compare_arrays(refpa, item['pa'], sorted(selected), context['p'])
    match = state['bore_observation']['crossing_ids'] == refstate['bore_observation']['crossing_ids']
    require(clean(glob) == screen['global_screen'] and clean(local) == screen['local_screen'] and
            screen['crossing_ids_match'] == match, 'saved screen does not replay with the actual fine-reference denominators')
    diagnostics = []
    for value in local['rows']:
        if value['stress_relative'] > .1 or value['plastic_relative'] > .2:
            pid = value['particle_id']; k = int(np.flatnonzero(refpa.particle_id == pid)[0])
            refnorm = float(np.linalg.norm(refpa.S.reshape(-1, 9)[k]))
            trialnorm = float(np.linalg.norm(item['pa'].S.reshape(-1, 9)[k]))
            stressdiff = float(np.linalg.norm((item['pa'].S-refpa.S).reshape(-1, 9)[k]))
            diagnostics.append(dict(**value, reference_stress_norm_Pa=refnorm,
                candidate_stress_norm_Pa=trialnorm, stress_difference_norm_Pa=stressdiff,
                stress_denominator_Pa=max(refnorm, 1e6), stress_floor_active=refnorm < 1e6))
    return dict(step=number, time_s=state['time_s'], passed=glob['pass'] and local['passed'] and match,
        screen_sha256=digest(directory/f'step_{number:03d}/temporal_screen.json'),
        global_screen=glob, local_screen=local, local_particle_count=len(selected),
        cumulative_event_and_blocked_ids=sorted(seeds), crossing_ids_match=match,
        local_failed_stress_or_plastic_rows=diagnostics,
        stress_and_plastic_denominator='ACTUAL_ACCEPTED_62_5US_FINE_REFERENCE'), screen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--multiplier', type=int, choices=(4, 8, 16), required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    folder, output = BASE/f'mult_{args.multiplier:02d}', project_path(args.output)
    require(folder.is_dir() and not output.exists() and not output.is_relative_to(folder),
            'existing trial and separate new readback output required')
    output.mkdir(parents=True, exist_ok=False)
    begin = time.perf_counter(); reader = Reader()
    report = dict(status='READBACK_RUNNING', command=[sys.executable, *sys.argv], multiplier=args.multiplier,
        physical_solver_evaluations=0, nnps_evaluations=0, geometry_certificate_evaluations=0,
        accepted_steps=0, accepted_time_advanced_s=0., restart_allowed=False, deposition_validated=False)
    try:
        for path in (Path(__file__), Path(reference.__file__), Path(old.__file__)):
            reader.bind(path)
        meta = reader.read(folder/'manifest.json')
        supervision = reader.read(project_path(meta['supervision_receipt']))
        require(meta['status'] in (PASS_PROCESS, FAIL_PROCESS) and
                meta['classification'] == CLASS and meta['diagnostic_only'] is True and
                meta['restart_allowed'] is False and meta['accepted_steps'] == 0 and
                meta['accepted_time_advanced_s'] == 0. and meta['multiplier'] == args.multiplier and
                supervision['exit_code'] == meta['exit_code'] and
                (meta['exit_code'] == 0) == (meta['status'] == PASS_PROCESS), 'trial terminal state or diagnostic scope differs')
        require(meta['wall_seconds'] <= 540.5 and supervision['limit_s'] == 600 and
                supervision['external_wall_s'] <= 605, 'trial exceeded registered runtime budget')
        require(not list(folder.rglob('accepted.npz')), 'diagnostic branch improperly published restart state')
        context = trial_context(reader, meta, args.multiplier)
        rows = reader.read(folder/'candidate_history.json') if (folder/'candidate_history.json').exists() else []
        require(len(rows) <= context['count'], 'diagnostic exceeded registered step count')
        previous, previous_state, direct = context['initial'], context['parent'], context['parent_sha256']
        items, screens, rawscreens, seeds = [], [], [], set(context['seeds'])
        for number, row in enumerate(rows, 1):
            require(row['step'] == number and row['parent_sha256'] == direct,
                    'diagnostic candidate chain differs from its true immediate parent')
            item = candidate_audit(reader, context['dt'], folder, row, previous, previous_state,
                                   context, context['frozen'])
            items.append(item)
            previous, previous_state, direct = item['pa'], item['state'], row['sha256']
            screenpath = folder/f'step_{number:03d}/temporal_screen.json'
            if screenpath.exists():
                require(number == len(screens)+1 and all(s['passed'] for s in screens),
                        'trial proceeded past a missing or failed endpoint screen')
                screen, raw = replay_screen(reader, folder, item, context, number, args.multiplier, seeds)
                screens.append(screen); rawscreens.append(raw)
            else:
                require(meta['status'] == FAIL_PROCESS and number == len(rows),
                        'only terminal failed candidate may lack a temporal screen')
        extras = []
        known = {f"step_{row['step']:03d}" for row in rows}
        for directory in sorted(folder.glob('step_*')):
            if directory.name not in known:
                files = [reader.bind(path).relative_to(ROOT).as_posix() for path in sorted(directory.rglob('*')) if path.is_file()]
                extras.append(dict(directory=directory.relative_to(ROOT).as_posix(), files=files,
                    qualification='FAILED_OR_PARTIAL_CANDIDATE_NOT_PHYSICALLY_QUALIFIED_NOT_RESTARTABLE'))
        complete = len(items) == len(screens) == context['count'] and all(s['passed'] for s in screens)
        if meta['status'] == PASS_PROCESS:
            require(complete and not extras and meta['completed_candidates'] == len(rows) and
                    clean(meta['screens']) == clean(rawscreens) and
                    clean(meta['interval_energy']) == clean(previous_state['interval_energy']),
                    'successful trial scope or saved ledger is incomplete')
            top_ids = sorted({e['particle_id'] for row in rows for e in row['events'] if e['feature'] == 'TOP'})
            require(top_ids and meta['top_particle_ids'] == top_ids and meta['source_unchanged'] is True and
                    meta['comparisons_are_against_62_5us'] is True, 'successful trial lacks its actual TOP event or source binding')
        else:
            require(not complete, 'failed manifest conflicts with a completely passing registered trial')
        if any(not s['passed'] for s in screens):
            require(meta['status'] == FAIL_PROCESS and
                    meta['error'] == "ValueError('larger physical timestep disagrees with reference')",
                    'failed endpoint comparison was not retained as a failed trial')
        reader.finish()
        report.update(status='PASS_SAVED_TIMESTEP_EVIDENCE_CONSISTENCY_NOT_PROCESS_VALIDATION',
            original_process_status=meta['status'], original_process_error=meta.get('error'),
            diagnostic_screen_passed=meta['status'] == PASS_PROCESS and complete,
            diagnostic_only=True, candidate_count=len(items), completed_screen_count=len(screens),
            actual_candidate_end_s=previous_state['time_s'], original_accepted_parent_time_s=context['parent']['time_s'],
            parent_sha256=context['parent_sha256'], reference_readback_sha256=context['reference_audit_sha256'],
            dt_s=context['dt'], registered_candidate_count=context['count'],
            candidates=[item['summary'] for item in items], screens=screens, unregistered_partial_evidence=extras,
            inherited_energy_prefix_records=len(context['prefix']),
            final_candidate_interval_energy={k:v for k,v in previous_state['interval_energy'].items() if k != 'records'},
            physical_step_wall_seconds=meta['wall_seconds'], forecast_wall_seconds=meta['forecast_wall_s'],
            contact_history_convergence='NOT_ESTABLISHED', evidence_sha256=reader.evidence,
            source_entry_checks_including_repeats=reader.source_entries)
        code = 0
    except Exception:
        report.update(status='FAIL_SAVED_TIMESTEP_EVIDENCE_READBACK', error=traceback.format_exc(), evidence_sha256=reader.evidence)
        code = 2
    report.update(exit_code=code, wall_seconds=time.perf_counter()-begin)
    (output/'summary.json').write_text(json.dumps(clean(report), indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({key:report.get(key) for key in ('status', 'exit_code', 'original_process_status',
        'diagnostic_screen_passed', 'candidate_count', 'dt_s', 'wall_seconds', 'error')}, indent=2))
    return code


if __name__ == '__main__':
    raise SystemExit(main())

