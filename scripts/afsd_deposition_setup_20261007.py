"""Independent rectangular fill-then-track setup on unchanged PySPH operators.

Only a new initial domain, finite inventory and machine program are supplied.
No checkpoint is imported or cropped, no free-material path is assigned, and
no command-plan audit or detached candidate is declared to be a deposition.
"""
from pathlib import Path
import copy
import json
import math
import numpy as np

import afsd_process_framework_20261007 as parent
from afsd_process_motion_20261007 import ProcessProgram
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.interval_acceptance import (
    initial_interval_energy_state, initial_bore_crossing_state,
    validate_interval_energy_state, validate_bore_crossing_state,
    require_interval_energy, observe_square_bore_crossings)

ROOT = parent.ROOT
CLASSIFICATION = 'INDEPENDENT_RECTANGULAR_FILL_THEN_TRACK_AFSD_APPROXIMATE_V1'
DEFAULT_CONFIG = ROOT/'configs/deposition_fill_then_track_20261007.json'
CONFIG_CLASSIFICATION = 'ENGINEERING_APPROXIMATE_RECTANGULAR_FILL_THEN_TRACK_AFSD'
GEOMETRY_KIND = 'INDEPENDENT_RECTANGULAR_SUBSTRATE_AND_FINITE_SQUARE_BAR'
require = parent.require
digest = parent.digest
identity_digest = parent.identity_digest
old = parent.old


def _finite(value, name):
    value = float(value)
    require(math.isfinite(value), 'finite '+name+' required')
    return value


def _cells(length, dx, name):
    length = _finite(length, name)
    count = round(length/dx)
    require(length > 0 and count >= 1 and abs(count*dx-length) <= 1e-12,
            'positive full particle cells required: '+name)
    return count


def geometry_and_inventory_audit(p, program, geometry, intent):
    """Explicit rectangular command budget; necessary allowances, not physics."""
    require(geometry.get('kind') == GEOMETRY_KIND and
            geometry.get('status') == 'DECLARED_ENGINEERING_ASSUMPTION' and
            bool(geometry.get('source')) and geometry.get('substrate_width_compatibility') ==
            'Y_EXTENT_ONLY_NOT_A_SQUARE_DOMAIN', 'explicit sourced rectangular geometry required')
    positive = ('dx', 'h_ratio', 'rho0', 'substrate_depth', 'substrate_width', 'width',
                'bar_length', 'bore_width', 'tool_radius', 'layer', 'actuator_depth', 'tool_top')
    require(all(_finite(p[key], key) > 0 for key in positive), 'positive dimensions required')
    bounds = np.array([geometry['substrate_'+axis+'_bounds_m'] for axis in 'xyz'], float)
    require(bounds.shape == (3, 2) and np.isfinite(bounds).all() and
            np.all(bounds[:, 1] > bounds[:, 0]), 'finite ordered rectangular bounds required')
    require(np.array_equal(np.asarray(geometry['bar_axis_xy_m'], float), [0., 0.]) and
            geometry['bar_bottom_z_m'] == 0., 'registered bar axis and bottom are world zero')
    require(np.array_equal(bounds[2], [-p['substrate_depth'], 0.]) and
            abs(np.diff(bounds[1])[0]-p['substrate_width']) <= 1e-12,
            'substrate z/y compatibility fields differ from actual rectangle')
    require(program.layer == p['layer'], 'motion and geometry shoulder heights differ')
    require(p['width'] <= p['bore_width'] and
            p['tool_radius'] > p['bore_width']/math.sqrt(2) and p['tool_top'] > p['layer'],
            'finite annular shoulder with an open square bore required')
    nx, ny, nz = [_cells(b-a, p['dx'], 'substrate_'+axis)
                  for axis, (a, b) in zip('xyz', bounds)]
    nbxy = _cells(p['width'], p['dx'], 'bar width')
    nbz = _cells(p['bar_length'], p['dx'], 'bar length')
    ncap = _cells(p['actuator_depth'], p['dx'], 'actuator depth')
    require(p['bar_length'] > p['layer']+p['actuator_depth']+p['dx']/2,
            'positive finite axial inventory with preserved actuator cap required')
    angle, feed, distance = program.integral(program.end_s)
    fill_end = _finite(intent['stationary_feed_end_s'], 'stationary feeding endpoint')
    require(0 < fill_end < program.end_s and
            any(fill_end == phase['end'] for phase in program.phases),
            'stationary feeding must end at an explicit phase boundary')
    fill_angle, fill_feed, fill_distance = program.integral(fill_end)
    require(fill_distance == 0. and all(phase['target'][2] == 0.
            for phase in program.phases if phase['end'] <= fill_end),
            'stationary feeding phase must not traverse')
    require(intent.get('traverse_requires_observed_spreading') is True and
            intent.get('success_is_not_command_completion') is True and
            bool(intent.get('required_observation')) and bool(intent.get('insufficient_spreading_action')),
            'explicit observation-dependent process intent required')
    require(abs(distance-_finite(intent['target_total_traverse_m'], 'target travel')) <= 1e-12 and
            abs(fill_feed-_finite(intent['planned_axial_feed_before_traverse_m'], 'fill feed')) <= 1e-12,
            'declared filling/travel budgets differ from integrated commands')
    radius, buffer = p['tool_radius'], 2*p['h_ratio']*p['dx']
    margins = dict(initial_x_minus=-bounds[0, 0]-radius,
                   final_x_plus=bounds[0, 1]-distance-radius,
                   y_minus=-bounds[1, 0]-radius, y_plus=bounds[1, 1]-radius)
    usable_feed = p['bar_length']-p['actuator_depth']-p['layer']-p['dx']/2
    failures = []
    if min(margins.values()) < buffer-1e-14:
        failures.append('RECTANGLE_DOES_NOT_COVER_FULL_COMMAND_PATH_PLUS_KERNEL_BUFFER')
    if feed > usable_feed+1e-14:
        failures.append('FINITE_FEEDSTOCK_ACTUATOR_RESERVE_EXHAUSTED')
    nsub, nbar = nx*ny*nz, nbxy*nbxy*nbz
    area = p['width']**2
    return dict(status='PASS_RECTANGULAR_DOMAIN_AND_FINITE_INVENTORY_PLAN' if not failures else
                'FAIL_RECTANGULAR_DOMAIN_OR_INVENTORY_PLAN', passed=not failures, failures=failures,
        geometry_kind=GEOMETRY_KIND, actual_substrate_bounds_m=bounds.tolist(),
        substrate_width_compatibility='Y_EXTENT_ONLY_NOT_A_SQUARE_DOMAIN',
        particle_count=nsub+nbar, substrate_count=nsub, bar_count=nbar,
        substrate_cell_counts=[nx, ny, nz], bar_cell_counts=[nbxy, nbxy, nbz],
        actuator_particle_count=nbxy*nbxy*ncap, cell_volume_m3=p['dx']**3,
        initial_total_mass_kg=(nsub+nbar)*p['rho0']*p['dx']**3,
        initial_substrate_volume_m3=nsub*p['dx']**3, initial_bar_volume_m3=nbar*p['dx']**3,
        programmed_duration_s=program.end_s, programmed_feed_m=float(feed),
        programmed_traverse_m=float(distance), programmed_revolutions=float(angle/(2*math.pi)),
        required_edge_buffer_m=buffer, edge_margins_m={k: float(v) for k, v in margins.items()},
        final_forward_edge_margin_m=float(margins['final_x_plus']),
        lateral_edge_margin_m=float(min(margins['y_minus'], margins['y_plus'])),
        usable_axial_command_m=usable_feed, remaining_axial_command_margin_m=float(usable_feed-feed),
        nominal_command_swept_volume_m3=float(area*feed),
        stationary_feed_end_s=fill_end, nominal_stationary_feed_m=float(fill_feed),
        nominal_stationary_swept_volume_m3=float(area*fill_feed),
        stationary_tool_revolutions=float(fill_angle/(2*math.pi)),
        total_rotation_source='INTEGRATED_MACHINE_COMMAND_NOT_MATERIAL_ROTATION',
        stationary_fill_status='PLANNED_COMMAND_ONLY_NOT_PROVEN_FILLED',
        replenishment='FINITE_BAR_NO_AUTOMATIC_INSERTION', domain_convergence='NOT_TESTED',
        deposition_status='NOT_PREDICTED_BY_GEOMETRY_OR_COMMAND_INTEGRALS',
        free_material_transport='NOT_PREDICTED_BY_COMMAND_INTEGRALS')


def _resolve_configuration(config):
    require(config['schema_version'] == 1 and config['classification'] == CONFIG_CLASSIFICATION,
            'unsupported independent deposition configuration')
    base = (ROOT/config['base_parameters']).resolve()
    require(base.is_relative_to(ROOT) and old.expanded.sha(base) == config['base_parameters_sha256'],
            'registered base parameters changed')
    p = json.loads(base.read_text(encoding='utf-8'))['parameters']
    for key, row in config['parameter_overrides'].items():
        require(key in ('substrate_width', 'bar_length', 'bore_width') and row['unit'] == 'm' and
                row['status'] == 'DECLARED_ENGINEERING_ASSUMPTION' and bool(row['source']),
                'unregistered deposition geometry override')
        p[key] = _finite(row['value'], key)
    require(config['thermal_initialization'] == 'INHERITED_450C_HOT_START_NOT_COLD_START' and
            config['interface_initialization'] == 'INHERITED_PREBONDED_INITIAL_LOWER_BAR_NOT_NEW_DEPOSITION' and
            config['feed_mode'] == 'FINITE_MATERIAL_LABELLED_TOP_AXIAL_ACTUATOR' and
            config['bar_rotation'] == 'FREE_XY_CONTACT_DRIVEN_BY_ROTATING_SQUARE_BORE',
            'unsupported initial thermomechanical state or material driver')
    program = ProcessProgram(config['phases'], p['layer'])
    audit = geometry_and_inventory_audit(p, program, config['geometry'], config['process_intent'])
    require(audit['passed'], ', '.join(audit['failures']))
    return p, program, audit


def load_configuration(path=DEFAULT_CONFIG):
    config = json.loads(Path(path).read_text(encoding='utf-8'))
    p, program, audit = _resolve_configuration(config)
    return config, p, program, audit


def source_hashes(config_path=DEFAULT_CONFIG):
    result = parent.source_hashes()
    paths = [Path(__file__).resolve(), Path(config_path).resolve(),
             ROOT/'tests/review/test_afsd_deposition_setup_20261007.py']
    for path in paths:
        require(path.is_relative_to(ROOT), 'project-local new deposition sources required')
        result[path.relative_to(ROOT).as_posix()] = old.expanded.sha(path)
    return result


def _rectangular_particles(p, geometry, audit):
    def centers(bounds, count):
        return bounds[0]+(np.arange(count)+.5)*p['dx']
    bounds = audit['actual_substrate_bounds_m']
    axes = [centers(b, n) for b, n in zip(bounds, audit['substrate_cell_counts'])]
    substrate = np.array(np.meshgrid(*axes, indexing='ij')).reshape(3, -1).T
    barxy = centers([-p['width']/2, p['width']/2], audit['bar_cell_counts'][0])
    barz = centers([0., p['bar_length']], audit['bar_cell_counts'][2])
    bar = np.array(np.meshgrid(barxy, barxy, barz, indexing='ij')).reshape(3, -1).T
    positions = np.vstack((substrate, bar))
    component = np.r_[np.zeros(len(substrate), dtype=int), np.ones(len(bar), dtype=int)]
    pa = make_array(positions, component, p)
    # The shared factory also serves old moving-actuator initializations. This
    # independent smooth-entry state starts every material particle at rest.
    pa.u[:] = 0.; pa.v[:] = 0.; pa.w[:] = 0.
    return pa


def _geometry_readback(pa, p, audit):
    points = old.positions(pa)
    expected_mass = p['rho0']*p['dx']**3
    require(len(points) == audit['particle_count'] and np.isfinite(points).all() and
            len(np.unique(points, axis=0)) == len(points), 'finite unique initial particle cells required')
    require(np.array_equal(points, np.column_stack((pa.x0, pa.y0, pa.z0))) and
            np.array_equal(pa.particle_id, np.arange(1, len(points)+1)) and
            np.all(pa.m == expected_mass), 'initial reference coordinates, IDs or physical mass differ')
    substrate, bar = pa.component_id == 0, pa.component_id == 1
    actual_bounds = np.column_stack((points[substrate].min(axis=0)-p['dx']/2,
                                    points[substrate].max(axis=0)+p['dx']/2))
    require(np.allclose(actual_bounds, audit['actual_substrate_bounds_m'], atol=1e-15, rtol=0.),
            'allocated substrate does not cover the explicit rectangle')
    require(int(substrate.sum()) == audit['substrate_count'] and int(bar.sum()) == audit['bar_count'] and
            int(pa.actuator.sum()) == audit['actuator_particle_count'], 'component/cap allocation differs')
    return dict(status='PASS_INITIAL_GEOMETRY_IDENTITY_ONLY_NOT_MECHANICAL_QUALIFICATION',
        actual_substrate_bounds_m=actual_bounds.tolist(), actual_bar_cell_bounds_m=
        np.column_stack((points[bar].min(axis=0)-p['dx']/2, points[bar].max(axis=0)+p['dx']/2)).tolist(),
        particle_count=len(points), substrate_count=int(substrate.sum()), bar_count=int(bar.sum()),
        actuator_particle_count=int(pa.actuator.sum()), total_mass_kg=float(pa.m.sum()),
        initially_below_shoulder_bar_count=int(np.count_nonzero(bar & (pa.z0 < p['layer']))),
        initial_material_at_rest=bool(np.all(pa.u == 0.) and np.all(pa.v == 0.) and np.all(pa.w == 0.)),
        old_checkpoint_read=False, material_trimmed_or_repositioned=False,
        immutable_material_identity_sha256=identity_digest(pa))


def independent_initial_state(config_path=DEFAULT_CONFIG):
    config, p, program, audit = load_configuration(config_path)
    backend = old.expanded.backend_provenance()
    pa = _rectangular_particles(p, config['geometry'], audit)
    geometry_audit = _geometry_readback(pa, p, audit)
    bottom = (pa.component_id == 0) & (pa.z0 < -p['substrate_depth']+p['dx'])
    free = np.ones((len(pa.x), 3), dtype=bool)
    free[bottom] = False; free[pa.actuator.astype(bool), 2] = False
    state = dict(time_s=0., step=0, classification=CLASSIFICATION,
        contact_history=np.zeros((len(pa.x), 8, 3)), free=free,
        interval_energy=initial_interval_energy_state(), parameters=copy.deepcopy(p),
        framework_config=copy.deepcopy(config), framework_config_sha256=digest(config),
        process_program=program.to_dict(), machine_history=[], domain_inventory_plan=audit,
        geometry_definition=copy.deepcopy(config['geometry']),
        geometry_definition_sha256=digest(config['geometry']), independent_geometry_audit=geometry_audit,
        independent_geometry_audit_sha256=digest(geometry_audit),
        independent_initialization=True, parent_checkpoint=None, accepted=False, restart_allowed=False,
        diagnostic_only=True, physical_qualification='NOT_RUN', accepted_time_advanced_s=0.,
        source_sha256=source_hashes(config_path), backend=backend,
        initial_mass_kg=float(pa.m.sum()), initial_particle_ids=pa.particle_id.astype(int).tolist(),
        immutable_material_identity_sha256=identity_digest(pa), crossing_ids=[], crossing_mass_kg=0.)
    state['bore_observation'] = initial_bore_crossing_state(old.positions(pa), pa.particle_id,
        pa.component_id, pa.m, 0., program.global_pose(0.), p['bore_width'])
    validate_state(state)
    return pa, state


def validate_state(state):
    require(state.get('classification') == CLASSIFICATION and state.get('accepted') is False and
            state.get('restart_allowed') is False and state.get('diagnostic_only') is True and
            state.get('independent_initialization') is True and state.get('parent_checkpoint') is None,
            'new independent rectangular state required; old checkpoints cannot be relabelled')
    config = state['framework_config']
    require(digest(config) == state['framework_config_sha256'], 'configuration history changed')
    p, program, audit = _resolve_configuration(config)
    require(p == state['parameters'] and audit == state['domain_inventory_plan'] and
            program.to_dict() == state['process_program'], 'parameters, geometry budget or machine program changed')
    require(state['geometry_definition'] == config['geometry'] and
            digest(state['geometry_definition']) == state['geometry_definition_sha256'],
            'explicit rectangular geometry history changed')
    ga = state['independent_geometry_audit']
    require(digest(ga) == state['independent_geometry_audit_sha256'] and
            ga['immutable_material_identity_sha256'] == state['immutable_material_identity_sha256'] and
            ga['particle_count'] == audit['particle_count'] and ga['total_mass_kg'] == state['initial_mass_kg'] and
            ga['old_checkpoint_read'] is False and ga['material_trimmed_or_repositioned'] is False,
            'independent geometry/identity audit changed')
    validate_interval_energy_state(state['interval_energy'])
    validate_bore_crossing_state(state['bore_observation'])
    require(state['bore_observation']['initial_time_s'] == 0. and
            state['bore_observation']['bore_width_m'] == p['bore_width'] and
            state['bore_observation']['particle_ids'] == sorted(state['initial_particle_ids']),
            'initial identity or bore-observer origin changed')
    rows = state['machine_history']; t = 0.
    require(len(rows) == state['step'] == state['interval_energy']['steps'], 'lost machine/energy prefix')
    for row in rows:
        require(row == program.step(t, row['dt_s']).to_dict(), 'machine interval history changed')
        t = row['t1_s']
    require(t == state['time_s'] == state['bore_observation']['last_time_s'], 'state times differ')
    require(state['crossing_ids'] == state['bore_observation']['crossing_ids'] and
            state['crossing_mass_kg'] == state['bore_observation']['crossing_mass_kg'],
            'passage summary differs from persistent history')
    return program


class ProcessTrial(parent.ProcessTrial):
    """Unchanged residual, tangent and material-cap prescription on a rectangle."""
    def __init__(self, pa, state, dt):
        program = validate_state(state)
        self.pa, self.state, self.p = pa, state, dict(state['parameters'])
        p = self.p
        require(identity_digest(pa) == state['immutable_material_identity_sha256'] and
                set(pa.particle_id.astype(int)) == set(state['initial_particle_ids']) and
                abs(pa.m.sum()/state['initial_mass_kg']-1.) <= 1e-12,
                'finite material identity or inventory changed')
        order = np.argsort(pa.particle_id); bore = state['bore_observation']
        require(pa.particle_id[order].astype(int).tolist() == bore['particle_ids'] and
                pa.component_id[order].astype(int).tolist() == bore['component_ids'] and
                pa.m[order].tolist() == bore['masses_kg'], 'actual per-ID inventory differs from observer')
        self.x0 = old.positions(pa); self.s0 = pa.S.reshape(-1, 3, 3).copy()
        self.bulk = old.bulk_modulus(p['shear'], p['poisson'])
        self.law = old.FKS(**{key: p[key] for key in old.FKS.__dataclass_fields__})
        self.motion = program.step(state['time_s'], dt)
        self.contact = old.expanded.BatchedChainRuleContact(tool_radius=p['tool_radius'],
            bore_width=p['bore_width'], tool_height=p['tool_top']-p['layer'],
            particle_radius=p['dx']/2, normal_stiffness=p['normal_stiffness'],
            tangential_stiffness=p['normal_stiffness'], friction=p['friction'], base_z=-p['substrate_depth'])
        bottom = (pa.component_id == 0) & (pa.z0 < -p['substrate_depth']+p['dx'])
        cap = pa.actuator.astype(bool)
        require(np.array_equal(cap, (pa.component_id == 1) &
                (pa.z0 > p['bar_length']-p['actuator_depth'])), 'material-labelled actuator cap changed')
        expected_free = np.ones((len(pa.x), 3), dtype=bool)
        expected_free[bottom] = False; expected_free[cap, 2] = False
        require(np.array_equal(state['free'], expected_free), 'only bottom xyz and original top-cap z may be fixed')


def advance_candidate(trial, result):
    """Retain fresh-map, original force/heat and full-prefix energy/passages gates."""
    require(isinstance(trial, ProcessTrial), 'explicit independent rectangular trial required')
    trial._require_interval(result['dt_s']); validate_state(trial.state)
    fresh = trial.evaluate(result['displacement'], result['dt_s'])
    require(np.array_equal(result['i'], fresh['i']) and np.array_equal(result['j'], fresh['j']) and
            result['graph_fingerprint'] == fresh['graph_fingerprint'], 'fresh common graph differs')
    energy = old.original_qs.energy_diagnostics(trial, fresh)
    interval = require_interval_energy(trial.state['interval_energy'], energy)
    pa, state = old.original_qs.accepted_copy(trial, fresh)
    observed = observe_square_bore_crossings(trial.state['bore_observation'], trial.x0,
        old.positions(pa), pa.particle_id, pa.component_id, pa.m,
        trial.motion.t0, trial.motion.t1, trial.motion.pose)
    state.update(interval_energy=interval, bore_observation=observed['next_state'],
        crossing_ids=observed['crossing_ids'], crossing_mass_kg=observed['crossing_mass_kg'],
        machine_history=trial.state['machine_history']+[trial.motion.to_dict()],
        accepted=False, restart_allowed=False, diagnostic_only=True, accepted_time_advanced_s=0.,
        physical_qualification='DETACHED_CANDIDATE_REQUIRES_PATH_TIME_AND_DEPOSITION_SCREEN')
    require(state['interval_energy']['records'][:-1] == trial.state['interval_energy']['records'],
            'full energy history changed')
    require(identity_digest(pa) == trial.state['immutable_material_identity_sha256'],
            'candidate altered immutable material identity')
    validate_state(state)
    return pa, state, fresh


def operator_reuse_audit():
    inherited = parent.operator_reuse_audit()
    require(ProcessTrial.evaluate is parent.ProcessTrial.evaluate and
            ProcessTrial.tangent is parent.ProcessTrial.tangent and
            ProcessTrial.prescribed is parent.ProcessTrial.prescribed and
            ProcessTrial._require_interval is parent.ProcessTrial._require_interval,
            'frozen residual, tangent or axial constraint implementation changed')
    return dict(status='PASS_UNCHANGED_INHERITED_THERMOMECHANICAL_OPERATORS',
        inherited_audit=inherited, changed='independent rectangular geometry, finite inventory and machine program',
        runtime_monkeypatch=False, material_trajectory_prescribed=False)
