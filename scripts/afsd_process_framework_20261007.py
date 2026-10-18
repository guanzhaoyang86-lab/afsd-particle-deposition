"""New hot-process framework on the frozen PySPH thermomechanical operators.

No old checkpoint migration or runtime patching. This module deliberately
inherits the reviewed ramp trial's evaluate/tangent unchanged and supplies an
independent phase program, domain, state protocol, and shared interval motion.
Detached candidates are not accepted physical process results.
"""
from pathlib import Path
import copy
import hashlib
import json
import numpy as np
import ramped_hot_trial as old
from afsd_process_motion_20261007 import ProcessProgram, ProcessInterval, domain_and_inventory_audit
from afsd_backend.pysph import discrete_process
from afsd_backend.pysph.interval_acceptance import (initial_bore_crossing_state,
    validate_interval_energy_state, validate_bore_crossing_state, require_interval_energy,
    observe_square_bore_crossings)

ROOT = old.ROOT
CLASSIFICATION = 'INDEPENDENT_PHASED_HOT_AFSD_FRAMEWORK_APPROXIMATE_V1'
DEFAULT_CONFIG = ROOT/'configs/process_framework_20261007.json'
require = old.require


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def identity_digest(pa):
    ids = np.asarray(pa.particle_id)
    require(ids.ndim == 1 and ids.dtype.kind in 'iu' and len(np.unique(ids)) == len(ids),
            'unique persistent particle IDs required')
    order = np.argsort(ids)
    keys = ('particle_id', 'component_id', 'material_id', 'tracer_type', 'wholebar_tag',
            'actuator', 'x0', 'y0', 'z0', 'm', 'h')
    return digest({key: np.asarray(getattr(pa, key))[order].tolist() for key in keys})


def load_configuration(path=DEFAULT_CONFIG):
    path = Path(path).resolve()
    config = json.loads(path.read_text(encoding='utf-8'))
    require(config['schema_version'] == 1 and config['classification'] ==
            'ENGINEERING_APPROXIMATE_HOT_PROCESS_FRAMEWORK', 'unsupported framework configuration')
    base = (ROOT/config['base_parameters']).resolve()
    require(base.is_relative_to(ROOT) and old.expanded.sha(base) == config['base_parameters_sha256'],
            'base parameters changed')
    p = json.loads(base.read_text(encoding='utf-8'))['parameters']
    for key, row in config['parameter_overrides'].items():
        require(key in ('substrate_width', 'bore_width') and row['unit'] == 'm' and
            row['status'] == 'DECLARED_ENGINEERING_ASSUMPTION' and row['source'],
            'unregistered geometric override')
        p[key] = float(row['value'])
    require(config['thermal_initialization'] == 'INHERITED_450C_HOT_START_NOT_COLD_START' and
        config['interface_initialization'] == 'INHERITED_PREBONDED_INITIAL_LOWER_BAR_NOT_NEW_DEPOSITION' and
        config['feed_mode'] == 'FINITE_MATERIAL_LABELLED_TOP_AXIAL_ACTUATOR' and
        config['bar_rotation'] == 'FREE_XY_CONTACT_DRIVEN_BY_ROTATING_SQUARE_BORE',
        'unsupported initialization or material driver')
    program = ProcessProgram(config['phases'], p['layer'])
    audit = domain_and_inventory_audit(p, program)
    require(audit['passed'], ', '.join(audit['failures']))
    return config, p, program, audit


def source_hashes(config_path=DEFAULT_CONFIG):
    paths = {Path(__file__).resolve(), Path(config_path).resolve(),
        ROOT/'scripts/afsd_process_motion_20261007.py'}
    result = old.source_hashes()
    result.update({p.relative_to(ROOT).as_posix(): old.expanded.sha(p) for p in paths})
    return result


def independent_initial_state(config_path=DEFAULT_CONFIG):
    config, p, program, audit = load_configuration(config_path)
    backend = old.expanded.backend_provenance()
    pa, state = discrete_process.independent_initial_state(p)
    state.update(classification=CLASSIFICATION, parameters=copy.deepcopy(p),
        framework_config=copy.deepcopy(config), framework_config_sha256=digest(config),
        process_program=program.to_dict(), machine_history=[], domain_inventory_plan=audit,
        independent_initialization=True, parent_checkpoint=None,
        accepted=False, restart_allowed=False, diagnostic_only=True,
        physical_qualification='NOT_RUN', accepted_time_advanced_s=0.,
        source_sha256=source_hashes(config_path), backend=backend,
        initial_mass_kg=float(pa.m.sum()), initial_particle_ids=pa.particle_id.astype(int).tolist(),
        immutable_material_identity_sha256=identity_digest(pa))
    state['bore_observation'] = initial_bore_crossing_state(old.positions(pa), pa.particle_id,
        pa.component_id, pa.m, 0., program.global_pose(0.), p['bore_width'])
    state['crossing_ids'] = []; state['crossing_mass_kg'] = 0.
    require(len(pa.x) == audit['particle_count'], 'domain allocation differs from plan')
    validate_state(state)
    return pa, state


def validate_state(state):
    require(state.get('classification') == CLASSIFICATION and state.get('accepted') is False and
        state.get('restart_allowed') is False and state.get('diagnostic_only') is True and
        state.get('independent_initialization') is True and state.get('parent_checkpoint') is None,
        'new independent framework state required; old or accepted state cannot be relabelled')
    config = state['framework_config']; p = state['parameters']
    require(digest(config) == state['framework_config_sha256'], 'configuration history changed')
    base = ROOT/config['base_parameters']
    require(old.expanded.sha(base) == config['base_parameters_sha256'], 'base source changed')
    expected = json.loads(base.read_text(encoding='utf-8'))['parameters']
    expected.update({k: float(v['value']) for k, v in config['parameter_overrides'].items()})
    require(p == expected, 'parameter history changed')
    program = ProcessProgram(config['phases'], p['layer'])
    require(program.to_dict() == state['process_program'] and
        domain_and_inventory_audit(p, program) == state['domain_inventory_plan'] and
        state['domain_inventory_plan']['passed'], 'program/domain plan differs')
    validate_interval_energy_state(state['interval_energy'])
    validate_bore_crossing_state(state['bore_observation'])
    rows = state['machine_history']; t = 0.
    require(len(rows) == state['step'] == state['interval_energy']['steps'], 'lost process/energy prefix')
    for row in rows:
        require(row == program.step(t, row['dt_s']).to_dict(), 'machine interval history changed')
        t = row['t1_s']
    require(t == state['time_s'] == state['bore_observation']['last_time_s'], 'state times differ')
    require(state['crossing_ids'] == state['bore_observation']['crossing_ids'] and
        state['crossing_mass_kg'] == state['bore_observation']['crossing_mass_kg'],
        'passage summaries differ from persistent observation history')
    return program


class ProcessTrial(old.RampedHotTrial):
    """Only initialization and commanded interval change, not physics operators."""
    def __init__(self, pa, state, dt):
        program = validate_state(state)
        self.pa, self.state, self.p = pa, state, dict(state['parameters'])
        p = self.p
        require(set(pa.particle_id.astype(int)) == set(state['initial_particle_ids']) and
                abs(pa.m.sum()/state['initial_mass_kg']-1) <= 1e-12 and
                identity_digest(pa) == state['immutable_material_identity_sha256'],
                'finite material inventory or persistent IDs changed')
        order = np.argsort(pa.particle_id); bore = state['bore_observation']
        require(pa.particle_id[order].astype(int).tolist() == bore['particle_ids'] and
            pa.component_id[order].astype(int).tolist() == bore['component_ids'] and
            pa.m[order].tolist() == bore['masses_kg'], 'observer and actual per-ID material inventory differ')
        self.x0 = old.positions(pa); self.s0 = pa.S.reshape(-1, 3, 3).copy()
        self.bulk = old.bulk_modulus(p['shear'], p['poisson'])
        self.law = old.FKS(**{k: p[k] for k in old.FKS.__dataclass_fields__})
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
        require(np.array_equal(state['free'], expected_free), 'only bottom xyz/top-cap z may be constrained')

    def _require_interval(self, dt):
        require(isinstance(self.motion, ProcessInterval) and dt == self.motion.dt_s and
            self.state['time_s'] == self.motion.t0 and
            self.state['process_program'] == self.motion.history.to_dict(),
            'all operators must use one exact commanded interval')


def advance_candidate(trial, result):
    """Fresh map plus unchanged force, heat, full-prefix energy and passage gates."""
    require(isinstance(trial, ProcessTrial), 'explicit process trial required')
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
        accepted=False, restart_allowed=False, diagnostic_only=True,
        accepted_time_advanced_s=0., physical_qualification='DETACHED_CANDIDATE_REQUIRES_PATH_AND_TIME_SCREEN')
    require(state['interval_energy']['records'][:-1] == trial.state['interval_energy']['records'],
            'full energy history changed')
    validate_state(state)
    return pa, state, fresh


def operator_reuse_audit():
    inherited = old.static_adaptation_audit()
    require(ProcessTrial.evaluate is old.RampedHotTrial.evaluate and
        ProcessTrial.tangent is old.RampedHotTrial.tangent and
        ProcessTrial.prescribed is old.RampedHotTrial.prescribed,
        'frozen residual, tangent or axial constraint implementation changed')
    return dict(status='PASS_UNCHANGED_INHERITED_THERMOMECHANICAL_OPERATORS',
        inherited_audit=inherited, changed='configuration, independent initialization and machine interval',
        runtime_monkeypatch=False, material_trajectory_prescribed=False)
