"""Prepare independent length controls; initialization only, no physical solve."""
from pathlib import Path
import copy
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / 'length_controls'
sys.path.insert(0, str(ROOT/'scripts'))
import afsd_deposition_setup_20261007 as model

OUT.mkdir(exist_ok=True)
original = ROOT/'configs/deposition_fill_then_track_20261007.json'
base = json.loads(original.read_text(encoding='utf-8'))
# The common unused future track is shortened so BOTH finite inventories pass.
# The actual comparison scope is stationary feeding, through at most 0.12044 s.
base['process_intent']['target_total_traverse_m'] = .001
base['phases'][-1]['name'] = 'diagnostic_unused_track_to_1mm'
base['phases'][-1]['duration_s'] = (.001-.5*.020*.00212)/.00212
base['phases'][-1]['source'] = ('Unused common tail for finite-inventory planning of both controls; '
    'all comparison solves must remain in stationary feeding before 0.220 s.')
base['diagnostic_control'] = dict(
    classification='INDEPENDENT_HOT_FREE_BAR_LENGTH_SENSITIVITY_NOT_PHYSICAL_VALIDATION',
    independent_zero_time_initialization=True, checkpoint_import=False,
    comparison_physical_end_limit_s=.12043920898437391,
    changed_factor='bar_length and resulting position of the same 0.5 mm material actuator cap',
    limitations=['Not a proposed calibrated production geometry or a representation of cold feed/support.',
        'Existing simulation state must not be cropped or relabelled as either control.',
        'A response difference would be sensitivity evidence, not proof that the shorter bar is correct.'])
cases, parameters, programs, states = [], [], [], []
for name, length in [('long_6p5mm',.0065),('short_4p0mm',.004)]:
    config = copy.deepcopy(base)
    config['parameter_overrides']['bar_length']['value'] = length
    config['parameter_overrides']['bar_length']['source'] = (
        'Declared independent control for free hot-bar length; not a modification of the saved continuation.')
    # Replace a now-inapplicable inherited prose limitation; retain the others.
    config['limitations'] = [v for v in config['limitations'] if not v.startswith('The 6.5 mm bar')]
    config['limitations'].append('The bar remains deformable above the 2.5 mm tool top; only the top 0.5 mm cap z is prescribed.')
    path = OUT/(name+'.json')
    path.write_text(json.dumps(config,indent=2)+'\n',encoding='utf-8')
    resolved,p,program,audit = model.load_configuration(path)
    pa,state = model.independent_initial_state(path)
    assert state['time_s']==0. and state['step']==0 and state['parent_checkpoint'] is None
    assert not state['accepted'] and not state['restart_allowed']
    assert np.count_nonzero(pa.actuator)==25
    assert np.all(pa.temperature==p['temperature'])
    assert np.all(pa.ep==0.) and np.all(pa.S==0.)
    cases.append(dict(name=name,configuration=str(path.relative_to(ROOT)),
        configuration_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        geometry_inventory_audit=audit,free_degrees_of_freedom=int(state['free'].sum()),
        initial_geometry_identity=state['independent_geometry_audit'],
        simulation_status='NOT_RUN',nonlinear_cost_benchmark_status='NOT_RUN'))
    parameters.append(p);programs.append(program);states.append((pa,state))

changed = [k for k in parameters[0] if parameters[0][k] != parameters[1][k]]
assert changed == ['bar_length']
assert programs[0].to_dict() == programs[1].to_dict()
original_config,original_p,original_program,original_audit = model.load_configuration(original)
for time_s in [0.,.010,.020,.050,.08958959960937454,.12043920898437391,.220]:
    np.testing.assert_array_equal(programs[0].integral(time_s), original_program.integral(time_s))
result = dict(classification='PASS_INITIAL_CONFIGURATION_AND_IDENTITY_ONLY_NOT_MECHANICAL_VALIDATION',
    original_configuration=str(original.relative_to(ROOT)),
    original_configuration_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
    changed_resolved_physics_parameters=changed,common_process_program=True,
    original_stationary_command_integrals_preserved=True,
    physical_steps_executed=0,simulation_status='NOT_RUN',cases=cases,
    next_execution='Measure nonlinear cost and check original physics/time/geometry gates for each independent setup before assigning a bounded comparison run.',
    comparison_metrics=['Actuator reaction and actual imposed displacement at matched physical time.',
        'Fixed lower-layer axial motion, upper-bar shortening and centroid-based radial size.',
        'Supplied crossing IDs/mass and shoulder occupancy/contact, keeping initial lower material separate.',
        'Original force/energy/geometry/time-screen outcomes; retain rejected attempts.'])
(OUT/'preparation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(dict(status=result['classification'],physical_steps_executed=0,
    cases=[dict(name=c['name'],particles=c['geometry_inventory_audit']['particle_count'],
        free_dof=c['free_degrees_of_freedom'],simulation_status=c['simulation_status']) for c in cases]),indent=2))
