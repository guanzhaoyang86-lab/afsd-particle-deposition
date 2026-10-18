"""Independent, material-free qualification of the declared startup motion.

No residual, NNPS, material update, thermal solve, tangent assembly or root is
executed. Geometric certificates apply to the declared per-step world chord
and shortest SO(3) rotation, not to the continuous global smoothstep curve.
The receipt qualifies this interface only; actual coupled startup remains a
separately supervised diagnostic with every original acceptance gate.
"""
from pathlib import Path
import argparse
import ast
import copy
import hashlib
import inspect
import json
import math
import sys
import textwrap
import time
import traceback
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT/'src') not in sys.path:
    sys.path.insert(0, str(ROOT/'src'))
from ramped_hot_motion import RampedHotMotion, RAMP_DURATION_S
import ramped_hot_trial as model
import expanded_hot_domain as domain
from afsd_backend.pysph.quasistatic_contact import _rotation_path
from afsd_backend.pysph.quasistatic_process import pose as old_pose
from afsd_backend.pysph.interval_acceptance import (
    initial_bore_crossing_state, observe_square_bore_crossings)
from batched_contact import BatchedChainRuleContact
from contact_shell_finite_enclosure import audit_contact_shell_finite_enclosure
from visible_support_path_tight import audit_visible_support_path

STATUS = 'PASS_RAMPED_MOTION_COMPONENT_QUALIFICATION_NOT_PROCESS_ACCEPTANCE'


def require(value, message):
    if not bool(value):
        raise AssertionError(message)


def rejected(call):
    try:
        call()
    except (ValueError, ArithmeticError, TypeError, AttributeError):
        return
    raise AssertionError('invalid input was accepted')


def equal_pose(a, b, tolerance=0.):
    for key in ('origin', 'rotation'):
        require(np.allclose(a[key], b[key], atol=tolerance, rtol=0.), 'pose differs: '+key)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    paths = list((ROOT/'src/afsd_backend/pysph').glob('*.py'))
    paths += [Path(__file__), domain.PARAMETER_SOURCE, domain.CONTRACT]
    paths += [ROOT/'scripts'/name for name in (
        'ramped_hot_motion.py', 'ramped_hot_trial.py', 'expanded_hot_domain.py',
        'visible_quasistatic.py', 'batched_contact.py', 'contact_shell_audit.py',
        'contact_shell_enclosure.py', 'contact_shell_finite_enclosure.py',
        'visible_support.py', 'visible_support_numerics.py',
        'visible_support_path.py', 'visible_support_path_tight.py',
        'audit_solid_occluded_neighbors.py')]
    return {p.resolve().relative_to(ROOT).as_posix(): sha(p) for p in sorted(set(paths))}


def integrated_history(p):
    motion = RampedHotMotion(p)
    T = RAMP_DURATION_S
    require(T == .020 and motion.velocity_factor(0.) == 0., 'registered rest startup changed')
    require(motion.velocity_factor(T) == motion.velocity_factor(2*T) == 1., 'steady limit differs')
    require(np.isclose(p['omega']*60/(2*math.pi), 300., atol=1e-12, rtol=0.), 'nominal rpm changed')
    # Four-node Gauss integration is exact for the cubic velocity on each side
    # of T, and is independent of the production factored quartic difference.
    nodes, weights = np.polynomial.legendre.leggauss(4)
    worst = 0.
    for a, b in ((0., .0005), (.0005, .001), (.005, .013), (.019, .021), (.020, .027)):
        intervals = [(a, min(b, T))] if a < T else []
        if b > T:
            intervals.append((max(a, T), b))
        value = math.fsum(.5*(hi-lo)*float(np.dot(weights,
            [motion.velocity_factor(.5*(lo+hi)+.5*(hi-lo)*z) for z in nodes]))
            for lo, hi in intervals)
        error = abs(motion.exposure_increment(a, b)-value)
        worst = max(worst, error)
        require(error <= 3e-18, 'smoothstep feed integral differs from cubic quadrature')
    require(motion.exposure(T) == T/2 and motion.exposure(2*T) == 1.5*T,
            'integrated startup offset differs')
    epsilon = 1e-8
    for point, expected in ((0., 0.), (T, 1.)):
        derivative = ((motion.exposure(point+epsilon)-motion.exposure(point))/epsilon
                      if point == 0. else
                      (motion.exposure(point+epsilon)-motion.exposure(point-epsilon))/(2*epsilon))
        require(abs(derivative-expected) < 1e-8, 'integrated velocity boundary is discontinuous')
    return dict(ramp_duration_s=T, steady_nominal_rpm=p['omega']*60/(2*math.pi),
                maximum_integral_difference_s=worst, no_time_or_modulus_scaling=True)


def refinement_endpoints(p):
    motion = RampedHotMotion(p)
    worst = 0.
    for t0, dt in ((0., .001), (.003, .001), (.0195, .001), (.023, .001)):
        whole = motion.step(t0, dt)
        for count in (2, 4, 8):
            grid = np.linspace(t0, t0+dt, count+1)
            pieces = [motion.step(float(a), float(b-a)) for a, b in zip(grid[:-1], grid[1:])]
            equal_pose(whole.pose(t0), pieces[0].pose(pieces[0].t0))
            equal_pose(whole.pose(t0+dt), pieces[-1].pose(pieces[-1].t1), 3e-16)
            error = abs(math.fsum(s.feed_displacement_m for s in pieces)-whole.feed_displacement_m)
            worst = max(worst, error)
            require(error <= 2e-20, 'refinement changed accumulated actuator displacement')
            for first, second in zip(pieces[:-1], pieces[1:]):
                equal_pose(first.pose(first.t1), second.pose(second.t0), 3e-16)
    return dict(maximum_partition_feed_difference_m=worst,
                endpoint_invariance='WITH_BINARY64_ROUNDOFF', intermediate_paths_may_differ=True)


def step_geodesic(p):
    motion = RampedHotMotion(p)
    error = 0.
    for t0, dt in ((0., .010), (.019, .002), (.021, .003)):
        step = motion.step(t0, dt)
        first, last = step.pose(step.t0), step.pose(step.t1)
        angle, rotate = _rotation_path(first['rotation'], last['rotation'])
        for fraction in (0., .125, .25, .5, .75, 1.):
            actual = step.pose(step.t0+fraction*(step.t1-step.t0))
            expected = dict(origin=(1-fraction)*first['origin']+fraction*last['origin'], rotation=rotate(fraction))
            equal_pose(actual, expected, 8e-16)
            error = max(error, float(np.max(np.abs(actual['rotation']-expected['rotation']))))
            require(np.isclose(step.prescribed_feed_at(step.t0+fraction*(step.t1-step.t0)),
                fraction*step.feed_displacement_m, atol=2e-20, rtol=0.), 'actuator path is not step linear')
        require(angle < math.pi, 'ambiguous step path')
    curved = motion.step(0., .010)
    mismatch = float(np.max(np.abs(curved.pose(.005)['rotation']-motion.global_pose(.005)['rotation'])))
    require(mismatch > .01, 'fixture does not distinguish global and declared step paths')
    return dict(maximum_rotation_matrix_difference=error, global_vs_step_midpoint_difference=mismatch,
                declared_path='LINEAR_WORLD_PARTICLES_AND_ORIGIN_SHORTEST_SO3_TOOL_ROTATION')


def contact_and_certificates(p):
    motion = RampedHotMotion(p)
    step = motion.step(0., .010)
    a, b = step.pose(step.t0), step.pose(step.t1)
    local = np.array([[-.0002, .0001, .00055], [.0002, -.0001, .00065]])
    x0 = local@a['rotation'].T+a['origin']
    x1 = local@b['rotation'].T+b['origin']
    saved = (x0.copy(), x1.copy())
    contact = BatchedChainRuleContact(tool_radius=p['tool_radius'], bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'], particle_radius=p['dx']/2,
        normal_stiffness=p['normal_stiffness'], tangential_stiffness=p['normal_stiffness'],
        friction=p['friction'], base_z=-p['substrate_depth'])
    first, last = contact.geometry(x0, a), contact.geometry(x1, b)
    sweep = contact._sweep(x0, x1, a['origin'], b['origin'], a['rotation'], b['rotation'], first, last)
    require(sweep['certified'] is True, 'contact solid sweep failed')
    shell = audit_contact_shell_finite_enclosure(contact, x0, x1, a, b)
    require(shell['status'] == 'CERTIFIED_NO_HIDDEN_SHELL_CONTACT', 'shared shell path not certified')
    ids = np.array([73, 11], dtype=np.int64)
    path = audit_visible_support_path(x0, x1, np.array([0]), np.array([1]), a, b, p, ids)
    require(path['status'] == 'CERTIFIED_VISIBLE_SUPPORT_PATH', 'shared support path not certified')
    angle, rotate = _rotation_path(a['rotation'], b['rotation'])
    max_local_error = 0.
    for fraction in (0., .25, .5, .75, 1.):
        x = (1-fraction)*x0+fraction*x1
        geometry = contact.geometry(x, step.pose(step.t0+fraction*step.dt_s))
        expected = (x-((1-fraction)*a['origin']+fraction*b['origin']))@rotate(fraction)
        max_local_error = max(max_local_error, float(np.max(np.abs(geometry['local_positions_m']-expected))))
    require(max_local_error < 2e-18, 'geometry does not use the certificate path')
    require(np.array_equal(x0, saved[0]) and np.array_equal(x1, saved[1]), 'geometry mutated endpoint arrays')
    return dict(contact_sweep=sweep, shell_status=shell['status'], support_path_status=path['status'],
        shell_query_count=shell['geometry_queries'], support_query_count=path['analytic_batch_queries'],
        maximum_body_coordinate_difference_m=max_local_error, force_evaluations=0)


def observer_identity(p):
    motion = RampedHotMotion(p)
    step = motion.step(0., .010)
    mid = step.pose(.005)
    half = .5*p['bore_width']
    # At this sharp corner, using the global smooth ramp at the crossing time
    # instead of the numerical step rotation changes the through-bore decision.
    locations = np.array([[.999*half, .999*half, 0.], [0., 0., 0.],
                          [.2*half, 0., 0.], [1.1*half, 0., 0.]])@mid['rotation'].T+mid['origin']
    x0, x1 = locations.copy(), locations.copy()
    x0[:, 2] += .0001
    x1[:, 2] -= .0001
    x0[2, 2] = locations[2, 2]-.00005
    ids = np.array([73, 11, 97, 43], dtype=np.int64)
    components = np.array([1, 2, 1, 1], dtype=np.int64)
    masses = np.array([2., 3., 4., 5.])*1e-7
    initial = initial_bore_crossing_state(x0, ids, components, masses, 0., step.pose(0.), p['bore_width'])
    before = copy.deepcopy(initial)
    result = observe_square_bore_crossings(initial, x0, x1, ids, components, masses, 0., .010, step.pose)
    wrong = observe_square_bore_crossings(initial, x0, x1, ids, components, masses, 0., .010, motion.global_pose)
    require(result['crossing_ids'] == [73] and wrong['crossing_ids'] == [], 'observer fixture failed to distinguish paths')
    order = np.array([3, 1, 0, 2])
    permuted = observe_square_bore_crossings(initial, x0[order], x1[order], ids[order], components[order], masses[order],
                                            0., .010, step.pose)
    require(result == permuted and initial == before, 'passive observer identity or history changed')
    record = result['new_events'][0]
    require(record['step_fraction'] == .5 and record['time_s'] == .005 and
            result['crossing_mass_kg'] == masses[0], 'crossing time or mass is wrong')
    return dict(declared_step_crossing_ids=result['crossing_ids'], wrong_global_curve_crossing_ids=wrong['crossing_ids'],
        id_permutation_equal=True, initial_below_and_substrate_ineligible=True, original_history_unchanged=True)


def invalid_intervals(p):
    motion = RampedHotMotion(p)
    for t in (-1e-12, float('nan'), float('inf')):
        rejected(lambda t=t: motion.global_pose(t))
    for dt in (0., -1e-6, float('nan'), float('inf'), .2):
        rejected(lambda dt=dt: motion.step(.02, dt))
    rejected(lambda: motion.exposure_increment(.002, .001))
    step = motion.step(.001, .0005)
    rejected(lambda: step.pose(np.nextafter(step.t0, -np.inf)))
    rejected(lambda: step.pose(np.nextafter(step.t1, np.inf)))
    rejected(lambda: setattr(step, 't0', 4.))
    return dict(reversed_nonfinite_unresolvable_geodesic_and_extrapolated_inputs='REJECTED', immutable_interval=True)


def initial_identity_and_trial(p, provenance):
    oldpa, oldstate = domain.independent_initial_state(p, provenance)
    pa, state = model.independent_initial_state()
    require(len(pa.x) == 918 and state['time_s'] == 0. and state['step'] == 0, 'new initial inventory/time differs')
    require(set(pa.properties) == set(oldpa.properties), 'initial particle schema differs')
    for name, prop in pa.properties.items():
        require(np.array_equal(prop.get_npy_array(), oldpa.properties[name].get_npy_array()), 'initial field changed: '+name)
    require(state['parameters'] == p, 'ramp changed material/geometry parameters')
    for name in ('free', 'contact_history'):
        require(np.array_equal(state[name], oldstate[name]), 'initial constraint/history changed: '+name)
    for name in ('bore_observation', 'interval_energy'):
        require(state[name] == oldstate[name], 'initial passive/energy state differs: '+name)
    require(state['classification'] != oldstate['classification'], 'new boundary model mislabeled as old')
    require(state['accepted'] is False and state['restart_allowed'] is False, 'initial state is accepted')
    saved = domain.clean(copy.deepcopy(state))
    trial = model.make_trial(pa, state, .0005)
    expected = np.zeros((len(pa.x), 3))
    expected[pa.actuator.astype(bool), 2] = -p['feed']*RampedHotMotion(p).exposure_increment(0., .0005)
    require(np.array_equal(trial.prescribed(.0005), expected), 'prescribed physical axes/feed integral differ')
    equal_pose(trial.motion.pose(0.), old_pose(p, 0.))
    require(domain.clean(state) == saved, 'trial construction reset/changed caller state')
    rejected(lambda: model.make_trial(oldpa, oldstate, .0005))
    return dict(particle_count=918, all_initial_particle_fields_bitwise_equal=True,
        material_parameters_unchanged=True, old_state_rejected=True, initial_rank_t0_pose_bitwise_equal=True,
        prescribed_particle_count=int(pa.actuator.sum()), prescribed_feed_displacement_m=float(expected[:,2].min()),
        constructor_evaluations_only=True, caller_state_unchanged=True)


def trial_source_audit():
    # This check is independent of any AST audit shipped inside the adapter.
    # An allow-list of the only intended differences is validated below after
    # normalization; extra numerical edits therefore cannot pass unnoticed.
    return independent_trial_ast_audit()


def independent_trial_ast_audit():
    from visible_quasistatic import VisibleSupportTrial
    from afsd_backend.pysph.quasistatic_process import QuasistaticTrial
    from afsd_backend.pysph import quasistatic_tangent
    def function(value):
        return ast.parse(textwrap.dedent(inspect.getsource(value))).body[0]
    def dump(value):
        return ast.dump(value, include_attributes=False)
    class RestoreNominalNames(ast.NodeTransformer):
        count = 0
        def visit_Call(self, node):
            self.generic_visit(node)
            if dump(node.func) == dump(ast.parse('self.motion.pose', mode='eval').body):
                require(len(node.args) == 1 and not node.keywords, 'unexpected interval pose signature')
                node.func = ast.Name(id='pose', ctx=ast.Load())
                node.args.insert(0, ast.Name(id='p', ctx=ast.Load()))
                self.count += 1
            return node
        def visit_Name(self, node):
            if node.id == 'SPATIAL_CLASSIFICATION':
                node.id = 'CLASSIFICATION'
            return node
    normalizer = RestoreNominalNames()
    normalized = normalizer.visit(function(model.RampedHotTrial.evaluate))
    require(normalizer.count == 4 and dump(normalized) == dump(function(VisibleSupportTrial.evaluate)),
            'evaluate contains changes beyond the four poses and classification alias')
    require(model.SPATIAL_CLASSIFICATION == 'HOT_QS_VISIBLE_SUPPORT_CHAIN_RULE_CONTACT_APPROXIMATE_V1',
            'spatial classification value changed')
    tangent = function(model.RampedHotTrial.tangent)
    guards = ast.parse("self._require_interval(result['dt_s'])\ngraph=result['visibility_graph']\n"
        "for key in ('i','j','g'):\n"
        "    if not np.array_equal(result[key],graph[key]):\n"
        "        raise ValueError('tangent must use the exact residual visible graph')\n").body
    require([dump(n) for n in tangent.body[2:5]] == [dump(n) for n in guards],
            'tangent interval/visible-graph guards differ from explicit allowance')
    del tangent.body[2:5]
    imports = [n for n in tangent.body if isinstance(n, ast.ImportFrom)]
    require(len(imports) == 1 and imports[0].module == 'afsd_backend.pysph.quasistatic_tangent'
            and imports[0].level == 0, 'standalone tangent must use the real absolute import')
    for item in imports[0].names:
        require(hasattr(quasistatic_tangent, item.name) and item.asname is None, 'tangent import cannot resolve')
    imports[0].module = 'quasistatic_tangent'
    imports[0].level = 1
    normalizer = RestoreNominalNames()
    normalized = normalizer.visit(tangent)
    require(normalizer.count == 2 and dump(normalized) == dump(function(QuasistaticTrial.tangent)),
            'tangent contains numerical edits outside the explicit import/pose/guard allowance')
    advance = function(model.advance_candidate)
    calls = [n for n in ast.walk(advance) if isinstance(n, ast.Call)]
    observers = [n for n in calls if isinstance(n.func, ast.Name) and n.func.id == 'observe_square_bore_crossings']
    require(len(observers) == 1 and [ast.unparse(n) for n in observers[0].args[-3:]] ==
            ['trial.motion.t0', 'trial.motion.t1', 'trial.motion.pose'], 'observer does not use the same step motion')
    accepted = [ast.unparse(n.func) for n in calls if isinstance(n.func, ast.Attribute) and n.func.attr == 'accepted_copy']
    require(accepted == ['original_qs.accepted_copy'], 'hardcoded legacy discrete observer was reintroduced')
    # Calling invalid dt validates the pre-operator guard without evaluating
    # the physics; source equivalence alone does not establish its execution.
    pa, state = model.independent_initial_state()
    trial = model.make_trial(pa, state, .0005)
    rejected(lambda: trial.prescribed(.00025))
    return dict(evaluate_pose_replacements=4, tangent_pose_replacements=2,
        unchanged_evaluate_ast=True, unchanged_tangent_numerical_ast=True,
        standalone_tangent_import_resolves=True, step_observer_wiring=True,
        contact_fd_epsilon_m=1e-9, no_constitutive_or_thermal_edits=True,
        allowed_non_pose_deltas=['original spatial classification alias', 'absolute tangent import',
            'existing visible graph guard plus explicit interval guard'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='outputs/ramped_hot_motion_20261005/v1')
    args = parser.parse_args()
    out = (ROOT/args.output).resolve()
    require(out.is_relative_to(ROOT) and not out.exists(), 'new output directory inside project required')
    out.mkdir(parents=True, exist_ok=False)
    begin = time.perf_counter()
    frozen = source_hashes()
    report = dict(status='RUNNING_MATERIAL_FREE_COMPONENT_CHECKS', source_and_input_sha256=frozen,
        command=[sys.executable, *sys.argv], checks=[], physical_evaluations=0, material_evaluations=0,
        thermal_evaluations=0, nnps_evaluations=0, root_evaluations=0, accepted_time_advanced_s=0.,
        accepted=False, restart_allowed=False, initial_particle_count=918,
        scope='GEOMETRY_KINEMATICS_IDENTITY_AND_STATIC_ADAPTATION_ONLY',
        limitations=['No coupled root, force response, thermal response, stability or process acceptance is qualified.',
            'Subdivision preserves endpoint motion but changes the declared numerical interior path.',
            'The step certificate is an engineering floating-point bound for the declared interpolation.'])
    try:
        p, provenance = domain.parameters(.008)
        checks = [(f.__name__, lambda f=f: f(p)) for f in (
            integrated_history, refinement_endpoints, step_geodesic, contact_and_certificates,
            observer_identity, invalid_intervals)]
        checks += [('initial_identity_and_trial', lambda: initial_identity_and_trial(p, provenance)),
                   ('independent_trial_ast_audit', trial_source_audit)]
        for name, function in checks:
            before = time.perf_counter()
            result = function()
            report['checks'].append(dict(name=name, status='PASS', metrics=result,
                                         wall_seconds=time.perf_counter()-before))
        require(source_hashes() == frozen, 'source changed during checks')
        report['status'] = STATUS
        report['exit_code'] = 0
    except Exception:
        report.update(status='FAIL_RAMPED_MOTION_COMPONENT_QUALIFICATION', exit_code=1, error=traceback.format_exc())
    report['wall_seconds'] = time.perf_counter()-begin
    (out/'receipt.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=report['status'], checks=len(report['checks']), wall_seconds=report['wall_seconds'])))
    if report['exit_code']:
        print(report['error'], file=sys.stderr)
        raise SystemExit(1)


if __name__ == '__main__':
    main()
