"""Independent physical-unit proposal tests; no real AFSD advancement."""
from pathlib import Path
import hashlib
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.stage_preconditioner import (
    CLASSIFICATION, MAXIMUM_BACKWARD_ERROR, StagePreconditioner)
from afsd_backend.pysph.structural_mobility import ConstrainedResponse
from afsd_backend.pysph.stage_iteration import FIELDS
from afsd_backend.pysph.particles import make_array


def fields():
    return dict(position_m=np.zeros((2, 3)), velocity_m_s=np.zeros((2, 3)),
        density_kg_m3=np.array([3., 4.]), temperature_K=np.array([400., 401.]),
        mean_velocity_m_s=np.zeros((2, 3)))


def spring(constrained=False, dt=1., mass_scale=1., time_scale=1.):
    k = np.zeros((6, 6))
    k[np.ix_([0, 3], [0, 3])] = np.array([[4., -4.], [-4., 4.]])*mass_scale/time_scale**2
    free = np.ones(6, bool)
    if constrained:
        free[3] = False
    op = ConstrainedResponse(np.repeat([2., 3.], 3)*mass_scale, k, dt, free)
    return StagePreconditioner.from_response(op, np.array([19, 4]))


def copied(value):
    return {key: a.copy() for key, a in value.items()}


def pa_hash(pa):
    result = hashlib.sha256()
    for key in sorted(pa.properties):
        result.update(key.encode())
        result.update(pa.properties[key].get_npy_array().tobytes())
    return result.hexdigest()


class CandidateAlgebraTests(unittest.TestCase):
    def assertClose(self, actual, expected):
        np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-13)

    def test_01_literal_two_dof_correction_and_full_linear_evidence(self):
        c = fields(); mapping = copied(c)
        mapping['mean_velocity_m_s'][:, 0] = [1., 2.]
        got = spring().apply(c, mapping)
        self.assertClose(got['correction']['mean_velocity_m_s'][:, 0], [14/11, 20/11])
        self.assertClose(got['correction']['position_m'][:, 0], [3/11, -2/11])
        self.assertClose(got['correction']['velocity_m_s'][:, 0], [6/11, -4/11])
        info = got['diagnostics']
        self.assertClose(info['rhs_Ns'], [2., 0., 0., 6., 0., 0.])
        self.assertClose(info['algebraic_reaction_impulse_Ns'], np.zeros(6))
        np.testing.assert_array_equal(info['particle_id'], [19, 4])
        self.assertLessEqual(info['backward_error_inf'], MAXIMUM_BACKWARD_ERROR)
        self.assertEqual(info['classification'], CLASSIFICATION)
        self.assertFalse(info['accepted']); self.assertFalse(info['physical_map'])

    def test_02_independent_linear_root_with_five_field_coupling(self):
        # Literal 2-mass midpoint oscillator, NOT generated with response/K/apply.
        # Independent root A=[[3,-1],[-1,4]], rhs=(2,-3).
        c = fields(); c['velocity_m_s'][:, 0] = [.7, -.4]
        c['position_m'][:, 0] = [.8, -.2]
        p = spring()
        def physical_map(y):
            result = fields()
            spring_force = -2.*(y['position_m'][0, 0]-y['position_m'][1, 0])
            result['velocity_m_s'][:, 0] = [1.+spring_force/2., -1.-spring_force/3.]
            result['mean_velocity_m_s'][:, 0] = .5*(np.array([1., -1.])+result['velocity_m_s'][:, 0])
            result['position_m'] = result['mean_velocity_m_s'].copy()
            # One-way manufactured scalar coupling remains in the real map.
            result['density_kg_m3'] = [3.+.25*y['mean_velocity_m_s'][0, 0], 4.]
            result['temperature_K'] = [400., 401.+.5*y['mean_velocity_m_s'][1, 0]]
            return result
        for _ in range(4):
            c = p.apply(c, physical_map(c))['proposal']
        expected = fields()
        expected['mean_velocity_m_s'][:, 0] = [5/11, -7/11]
        expected['position_m'][:, 0] = [5/11, -7/11]
        expected['velocity_m_s'][:, 0] = [-1/11, -3/11]
        expected['density_kg_m3'] = np.array([3.+5/44, 4.])
        expected['temperature_K'] = np.array([400., 401.-7/22])
        for key in FIELDS:
            self.assertClose(c[key], expected[key])
            self.assertClose(physical_map(c)[key]-c[key], np.zeros_like(c[key]))

    def test_03_nonzero_prescribed_correction_not_full_actuator_command(self):
        c = fields(); c['mean_velocity_m_s'][:, 0] = [-3., -8.]
        mapping = copied(c); mapping['mean_velocity_m_s'][:, 0] += [1., .2]
        got = spring(constrained=True).apply(c, mapping)
        self.assertClose(got['correction']['mean_velocity_m_s'][:, 0], [11/15, 1/5])
        self.assertClose(got['diagnostics']['algebraic_reaction_impulse_Ns'], [0., 0., 0., -8/15, 0., 0.])
        self.assertClose(got['diagnostics']['shifted_free_rhs_Ns'], [2.2, 0., 0., 0., 0.])
        self.assertClose(got['proposal']['mean_velocity_m_s'][1, 0], -7.8)

    def test_04_kinematic_defects_and_direct_scalar_corrections_retained(self):
        rng = np.random.default_rng(912)
        c = fields(); mapping = {key: c[key]+rng.normal(size=c[key].shape) for key in FIELDS}
        dt = .3; got = spring(dt=dt).apply(c, mapping)
        d = {key: mapping[key]-c[key] for key in FIELDS}; delta = got['correction']
        self.assertClose(delta['position_m']-dt*delta['mean_velocity_m_s'],
                         d['position_m']-dt*d['mean_velocity_m_s'])
        self.assertClose(delta['velocity_m_s']-2*delta['mean_velocity_m_s'],
                         d['velocity_m_s']-2*d['mean_velocity_m_s'])
        for key in ('density_kg_m3', 'temperature_K'):
            np.testing.assert_array_equal(delta[key], d[key])

    def test_05_event_average_root_not_forced_to_endpoint_average(self):
        c = fields()
        # Independent slip reversal reference: displacement=-1/18 m over 1 s;
        # endpoint mean=(1-2/3)/2=1/6 differs from real path average -1/18.
        c['position_m'][:, 0] = -1/18
        c['velocity_m_s'][:, 0] = -2/3
        c['mean_velocity_m_s'][:, 0] = -1/18
        got = spring().apply(c, copied(c))
        for key in FIELDS:
            np.testing.assert_array_equal(got['proposal'][key], c[key])
            np.testing.assert_array_equal(got['correction'][key], np.zeros_like(c[key]))
        self.assertNotEqual(got['proposal']['mean_velocity_m_s'][0, 0], 1/6)

    def test_06_mass_only_is_direct_map_including_independent_mean_velocity(self):
        op = ConstrainedResponse(np.repeat([2., 3.], 3), np.zeros((6, 6)), 5., np.ones(6, bool))
        p = StagePreconditioner.from_response(op, np.array([19, 4]))
        c = fields(); mapping = {key: c[key]+.5 for key in FIELDS}
        got = p.apply(c, mapping)
        for key in FIELDS:
            np.testing.assert_array_equal(got['proposal'][key], mapping[key])

    def test_07_dimensional_covariance_manufactured_only(self):
        # This is a unit-coordinate transformation, NOT production mass/time scaling.
        rng = np.random.default_rng(814); c = fields()
        mapping = {key: c[key]+rng.normal(size=c[key].shape) for key in FIELDS}
        scales = dict(position_m=5., velocity_m_s=5/3, mean_velocity_m_s=5/3,
                      density_kg_m3=7/5**3, temperature_K=2.)
        first = spring().apply(c, mapping)
        transformed = spring(dt=3., mass_scale=7., time_scale=3.).apply(
            {k: v*scales[k] for k, v in c.items()}, {k: v*scales[k] for k, v in mapping.items()})
        for key in FIELDS:
            self.assertClose(transformed['correction'][key], first['correction'][key]*scales[key])
        self.assertClose(transformed['diagnostics']['rhs_Ns'], first['diagnostics']['rhs_Ns']*7*5/3)

    def test_08_no_false_zero_correction_and_independent_inverse(self):
        rng = np.random.default_rng(481)
        for constrained in (False, True):
            p = spring(constrained); c = fields()
            mapping = {key: c[key]+rng.normal(size=c[key].shape) for key in FIELDS}
            got = p.apply(c, mapping); delta = got['correction']
            bar = delta['mean_velocity_m_s'][:, 0]
            recovered = np.array([(3*bar[0]-bar[1])/2, (-bar[0]+4*bar[1])/3])
            if constrained:
                recovered[1] = bar[1]
            self.assertClose(recovered, (mapping['mean_velocity_m_s']-c['mean_velocity_m_s'])[:, 0])
            self.assertGreater(np.linalg.norm(delta['mean_velocity_m_s']), 0.)

    def test_09_backward_error_independent_formula(self):
        c = fields(); mapping = copied(c)
        mapping['mean_velocity_m_s'][:, 0] = [.31, -.71]
        result = spring().apply(c, mapping); info = result['diagnostics']
        x = result['correction']['mean_velocity_m_s'][:, 0]
        a = np.array([[3., -1.], [-1., 4.]])
        b = np.array([.62, -2.13])
        numerator = np.max(np.abs(a@x-b))
        denominator = 5.*np.max(np.abs(x))+np.max(np.abs(b))
        self.assertClose(info['backward_error_numerator_Ns'], numerator)
        self.assertClose(info['backward_error_denominator_Ns'], denominator)
        self.assertClose(info['backward_error_inf'], numerator/denominator)

    def test_10_inputs_and_outputs_do_not_alias(self):
        c = fields(); mapping = {key: c[key]+1. for key in FIELDS}
        c0, map0 = copied(c), copied(mapping)
        p = spring(); result = p.apply(c, mapping)
        result['proposal']['position_m'][:] = 99.
        result['diagnostics']['direct_defect']['position_m'][:] = 20.
        result['diagnostics']['particle_id'][:] = 0
        for key in FIELDS:
            np.testing.assert_array_equal(c[key], c0[key]); np.testing.assert_array_equal(mapping[key], map0[key])
        np.testing.assert_array_equal(p.particle_ids, [19, 4])
        self.assertFalse(np.any(result['correction']['position_m'] == 20.))

    def test_11_invalid_nonfinite_shape_and_identity_reject(self):
        p = spring(); c = fields()
        with self.assertRaises(ValueError): p.apply({'position_m': c['position_m']}, c)
        for key in FIELDS:
            bad = copied(c); bad[key] = bad[key][1:]
            with self.assertRaises(ValueError): p.apply(c, bad)
            bad = copied(c); bad[key].flat[0] = np.nan
            with self.assertRaises(ValueError): p.apply(bad, c)
        with self.assertRaises(ValueError): StagePreconditioner.from_response(p.response, np.array([3, 3]))
        with self.assertRaises(ValueError): StagePreconditioner.from_response(p.response, np.array([3., 4.]))

    def test_12_no_temperature_density_clipping_in_proposal(self):
        c = fields(); mapping = copied(c)
        mapping['temperature_K'][:] = -12.; mapping['density_kg_m3'][:] = -30.
        result = spring().apply(c, mapping)
        np.testing.assert_array_equal(result['proposal']['temperature_K'], [-12., -12.])
        np.testing.assert_array_equal(result['proposal']['density_kg_m3'], [-30., -30.])
        # Caller MUST reject these via physical evaluate; proposal is unaccepted.
        self.assertFalse(result['diagnostics']['accepted'])


class ActualThreeDimensionalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        xyz = np.array(np.meshgrid([-.0005, 0., .0005], [-.0005, 0., .0005],
                                   [-.0005, 0., .0005], indexing='ij')).reshape(3, -1).T
        cls.p = dict(dx=.0005, h_ratio=1.3, rho0=2700., width=.0015,
            marker_low=-.001, marker_high=.001, bar_length=.01,
            actuator_depth=.0005, temperature=400., traverse=.00212, feed=.00212,
            shear=26.3e9, poisson=.33, release_time=1.)
        cls.pa = make_array(xyz, np.zeros(27, int), cls.p)
        cls.before = pa_hash(cls.pa)
        cls.pc = StagePreconditioner(cls.pa, cls.p, 1e-7, 0.)
        cls.c = dict(position_m=xyz.copy(), velocity_m_s=np.zeros_like(xyz),
            mean_velocity_m_s=np.zeros_like(xyz), density_kg_m3=cls.pa.rho.copy(),
            temperature_K=cls.pa.temperature.copy())

    def test_13_three_rigid_translations_and_rotations_unaltered(self):
        xyz = self.c['position_m']
        for axis in np.eye(3):
            for mode in (np.broadcast_to(axis, xyz.shape).copy(), np.cross(axis, xyz)):
                mapping = copied(self.c); mapping['mean_velocity_m_s'] = mode.copy()
                result = self.pc.apply(self.c, mapping)
                delta = result['correction']['mean_velocity_m_s']
                self.assertLessEqual(np.linalg.norm(delta-mode), 2e-11*np.linalg.norm(mode))
                self.assertLessEqual(result['diagnostics']['backward_error_inf'], 1e-12)

    def test_14_passive_labels_do_not_change_operator_or_correction(self):
        q = self.pa.extract_particles(np.arange(27, dtype=np.uint32))
        q.tracer_type[:] = np.arange(27) % 3; q.wholebar_tag[:] = np.arange(27) % 4
        other = StagePreconditioner(q, self.p, 1e-7, 0.)
        np.testing.assert_array_equal(self.pc.response.matrix, other.response.matrix)
        mapping = copied(self.c); mapping['mean_velocity_m_s'] = np.sin(np.arange(81)).reshape(27, 3)
        first, second = self.pc.apply(self.c, mapping), other.apply(self.c, mapping)
        for key in FIELDS:
            np.testing.assert_array_equal(first['correction'][key], second['correction'][key])
        self.assertEqual(self.before, pa_hash(self.pa))

    def test_15_original_actuator_constraint_and_release_interval_guard(self):
        q = self.pa.extract_particles(np.arange(27, dtype=np.uint32))
        q.z0[-1] = .0098
        p = StagePreconditioner(q, self.p, 1e-7, 0.)
        self.assertEqual(np.count_nonzero(~p.response.free), 1)
        self.assertFalse(p.response.free[-1])
        mapping = copied(self.c); mapping['mean_velocity_m_s'][-1, 2] = -.004
        result = p.apply(self.c, mapping)
        self.assertEqual(result['correction']['mean_velocity_m_s'][-1, 2], -.004)
        self.assertEqual(result['diagnostics']['constraint_phase'], 'DRIVEN')
        with self.assertRaisesRegex(ValueError, 'release must lie'):
            StagePreconditioner(q, self.p, .2, .9)
        released = StagePreconditioner(q, self.p, 1e-7, 1.)
        self.assertTrue(np.all(released.response.free))

    def test_16_source_particle_history_and_physical_masses_unchanged(self):
        self.assertEqual(self.before, pa_hash(self.pa))
        np.testing.assert_array_equal(self.pc.response.mass_per_dof, np.repeat(self.pa.m, 3))
        self.assertEqual(self.pc.metadata['frozen_geometry'], 'ACCEPTED_STATE')
        self.assertFalse(self.pc.metadata['full_tangent'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
