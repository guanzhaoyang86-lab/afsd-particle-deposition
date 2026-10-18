"""Independent analytic contact paths and literal coupled block references."""
from pathlib import Path
import copy
import sys
import unittest
from unittest.mock import patch
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.contact_force_secant import response
from afsd_backend.pysph.contact_kick import kick
from afsd_backend.pysph.contact_material_response import ContactMaterialResponse, BoundContactMaterialSearch
from afsd_backend.pysph.material_stage_response import MaterialStageResponse
from test_material_stage_response import coupled, coupled_fields, fields


class ContactSecantTests(unittest.TestCase):
    def close(self, a, b):
        np.testing.assert_allclose(a, b, rtol=2e-5, atol=2e-8)

    def one(self, f=0., v=1., dt=2.):
        return response([v, 0, 0], 1., [f, 0, -1], [[0, 0, 1]], [1.], [[0, 0, 0]], 1., dt)

    def test_no_contact_with_prescribed_axis(self):
        r = response([1, 2, 3], 2, [4, 5, 6], [], [], [], .3, .5, [1, 0, 1], [1, 2, 3])
        self.close(r['B'], np.diag([.125, 0, .125]))
        self.close(r['C'], np.diag([.25, 0, .25]))
        self.assertEqual(r['kick_calls'], 0)

    def test_stopping_distance_mean_is_not_half_endpoint(self):
        r = self.one()
        self.close(r['baseline']['mean_velocity_m_s'], [.25, 0, 0])
        self.close(r['B'][0, 0], .25)
        self.close(r['C'][0, 0], 0.)
        self.assertFalse(r['full_tangent'])

    def test_reverse_sliding_literal_piecewise_path(self):
        # a+=-3, stop=1/3; a-=-1, end=-5/3, distance=-11/9.
        r = self.one(-2.)
        self.close(r['baseline']['velocity_m_s'][0], -5/3)
        self.close(r['baseline']['mean_velocity_m_s'][0], -11/18)
        self.close(r['B'][0, 0], 22/27)
        self.close(r['C'][0, 0], 16/9)
        self.assertTrue(all(q['raw_heat_certificate_pass'] for q in r['probes']))

    def test_two_faces_shared_mass_and_actuator_reaction(self):
        args = ([0, 0, 1], 2., [0, 0, 0], [[1, 0, 0], [0, 1, 0]],
                [1., 2.], [[0, 0, 0], [0, 0, 0]], 1., 2.)
        physical = kick(*args, free_axes=[0, 0, 1], prescribed_velocity=[0, 0, 1])
        r = response(*args, free_axes=[0, 0, 1], prescribed_velocity=[0, 0, 1], baseline_result=physical)
        self.close(r['baseline']['mean_velocity_m_s'], [0, 0, 1/6])
        self.close(r['B'], np.diag([0, 0, 1/18]))
        self.close(r['C'], np.zeros((3, 3)))
        self.close(physical['constraint_impulse'], [-2, -4, 0])
        self.close(physical['constraint_work_J'], 0.)

    def test_cone_boundary_secant_is_named_and_branch_change_preserved(self):
        r = self.one(1., v=0., dt=1.)
        self.close(r['B'][0, 0], .25)
        self.close(r['C'][0, 0], .5)
        self.assertTrue(r['branch_changes'])
        self.assertEqual(r['classification'], 'CONTACT_FORCE_SECANT_APPROXIMATE_NOT_JACOBIAN')

    def test_baseline_mismatch_rejected(self):
        physical = kick([1, 0, 0], 1., [0, 0, -1], [[0, 0, 1]], [1.], [[0, 0, 0]], 1., 2.)
        physical['velocity'][0] += .1
        with self.assertRaisesRegex(ValueError, 'baseline'):
            response([1, 0, 0], 1., [0, 0, -1], [[0, 0, 1]], [1.], [[0, 0, 0]], 1., 2., baseline_result=physical)


class ContactMaterialAlgebraTests(unittest.TestCase):
    def close(self, a, b):
        np.testing.assert_allclose(a, b, rtol=2e-12, atol=2e-13)

    def block(self, prescribed=False):
        kw = coupled()
        if prescribed:
            kw['free_axes'][1, 0] = False
        material = MaterialStageResponse(**kw)
        B = np.tile(np.eye(3), (2, 1, 1))
        B[0] *= .2; B[1] *= .3
        C = B*3  # Deliberately different from 2B.
        if prescribed:
            B[1, 0, :] = 0; B[1, :, 0] = 0
            C[1, 0, :] = 0; C[1, :, 0] = 0
        return ContactMaterialResponse(material, B, C)

    def test_literal_uneliminated_contact_material_density_temperature(self):
        p = self.block(); c, g = coupled_fields(); d = p.apply(c, g)['correction']
        # Unknowns (mean1,mean2,rho1,T1). Jv=[[-1,1],[0,-2]],
        # Jr=(1,-1), JT=(1,-1), kappa=-1, E=(1,-1).
        ref = np.linalg.solve([[1.2, -.2, -.2, -.2], [0, 1.6, .3, .3],
                               [1, -1, 1, 0], [0, 0, 0, 1]], [1, 2, .5, .25])
        self.close([d['mean_velocity_m_s'][0, 0], d['mean_velocity_m_s'][1, 0],
                    d['density_kg_m3'][0], d['temperature_K'][0]], ref)
        force = np.array([-ref[0]+ref[1]+ref[2]+ref[3], -2*ref[1]-ref[2]-ref[3]])
        self.close(d['position_m'][:, 0], 2*np.array([.2, .3])*force)
        self.close(d['velocity_m_s'][:, 0], np.array([.6, .9])*force)

    def test_full_prescribed_columns_survive_density_backsolve(self):
        p = self.block(True); c, g = coupled_fields(); d = p.apply(c, g)['correction']
        # delta2=2; 1.4 delta1 = 1+.4*2+.2*(.5+.25) = 1.95.
        self.close(d['mean_velocity_m_s'][:, 0], [39/28, 2])
        self.close(d['density_kg_m3'][0], 31/28)
        self.close(d['velocity_m_s'][1], [0, 0, 0])

    def test_analytic_free_contact_reduces_to_material_response(self):
        m = MaterialStageResponse(**coupled()); c, g = coupled_fields()
        B = np.array([np.eye(3)*m.dt/(2*mass) for mass in [2, 3]])
        p = ContactMaterialResponse(m, B, 2*B)
        a = p.apply(c, g); b = m.apply(c, g)
        for k in c:
            self.close(a['correction'][k], b['correction'][k])

    def test_zero_and_pure_state(self):
        p = self.block(); c, g = coupled_fields(); before = copy.deepcopy((c, g))
        p.apply(c, g)
        for i, q in enumerate((c, g)):
            for k in c:
                np.testing.assert_array_equal(q[k], before[i][k])
        for v in p.apply(c, c)['correction'].values():
            self.close(v, np.zeros_like(v))
        self.assertFalse(p.matrix.flags.writeable)

    def test_prescribed_nonzero_response_rejected(self):
        kw = coupled(); kw['free_axes'][1, 0] = False
        with self.assertRaisesRegex(ValueError, 'prescribed'):
            ContactMaterialResponse(MaterialStageResponse(**kw), np.ones((2, 3, 3)), np.ones((2, 3, 3)))

    def test_singular_search_matrix_no_regularization(self):
        m = MaterialStageResponse(**coupled()); c, g = coupled_fields()
        # Exposing a singular search system must fail; no pseudo-inverse.
        p = ContactMaterialResponse(m, np.zeros((2, 3, 3)), np.zeros((2, 3, 3)))
        p.matrix = np.zeros((6, 6))
        with patch('numpy.linalg.pinv', side_effect=AssertionError('forbidden')):
            with self.assertRaises(np.linalg.LinAlgError):
                p.apply(c, g)


class BindingTests(unittest.TestCase):
    class Trial:
        def __init__(self):
            self.hash = 'parent'; self.calls = 0; self.fail = False
        def accepted_hash(self):
            return self.hash
        def candidate(self):
            return fields(2)
        def evaluate(self, c, dt):
            self.calls += 1
            if self.fail:
                raise ValueError('trial failed')
            return dict(map=copy.deepcopy(c), accepted_state_hash=self.hash,
                        dt_s=dt, time_before_s=0.)

    def bound(self):
        t = self.Trial()
        return BoundContactMaterialSearch(t, {}, {'t': 0.}, {},
            np.zeros((2, 3, 3)), np.zeros((2, 3, 3)))

    def test_single_delegation_lazy_assembly_and_consumed_context(self):
        b = self.bound(); c = b.candidate(); e = b.evaluate(c, 2.)
        with patch.object(MaterialStageResponse, 'from_journal', return_value=MaterialStageResponse(**coupled())) as factory:
            r = b.apply(c, e['map'])
            self.assertEqual(factory.call_count, 1)
            self.assertEqual(b.trial.calls, 1)
            self.assertEqual(r['diagnostics']['classification'], 'REFRESHED_MATERIAL_FROZEN_CONTACT_SECANT_SEARCH_APPROXIMATE')
        with self.assertRaises(ValueError):
            b.apply(c, e['map'])

    def test_stale_candidate_map_journal_parent_all_rejected(self):
        for kind in ('candidate', 'map', 'journal', 'parent'):
            b = self.bound(); c = b.candidate(); e = b.evaluate(c, 2.)
            m = copy.deepcopy(e['map'])
            if kind == 'candidate': c['temperature_K'][0] += 1
            if kind == 'map': m['temperature_K'][0] += 1
            if kind == 'journal': e['dt_s'] += 1
            if kind == 'parent': b.trial.hash = 'other'
            with self.subTest(kind=kind), self.assertRaisesRegex(ValueError, 'stale'):
                b.apply(c, m)

    def test_failed_trial_invalidates_old_context(self):
        b = self.bound(); c = b.candidate(); e = b.evaluate(c, 2.)
        b.trial.fail = True
        with self.assertRaisesRegex(ValueError, 'trial failed'): b.evaluate(c, 2.)
        with self.assertRaisesRegex(ValueError, 'no completed'): b.apply(c, e['map'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
