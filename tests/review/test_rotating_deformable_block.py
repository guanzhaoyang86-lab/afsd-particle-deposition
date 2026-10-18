"""Independent short regression gates; not a continuum trajectory reference."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import probe_rotating_deformable_block as c
from afsd_backend.pysph.mechanics import corrected_operators, internal_force, positions, velocities
from afsd_backend.pysph.material import update, FKS
from afsd_backend.pysph.output import restore, checkpoint
from afsd_scaffold.vtk_legacy import read_point_cloud


class RotatingDeformableShort(unittest.TestCase):
    def test_01_si_sources_mass_and_true_3d_initial_geometry(self):
        config = json.loads((ROOT/'configs/benchmarks/rotating_deformable_block.json').read_text(encoding='utf-8'))
        for item in config['parameters'].values():
            self.assertTrue(all(k in item for k in ('value', 'unit', 'source', 'status')))
            self.assertIsNotNone(item['value'])
        b = c.Block(5e-9); x = positions(b.pa)
        self.assertEqual(len(x), 64)
        for a in range(3): self.assertEqual(len(np.unique(x[:, a])), 4)
        self.assertEqual(np.linalg.matrix_rank(x-x.mean(axis=0)), 3)
        np.testing.assert_allclose(b.pa.m, 3.375e-7, rtol=0, atol=1e-21)
        self.assertAlmostEqual(b.pa.m.sum(), 2.16e-5, delta=1e-18)
        self.assertEqual(np.count_nonzero(b.pa.boundary_role == 0), 32)
        self.assertEqual(np.count_nonzero(b.pa.boundary_role == 1), 16)
        self.assertEqual(np.count_nonzero(b.pa.boundary_role == 2), 16)
        self.assertFalse(config['production_authorized'])

    def test_02_first_kick_independent_literal_force_solution(self):
        b = c.Block(5e-9); x0 = positions(b.pa).copy(); dt = 5e-9
        top = b.pa.boundary_role == 2
        wall = np.column_stack((-20*x0[top, 1], 20*x0[top, 0], np.zeros(16)))
        impulse = .2*dt*wall/np.linalg.norm(wall, axis=1)[:, None]
        expected_v = np.zeros((64, 3)); expected_v[top] = impulse/3.375e-7
        expected_x = x0+.5*dt*expected_v
        b.advance()
        np.testing.assert_allclose(velocities(b.pa), expected_v, rtol=0, atol=1e-12)
        np.testing.assert_allclose(positions(b.pa), expected_x, rtol=0, atol=2e-17)
        np.testing.assert_array_equal(b.pa.S, np.zeros(64*9))
        np.testing.assert_array_equal(b.pa.ep, np.zeros(64))
        np.testing.assert_array_equal(b.pa.rho, np.full(64, 2700.))
        traces = [a for a in b.last_traces if int(a['particle_id']) in set(b.pa.particle_id[top])]
        self.assertEqual(len(traces), 16)
        q = sum(float(np.sum(t['result']['friction_heat_J'])) for t in traces)
        wall_work = sum(float(t['result']['wall_work_J']) for t in traces)
        expected_k = float(np.sum(impulse*impulse)/(2*3.375e-7))
        expected_wall = float(np.sum(impulse*wall))
        self.assertAlmostEqual(q, expected_wall-expected_k, delta=1e-18)
        self.assertAlmostEqual(wall_work, expected_wall, delta=1e-18)
        for t in traces:
            np.testing.assert_allclose(t['result']['normal_impulse'], [[0, 0, -dt]], atol=1e-20, rtol=0)
            np.testing.assert_allclose(t['result']['constraint_impulse'], [0, 0, dt], atol=1e-20, rtol=0)

    def test_03_affine_three_dimensional_gradient_literal(self):
        b = c.Block(5e-9); x = positions(b.pa)
        a = np.array([[.3, -.7, .4], [1.2, .1, -.6], [-.2, .9, .8]])
        v = x@a.T+[.02, -.01, .03]
        b.pa.u[:], b.pa.v[:], b.pa.w[:] = v.T
        i, j, g, grad, rate, correction = corrected_operators(b.pa)
        np.testing.assert_allclose(grad, np.broadcast_to(a, grad.shape), rtol=0, atol=1e-10)
        np.testing.assert_allclose(rate, -2700*np.trace(a), rtol=0, atol=1e-7)
        self.assertLess(float(np.linalg.cond(np.linalg.inv(correction)).max()), 1e6)
        for k in range(64):
            neighbors = np.r_[j[i == k], i[j == k]]
            self.assertEqual(np.linalg.matrix_rank(x[neighbors]-x[k]), 3)

    def test_04_nonaffine_adjoint_power_total_force_and_torque(self):
        b = c.Block(5e-9); x = positions(b.pa)
        rng = np.random.default_rng(704)
        v = rng.normal(size=(64, 3))
        b.pa.u[:], b.pa.v[:], b.pa.w[:] = v.T
        sigma = rng.normal(size=(64, 3, 3)); sigma = (sigma+sigma.swapaxes(1, 2))*1e4
        i, j, g, grad, _, correction = corrected_operators(b.pa)
        f = internal_force(b.pa, sigma, i, j, g, correction=correction)
        # Independent contraction identity, not a force-generated reference trajectory.
        power = float(np.sum(f*v))
        continuum = float(np.einsum('i,iab,iab->', b.pa.m/b.pa.rho, sigma, grad))
        scale = float(np.sum(np.linalg.norm(f, axis=1)*np.linalg.norm(v, axis=1)))
        self.assertLess(abs(power+continuum), 1e-12+1e-10*scale)
        self.assertLess(np.linalg.norm(f.sum(axis=0)), 1e-12+1e-10*np.linalg.norm(f, axis=1).sum())
        self.assertLess(np.linalg.norm(np.cross(x, f).sum(axis=0)), 1e-14+1e-10*np.linalg.norm(np.cross(x, f), axis=1).sum())

    def test_05_nonzero_initial_stress_rigid_rotation_independent_matrix(self):
        s0 = np.array([[[11e6, 2e6, 0], [2e6, -7e6, 1e6], [0, 1e6, -4e6]]])
        spin = np.array([[[0., -3., 0.], [3., 0., 0.], [0., 0., 0.]]])
        angle = .003
        r = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        s, ep, plastic, _, _ = update(s0, np.zeros(1), spin, np.array([293.15]), .001, 26.3e9, FKS())
        np.testing.assert_allclose(s[0], r@s0[0]@r.T, rtol=0, atol=1e-6)
        np.testing.assert_array_equal(ep, [0.]); np.testing.assert_array_equal(plastic, [0.])

    def test_06_force_transmits_to_unloaded_interior_in_three_dimensions(self):
        b = c.Block(5e-9); b.advance(); b.advance(); b.advance()
        free = b.pa.boundary_role == 0
        self.assertGreater(np.linalg.norm(velocities(b.pa)[free]), 1e-10)
        self.assertGreater(float(np.max(np.abs(b.pa.w[free]))), 1e-12)
        s = b.pa.S.reshape(-1, 3, 3)
        self.assertGreater(float(np.max(np.abs(s[:, 1, 2]))), 1.)
        self.assertGreater(float(np.max(np.abs(s[:, 0, 2]))), 1.)

    def test_07_observation_label_swap_and_zero_no_physics_change(self):
        a, b, d = [c.Block(5e-9) for _ in range(3)]
        b.pa.tracer_type[:] = 99-b.pa.tracer_type; b.pa.wholebar_tag[:] = 17
        d.pa.tracer_type[:] = 0; d.pa.wholebar_tag[:] = 0
        for _ in range(4): a.advance(); b.advance(); d.advance()
        for candidate in (b, d):
            for k in ('x', 'y', 'z', 'u', 'v', 'w', 'rho', 'S', 'ep', 'plastic_work'):
                np.testing.assert_array_equal(getattr(a.pa, k), getattr(candidate.pa, k))
            self.assertEqual(a.cumulative, candidate.cumulative)

    def test_08_nonzero_state_checkpoint_two_step_continuation(self):
        a = c.Block(5e-9); a.advance(); a.advance(); a.advance()
        self.assertGreater(np.linalg.norm(a.pa.S), 1.)
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; a.checkpoint(path)
            b = c.Block.from_checkpoint(path)
            for _ in range(2): a.advance(); b.advance()
        for k in a.pa.properties:
            np.testing.assert_array_equal(getattr(a.pa, k), getattr(b.pa, k))
        self.assertEqual(a.cumulative, b.cumulative)
        self.assertEqual(a.time, b.time)

    def test_09_nonzero_state_permutation_tracks_roles_and_ids(self):
        a = c.Block(5e-9); a.advance(); a.advance(); a.advance()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; a.checkpoint(path)
            pa, state = restore(path)
            order = np.random.default_rng(71).permutation(64).astype(np.int32)
            b = c.Block(5e-9, pa=pa.extract_particles(order), state=state)
            for _ in range(2): a.advance(); b.advance()
        ids = np.argsort(b.pa.particle_id)
        for k in ('particle_id', 'tracer_type', 'boundary_role', 'm'):
            np.testing.assert_array_equal(getattr(a.pa, k), getattr(b.pa, k)[ids])
        for k in ('x', 'y', 'z', 'u', 'v', 'w', 'rho', 'ep'):
            np.testing.assert_allclose(getattr(a.pa, k), getattr(b.pa, k)[ids], rtol=0, atol=1e-12)
        np.testing.assert_allclose(a.pa.S.reshape(64, 9), b.pa.S.reshape(64, 9)[ids], rtol=0, atol=1e-6)
        for k in a.cumulative: self.assertAlmostEqual(a.cumulative[k], b.cumulative[k], delta=1e-12)

    def test_10_rejected_material_and_kick_candidate_rolls_back(self):
        b = c.Block(5e-9); b.advance(); b.advance()
        before = {k: getattr(b.pa, k).copy() for k in b.pa.properties}
        ledger = copy.deepcopy(b.cumulative); time = b.time
        accepted_journal = b.last_journal; accepted_traces = b.last_traces
        with patch.object(c, 'kick', side_effect=ValueError('INJECTED_CONTACT_FAILURE')):
            with self.assertRaisesRegex(ValueError, 'INJECTED_CONTACT_FAILURE'): b.advance()
        for k, v in before.items(): np.testing.assert_array_equal(getattr(b.pa, k), v)
        self.assertEqual(b.cumulative, ledger); self.assertEqual(b.time, time)
        self.assertEqual(b.rejected_candidate['classification'], 'REJECTED')
        self.assertFalse(b.rejected_candidate['restart_allowed'])
        self.assertIs(b.last_journal, accepted_journal)
        self.assertIs(b.last_traces, accepted_traces)
        self.assertEqual(b.rejected_candidate['journal']['step'], b.step+1)
        self.assertEqual(b.rejected_candidate['journal']['stage'], 'CANDIDATE_MATERIAL_BEFORE_CONTACT')
        self.assertEqual(c.rejected_jsonable({'v': np.array([np.nan, np.inf])}),
            {'v': [{'nonfinite_float': 'NaN'}, {'nonfinite_float': 'Infinity'}]})

    def test_11_vtk_contains_solved_material_not_fake_thermal(self):
        b = c.Block(5e-9); b.advance(); b.advance()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'frame.vtk'; b.vtk(path); v = read_point_cloud(path)
        np.testing.assert_array_equal(v['points'], positions(b.pa))
        np.testing.assert_array_equal(v['scalars']['particle_id'], b.pa.particle_id)
        self.assertIn('density_kg_per_m3', v['scalars'])
        self.assertIn('stress_yz_Pa', v['scalars'])
        self.assertFalse(any('temperature' in k.lower() for k in v['scalars']))
        self.assertNotIn('internal_energy_J_per_kg', v['scalars'])

    def test_12_acoustic_bound_and_no_unregistered_dt(self):
        b = c.Block(5e-9)
        bulk = 2*26.3e9*1.33/(3*.34)
        bound = .1*.0005*1.3/np.sqrt((bulk+4*26.3e9/3)/2700.)
        self.assertAlmostEqual(b.stable_dt(), bound, delta=1e-20)
        self.assertLessEqual(5e-9, bound)
        with self.assertRaises(ValueError): c.Block(1e-6)

    def test_13_rejected_checkpoint_cannot_restart(self):
        b = c.Block(5e-9); b.advance()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; b.checkpoint(path)
            pa, state = restore(path)
            self.assertFalse(state['production_authorized'])
            state['classification'] = 'REJECTED'; state['restart_allowed'] = False
            checkpoint(Path(td)/'rejected.npz', pa, state)
            with self.assertRaises(ValueError): c.Block.from_checkpoint(Path(td)/'rejected.npz')

    def test_14_post_material_acceptance_failure_restores_nonzero_stress(self):
        # Failure after S/rho/x/v updates tests full atomicity, not only a
        # contact exception before material properties have been assigned.
        b = c.Block(5e-9); b.advance(); b.advance()
        before = {k: getattr(b.pa, k).copy() for k in b.pa.properties}
        old_time = b.time; old_step = b.step; ledger = dict(b.cumulative)
        with patch.object(b, '_validate_candidate', side_effect=ValueError('INJECTED_ACCEPTANCE_FAILURE')):
            with self.assertRaisesRegex(ValueError, 'INJECTED_ACCEPTANCE_FAILURE'): b.advance()
        self.assertFalse(np.array_equal(b.rejected_candidate['properties']['S'], before['S']))
        self.assertFalse(np.array_equal(b.rejected_candidate['properties']['u'], before['u']))
        for key, values in before.items(): np.testing.assert_array_equal(getattr(b.pa, key), values)
        self.assertEqual(b.time, old_time); self.assertEqual(b.step, old_step)
        self.assertEqual(b.cumulative, ledger)
        self.assertFalse(b.rejected_candidate['restart_allowed'])


if __name__ == '__main__': unittest.main(verbosity=2)
