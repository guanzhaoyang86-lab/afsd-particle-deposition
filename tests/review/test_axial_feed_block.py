"""Independent short gates for the local axial-feed SPH companion.

At most four actual steps per instance. Release tests explicitly transplant a
nonzero short state to the release boundary; they are local event tests, NOT
evidence that an 800 ns trajectory or a deposition experiment completed.
Literal geometry, initial energy/momentum and impulse formulas do not call the
contact or actuator implementations to manufacture their expected answers.
"""
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
import probe_axial_feed as c
from afsd_backend.pysph.mechanics import positions, velocities, corrected_operators
from afsd_backend.pysph.output import restore, checkpoint
from afsd_scaffold.vtk_legacy import read_point_cloud

MASS = 3.375e-7
KN = 657500.
DT = 5e-9
FEED = .00212
RELEASE = 4e-7
K0 = 1.213488e-11
P0 = np.array([0., 0., -1.1448e-8])


def top_from_initial(pa):
    return np.asarray(pa.z0) > -.0085


class AxialFeedShort(unittest.TestCase):
    def test_01_unscaled_three_dimensional_cube_and_actual_free_axes(self):
        b = c.Block(DT); x = positions(b.pa)
        self.assertEqual(x.shape, (64, 3))
        self.assertEqual(np.linalg.matrix_rank(x-x.mean(axis=0)), 3)
        for axis in range(3): self.assertEqual(len(np.unique(x[:, axis])), 4)
        np.testing.assert_allclose(np.unique(x[:, 0]), [-.00075, -.00025, .00025, .00075], rtol=0., atol=1e-18)
        np.testing.assert_allclose(np.unique(x[:, 2]), [-.00975, -.00925, -.00875, -.00825], rtol=0., atol=2e-18)
        np.testing.assert_allclose(b.pa.m, MASS, rtol=0., atol=1e-21)
        self.assertAlmostEqual(float(b.pa.m.sum()), 2.16e-5, delta=1e-18)
        top = top_from_initial(b.pa)
        self.assertEqual(np.count_nonzero(top), 16)
        b.advance()
        for trace in b.last_traces:
            expected = [1., 1., 0.] if top[trace['particle_index']] else [1., 1., 1.]
            self.assertEqual(list(trace['free_axes']), expected)

    def test_02_literal_initial_energy_momentum_and_unstrained_state(self):
        b = c.Block(DT); top = top_from_initial(b.pa); v = velocities(b.pa)
        expected = np.zeros((64, 3)); expected[top, 2] = -FEED
        np.testing.assert_array_equal(v, expected)
        np.testing.assert_array_equal(b.pa.S, np.zeros(64*9))
        np.testing.assert_array_equal(b.pa.rho, np.full(64, 2700.))
        self.assertAlmostEqual(float(.5*np.sum(b.pa.m[:, None]*v*v)), K0, delta=1e-25)
        np.testing.assert_allclose(np.sum(b.pa.m[:, None]*v, axis=0), P0, rtol=0., atol=1e-22)
        self.assertAlmostEqual(b.record()['kinetic_J'], K0, delta=1e-25)
        self.assertAlmostEqual(b.record()['initial_kinetic_J'], K0, delta=1e-25)
        np.testing.assert_allclose(b.record()['initial_linear_momentum_kg_m_s'], P0, rtol=0., atol=1e-22)
        np.testing.assert_array_equal(b.record()['control']['eligible'], top)
        self.assertAlmostEqual(b.record()['energy_diagnostic_residual_J'], 0., delta=1e-25)
        self.assertLess(abs(b.initial_normal_potential_J), 1e-25)
        np.testing.assert_allclose(b.record()['normal_contact']['normal_force_N'], 0., rtol=0., atol=1e-12)

    def test_03_sources_units_removed_prescribed_load_and_open_tool(self):
        b = c.Block(DT)
        cfg = json.loads((ROOT/'configs/benchmarks/axial_feed.json').read_text(encoding='utf-8'))
        self.assertFalse(cfg['production_authorized'])
        self.assertIn('normal_force', cfg['remove_inherited_parameters'])
        self.assertNotIn('normal_force', b.p)
        for item in cfg['parameters'].values():
            self.assertTrue(all(item.get(k) is not None for k in ('value', 'unit', 'source', 'status')))
        self.assertEqual(b.p['friction'], 0.)
        self.assertEqual(b.p['omega'], 0.)
        self.assertEqual(b.p['bore_width'], .008)
        self.assertEqual(b.p['width'], .002)
        self.assertEqual(b.p['substrate_depth'], .01)
        b.advance()
        for trace in b.last_traces:
            for face in trace['contact_rows']:
                self.assertEqual(int(face['feature_id']), 0)
                np.testing.assert_array_equal(face['normal'], [0., 0., 1.])

    def test_04_only_initial_top_z_is_prescribed_and_reaction_work_is_saved(self):
        b = c.Block(DT); top = top_from_initial(b.pa); b.advance(); b.advance()
        np.testing.assert_allclose(b.pa.w[top], -FEED, rtol=0., atol=1e-12)
        np.testing.assert_allclose(b.pa.z[top]-b.pa.z0[top], -FEED*2*DT, rtol=0., atol=2e-17)
        reactions = []
        for trace in b.last_traces:
            r = trace['result']; ja = np.asarray(r['constraint_impulse'])
            free = np.asarray(trace['free_axes'], dtype=bool)
            np.testing.assert_allclose(ja[free], 0., rtol=0., atol=1e-17)
            work = sum(float(np.dot(a['constraint_impulse'], .5*(a['v_before']+a['v_after'])))
                       for a in r['intervals'])
            self.assertAlmostEqual(r['constraint_work_J'], work, delta=1e-22)
            reactions.append(ja)
        self.assertGreater(float(np.linalg.norm(reactions)), 1e-15)
        self.assertGreater(abs(b.cumulative['constraint_work_J']), 1e-20)

    def test_05_first_kick_literal_impulse_and_projected_mobility_identity(self):
        b = c.Block(DT); b.advance()
        for trace in b.last_traces:
            r = trace['result']; v0 = np.asarray(trace['material_velocity_before_m_s'])
            force = np.asarray(trace['internal_force_N'])
            jn = np.zeros(3)
            for face in trace['contact_rows']:
                jn += DT*float(face['normal_force_N'])*np.asarray(face['normal'])
            free = np.asarray(trace['free_axes'])
            # Independent mu=0 explicit projected impulse formula; force is
            # measured input, not an independently validated constitutive law.
            expected = v0+free*(DT*force+jn)/MASS
            np.testing.assert_allclose(r['velocity'], expected, rtol=0., atol=1e-12)
            expected_reaction = (free-1.)*(DT*force+jn)
            np.testing.assert_allclose(r['constraint_impulse'], expected_reaction, rtol=0., atol=1e-20)
            np.testing.assert_allclose(np.asarray(r['tangential_impulse']), 0., rtol=0., atol=1e-20)
            self.assertEqual(float(np.sum(r['friction_heat_J'])), 0.)

    def test_06_two_steps_transmit_stress_and_motion_to_nonactuated_material(self):
        b = c.Block(DT); top = top_from_initial(b.pa)
        b.advance(); b.advance()
        self.assertGreater(float(np.linalg.norm(velocities(b.pa)[~top])), 1e-9)
        self.assertGreater(float(np.max(np.abs(b.pa.w[~top]))), 1e-9)
        s = b.pa.S.reshape(-1, 3, 3)
        self.assertGreater(float(np.max(np.abs(s[:, 2, 2]))), 1.)
        self.assertGreater(float(np.max(np.abs(s[:, 0, 2]))), 1.)
        self.assertGreater(float(np.max(np.abs(s[:, 1, 2]))), 1.)
        self.assertGreater(float(np.max(np.abs(b.pa.rho-2700.))), 1e-8)
        i, j, g, L, drho, B = corrected_operators(b.pa); x = positions(b.pa)
        for k in range(64):
            neighbors = np.r_[j[i == k], i[j == k]]
            self.assertEqual(np.linalg.matrix_rank(x[neighbors]-x[k]), 3)
        self.assertLess(float(np.max(np.linalg.cond(np.linalg.inv(B)))), 1e6)

    def test_07_kinetic_balance_subtracts_initial_energy_and_keeps_normal_work(self):
        b = c.Block(DT); b.advance(); b.advance(); r = b.record(); d = b.cumulative
        residual = (r['kinetic_J']-K0+d['friction_raw_J']-d['wall_work_J']-
                    d['constraint_work_J']-d['internal_work_J']-d['normal_relative_work_J'])
        self.assertAlmostEqual(residual, 0., delta=1e-18)
        self.assertEqual(d['friction_raw_J'], 0.)
        self.assertEqual(d['wall_work_J'], 0.)
        self.assertEqual(d['normal_wall_work_J'], 0.)
        normal = r['normal_contact']
        du = float(np.sum(normal['potential_J']))-b.initial_normal_potential_J
        self.assertAlmostEqual(d['normal_temporal_residual_J'], du+d['normal_relative_work_J'], delta=1e-18)
        momentum = np.sum(b.pa.m[:, None]*velocities(b.pa), axis=0)
        external = np.array([d[f'base_impulse_{axis}_Ns']+d[f'actuator_impulse_{axis}_Ns']
                             for axis in 'xyz'])
        np.testing.assert_allclose(momentum-P0, external, rtol=0., atol=1e-18)
        angular = np.sum(b.pa.m[:, None]*np.cross(positions(b.pa), velocities(b.pa)), axis=0)
        external_l = np.array([d[f'base_angular_impulse_{axis}_Nms']+d[f'actuator_angular_impulse_{axis}_Nms']
                               for axis in 'xyz'])
        np.testing.assert_allclose(angular, external_l, rtol=0., atol=1e-21)

    def test_08_tags_roles_and_cached_actuator_flag_are_observational(self):
        a, b, d = [c.Block(DT) for _ in range(3)]
        b.pa.tracer_type[:] = 99-b.pa.tracer_type; b.pa.wholebar_tag[:] = 17; b.pa.as_rs_tag[:] = -5
        b.pa.boundary_role[:] = 1; b.pa.actuator[:] = 0
        d.pa.tracer_type[:] = 0; d.pa.wholebar_tag[:] = 0; d.pa.as_rs_tag[:] = 0
        d.pa.boundary_role[:] = 0; d.pa.actuator[:] = 1
        for _ in range(2): a.advance(); b.advance(); d.advance()
        for other in (b, d):
            for key in ('x', 'y', 'z', 'u', 'v', 'w', 'rho', 'S', 'ep', 'plastic_work'):
                np.testing.assert_array_equal(getattr(a.pa, key), getattr(other.pa, key))
            self.assertEqual(a.cumulative, other.cumulative)

    def test_09_nonzero_restart_preserves_fields_initial_energy_and_ledger(self):
        a = c.Block(DT); a.advance(); a.advance()
        self.assertGreater(float(np.linalg.norm(a.pa.S)), 1.)
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; a.checkpoint(path); b = c.Block.from_checkpoint(path)
            for _ in range(2): a.advance(); b.advance()
        for key in a.pa.properties: np.testing.assert_array_equal(getattr(a.pa, key), getattr(b.pa, key))
        self.assertEqual(a.cumulative, b.cumulative)
        self.assertEqual(a.time, b.time); self.assertEqual(a.step, b.step)
        self.assertEqual(a.record()['energy_diagnostic_residual_J'], b.record()['energy_diagnostic_residual_J'])

    def test_10_permutation_follows_original_id_and_material_subset(self):
        a = c.Block(DT); a.advance(); a.advance()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; a.checkpoint(path); pa, state = restore(path)
            permutation = np.random.default_rng(73).permutation(64).astype(np.int32)
            b = c.Block(DT, pa=pa.extract_particles(permutation), state=state)
            for _ in range(2): a.advance(); b.advance()
        order = np.argsort(b.pa.particle_id)
        for key in ('particle_id', 'tracer_type', 'm', 'x0', 'y0', 'z0'):
            np.testing.assert_array_equal(getattr(a.pa, key), getattr(b.pa, key)[order])
        for key in ('x', 'y', 'z', 'u', 'v', 'w', 'rho', 'ep'):
            np.testing.assert_allclose(getattr(a.pa, key), getattr(b.pa, key)[order], rtol=0., atol=1e-12)
        np.testing.assert_allclose(a.pa.S.reshape(64, 9), b.pa.S.reshape(64, 9)[order], rtol=0., atol=1e-6)
        for key in a.cumulative:
            np.testing.assert_allclose(a.cumulative[key], b.cumulative[key], rtol=0., atol=1e-12)

    def test_11_rejected_candidate_rolls_back_all_fields_and_accepted_evidence(self):
        b = c.Block(DT); b.advance(); b.advance()
        before = {key: getattr(b.pa, key).copy() for key in b.pa.properties}
        ledger = copy.deepcopy(b.cumulative); old_time = b.time; old_step = b.step
        journal, traces = b.last_journal, b.last_traces; validate = b._validate_candidate
        def reject_after_accounting(increment):
            validate(increment)
            raise ValueError('INJECTED_AFTER_AXIAL_LEDGER')
        with patch.object(b, '_validate_candidate', side_effect=reject_after_accounting):
            with self.assertRaisesRegex(ValueError, 'INJECTED_AFTER_AXIAL_LEDGER'): b.advance()
        self.assertFalse(np.array_equal(b.rejected_candidate['properties']['S'], before['S']))
        for key, value in before.items(): np.testing.assert_array_equal(getattr(b.pa, key), value)
        self.assertEqual(b.cumulative, ledger); self.assertEqual(b.time, old_time); self.assertEqual(b.step, old_step)
        self.assertIs(b.last_journal, journal); self.assertIs(b.last_traces, traces)
        self.assertEqual(b.rejected_candidate['classification'], 'REJECTED')
        self.assertFalse(b.rejected_candidate['restart_allowed'])
        self.assertIn('normal_energy', b.rejected_candidate['journal'])

    def test_12_rejected_checkpoint_cannot_restart(self):
        b = c.Block(DT); b.advance()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'accepted.npz'; b.checkpoint(path); pa, state = restore(path)
            self.assertFalse(state['production_authorized'])
            state['classification'] = 'REJECTED'; state['restart_allowed'] = False
            rejected = Path(td)/'rejected.npz'; checkpoint(rejected, pa, state)
            with self.assertRaises(ValueError): c.Block.from_checkpoint(rejected)

    def test_13_acoustic_and_penalty_timestep_bounds_use_unscaled_parameters(self):
        b = c.Block(DT); bulk = 2*26.3e9*1.33/(3*.34)
        acoustic = .1*.0005*1.3/(np.sqrt((bulk+4*26.3e9/3)/2700.)+FEED)
        normal = .1*np.sqrt(MASS/KN)
        self.assertAlmostEqual(b.stable_dt(), min(acoustic, normal), delta=1e-20)
        self.assertLessEqual(DT, b.stable_dt())
        with self.assertRaises(ValueError): c.Block(1e-6)

    def test_14_vtk_has_solved_three_dimensional_fields_not_fabricated_heat(self):
        b = c.Block(DT); b.advance(); b.advance()
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'frame.vtk'; b.vtk(path); vtk = read_point_cloud(path)
        np.testing.assert_array_equal(vtk['points'], positions(b.pa))
        np.testing.assert_array_equal(vtk['scalars']['particle_id'], b.pa.particle_id)
        for key in ('density_kg_per_m3', 'stress_xz_Pa', 'stress_yz_Pa',
                    'normal_force_N', 'normal_potential_J', 'contact_active', 'signed_gap_m'):
            self.assertIn(key, vtk['scalars'])
        self.assertFalse(any('temperature' in key.lower() for key in vtk['scalars']))
        self.assertNotIn('internal_energy_J_per_kg', vtk['scalars'])

    def test_15_local_release_event_never_resets_velocity_or_keeps_reaction(self):
        # Synthetic time transplantation only. This is not a completed 400 ns run.
        b = c.Block(DT); b.advance(); b.advance(); v = velocities(b.pa).copy()
        b.time = RELEASE; b.step = int(round(RELEASE/DT))
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'local_release.npz'; b.checkpoint(path); restored = c.Block.from_checkpoint(path)
            np.testing.assert_array_equal(velocities(restored.pa), v)
            b.advance(); restored.advance()
        np.testing.assert_array_equal(velocities(b.pa), velocities(restored.pa))
        self.assertGreater(float(np.linalg.norm(velocities(b.pa))), 0.)
        for trace in b.last_traces:
            self.assertEqual(list(trace['free_axes']), [1., 1., 1.])
            np.testing.assert_array_equal(trace['material_velocity_before_m_s'], v[trace['particle_index']])
            np.testing.assert_allclose(trace['result']['constraint_impulse'], 0., rtol=0., atol=1e-17)
            self.assertAlmostEqual(trace['result']['constraint_work_J'], 0., delta=1e-22)

    def test_16_interval_straddling_release_rejected_before_any_mutation(self):
        b = c.Block(DT); b.advance(); b.time = RELEASE-.5*DT
        before = {key: getattr(b.pa, key).copy() for key in b.pa.properties}
        ledger = copy.deepcopy(b.cumulative); old_time = b.time; old_step = b.step
        with self.assertRaisesRegex(ValueError, 'release.*boundary'): b.advance()
        for key, value in before.items(): np.testing.assert_array_equal(getattr(b.pa, key), value)
        self.assertEqual(b.cumulative, ledger); self.assertEqual(b.time, old_time); self.assertEqual(b.step, old_step)


if __name__ == '__main__': unittest.main(verbosity=2)
