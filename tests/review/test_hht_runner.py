"""Read-only runner wiring checks; no HHT advance, trajectory or output stage.

Bounds registered in this file before its first run: literal initial arrays,
tagged physical results and reordered contact rows are bitwise checks; total forces
allow 5e-12 N absolute error, accelerations 2e-8 m/s^2, and stresses 1e-9 Pa.
These roundoff allowances cover O(10 N) contact forces, O(1e5 m/s^2)
accelerations and O(1e4 Pa) stress; they are not temporal accuracy gates.
All cases initialize the same literal physical fields at the stated wall-time
layers. They do not claim an evolved physical state at those later times.
"""
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))

import probe_hht_coupled as driver
import reference_hht_block as independent
from afsd_backend.pysph.hht_coupled import initialize
from afsd_backend.pysph.mechanics import positions, velocities
from afsd_backend.pysph.updated_trial import particle_hash


def record(name, **values):
    print('HHT_RUNNER_EVIDENCE '+json.dumps(dict(name=name, trajectory_advanced=False,
          **values), sort_keys=True, allow_nan=False), flush=True)


class RunnerWiring(unittest.TestCase):
    def test_01_literal_fixture_initial_physical_fields(self):
        pa, parameters = driver.setup()
        indices = np.array([(i, j, k) for i in range(3)
                            for j in range(3) for k in range(3)])
        points = .005*indices
        volume = .005**3*np.prod(np.where(indices == 1, 1., .5), axis=1)
        mass = 2700.*volume
        velocity = np.zeros((27, 3)); velocity[:, 2] = .02
        stress = np.tile(np.diag([20000., -10000., -10000.]), (27, 1, 1))
        temperature = 400.+2.*points[:, 0]/.01
        np.testing.assert_array_equal(positions(pa), points)
        np.testing.assert_array_equal(velocities(pa), velocity)
        np.testing.assert_array_equal(pa.m, mass)
        np.testing.assert_array_equal(pa.rho, np.full(27, 2700.))
        np.testing.assert_array_equal(pa.h, np.full(27, 1.3*.005))
        np.testing.assert_array_equal(pa.S.reshape(27, 3, 3), stress)
        np.testing.assert_array_equal(pa.temperature, temperature)
        np.testing.assert_array_equal(pa.ep, np.zeros(27))
        np.testing.assert_array_equal(pa.particle_id, np.arange(1, 28))
        np.testing.assert_array_equal(np.column_stack((pa.x0, pa.y0, pa.z0)), points)
        np.testing.assert_array_equal(pa.tracer_type, np.arange(27) % 3)
        np.testing.assert_array_equal(pa.wholebar_tag, np.arange(27) % 3)
        np.testing.assert_array_equal(pa.actuator, np.zeros(27, dtype=int))
        for key, value in dict(hht_weight=.9, friction=.2, normal_preload=2e-6,
                normal_stroke=2e-7, normal_stiffness=7e6, wall_speed=.2,
                end_time=2e-6, conductivity=167., cp=900., heat_partition=.5,
                contact_conductance=1000., plastic_heat_fraction=1.).items():
            self.assertEqual(parameters[key], value, key)
        self.assertAlmostEqual(float(pa.m.sum()), .0027, delta=2e-18)
        record('literal_inputs', particle_count=27, mass_kg=float(pa.m.sum()),
               initial_nonzero_velocity=True, initial_nonzero_stress=True)

    def test_02_independent_force_and_stress_at_both_control_limits(self):
        pa, parameters = driver.setup(); before = particle_hash(pa)
        data = independent.fixture()
        x, v, S, rho, T, _ = independent.unpack(independent.packed_initial(data))
        fields = independent.pair_fields(x, v, S, rho, data)
        face_limits = {}
        for at, right in ((0., False), (1e-6, False), (1e-6, True), (2e-6, True)):
            with self.subTest(time_s=at, right=right):
                state = initialize(pa, parameters, driver.boundary_for(parameters, right), time=at)
                boundary = independent.boundary_fields(at, x, v, T, data, right)
                expected_force = fields['force']+boundary['force']
                frame = driver.frame(pa, state, {})
                np.testing.assert_allclose(state['force_N'], expected_force, atol=5e-12, rtol=0.)
                np.testing.assert_allclose(state['acceleration_m_s2'],
                    expected_force/data['mass_kg'][:, None], atol=2e-8, rtol=0.)
                np.testing.assert_allclose(frame['face_force_N'], boundary['face_force'],
                                           atol=5e-12, rtol=0.)
                np.testing.assert_allclose(frame['stress_Pa'], fields['stress'], atol=1e-9, rtol=0.)
                self.assertGreater(float(np.linalg.norm(state['acceleration_m_s2'])), 0.)
                self.assertEqual(particle_hash(pa), before)
                self.assertEqual(len(state['boundary_snapshot']['contacts']), 18)
                for contact in state['boundary_snapshot']['contacts']:
                    self.assertEqual(contact['conductance_W_K'], .025)
                    side = int(contact['contact_id'].split(':')[0])
                    self.assertEqual(contact['wall_temperature_K'], [410., 420.][side])
                    self.assertEqual(contact['wall_velocity_m_s'][2], -.2 if right else .2)
                if at == 1e-6:
                    face_limits[right] = frame['face_force_N'].copy()
                record('independent_force_layer', time_s=at, side='RIGHT' if right else 'LEFT',
                    force_error_N=float(np.max(np.abs(state['force_N']-expected_force))),
                    face_force_N=frame['face_force_N'].tolist())
        np.testing.assert_array_equal(face_limits[False][:, 0], face_limits[True][:, 0])
        np.testing.assert_allclose(face_limits[False][:, 2], -face_limits[True][:, 2],
                                   atol=5e-12, rtol=0.)

    def test_03_passive_labels_do_not_enter_initial_physics(self):
        pa, parameters = driver.setup()
        tagged = pa.extract_particles(np.arange(27, dtype=np.uint32))
        tagged.tracer_type[:] = (2*np.arange(27)+1) % 3
        tagged.wholebar_tag[:] = np.arange(27)[::-1] % 5
        changed_hash = particle_hash(tagged)
        for right in (False, True):
            boundary = driver.boundary_for(parameters, right)
            ordinary = initialize(pa, parameters, boundary, time=1e-6)
            changed = initialize(tagged, parameters, boundary, time=1e-6)
            np.testing.assert_array_equal(ordinary['force_N'], changed['force_N'])
            np.testing.assert_array_equal(ordinary['acceleration_m_s2'], changed['acceleration_m_s2'])
            first, second = driver.frame(pa, ordinary, {}), driver.frame(tagged, changed, {})
            for key in ('position_m', 'velocity_m_s', 'stress_Pa', 'density_kg_m3',
                        'temperature_K', 'ep', 'face_force_N'):
                np.testing.assert_array_equal(first[key], second[key], err_msg=key)
        self.assertEqual(particle_hash(tagged), changed_hash)
        self.assertNotEqual(particle_hash(pa), changed_hash)
        record('passive_labels', physical_initialization='BITWISE_EQUAL',
               labels_preserved=True)

    def test_04_particle_reorder_preserves_identity_bound_contact_cohorts(self):
        pa, parameters = driver.setup()
        # This permutation mixes both faces and interior nodes, not just reversal.
        order = (7*np.arange(27)+4) % 27
        reordered = pa.extract_particles(order.astype(np.uint32))
        pristine = particle_hash(reordered)
        np.testing.assert_array_equal(reordered.particle_id, pa.particle_id[order])
        for at, right in ((0., False), (1e-6, True)):
            boundary = driver.boundary_for(parameters, right)
            before = initialize(pa, parameters, boundary, time=at)
            after = initialize(reordered, parameters, boundary, time=at)
            # NNPS accumulation can change its last bits after a permutation;
            # the independent force roundoff allowance applies to total forces.
            np.testing.assert_allclose(after['force_N'], before['force_N'][order],
                                       atol=5e-12, rtol=0.)
            np.testing.assert_allclose(after['acceleration_m_s2'],
                before['acceleration_m_s2'][order], atol=2e-8, rtol=0.)
            rows0 = {row['contact_id']: row for row in before['boundary_snapshot']['contacts']}
            rows1 = {row['contact_id']: row for row in after['boundary_snapshot']['contacts']}
            self.assertEqual(set(rows0), set(rows1))
            for cid, row0 in rows0.items():
                row1 = rows1[cid]; index = row1['particle_index']
                self.assertEqual(row1['particle_id'], row0['particle_id'])
                self.assertEqual(int(reordered.particle_id[index]), row1['particle_id'])
                face = int(cid.split(':')[0])
                self.assertEqual(reordered.x0[index], [0., .01][face])
                for key in ('normal', 'tangential_force_N', 'wall_velocity_m_s',
                            'normal_force_N', 'penetration_m', 'normal_potential_J',
                            'conductance_W_K', 'wall_temperature_K'):
                    np.testing.assert_array_equal(row1[key], row0[key], err_msg=key)
            np.testing.assert_array_equal(driver.frame(pa, before, {})['face_force_N'],
                driver.frame(reordered, after, {})['face_force_N'])
        self.assertEqual(particle_hash(reordered), pristine)
        record('reordered_cohorts', preserved_contacts=18, persistent_ids=True,
               no_position_based_current_cohort_reassignment=True)


if __name__ == '__main__':
    unittest.main(verbosity=2)
