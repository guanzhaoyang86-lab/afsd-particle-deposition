"""Independent synthetic observation references; these are not AFSD trials."""
import copy
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.transport_observer import observe


def fixture():
    p = dict(layer=1., traverse=.2, dx=.1, bar_length=3., actuator_depth=.5,
             feed=.1, release_time=10.)
    xyz = np.array([[0., 0., 1.2], [.2, .1, .8], [.4, -.1, 2.8], [0., 0., -.2]])
    a = SimpleNamespace(particle_id=np.array([8, 2, 6, 4]), component_id=np.array([1, 1, 1, 0]),
        material_id=np.ones(4, int), m=np.array([2., 3., 5., 7.]), rho=np.array([4., 6., 10., 14.]),
        ep=np.array([.1, 0., 0., 0.]), tracer_type=np.array([1, 2, 0, 0]),
        wholebar_tag=np.array([1, 2, 2, 0]), as_rs_tag=np.array([0, -1, 1, 0]),
        actuator=np.array([0, 0, 1, 0]))
    for i, axis in enumerate('xyz'):
        setattr(a, axis, xyz[:, i].copy()); setattr(a, axis+'0', xyz[:, i].copy())
    return a, copy.deepcopy(a), p


class TransportTests(unittest.TestCase):
    def test_linear_displacement_exact_mass_volume_reference(self):
        a, b, p = fixture(); b.x += .1; b.z[:2] += [-.4, .4]; b.ep[0] = .2
        r = observe(a, b, 0., 1., p)
        self.assertEqual(r['per_particle']['particle_id'], [2, 4, 6, 8])
        self.assertEqual(r['crossing']['downward']['particle_ids'], [8])
        self.assertEqual(r['crossing']['upward']['particle_ids'], [2])
        self.assertEqual(r['crossing']['signed_downward_event_mass_kg'], -1.)
        self.assertEqual(r['groups']['bar']['mass_kg'], 10.)
        self.assertEqual(r['groups']['bar']['volume_before_m3'], 1.5)
        np.testing.assert_allclose(r['groups']['free_bar']['lab_mass_weighted_displacement_m'], [.1, 0., .08], atol=1e-15)
        self.assertEqual(r['plasticity']['continued_increase']['particle_ids'], [8])
        self.assertEqual(r['plasticity']['continued_increase']['mass_kg'], 2.)
        self.assertEqual(r['deposition_status'], 'NOT_MEASURED')

    def test_independent_permutations_preserve_all_observation(self):
        a, b, p = fixture(); b.x += .03; expected = observe(a, b, 0., 1., p)
        for pa, permutation in ((a, [3, 0, 1, 2]), (b, [2, 1, 3, 0])):
            for key, value in vars(pa).items(): setattr(pa, key, value[permutation])
        self.assertEqual(observe(a, b, 0., 1., p), expected)

    def test_pure_reference_translation_not_lab_motion(self):
        a, b, p = fixture(); r = observe(a, b, 2., 3., p)
        np.testing.assert_array_equal(r['per_particle']['lab_displacement_m'], np.zeros((4, 3)))
        np.testing.assert_allclose(r['per_particle']['tool_translation_displacement_m'], np.tile([-.2, 0., 0.], (4, 1)), atol=1e-15)
        self.assertEqual(r['crossing']['status'], 'NO_CENTER_CROSSING_NOT_DEPOSITION')
        self.assertEqual(r['macro_transport_status'], 'NOT_ESTABLISHED_BY_ENDPOINT_MOTION')

    def test_micro_motion_has_no_macro_or_deposition_promotion(self):
        a, b, p = fixture(); b.x += 1e-12; r = observe(a, b, 0., 1., p)
        self.assertLess(r['groups']['free_bar']['max_lab_displacement_over_dx'], 1.1e-11)
        self.assertEqual(r['macro_transport_status'], 'NOT_ESTABLISHED_BY_ENDPOINT_MOTION')
        self.assertEqual(r['crossing']['status'], 'NO_CENTER_CROSSING_NOT_DEPOSITION')

    def test_repeated_opposite_crossings_deduplicate_unique_mass_not_events(self):
        a, b, p = fixture(); b.z[0] = .9
        first = observe(a, b, 0., 1., p); c = copy.deepcopy(b); c.z[0] = 1.1
        second = observe(b, c, 1., 2., p, first['cumulative']); d = copy.deepcopy(c); d.z[0] = .8
        last = observe(c, d, 2., 3., p, second['cumulative'])
        self.assertEqual(last['cumulative']['downward_unique_mass_kg'], 2.)
        self.assertEqual(last['cumulative']['downward_event_mass_kg'], 4.)
        self.assertEqual(last['cumulative']['upward_event_mass_kg'], 2.)
        self.assertEqual(last['cumulative']['downward_event_count'], 2)
        with self.assertRaisesRegex(ValueError, 'endpoint time'): observe(c, d, 2., 3., p, first['cumulative'])

    def test_plane_endpoint_half_open_convention(self):
        a, b, p = fixture(); b.z[0] = 1.; r = observe(a, b, 0., 1., p)
        self.assertEqual(r['crossing']['downward']['particle_count'], 0)
        c = copy.deepcopy(b); c.z[0] = .99; r = observe(b, c, 1., 2., p)
        self.assertEqual(r['crossing']['downward']['particle_ids'], [8])
        r = observe(c, b, 2., 3., p)
        self.assertEqual(r['crossing']['upward']['particle_ids'], [8])

    def test_tags_do_not_change_group_crossing_or_plastic_measurements(self):
        a, b, p = fixture(); b.z[0] -= .4; expected = observe(a, b, 0., 1., p)
        for pa in (a, b):
            pa.tracer_type[:] = -9; pa.as_rs_tag *= -1; pa.actuator[:] = 1
        changed = observe(a, b, 0., 1., p)
        for key in ('groups', 'crossing', 'plasticity', 'cumulative'): self.assertEqual(changed[key], expected[key])
        self.assertEqual(changed['groups']['free_bar']['particle_ids'], [2, 8])

    def test_does_not_mutate_or_alias_inputs_history(self):
        a, b, p = fixture(); snapshots = [copy.deepcopy(vars(q)) for q in (a, b)]
        r = observe(a, b, 0., 1., p); history = copy.deepcopy(r['cumulative'])
        observe(b, b, 1., 2., p, history); self.assertEqual(history, r['cumulative'])
        r['per_particle']['position_after_m'][0][0] = 999.
        for pa, snapshot in zip((a, b), snapshots):
            for key in snapshot: np.testing.assert_array_equal(getattr(pa, key), snapshot[key])

    def test_mass_id_component_origin_changes_fail_not_invent_inflow(self):
        for key in ('m', 'particle_id', 'component_id', 'material_id', 'x0', 'y0', 'z0'):
            with self.subTest(key=key):
                a, b, p = fixture(); getattr(b, key)[0] += 1
                with self.assertRaises(ValueError): observe(a, b, 0., 1., p)
        a, b, p = fixture()
        for key, value in vars(b).items(): setattr(b, key, value[:-1])
        with self.assertRaisesRegex(ValueError, 'identity/mass'): observe(a, b, 0., 1., p)

    def test_invalid_fields_and_plastic_history_fail(self):
        for key, value in (('rho', 0.), ('x', np.nan), ('ep', -.1), ('particle_id', 2)):
            a, b, p = fixture(); getattr(b, key)[0] = value
            with self.assertRaises(ValueError): observe(a, b, 0., 1., p)
        a, b, p = fixture(); b.ep[0] = .09
        with self.assertRaisesRegex(ValueError, 'decreased'): observe(a, b, 0., 1., p)


if __name__ == '__main__': unittest.main(verbosity=2)
