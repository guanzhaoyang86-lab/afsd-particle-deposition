"""Geometry-only independent protocol tests; contains no experimental data."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_validation.measurement import MeasurementError, measure_proxies


def fixture():
    # Synthetic SI positions; rows 0,1 AS and row2 RS are in the selected slab.
    points = [[.002, -.001, 0.], [.003, -.002, 0.], [.004, .001, 0.],
              [.005, 0., 0.], [.006, .016, 0.], [.007, -.001, .002],
              [.020, -.001, 0.]]
    scalars = dict(particle_id=list(range(10, 17)), material_id=[1]*7,
                   tracer_type=[1]*7, wholebar_tag=[2]*7,
                   mass_kg=[2., 2., 1., 1., 1., 1., 1.], density_kg_per_m3=[2.]*7)
    spec = dict(frame='LAB', label_scope='LOCAL', tracer_value=1, material_id=1,
                roi_x_m=[.001, .010], section_origin_m=[0., 0., 0.],
                section_normal=[0., 0., 1.], section_thickness_m=.001,
                startup_x_m=[0., .001], termination_x_m=[.01, .02],
                region='STEADY_ASSUMED', roi_source='independent geometry fixture, not paper ROI',
                roi_status='ASSUMED', time_s=.1, omega_rad_s=10., traverse_m_s=.00212,
                tool_axis_y_m=0., side_band_m=.01588)
    return points, scalars, spec


class MeasurementReview(unittest.TestCase):
    def test_independent_geometric_count_and_volume_are_not_area(self):
        points, fields, spec = fixture()
        result = measure_proxies(points, fields, spec)
        self.assertEqual(result['R_count']['value'], 2.)
        self.assertEqual(result['R_volume']['value'], 4.)
        self.assertEqual(result['selected_particle_ids'], {'AS': [10, 11], 'RS': [12], 'CENTERLINE': [13]})
        self.assertEqual(result['experimental_area_ratio']['status'], 'BLOCKED')
        self.assertIsNone(result['experimental_area_ratio']['value'])
        self.assertFalse(result['steady_state_verified'])

    def test_reverse_rotation_swaps_by_cross_product_not_names(self):
        points, fields, spec = fixture(); spec['omega_rad_s'] *= -1
        result = measure_proxies(points, fields, spec)
        self.assertEqual(result['R_count']['value'], .5)
        self.assertEqual(result['R_volume']['value'], .25)
        self.assertEqual(result['selected_particle_ids']['AS'], [12])

    def test_zero_denominator_and_absent_marker_are_undefined_not_zero(self):
        points, fields, spec = fixture(); fields['tracer_type'][2] = 0
        result = measure_proxies(points, fields, spec)
        self.assertIsNone(result['R_count']['value']); self.assertIsNone(result['R_volume']['value'])
        self.assertEqual(result['R_count']['status'], 'UNDEFINED')
        fields['tracer_type'] = [0]*7
        result = measure_proxies(points, fields, spec)
        self.assertIsNone(result['R_count']['value'])
        self.assertEqual(result['R_count']['numerator_AS'], 0)

    def test_local_and_wholebar_are_separate_explicit_scopes(self):
        points, fields, spec = fixture()
        fields['wholebar_tag'][0] = 0
        spec['label_scope'] = 'WHOLEBAR_DIAGNOSTIC'; spec['tracer_value'] = 2
        result = measure_proxies(points, fields, spec)
        self.assertEqual(result['R_count']['value'], 1.)
        self.assertEqual(result['label_field'], 'wholebar_tag')
        spec['label_scope'] = 'LOCAL'; spec['tracer_value'] = 1
        self.assertEqual(measure_proxies(points, fields, spec)['R_count']['value'], 2.)

    def test_mixed_scopes_values_and_missing_tracer_fail(self):
        for value in ('LOCAL+WHOLEBAR', ['LOCAL', 'WHOLEBAR_DIAGNOSTIC']):
            points, fields, spec = fixture(); spec['label_scope'] = value
            with self.assertRaises(MeasurementError): measure_proxies(points, fields, spec)
        points, fields, spec = fixture(); spec['tracer_value'] = [1, 2]
        with self.assertRaises(MeasurementError): measure_proxies(points, fields, spec)
        points, fields, spec = fixture(); del fields['tracer_type']
        with self.assertRaises(MeasurementError): measure_proxies(points, fields, spec)

    def test_missing_or_unknown_geometry_never_guessed(self):
        for field in ('roi_x_m', 'section_origin_m', 'section_normal', 'section_thickness_m',
                      'startup_x_m', 'termination_x_m', 'roi_source', 'time_s'):
            points, fields, spec = fixture(); del spec[field]
            with self.assertRaises(MeasurementError): measure_proxies(points, fields, spec)
        points, fields, spec = fixture(); spec['roi_x_m'] = None
        with self.assertRaises(MeasurementError): measure_proxies(points, fields, spec)

    def test_frozen_axis_band_and_explicit_regions(self):
        for field, bad in (('tool_axis_y_m', .001), ('side_band_m', .01), ('omega_rad_s', 0),
                           ('frame', 'CO_ROTATING'), ('section_normal', [0., 0., 2.]),
                           ('roi_x_m', [0., .005])):
            points, fields, spec = fixture(); spec[field] = bad
            with self.assertRaises(MeasurementError): measure_proxies(points, fields, spec)
        points, fields, spec = fixture(); spec['region'] = 'STARTUP'; spec['roi_x_m'] = [0., .001]
        self.assertEqual(measure_proxies(points, fields, spec)['R_count']['status'], 'UNDEFINED')

    def test_no_input_mutation_and_permutation_invariance(self):
        points, fields, spec = fixture(); before = deepcopy((points, fields, spec))
        first = measure_proxies(points, fields, spec)
        self.assertEqual((points, fields, spec), before)
        second = measure_proxies(list(reversed(points)), {k: list(reversed(v)) for k, v in fields.items()}, spec)
        self.assertEqual(first, second)

    def test_invalid_state_or_duplicate_ids_rejected(self):
        for field, bad in (('mass_kg', 0.), ('density_kg_per_m3', float('nan')), ('particle_id', 11)):
            points, fields, spec = fixture(); fields[field][0] = bad
            with self.assertRaises(MeasurementError): measure_proxies(points, fields, spec)


if __name__ == '__main__':
    unittest.main(verbosity=2)
