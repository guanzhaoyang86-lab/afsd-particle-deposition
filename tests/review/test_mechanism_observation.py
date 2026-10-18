"""Independent SI geometry/volume/time arithmetic; never generated physics data."""
import copy
import math
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from measure_afsd_mechanisms import observe, process_costs


class MechanismObservations(unittest.TestCase):
    def fixture(self):
        # Invented observation-unit fixture only, NOT AFSD simulation/reference.
        p = dict(layer=1., bore_width=2., omega=0., traverse=0., feed=.05, dx=.1)
        data = dict(points=[[1.5, 0., .5], [.2, 0., .5], [1.5, 0., .5]],
                    vectors={'initial_position_m': [[0., 0., 1.2], [.2, 0., .5], [0., 0., 1.2]],
                             'velocity_m_per_s': [[0., 0., 0.]]*3},
                    scalars={'particle_id': [7, 11, 19], 'component_id': [1, 0, 1],
                             'actuator': [0, 0, 1], 'mass_kg': [2., 3., 2.],
                             'density_kg_per_m3': [4., 6., 4.],
                             'equivalent_plastic_strain': [1e-5, 1e-11, 1e-2],
                             'tracer_type': [1, 2, 9], 'wholebar_tag': [3, 4, 5]})
        return data, p

    def test_mass_over_density_is_volume_not_count_or_area(self):
        data, p = self.fixture()
        r = observe(data, p, 0.)
        self.assertEqual(r['outside_bore_layer_candidate_particles'], 1)
        self.assertEqual(r['outside_bore_layer_candidate_mass_kg'], 2.)
        self.assertEqual(r['outside_bore_layer_candidate_volume_m3'], .5)  #2kg/(4kg/m3)
        self.assertEqual(r['net_initially_above_layer_currently_below_mass_kg'], 4.)
        self.assertFalse(r['sustained_deposition_validated'])

    def test_free_plastic_observation_keeps_raw_below_floor_values(self):
        data, p = self.fixture()
        r = observe(data, p, 0.)
        self.assertEqual(r['free_plastic_raw_positive_count'], 2)
        self.assertEqual(r['free_plastic_above_floor_count'], 1)
        self.assertEqual(r['free_plastic_particle_ids'], [7])
        self.assertEqual(r['maximum_ep'], .01)  #actuator still part of raw field
        data['scalars']['equivalent_plastic_strain'][1] = -1e-20
        self.assertEqual(observe(data, p, 0.)['negative_ep_count'], 1)  #never clipped

    def test_labels_and_reordering_do_not_change_observables(self):
        data, p = self.fixture(); baseline = observe(data, p, 0.)
        changed = copy.deepcopy(data)
        changed['scalars']['tracer_type'] = [0, 99, -1]
        del changed['scalars']['wholebar_tag']
        self.assertEqual(observe(changed, p, 0.), baseline)
        order = [2, 0, 1]
        changed['points'] = [changed['points'][i] for i in order]
        for group in ('scalars', 'vectors'):
            for k, values in changed[group].items():
                changed[group][k] = [values[i] for i in order]
        self.assertEqual(observe(changed, p, 0.), baseline)

    def test_translating_rotating_bore_query_has_known_geometry(self):
        data, p = self.fixture()
        p.update(omega=math.pi/4, traverse=2.)
        # At1s rotate square45deg, translate2m. Point1.2m along LAB x from
        # toolcentre is INSIDE: its local components are1.2/sqrt2<1.
        data['points'][0] = [3.2, 0., .5]
        r = observe(data, p, 1.)
        self.assertEqual(r['outside_bore_layer_candidate_particles'], 0)
        self.assertEqual(r['prescribed_feed_m'], .05)

    def test_multiple_plastic_ids_are_canonical_under_reordering(self):
        data, p = self.fixture()
        data['scalars']['equivalent_plastic_strain'][1] = 2e-5
        data['scalars']['particle_id'] = [11, 7, 19]
        reference = observe(data, p, 0.)
        self.assertEqual(reference['free_plastic_particle_ids'], [7, 11])
        data['points'].reverse()
        for group in ('scalars', 'vectors'):
            for value in data[group].values():
                value.reverse()
        self.assertEqual(observe(data, p, 0.), reference)

    def test_closed_si_cost_references_and_zero_rate(self):
        p = dict(dx=.1, layer=.2, omega=2*math.pi, feed=.05, traverse=.1)
        costs = {r['event']: r for r in process_costs(p, .01, .2)}
        # One revolution1s/0.01s=100steps; each.2s gives20wallseconds.
        self.assertEqual(costs['tool_1revolution']['estimated_steps'], 100)
        self.assertEqual(costs['tool_1revolution']['estimated_wall_seconds'], 20.)
        self.assertEqual(costs['feed_1dx']['physical_seconds'], 2.)
        self.assertEqual(costs['feed_1dx']['estimated_wall_seconds'], 40.)
        self.assertEqual(costs['feed_1layer']['estimated_wall_seconds'], 80.)
        p['feed'] = 0.
        unavailable = {r['event']: r for r in process_costs(p, .01, .2)}['feed_1dx']
        self.assertEqual(unavailable['status'], 'NOT_AVAILABLE_ZERO_PRESCRIBED_RATE')
        self.assertNotIn('estimated_wall_seconds', unavailable)

    def test_invalid_fields_and_costs_fail_not_clip(self):
        data, p = self.fixture()
        data['scalars']['density_kg_per_m3'][0] = 0.
        with self.assertRaises(ValueError):
            observe(data, p, 0.)
        with self.assertRaises(ValueError):
            process_costs(p, 0., 1.)
        with self.assertRaises(ValueError):
            process_costs(p, .1, float('nan'))


if __name__ == '__main__':
    unittest.main(verbosity=2)
