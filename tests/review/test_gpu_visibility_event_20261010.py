"""Independent coupled-operator fixtures and bounded event routing regressions.

Manufactured components only; these are not loaded AFSD process validation.
"""
from pathlib import Path
import copy
import json
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import run_gpu_visibility_event_20261010 as event
import afsd_gpu_event_trial_20261010 as trial
from afsd_event_graph_20261010 import restrict_activation_graph
import check_visible_support_operators as manufactured


class CoupledOperators(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.x, cls.volume, h, cls.block, _, _ = manufactured.fixture(3)
        cls.pa = manufactured.particle_array(cls.x, cls.volume, h)
        cls.raw, cls.start = manufactured.graph(cls.pa)
        cls.ids = np.arange(len(cls.x), dtype=np.int64)
        i, j, g = cls.raw
        cls.incoming = [tuple(sorted((int(i[k]), int(j[k]))))
            for k in cls.start['blocked_raw_edge_indices']]
        cls.end_pose = dict(origin=np.array([0., 0., .01]), rotation=np.eye(3))
        cls.delayed = restrict_activation_graph(cls.x, cls.x, i, j, g,
            manufactured.POSE, cls.end_pose, manufactured.PARAMETERS, cls.ids, cls.incoming)
        cls.restored = restrict_activation_graph(cls.x, cls.x, i, j, g,
            cls.end_pose, cls.end_pose, manufactured.PARAMETERS, cls.ids, ())

    def test_single_graph_exactly_recovers_old_side_and_then_full_support(self):
        self.assertGreater(len(self.incoming), 0)
        for key, raw in zip(('i', 'j', 'g'), self.raw):
            np.testing.assert_array_equal(self.delayed[key], self.start[key])
            np.testing.assert_array_equal(self.restored[key], raw)

    def test_mechanical_coupling_is_absent_before_activation_and_restored_after(self):
        velocity = self.x@manufactured.AFFINE.T
        velocity[self.block == 1] = 0.
        stress = manufactured.stress_field(self.x)
        stress[self.block == 1] = 0.
        old_L, old_B = manufactured.geometry(self.x, velocity, self.pa, self.start)
        left_L, left_B = manufactured.geometry(self.x, velocity, self.pa, self.delayed)
        right_L, right_B = manufactured.geometry(self.x, velocity, self.pa, self.restored)
        np.testing.assert_array_equal(left_L, old_L)
        np.testing.assert_array_equal(left_B, old_B)
        left_F = manufactured.force(self.pa, stress, self.delayed, left_B)
        right_F = manufactured.force(self.pa, stress, self.restored, right_B)
        self.assertEqual(float(np.abs(left_L[self.block == 1]).max()), 0.)
        self.assertEqual(float(np.abs(left_F[self.block == 1]).max()), 0.)
        self.assertGreater(float(np.abs(right_L[self.block == 1]).max()), 1e-5)
        self.assertGreater(float(np.abs(right_F[self.block == 1]).max()), 1e-5)
        for graph, gradient, force in ((self.delayed, left_L, left_F), (self.restored, right_L, right_F)):
            self.assertLess(np.linalg.norm(force.sum(axis=0)), 1e-10)
            power = np.sum(velocity*force)+np.sum(self.volume[:, None, None]*stress*gradient)
            self.assertLess(abs(power), 1e-10)

    def test_heat_coupling_and_capacity_balance_across_activation(self):
        temperature = np.where(self.block == 0, 600., 300.)
        left = manufactured.heat(self.pa, self.x, self.delayed, temperature)
        right = manufactured.heat(self.pa, self.x, self.restored, temperature)
        np.testing.assert_allclose(left['temperature'], temperature, rtol=0., atol=1e-10)
        self.assertGreater(float(np.max(temperature[self.block == 0]-right['temperature'][self.block == 0])), 1e-5)
        self.assertGreater(float(np.max(right['temperature'][self.block == 1]-temperature[self.block == 1])), 1e-5)
        for result in (left, right):
            residual = np.sum(self.pa.m*manufactured.CP*(result['temperature']-temperature))
            self.assertLess(abs(residual), 1e-12)
            self.assertLess(float(np.abs(result['linear_system_residual_J']).max()), 1e-12)

    def test_original_gpu_equations_and_tangent_are_inherited(self):
        self.assertEqual(trial.operator_reuse_audit()['status'], 'PASS_EXPLICIT_ONE_SIDED_EVENT_GRAPH_ADAPTATION')
        self.assertIs(trial.ProcessTrial.tangent, trial.parent.ProcessTrial.tangent)
        self.assertIs(trial.advance_candidate, trial.parent.advance_candidate)


class OriginalFailedInterval(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        folder = ROOT/'outputs/gpu_acceleration_20261007/campaign_10h_20261010_v1/attempt_0044/coarse'
        cls.fixture = {key:json.loads((folder/name).read_text(encoding='utf-8')) for key,name in (
            ('routing','visibility_event_routing.json'), ('path','visible_path_audit.json'),
            ('shell','continuous_shell_audit.json'), ('geometry','geometry_audit.json'))}

    def test_actual_old44_has_exactly_the_recognized_incoming_edge(self):
        pairs, guess = event.classify_activation(copy.deepcopy(self.fixture))
        self.assertEqual(pairs, ((1334, 1371),))
        self.assertGreater(guess, .3593)
        self.assertLess(guess, .3594)

    def test_outgoing_event_is_rejected(self):
        item = copy.deepcopy(self.fixture)
        item['routing']['transition_brackets'][0]['direction'] = 'VISIBLE_TO_BLOCKED'
        with self.assertRaises(ValueError):
            event.classify_activation(item)

    def test_query_budget_exhaustion_is_not_activation(self):
        item = copy.deepcopy(self.fixture)
        item['path']['failure_kind'] = 'VISIBLE_PATH_QUERY_BUDGET_EXHAUSTED'
        with self.assertRaises(ValueError):
            event.classify_activation(item)

    def test_additional_intersection_at_another_edge_cannot_be_removed(self):
        item = copy.deepcopy(self.fixture)
        item['path']['found_edges'][0]['particle_id_i'] = 1333
        with self.assertRaises(ValueError):
            event.classify_activation(item)

    def test_unresolved_router_bracket_cannot_authorize_a_probe(self):
        route = copy.deepcopy(self.fixture['routing'])
        self.assertTrue(event.resolved_routing(route))
        route['transition_brackets'][0]['status'] = 'QUERY_BUDGET_EXHAUSTED'
        self.assertFalse(event.resolved_routing(route))

    def test_too_wide_router_bracket_cannot_authorize_a_probe(self):
        route = copy.deepcopy(self.fixture['routing'])
        route['transition_brackets'][0]['bracket_width_fraction'] = 2*route['fraction_tolerance']
        self.assertFalse(event.resolved_routing(route))

    def test_multiple_events_cannot_be_treated_as_one(self):
        item = copy.deepcopy(self.fixture)
        item['routing']['transition_brackets'] *= 2
        with self.assertRaises(ValueError):
            event.classify_activation(item)

    def test_noninitial_penetration_is_not_this_activation_case(self):
        item = copy.deepcopy(self.fixture)
        item['path']['found_edges'][0]['fraction'] = .1
        with self.assertRaises(ValueError):
            event.classify_activation(item)

    def test_contact_miss_cannot_be_fixed_by_graph_activation(self):
        item = copy.deepcopy(self.fixture)
        item['geometry']['transient_contacts_missed_by_endpoints'] = [1]
        with self.assertRaises(ValueError):
            event.classify_activation(item)


if __name__ == '__main__':
    unittest.main(verbosity=2)
