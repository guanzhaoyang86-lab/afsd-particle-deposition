"""Geometry-only graph-policy checks; no force solve, GPU or state advance."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import afsd_event_graph_20261010 as event_graph
from visible_support import filter_visible_support, _hash_array, _topology_hash
from visible_support_remote_domain_20261007 import VisibleSupportNumericalDomainError


class EventGraphTests(unittest.TestCase):
    def setUp(self):
        self.p = dict(tool_radius=2.9, bore_width=2., layer=0., tool_top=2.)
        self.pose = dict(origin=np.zeros(3), rotation=np.eye(3))
        # Edge 0 crosses the solid at the start and is below it at the end;
        # edge 1 is always below it; edge 2 stays geometrically blocked.
        self.x0 = np.array([[.5,0,1],[3.2,0,1],[.3,-.5,-1],[.8,-.5,-1],[.4,0,1],[3.1,0,1.]])
        self.x1 = self.x0.copy(); self.x1[:2,2] = -1
        self.i = np.array([0,2,4]); self.j = np.array([1,3,5])
        self.g = np.array([[.25,-.75,1.5],[1.,2.,3.],[4.,5.,6.]])
        # IDs deliberately differ from index ordering.
        self.ids = np.array([31,17,91,73,121,113])
        self.incoming = [[17,31]]

    def graph(self, incoming=None, **replace):
        args = dict(x0=self.x0, x1=self.x1, raw_i=self.i, raw_j=self.j,
            raw_g=self.g, pose0=self.pose, pose1=self.pose, p=self.p,
            particle_ids=self.ids, incoming_pid_pairs=self.incoming if incoming is None else incoming)
        args.update(replace)
        return event_graph.restrict_activation_graph(**args)

    def test_only_requested_endpoint_activation_is_delayed(self):
        graph = self.graph()
        np.testing.assert_array_equal(graph['endpoint_graph']['kept_raw_edge_indices'], [0,1])
        np.testing.assert_array_equal(graph['kept_raw_edge_indices'], [1])
        np.testing.assert_array_equal(graph['blocked_raw_edge_indices'], [0,2])
        np.testing.assert_array_equal(graph['i'], [2])
        np.testing.assert_array_equal(graph['j'], [3])
        np.testing.assert_array_equal(graph['g'], self.g[[1]])
        self.assertEqual(graph['activation_delay']['delayed_pid_pairs'], [[17,31]])
        self.assertEqual(graph['diagnostics']['activation_delayed_edge_count'], 1)
        self.assertEqual(graph['diagnostics']['endpoint_blocked_edge_count'], 1)

    def test_no_mask_preserves_all_original_physical_arrays_and_fingerprints(self):
        reference = filter_visible_support(self.x1,self.i,self.j,self.g,self.pose,self.p,self.ids)
        actual = self.graph([])
        for key, value in reference.items():
            if isinstance(value, np.ndarray):
                self.assertEqual(actual[key].dtype, value.dtype)
                self.assertEqual(actual[key].tobytes(), value.tobytes(), key)
        for key, value in reference['fingerprints'].items():
            self.assertEqual(actual['fingerprints'][key], value, key)
        self.assertEqual(actual['activation_delay']['delayed_edge_count'], 0)
        self.assertEqual(actual['diagnostics']['classification'], event_graph.CLASSIFICATION)

    def test_next_step_empty_mask_fully_restores_incoming_support(self):
        delayed = self.graph()
        restored = self.graph([], x0=self.x1)
        np.testing.assert_array_equal(restored['i'], delayed['endpoint_graph']['i'])
        np.testing.assert_array_equal(restored['g'], delayed['endpoint_graph']['g'])
        self.assertIn(0, restored['kept_raw_edge_indices'])
        self.assertFalse(delayed['activation_delay']['persistent_disabled_edges'])

    def test_still_blocked_incoming_does_not_change_endpoint_occlusion(self):
        graph = self.graph([[17,31],[113,121]])
        self.assertEqual(graph['activation_delay']['delayed_edge_count'], 1)
        np.testing.assert_array_equal(graph['strict_blocked_raw_edge_indices'], [2])
        self.assertEqual(graph['diagnostics']['strict_blocked_edge_count'], 1)
        self.assertEqual(graph['diagnostics']['blocked_edge_count'], 2)

    def test_unknown_self_duplicate_noncanonical_and_noninteger_masks_rejected(self):
        for mask in ([[17,999]], [[17,17]], [[31,17]], [[17,31],[17,31]],
                [[False,31]], [[17.,31]], [[17]], None, np.array([[17,31]])):
            with self.subTest(mask=mask), self.assertRaises(ValueError):
                self.graph(incoming_pid_pairs=mask)

    def test_start_visible_pair_cannot_be_removed(self):
        with self.assertRaisesRegex(ValueError, 'originally blocked'):
            self.graph([[73,91]])

    def test_known_ids_not_in_raw_support_cannot_make_stale_mask(self):
        with self.assertRaisesRegex(ValueError, 'absent from supplied raw support'):
            self.graph([[17,73]])

    def test_raw_input_geometry_and_gradients_remain_unchanged(self):
        values = [self.x0,self.x1,self.i,self.j,self.g,self.ids,self.pose['origin'],self.pose['rotation']]
        before = [value.copy() for value in values]
        p, incoming = copy.deepcopy(self.p), copy.deepcopy(self.incoming)
        graph = self.graph()
        for actual, expected in zip(values,before):
            np.testing.assert_array_equal(actual, expected)
        self.assertEqual(self.p,p); self.assertEqual(self.incoming,incoming)
        graph['i'][0] = 0; graph['g'][0] = 99
        np.testing.assert_array_equal(self.g,before[4])
        np.testing.assert_array_equal(graph['endpoint_graph']['g'],self.g[[0,1]])

    def test_fingerprints_describe_selected_graph_and_keep_original_geometry(self):
        graph = self.graph(); fingerprints = graph['fingerprints']
        self.assertEqual(fingerprints['visible_topology_sha256'], _topology_hash([[91,73]]))
        self.assertEqual(fingerprints['blocked_topology_sha256'], _topology_hash([[31,17],[121,113]]))
        self.assertEqual(fingerprints['visible_gradient_sha256'], _hash_array(self.g[[1]]))
        self.assertEqual(fingerprints['raw_gradient_sha256'], _hash_array(self.g))
        self.assertEqual(fingerprints['geometry_policy_sha256'],
            graph['endpoint_graph']['fingerprints']['geometry_policy_sha256'])
        self.assertNotEqual(fingerprints['activation_policy_sha256'], self.graph([])['fingerprints']['activation_policy_sha256'])

    def test_original_start_and_endpoint_numerical_domain_checks_both_apply(self):
        ambiguous = self.x0.copy()
        ambiguous[:2] = [[-2.,2.9,1.],[2.,2.9,1.]]
        for key in ('x0','x1'):
            with self.subTest(endpoint=key), self.assertRaises(VisibleSupportNumericalDomainError):
                self.graph([], **{key: ambiguous})

    def test_invalid_raw_gradients_shapes_and_particle_ids_rejected(self):
        for replacement in (dict(raw_g=self.g[:-1]), dict(raw_i=np.array([1,2,4])),
                dict(particle_ids=np.array([31,31,91,73,121,113])), dict(x0=self.x0[:-1])):
            with self.subTest(replacement=list(replacement)), self.assertRaises(ValueError):
                self.graph(**replacement)

    def test_old_save_graph_fields_and_partition_are_complete_and_metadata_serializable(self):
        graph = self.graph()
        kept, blocked = graph['kept_raw_edge_indices'], graph['blocked_raw_edge_indices']
        np.testing.assert_array_equal(np.sort(np.r_[kept,blocked]),np.arange(len(self.i)))
        self.assertEqual([r['raw_edge_index'] for r in graph['diagnostics']['blocked_edges']], blocked.tolist())
        metadata = json.loads(json.dumps(dict(diagnostics=graph['diagnostics'],fingerprints=graph['fingerprints'])))
        self.assertEqual(metadata['diagnostics']['activation_delay']['incoming_pid_pairs'],[[17,31]])
        delayed = next(row for row in metadata['diagnostics']['blocked_edges'] if row['raw_edge_index']==0)
        self.assertEqual(delayed['strict_total_length_m'],0.)
        self.assertGreater(delayed['start_filter_evidence']['strict_total_length_m'],0.)
        self.assertFalse(metadata['diagnostics']['activation_delay']['path_certified'])


class SavedEventGraphTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        stage = ROOT/'outputs/gpu_acceleration_20261007/campaign_10h_20261010_v1/attempt_0044'
        raw = (stage/'manifest.json').read_bytes()
        if hashlib.sha256(raw).hexdigest() != '0db409319a1ff9de2f8d3e7fc3a9e8615e2d9fcedef72f181132da65071ecd6f':
            raise AssertionError('Frozen failed-event manifest changed')
        manifest = json.loads(raw)
        for path in (stage/'coarse/visible_graph.npz',stage/'coarse/motion.json'):
            if hashlib.sha256(path.read_bytes()).hexdigest() != manifest['artifacts_sha256'][path.relative_to(ROOT).as_posix()]:
                raise AssertionError('Frozen failed-event evidence changed')
        with np.load(stage/'coarse/visible_graph.npz',allow_pickle=False) as saved:
            cls.saved = {name: saved[name] for name in saved.files}
        motion = json.loads((stage/'coarse/motion.json').read_text())
        cls.pose0 = dict(origin=motion['origin0_m'],rotation=motion['rotation0'])
        cls.pose1 = dict(origin=motion['origin1_m'],rotation=motion['rotation1'])
        config = json.loads((ROOT/'configs/deposition_fill_then_track_20261007.json').read_text())
        cls.p = json.loads((ROOT/config['base_parameters']).read_text())['parameters']
        cls.p.update({key:row['value'] for key,row in config['parameter_overrides'].items()})

    def graph(self, incoming, start=None):
        d = self.saved
        return event_graph.restrict_activation_graph(d['position0'] if start is None else start,
            d['position1'],d['raw_i'],d['raw_j'],d['raw_g'],self.pose0 if start is None else self.pose1,
            self.pose1,self.p,d['particle_id'],incoming)

    def test_real_failed44_removes_exactly_declared_incoming_pair_during_new_solve(self):
        graph = self.graph([[1334,1371]])
        self.assertEqual(graph['activation_delay']['delayed_pid_pairs'],[[1334,1371]])
        self.assertEqual(graph['activation_delay']['delayed_edge_count'],1)
        self.assertEqual(len(graph['i']),len(self.saved['i'])-1)
        np.testing.assert_array_equal(graph['endpoint_graph']['i'],self.saved['i'])
        np.testing.assert_array_equal(graph['endpoint_graph']['j'],self.saved['j'])
        np.testing.assert_array_equal(graph['endpoint_graph']['g'],self.saved['g'])
        np.testing.assert_array_equal(graph['g'],self.saved['raw_g'][graph['kept_raw_edge_indices']])
        self.assertFalse(graph['activation_delay']['event_time_resolved'])

    def test_real_failed44_no_mask_and_postevent_restore_complete_original_endpoint_graph(self):
        for graph in (self.graph([]),self.graph([],start=self.saved['position1'])):
            for name in ('i','j','g','kept_raw_edge_indices','blocked_raw_edge_indices'):
                np.testing.assert_array_equal(graph[name],self.saved[name])
            self.assertEqual(graph['activation_delay']['delayed_edge_count'],0)


if __name__ == '__main__':
    unittest.main()
