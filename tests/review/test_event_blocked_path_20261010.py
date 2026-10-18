"""Independent geometric positives/negatives and frozen saved-event evidence."""
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import afsd_event_blocked_path_20261010 as proof


class BlockedPrefixTests(unittest.TestCase):
    def setUp(self):
        self.p=dict(tool_radius=2.9,bore_width=2.,layer=0.,tool_top=2.)
        self.pose=dict(origin=np.zeros(3),rotation=np.eye(3))
        self.x0=np.array([[.5,0.,1.],[3.2,0.,1.]])
        self.x1=self.x0.copy();self.x1[:,2]=-1.
        self.i=np.array([0]);self.j=np.array([1]);self.ids=np.array([11,13])

    def audit(self,end=.49,**replace):
        args=dict(x0=self.x0,x1=self.x1,i=self.i,j=self.j,pose0=self.pose,pose1=self.pose,
            p=self.p,particle_ids=self.ids,end_fraction=end)
        args.update(replace)
        return proof.audit_blocked_prefix(**args)

    def test_exact_translation_before_known_exit_has_continuous_certificate(self):
        report=self.audit()
        self.assertEqual(report['status'],proof.CERTIFIED)
        self.assertTrue(report['certificates'])
        self.assertEqual(report['failure_kind'],None)
        self.assertEqual(report['unresolved_intervals'],[])
        self.assertLessEqual(report['analytic_batch_queries'],256)
        intervals=sorted((row['start_fraction'],row['end_fraction']) for row in report['certificates'])
        self.assertEqual(intervals[0][0],0.);self.assertEqual(intervals[-1][1],.49)
        for left,right in zip(intervals,intervals[1:]):
            self.assertEqual(left[1],right[0])

    def test_fixed_segment_witnesses_stay_inside_original_solid(self):
        report=self.audit()
        for row in report['certificates']:
            for witness in row['witnesses']:
                s=witness['witness_segment_parameter']
                for t in np.linspace(row['start_fraction'],row['end_fraction'],17):
                    x=self.x0+t*(self.x1-self.x0)
                    point=(1.-s)*x[0]+s*x[1]
                    self.assertLess(np.linalg.norm(point[:2]),self.p['tool_radius'])
                    self.assertTrue(abs(point[0])>1. or abs(point[1])>1.)
                    self.assertGreater(point[2],0.);self.assertLess(point[2],2.)

    def test_touching_and_beyond_known_exit_are_never_certified(self):
        for end in (.5,.51,1.):
            report=self.audit(end)
            self.assertEqual(report['status'],proof.UNRESOLVED)
            self.assertTrue(report['visible_samples'])

    def test_stationary_strictly_blocked_segment_certifies(self):
        report=self.audit(1.,x1=self.x0)
        self.assertEqual(report['status'],proof.CERTIFIED)
        self.assertEqual(report['physical_solver_evaluations'],0)
        self.assertEqual(report['edge_deletions'],0)

    def test_blocked_endpoints_with_hidden_open_middle_are_rejected(self):
        first=np.array([[-1.2,-.05,.8],[-1.1,.05,1.2]])
        last=first.copy();last[:,0]+=2.3
        report=self.audit(1.,x0=first,x1=last)
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertTrue(any(row['fraction']==.5 for row in report['visible_samples']))

    def test_rotating_tool_opens_and_reblocks_between_blocked_endpoints(self):
        x=np.array([[1.05,-.05,.8],[1.15,.05,1.2]])
        pose1=dict(origin=np.zeros(3),rotation=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]]))
        report=self.audit(1.,x0=x,x1=x,pose1=pose1)
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertTrue(report['visible_samples'])

    def test_curvature_box_bounds_analytic_circular_arc_beyond_endpoint_extrema(self):
        x=np.array([[1.3,-.1,.8],[1.5,.1,1.2]])
        def rotation(angle):
            c,s=np.cos(angle),np.sin(angle)
            return np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])
        report=self.audit(1.,x0=x,x1=x,
            pose0=dict(origin=np.zeros(3),rotation=rotation(-.2)),
            pose1=dict(origin=np.zeros(3),rotation=rotation(.2)))
        self.assertEqual(report['status'],proof.CERTIFIED)
        boxes=[r for r in report['certificates'] if r['method']=='FIXED_SEGMENT_PARAMETER_CURVATURE_BOX']
        self.assertTrue(boxes)
        beyond_unpadded_endpoints=False
        for row in boxes:
            for witness in row['witnesses']:
                s=witness['witness_segment_parameter'];world=(1-s)*x[0]+s*x[1]
                lower=np.array(witness['body_box_lower_m']);upper=np.array(witness['body_box_upper_m'])
                endpoints=np.array(witness['endpoint_body_positions_m'])
                for fraction in np.linspace(row['start_fraction'],row['end_fraction'],33):
                    angle=-.2+.4*fraction
                    actual=np.array([world[0]*np.cos(angle)+world[1]*np.sin(angle),
                        -world[0]*np.sin(angle)+world[1]*np.cos(angle),world[2]])
                    self.assertTrue(np.all(actual>=lower) and np.all(actual<=upper))
                    beyond_unpadded_endpoints |= bool(np.any(actual>endpoints.max(axis=0)+1e-10))
        self.assertTrue(beyond_unpadded_endpoints)
        self.assertEqual(report['analytic_batch_queries'],report['partition_queries']+report['fixed_witness_box_queries'])

    def test_samples_without_enclosure_budget_cannot_become_pass(self):
        report=self.audit(maximum_queries=2)
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertEqual(report['failure_kind'],'MAXIMUM_ANALYTIC_BATCH_QUERIES')
        self.assertEqual(report['certificates'],[])
        self.assertEqual(len(report['strict_samples']),2)

    def test_invalid_bounds_nonfinite_and_degenerate_input_rejected(self):
        for kwargs in (dict(end_fraction=0.),dict(end_fraction=float('nan')),
                dict(maximum_queries=257),dict(maximum_depth=33),dict(maximum_depth=True),
                dict(x1=np.full((2,3),float('inf'))),dict(x0=np.zeros((2,3)))):
            with self.subTest(kwargs=list(kwargs)),self.assertRaises((ValueError,ArithmeticError)):
                self.audit(**kwargs)

    def test_interior_collapsed_segment_is_unresolved(self):
        first=np.array([[1.2,0.,.8],[1.5,0.,1.2]])
        last=first[::-1].copy()
        report=self.audit(1.,x0=first,x1=last)
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertEqual(report['failure_kind'],'DEGENERATE_INCOMING_EDGE')

    def test_hidden_quarter_interval_collapse_cannot_pass_inside_point_box(self):
        # The midpoint of the edge stays strictly inside the tool at all
        # times, but edge length is zero at f=1/4 rather than a sampled midpoint.
        first=np.array([[1.3,0.,1.],[1.5,0.,1.]])
        last=np.array([[1.7,0.,1.],[1.1,0.,1.]])
        report=self.audit(1.,x0=first,x1=last)
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertEqual(report['failure_kind'],'DEGENERATE_INCOMING_EDGE')
        self.assertAlmostEqual(report['degenerate_edges'][0]['fraction'],.25)

    def test_original_numerical_domain_ambiguity_is_not_relaxed(self):
        tangent=np.array([[-2.,2.9,1.],[2.,2.9,1.]])
        report=self.audit(1.,x0=tangent,x1=tangent)
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertEqual(report['failure_kind'],'PARTITION_NUMERICAL_DOMAIN_UNRESOLVED')
        self.assertTrue(report['numerical_domain_failures'])

    def test_inputs_unchanged_and_geometry_queries_use_actual_contraction(self):
        a,b=self.x0.copy(),self.x1.copy()
        report=self.audit()
        np.testing.assert_array_equal(self.x0,a);np.testing.assert_array_equal(self.x1,b)
        for row in report['query_log']:
            if row['kind']=='contracted':
                self.assertLess(row['queried_radius_m'],2.9)
                self.assertGreater(row['queried_half_bore_m'],1.)
                self.assertGreater(row['queried_zlo_m'],0.)
                self.assertLess(row['queried_zhi_m'],2.)
        self.assertNotEqual(report['input_sha256']['prefix_endpoint'],self.audit(.4)['input_sha256']['prefix_endpoint'])
        json.dumps(report)


class SavedBlockedPrefixTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        stage=ROOT/'outputs/gpu_acceleration_20261007/campaign_10h_20261010_v1/attempt_0044'
        blob=(stage/'manifest.json').read_bytes()
        if hashlib.sha256(blob).hexdigest()!='0db409319a1ff9de2f8d3e7fc3a9e8615e2d9fcedef72f181132da65071ecd6f':
            raise AssertionError('Frozen failed-event manifest changed')
        manifest=json.loads(blob)
        for f in (stage/'coarse/visible_graph.npz',stage/'coarse/motion.json'):
            if hashlib.sha256(f.read_bytes()).hexdigest()!=manifest['artifacts_sha256'][f.relative_to(ROOT).as_posix()]:
                raise AssertionError('Frozen failed-event evidence changed')
        with np.load(stage/'coarse/visible_graph.npz',allow_pickle=False) as z:
            cls.x0,cls.x1,cls.ids=z['position0'],z['position1'],z['particle_id']
        motion=json.loads((stage/'coarse/motion.json').read_text())
        cls.pose0=dict(origin=motion['origin0_m'],rotation=motion['rotation0'])
        cls.pose1=dict(origin=motion['origin1_m'],rotation=motion['rotation1'])
        cfg=json.loads((ROOT/'configs/deposition_fill_then_track_20261007.json').read_text())
        cls.p=json.loads((ROOT/cfg['base_parameters']).read_text())['parameters']
        cls.p.update({k:v['value'] for k,v in cfg['parameter_overrides'].items()})

    def audit(self,end,**kw):
        return proof.audit_blocked_prefix(self.x0,self.x1,np.array([1333]),np.array([1370]),
            self.pose0,self.pose1,self.p,self.ids,end,**kw)

    def test_real44_strictly_blocked_early_prefix_is_continuously_certified(self):
        report=self.audit(.01)
        self.assertEqual(report['status'],proof.CERTIFIED)
        self.assertTrue(report['certificates'])

    def test_real44_prefix_crossing_saved_event_is_unresolved(self):
        report=self.audit(.36)
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertTrue(report['visible_samples'])
        self.assertEqual(report['visible_samples'][0]['particle_id_i'],1334)
        self.assertEqual(report['visible_samples'][0]['particle_id_j'],1371)

    def test_real44_prefix_immediately_before_saved_event_has_conservative_box_proof(self):
        # Saved event fraction ~0.359310138: this ends only 2.38e-7 fraction
        # before it. No candidate motion or physical result is changed.
        report=self.audit(.3593099)
        self.assertEqual(report['status'],proof.CERTIFIED)
        self.assertLessEqual(report['analytic_batch_queries'],256)
        self.assertTrue(any(r['method']=='FIXED_SEGMENT_PARAMETER_CURVATURE_BOX' for r in report['certificates']))
        intervals=sorted((r['start_fraction'],r['end_fraction']) for r in report['certificates'])
        self.assertEqual(intervals[0][0],0.)
        self.assertEqual(intervals[-1][1],.3593099)
        self.assertTrue(all(a[1]==b[0] for a,b in zip(intervals,intervals[1:])))

    def test_real44_near_event_exhaustion_fails_closed(self):
        report=self.audit(.359,maximum_queries=4)
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertLessEqual(report['analytic_batch_queries'],4)


if __name__=='__main__':
    unittest.main()
