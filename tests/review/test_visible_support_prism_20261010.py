"""Independent near-lip geometry, rotation, failure, and saved-path evidence."""
import ast
import hashlib
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import visible_support_prism_20261010 as proof
import visible_support_remote_domain_20261007 as frozen


def rotation(angle):
    c,s=np.cos(angle),np.sin(angle)
    return np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])


class PrismGeometryTests(unittest.TestCase):
    def setUp(self):
        self.p=dict(tool_radius=3.,bore_width=2.,layer=0.,tool_top=2.)
        self.pose=dict(origin=np.zeros(3),rotation=np.eye(3))
        self.i,self.j,self.ids=np.array([0]),np.array([1]),np.array([11,13])
        self.corner=np.array([[.5,0.,.5],[1.5,0.,-.50000002]])

    def audit(self,x0,x1=None,pose0=None,pose1=None,**kw):
        return proof.audit_visible_support_path(np.asarray(x0,float),
            np.asarray(x0 if x1 is None else x1,float),self.i,self.j,
            self.pose if pose0 is None else pose0,self.pose if pose1 is None else pose1,
            self.p,self.ids,**kw)

    def test_near_corner_translation_has_exact_separating_line(self):
        # The z=0 crossing is just inside x=1. Tangential y motion makes
        # the old isotropic enclosure inconclusive, while x/z is constant.
        end=self.corner+[0.,.2,0.]
        old=frozen.audit_visible_support_path(self.corner,end,self.i,self.j,
            self.pose,self.pose,self.p,self.ids,maximum_depth=0)
        self.assertEqual(old['status'],proof.UNRESOLVED)
        report=self.audit(self.corner,end,maximum_depth=0)
        self.assertEqual(report['status'],proof.CERTIFIED)
        self.assertEqual(report['analytic_batch_queries'],5)
        self.assertEqual(report['prism_batch_queries'],1)
        self.assertEqual(report['prism_curvature_method'],'CONSTANT_ROTATION_AFFINE_BODY_PATH')
        witnesses=[w for r in report['certificates'] if r.get('method') for w in r['witnesses']]
        self.assertEqual(len(witnesses),1)
        self.assertEqual(len(witnesses[0]['wall_separations']),4)
        self.assertTrue(any(all(v!=0 for v in w['normal']) for w in witnesses[0]['wall_separations']))
        self.assertTrue(all(w['strict_gap_m']>0. for w in witnesses[0]['wall_separations']))

    def test_all_four_signed_walls_are_required(self):
        for axis in (0,1):
            for sign in (-1,1):
                x=np.array([[0.,0.,.8],[0.,0.,1.2]])
                x[:,axis]=sign*np.array([1.2,1.4])
                report=self.audit(x)
                self.assertEqual(report['status'],proof.FOUND)
                self.assertEqual(report['prism_batch_queries'],0)

    def test_touching_corner_has_no_strict_prism_separation(self):
        x=self.corner.copy();x[1,2]=-.5
        report=self.audit(x,x+[0.,.2,0.],maximum_depth=0)
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertFalse(report['found_edges'])

    def test_rotation_visible_arc_uses_curvature_and_strict_projection_bounds(self):
        x=np.array([[0.,0.,1.],[1.96,0.,-1.]])
        report=self.audit(x,pose0=dict(origin=np.zeros(3),rotation=rotation(-.3)),
            pose1=dict(origin=np.zeros(3),rotation=rotation(.1)))
        self.assertEqual(report['status'],proof.CERTIFIED)
        prisms=[r for r in report['certificates'] if r.get('method')]
        self.assertTrue(prisms)
        for certificate in prisms:
            for witness in certificate['witnesses']:
                self.assertGreater(witness['second_derivative_bound_m'],0.)
                for t in np.linspace(certificate['start_fraction'],certificate['end_fraction'],17):
                    angle=-.3+.4*t
                    for s in np.linspace(0.,1.,9):
                        point=(1.-s)*x[0]+s*x[1]
                        body=np.array([point[0]*np.cos(angle)+point[1]*np.sin(angle),
                            -point[0]*np.sin(angle)+point[1]*np.cos(angle),point[2]])
                        for wall in witness['wall_separations']:
                            axis='xy'.index(wall['horizontal_axis'])
                            projected=np.dot([wall['horizontal_sign']*body[axis],body[2]],wall['normal'])
                            lo,hi=wall['hull_projection_interval_m']
                            self.assertGreaterEqual(projected,lo)
                            self.assertLessEqual(projected,hi)

    def test_rotation_clear_endpoints_and_midpoint_can_still_cross_at_three_quarters(self):
        # z=0 has body x=1.003*cos(theta). It is inside the bore at
        # theta=-.3,-.1,+.1 (endpoints/midpoint), but outside at theta=0.
        x=np.array([[0.,0.,1.],[2.006,0.,-1.]])
        report=self.audit(x,pose0=dict(origin=np.zeros(3),rotation=rotation(-.3)),
            pose1=dict(origin=np.zeros(3),rotation=rotation(.1)))
        self.assertEqual(report['status'],proof.FOUND)
        self.assertTrue(any(row['fraction']==.75 for row in report['found_edges']))

    def test_budget_counts_prism_and_fails_before_unfunded_proof(self):
        for budget in (2,4):
            report=self.audit(self.corner,self.corner+[0.,.2,0.],maximum_queries=budget)
            self.assertEqual(report['status'],proof.UNRESOLVED)
            self.assertEqual(report['failure_kind'],'MAXIMUM_ANALYTIC_BATCH_QUERIES')
            self.assertEqual(report['analytic_batch_queries'],budget)
            self.assertEqual(report['prism_batch_queries'],0)

    def test_hidden_quarter_collapse_cannot_become_hull_certificate(self):
        x0=np.array([[.2,0.,.5],[.4,0.,.5]])
        x1=np.array([[.6,0.,.5],[0.,0.,.5]])
        report=self.audit(x0,x1)
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertEqual(report['failure_kind'],'DEGENERATE_INTERIOR_SUPPORT_EDGE')
        self.assertAlmostEqual(report['degenerate_edges'][0]['fraction'],.25)

    def test_original_ambiguous_cylinder_domain_still_rejects(self):
        report=self.audit([[-2.,3.,1.],[2.,3.,1.]])
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertEqual(report['failure_kind'],'PARTITION_NUMERICAL_DOMAIN_UNRESOLVED')
        self.assertEqual(report['prism_batch_queries'],0)

    def test_missing_nonzero_rotation_curvature_does_not_become_zero(self):
        report=self.audit(self.corner,self.corner+[0.,-.2,0.],
            pose1=dict(origin=np.zeros(3),rotation=rotation(5e-7)),maximum_depth=0)
        self.assertEqual(report['status'],proof.UNRESOLVED)
        self.assertEqual(report['prism_curvature_method'],'SAT_UNAVAILABLE_OUTSIDE_INHERITED_CURVATURE_DOMAIN')
        self.assertEqual(report['prism_batch_queries'],0)

    def test_invalid_input_and_query_limits_reject(self):
        for kwargs in (dict(maximum_queries=257),dict(maximum_depth=33),dict(maximum_queries=True)):
            with self.assertRaises(ValueError):self.audit(self.corner,**kwargs)
        with self.assertRaises(ValueError):self.audit(np.full((2,3),float('nan')))

    def test_inputs_are_unchanged_and_ids_do_not_drive_proof(self):
        first,last=self.corner.copy(),self.corner+[0.,.2,0.]
        copies=[a.copy() for a in (first,last,self.i,self.j,self.ids)]
        report=self.audit(first,last)
        for before,after in zip(copies,(first,last,self.i,self.j,self.ids)):
            np.testing.assert_array_equal(before,after)
        self.ids=np.array([2,1]);renamed=self.audit(first,last)
        self.assertEqual(report['status'],renamed['status'])
        self.assertEqual(report['certificates'],renamed['certificates'])
        self.assertNotEqual(report['input_sha256']['retained_topology'],renamed['input_sha256']['retained_topology'])
        json.dumps(report)

    def test_ordinary_original_certificate_and_empty_graph_still_work(self):
        report=self.audit([[.2,0.,.5],[.8,0.,1.5]])
        self.assertEqual(report['status'],proof.CERTIFIED)
        self.assertEqual(report['analytic_batch_queries'],4)
        self.assertEqual(report['prism_batch_queries'],0)
        empty=np.array([],int)
        report=proof.audit_visible_support_path(self.corner,self.corner,empty,empty,
            self.pose,self.pose,self.p,self.ids)
        self.assertEqual(report['status'],proof.CERTIFIED)
        self.assertEqual(report['edge_count'],0)


class FrozenGuardPreservationTests(unittest.TestCase):
    def test_original_strict_domain_enclosure_and_binary_subdivision_are_preserved(self):
        old_path=ROOT/'scripts/visible_support_remote_domain_20261007.py'
        self.assertEqual(hashlib.sha256(old_path.read_bytes()).hexdigest(),proof.FROZEN_PARENT_SHA256)
        def function(path):
            return next(n for n in ast.parse(path.read_text()).body
                if isinstance(n,ast.FunctionDef) and n.name=='audit_visible_support_path')
        old=function(old_path);new=function(ROOT/'scripts/visible_support_prism_20261010.py')
        for name in ('query','local_at','strict_at','unresolved'):
            a=next(n for n in old.body if isinstance(n,ast.FunctionDef) and n.name==name)
            b=next(n for n in new.body if isinstance(n,ast.FunctionDef) and n.name==name)
            self.assertEqual(ast.dump(a),ast.dump(b),name)
        old_loop=next(n for n in ast.walk(old) if isinstance(n,ast.While))
        new_loop=next(n for n in ast.walk(new) if isinstance(n,ast.While))
        class RemoveNewProof(ast.NodeTransformer):
            def visit_Assign(self,node):
                if any(isinstance(t,ast.Name) and t.id=='prisms' for t in node.targets):return None
                return self.generic_visit(node)
            def visit_If(self,node):
                if isinstance(node.test,ast.Name) and node.test.id=='prisms':return None
                return self.generic_visit(node)
        self.assertEqual(ast.dump(old_loop),ast.dump(RemoveNewProof().visit(new_loop)))


class SavedPathTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cfg=json.loads((ROOT/'configs/deposition_fill_then_track_20261007.json').read_text())
        cls.p=json.loads((ROOT/cfg['base_parameters']).read_text())['parameters']
        cls.p.update({k:v['value'] for k,v in cfg['parameter_overrides'].items()})

    def saved(self,relative,expected):
        stage=ROOT/'outputs/gpu_acceleration_20261007'/relative
        for name,digest in expected.items():
            self.assertEqual(hashlib.sha256((stage/name).read_bytes()).hexdigest(),digest)
        with np.load(stage/'visible_graph.npz',allow_pickle=False) as z:
            x0,x1,ids,i,j=z['position0'],z['position1'],z['particle_id'],z['i'],z['j']
        m=json.loads((stage/'motion.json').read_text())
        report=proof.audit_visible_support_path(x0,x1,i,j,
            dict(origin=m['origin0_m'],rotation=m['rotation0']),
            dict(origin=m['origin1_m'],rotation=m['rotation1']),self.p,ids)
        previous=json.loads((stage/'visible_path_audit.json').read_text())
        self.assertEqual(report['input_sha256'],previous['input_sha256'])
        return report,previous

    def test_actual_post_event_full_graph_certified_with_original_budget(self):
        report,previous=self.saved('event_box_campaign_20261010_v1/sensitivity_one/coarse/after_event/proposal',{
            'visible_graph.npz':'fbf9310961b0ccf42361c9c7484d3d6ec8522c87b8ae4c71547d804e697ee1fc',
            'motion.json':'56bc854d21e20e92b49abdbb841f850974bc6deb22d42ce6ce922ac52043819d',
            'visible_path_audit.json':'fa34669f5e95d1087e6bb094fbab723641d784bf4d9c437723eae594844a12ae'})
        self.assertEqual(previous['failure_kind'],'MAXIMUM_ANALYTIC_BATCH_QUERIES')
        self.assertEqual(report['status'],proof.CERTIFIED)
        self.assertEqual(report['edge_count'],33720)
        self.assertLessEqual(report['analytic_batch_queries'],256)
        self.assertLessEqual(report['deepest_subdivision'],32)
        self.assertGreater(report['prism_batch_queries'],0)
        self.assertEqual(report['analytic_batch_queries'],report['strict_geometry_batch_queries']+
            report['enclosure_batch_queries']+report['prism_batch_queries'])
        self.assertFalse(report['found_edges'] or report['degenerate_edges'] or
            report['numerical_domain_failures'] or report['unresolved_intervals'])
        # Every retained edge must have a gap-free interval cover.
        coverage=[[] for _ in range(report['edge_count'])]
        for row in report['certificates']:
            for edge in row['edge_indices']:
                coverage[edge].append((row['start_fraction'],row['end_fraction']))
        for intervals in coverage:
            intervals.sort()
            self.assertEqual(intervals[0][0],0.);self.assertEqual(intervals[-1][1],1.)
            self.assertTrue(all(a[1]==b[0] for a,b in zip(intervals,intervals[1:])))
        self.assertEqual(report['physical_solver_evaluations'],0)
        self.assertEqual(report['edge_deletions'],0)

    def test_original_real_activation_remains_strictly_rejected_at_start(self):
        stage=ROOT/'outputs/gpu_acceleration_20261007/campaign_10h_20261010_v1/attempt_0044'
        blob=(stage/'manifest.json').read_bytes()
        self.assertEqual(hashlib.sha256(blob).hexdigest(),'0db409319a1ff9de2f8d3e7fc3a9e8615e2d9fcedef72f181132da65071ecd6f')
        manifest=json.loads(blob)
        expected={name:manifest['artifacts_sha256'][(stage/'coarse'/name).relative_to(ROOT).as_posix()]
            for name in ('visible_graph.npz','motion.json','visible_path_audit.json')}
        report,_=self.saved('campaign_10h_20261010_v1/attempt_0044/coarse',expected)
        self.assertEqual(report['status'],proof.FOUND)
        self.assertEqual(report['analytic_batch_queries'],2)
        self.assertEqual(report['prism_batch_queries'],0)
        self.assertTrue(any(row['fraction']==0. and {row['particle_id_i'],row['particle_id_j']}=={1334,1371}
            for row in report['found_edges']))


if __name__=='__main__':
    unittest.main()
