"""Adversarial components for the fixed-graph continuous visibility audit."""
import unittest
import numpy as np
from visible_support import filter_visible_support
from visible_support_path import audit_visible_support_path,CERTIFIED,FOUND,UNRESOLVED


def zrotation(angle):
    c,s=np.cos(angle),np.sin(angle)
    return np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])


class VisibleSupportPathTests(unittest.TestCase):
    def setUp(self):
        self.p=dict(tool_radius=3.,bore_width=2.,layer=0.,tool_top=2.)
        self.pose=dict(origin=np.zeros(3),rotation=np.eye(3))
        self.i=np.array([0]);self.j=np.array([1]);self.ids=np.array([91,13])

    def audit(self,x0,x1=None,pose1=None,**kwargs):
        return audit_visible_support_path(np.asarray(x0,float),np.asarray(x0 if x1 is None else x1,float),
            self.i,self.j,self.pose,pose1 or self.pose,self.p,self.ids,**kwargs)

    def test_stationary_bore_visible_with_full_certificate(self):
        r=self.audit([[.2,0,.5],[.8,0,1.5]])
        self.assertEqual(r['status'],CERTIFIED)
        self.assertEqual(r['analytic_batch_queries'],4)
        self.assertEqual(r['edge_deletions'],0)

    def test_rotating_lip_has_hidden_intersection_with_valid_centres(self):
        # At both endpoint rotations the edge crosses z=0 inside the square.
        # At 45 degrees its crossing lies at body x=1.1, outside the bore.
        x=np.array([[1.1*np.sqrt(2),1.1*np.sqrt(2),-1.],[0.,0.,1.]])
        pose1=dict(origin=np.zeros(3),rotation=zrotation(np.pi/2))
        for pose in (self.pose,pose1):
            g=filter_visible_support(x,self.i,self.j,np.ones((1,3)),pose,self.p,self.ids)
            self.assertEqual(len(g['i']),1)
        r=self.audit(x,pose1=pose1)
        self.assertEqual(r['status'],FOUND)
        self.assertTrue(any(row['fraction']==.5 for row in r['found_edges']))
        self.assertGreater(r['found_edges'][0]['total_inside_length_m'],0.)

    def test_common_rigid_endpoints_still_follow_actual_world_chord(self):
        x0=np.array([[3.1,0,.5],[3.1,0,1.5]])
        R=zrotation(np.pi/2);x1=x0@R.T
        r=self.audit(x0,x1,pose1=dict(origin=np.zeros(3),rotation=R))
        self.assertEqual(r['status'],FOUND)
        self.assertTrue(any(row['fraction']==.5 for row in r['found_edges']))

    def test_exchanged_endpoints_collapsing_in_solid_fail_closed(self):
        # Both endpoint segments are outside R; at mid-time both particles
        # coincide inside the wall. A length-only intersection test returns 0.
        x0=np.array([[1.5,-4.,1.],[1.5,-3.,1.]])
        x1=np.array([[1.5,4.,1.],[1.5,3.,1.]])
        for x in (x0,x1):
            self.assertEqual(len(filter_visible_support(x,self.i,self.j,np.ones((1,3)),self.pose,self.p,self.ids)['i']),1)
        r=self.audit(x0,x1)
        self.assertEqual(r['status'],UNRESOLVED)
        self.assertEqual(r['failure_kind'],'DEGENERATE_INTERIOR_SUPPORT_EDGE')
        self.assertEqual(r['degenerate_edges'][0]['fraction'],.5)
        self.assertEqual(r['degenerate_edges'][0]['reference_edge_length_m'],0.)

    def test_budget_exhaustion_never_becomes_sampled_pass(self):
        r=self.audit([[.2,0,.5],[.8,0,1.5]],maximum_queries=2)
        self.assertEqual(r['status'],UNRESOLVED)
        self.assertEqual(r['analytic_batch_queries'],2)
        self.assertTrue(r['unresolved_intervals'])

    def test_common_translation_body_constant(self):
        x=np.array([[.2,0,.5],[.8,0,1.5]]);shift=np.array([.4,-.3,.6])
        r=self.audit(x,x+shift,pose1=dict(origin=shift,rotation=np.eye(3)))
        self.assertEqual(r['status'],CERTIFIED)

    def test_released_endpoint_edge_rejected_at_start(self):
        x0=np.array([[0,0,1.],[4,0,1.]]);x1=np.array([[0,0,1.],[.5,0,1.]])
        endpoint=filter_visible_support(x1,self.i,self.j,np.ones((1,3)),self.pose,self.p,self.ids)
        self.assertEqual(len(endpoint['i']),1)
        r=self.audit(x0,x1)
        self.assertEqual(r['status'],FOUND)
        self.assertTrue(any(row['fraction']==0. for row in r['found_edges']))
        self.assertEqual(r['edge_deletions'],0)

    def test_depth_limit_and_subdivision_are_distinct(self):
        x0=np.array([[.99,0,.5],[.99,0,1.5]]);x1=x0+np.array([0,.1,0])
        unresolved=self.audit(x0,x1,maximum_depth=0)
        self.assertEqual(unresolved['status'],UNRESOLVED)
        resolved=self.audit(x0,x1)
        self.assertEqual(resolved['status'],CERTIFIED)
        self.assertGreater(resolved['deepest_subdivision'],0)

    def test_exact_grazing_never_certified(self):
        r=self.audit([[1.,-.5,.5],[1.,.5,1.5]],maximum_depth=0)
        self.assertEqual(r['status'],UNRESOLVED)
        self.assertFalse(r['found_edges'])

    def test_empty_graph_has_explicit_scope(self):
        x=np.array([[.2,0,.5],[.8,0,1.5]])
        r=audit_visible_support_path(x,x,np.array([],int),np.array([],int),self.pose,self.pose,self.p,self.ids)
        self.assertEqual(r['status'],CERTIFIED)
        self.assertEqual(r['edge_count'],0)

    def test_inputs_pure_and_ids_only_observational(self):
        x0=np.array([[.99,0,.5],[.99,0,1.5]]);x1=x0+np.array([0,.1,0])
        original=[a.copy() for a in (x0,x1,self.i,self.j,self.ids)]
        a=self.audit(x0,x1)
        self.ids=np.array([4,2]);b=self.audit(x0,x1)
        for old,now in zip(original[:4],(x0,x1,self.i,self.j)):np.testing.assert_array_equal(old,now)
        self.assertEqual(a['status'],b['status'])
        self.assertEqual(a['analytic_batch_queries'],b['analytic_batch_queries'])
        self.assertNotEqual(a['input_sha256']['retained_topology'],b['input_sha256']['retained_topology'])

    def test_bounds_and_nonfinite_inputs_rejected(self):
        for options in ({'maximum_queries':257},{'maximum_queries':1},{'maximum_depth':33}):
            with self.assertRaises(ValueError):self.audit([[.2,0,.5],[.8,0,1.5]],**options)
        with self.assertRaises(ValueError):self.audit([[np.nan,0,.5],[.8,0,1.5]])


if __name__=='__main__':unittest.main(verbosity=2)
