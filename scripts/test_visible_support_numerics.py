"""Limited-domain and numerical ambiguity rejection components."""
import unittest
import numpy as np
from visible_support_numerics import (check_partition_numerical_domain,
    require_visible_support_numerical_domain,VisibleSupportNumericalDomainError,PASS,UNRESOLVED)
from visible_support_path import audit_visible_support_path,UNRESOLVED as PATH_UNRESOLVED


class VisibleSupportNumericalDomainTests(unittest.TestCase):
    def setUp(self):
        self.p=dict(tool_radius=2.9,bore_width=2.,layer=0.,tool_top=2.)
        self.pose=dict(origin=np.zeros(3),rotation=np.eye(3))
        self.i=np.array([0]);self.j=np.array([1]);self.ids=np.array([1,2])

    def test_long_edge_cancellation_counterexample_rejected(self):
        # The true wall intersection is about 1.51 units long, but the old
        # unbounded discriminant calculation can lose it to cancellation.
        x=np.array([[-1e8,2.8,1.],[1e8,2.8,1.]])
        with self.assertRaises(VisibleSupportNumericalDomainError) as caught:
            require_visible_support_numerical_domain(x,self.i,self.j,self.pose,self.p,self.ids)
        self.assertIn('COORDINATE_EDGE_OR_RADIUS_SCALE_RATIO_EXCEEDS_EIGHT',caught.exception.report['reasons'])
        path=audit_visible_support_path(x,x,self.i,self.j,self.pose,self.pose,self.p,self.ids)
        self.assertEqual(path['status'],PATH_UNRESOLVED)
        self.assertEqual(path['failure_kind'],'PARTITION_NUMERICAL_DOMAIN_UNRESOLVED')

    def test_tangent_discriminant_uncertainty_rejected(self):
        start=np.array([[-2.,2.9,1.]])
        r=check_partition_numerical_domain(start,np.array([[4.,0.,0.]]),2.9,1.)
        self.assertEqual(r['status'],UNRESOLVED)
        self.assertEqual(r['uncertain_edge_indices'],[0])

    def test_resolved_inside_and_outside_cylinder_queries(self):
        for y in (2.8,3.):
            r=check_partition_numerical_domain(np.array([[-2.,y,1.]]),np.array([[4.,0.,0.]]),2.9,1.)
            self.assertEqual(r['status'],PASS)

    def test_units_translation_rotation_and_purity(self):
        for scale in (1e-4,1.,1e3):
            body=np.array([[.2,0,.5],[.8,0,1.5]])*scale
            angle=.61;c,s=np.cos(angle),np.sin(angle);R=np.array([[c,-s,0.],[s,c,0.],[0,0,1.]])
            origin=np.array([.2,-.1,.3])*scale;x=body@R.T+origin;saved=x.copy()
            p=dict(tool_radius=2.9*scale,bore_width=2*scale,layer=0.,tool_top=2*scale)
            r=require_visible_support_numerical_domain(x,self.i,self.j,dict(origin=origin,rotation=R),p,self.ids)
            self.assertEqual(r['status'],PASS);np.testing.assert_array_equal(x,saved)

    def test_norm_limit_and_axis_parallel_boundary_ambiguity(self):
        r=check_partition_numerical_domain(np.array([[6.,6.,0.]]),np.array([[0.,0.,1.]]),2.9,1.)
        self.assertEqual(r['status'],UNRESOLVED)
        r=check_partition_numerical_domain(np.array([[2.9,0.,.5]]),np.array([[0.,0.,1.]]),2.9,1.)
        self.assertEqual(r['status'],UNRESOLVED)

    def test_enclosure_query_itself_is_domain_checked(self):
        x0=np.array([[7.,0.,-1.],[7.,0.,-2.]])
        x1=np.array([[-7.,0.,-1.],[-7.,0.,-2.]])
        r=audit_visible_support_path(x0,x1,self.i,self.j,self.pose,self.pose,self.p,self.ids)
        self.assertEqual(r['status'],PATH_UNRESOLVED)
        self.assertEqual(r['numerical_domain_failures'][0]['kind'],'enclosure')


if __name__=='__main__':unittest.main(verbosity=2)
