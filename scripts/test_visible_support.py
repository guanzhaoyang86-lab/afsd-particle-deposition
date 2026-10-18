"""Independent geometry contracts for the new visible support policy."""
import unittest
import numpy as np
from visible_support import filter_visible_support
from audit_solid_occluded_neighbors import audit_solid_edge_visibility


class VisibleSupportGeometryTests(unittest.TestCase):
    def setUp(self):
        self.p=dict(tool_radius=3.,bore_width=2.,layer=0.,tool_top=2.)
        self.pose=dict(origin=np.zeros(3),rotation=np.eye(3))

    def edge(self,points,p=None,pose=None):
        return filter_visible_support(np.asarray(points,float),np.array([0]),np.array([1]),
            np.array([[.4,-.3,.7]]),pose or self.pose,p or self.p,np.array([91,13]))

    def test_open_bore_below_above_and_outside_are_visible(self):
        for points in ([[.2,0,.5],[.8,0,1.5]],[[2,0,-.1],[2.5,0,-.2]],
                       [[2,0,2.1],[2.5,0,2.2]],[[3.1,0,.2],[3.2,0,1.8]]):
            with self.subTest(points=points):
                self.assertEqual(self.edge(points)['diagnostics']['visible_edge_count'],1)

    def test_wall_and_aperture_lip_reject_valid_centres(self):
        for points in ([[0,0,1],[4,0,1]],[[1.5,0,-1],[.75,0,1]],
                       [[0,1.5,-1],[0,.75,1]],[[-4,0,1],[4,0,1]]):
            result=self.edge(points)
            self.assertEqual(result['diagnostics']['strict_blocked_edge_count'],1)
            self.assertEqual(len(result['i']),0)

    def test_subnanometre_strict_intersection_is_not_old_guard_tolerance(self):
        x=np.array([[1.-.25e-9,0,1.],[1.+.25e-9,0,1.]])
        old=audit_solid_edge_visibility(x,np.array([0]),np.array([1]),self.pose,self.p)
        self.assertEqual(old['occluded_edge_count'],0)
        new=self.edge(x)
        self.assertEqual(new['diagnostics']['strict_blocked_edge_count'],1)
        self.assertGreater(new['diagnostics']['maximum_strict_intersection_length_m'],0.)
        self.assertLess(new['diagnostics']['maximum_strict_intersection_length_m'],1e-9)

    def test_grazing_and_ulp_neighborhood_fail_conservatively(self):
        for points in ([[1,-.5,.5],[1,.5,1.5]], [[-2,3,1],[2,3,1]],
                       [[1.5,0,-1],[.5,0,1]],[[1.5,0,0],[2.5,0,0]]):
            result=self.edge(points)
            self.assertEqual(result['diagnostics']['blocked_edge_count'],1)
        for y in (np.nextafter(1.,0.),1.,np.nextafter(1.,2.)):
            self.assertEqual(self.edge([[0,y,.5],[.5,y,1.5]])['diagnostics']['blocked_edge_count'],1)
        self.assertEqual(self.edge([[0,1.-1e-6,.5],[.5,1.-1e-6,1.5]])['diagnostics']['visible_edge_count'],1)

    def test_thin_wall_two_sides_block(self):
        for thickness in (1e-2,1e-5,1e-7):
            p=dict(self.p,tool_radius=1.+thickness)
            result=self.edge([[1.-thickness,0,1.],[1.+2*thickness,0,1.]],p=p)
            self.assertEqual(result['diagnostics']['strict_blocked_edge_count'],1)
            self.assertAlmostEqual(result['diagnostics']['maximum_strict_intersection_length_m'],thickness,delta=2e-14)

    def test_rotated_translated_pose_covariance(self):
        angle=.71;c,s=np.cos(angle),np.sin(angle)
        rotation=np.array([[c,-s,0],[s,c,0],[0,0,1.]])
        origin=np.array([.8,-.3,.19]);body=np.array([[1.5,0,-1],[.75,0,1]])
        before=self.edge(body);after=self.edge(body@rotation.T+origin,pose=dict(origin=origin,rotation=rotation))
        self.assertEqual(before['diagnostics']['strict_blocked_edge_count'],after['diagnostics']['strict_blocked_edge_count'])
        self.assertAlmostEqual(before['diagnostics']['maximum_strict_intersection_length_m'],
            after['diagnostics']['maximum_strict_intersection_length_m'],delta=2e-14)

    def test_density_and_units_do_not_change_geometric_policy(self):
        for count in (3,11,37):
            for scale in (1e-4,1.,1e3):
                z=np.linspace(.1,1.9,count)*scale;x=np.column_stack((np.full(count,.7*scale),np.zeros(count),z))
                i=np.arange(count-1);j=i+1;g=np.tile([1.,2.,3.],(len(i),1))
                p=dict(tool_radius=3*scale,bore_width=2*scale,layer=0.,tool_top=2*scale)
                result=filter_visible_support(x,i,j,g,self.pose,p,np.arange(count)+10)
                self.assertEqual(result['diagnostics']['visible_edge_count'],count-1)

    def test_no_filter_exact_identity_raw_order_and_pure_inputs(self):
        x=np.array([[.1,0,.2],[.2,0,.4],[.3,0,.6],[.4,0,.8]])
        i=np.array([2,0,1],np.int32);j=np.array([3,2,3],np.int32)
        g=np.arange(9,dtype=np.float64).reshape(3,3);ids=np.array([8,2,7,4])
        original=[a.copy() for a in (x,i,j,g,ids)]
        result=filter_visible_support(x,i,j,g,self.pose,self.p,ids)
        for name,value in zip(('i','j','g'),(i,j,g)):
            self.assertEqual(result[name].tobytes(),value.tobytes());self.assertEqual(result[name].dtype,value.dtype)
        for before,after in zip(original,(x,i,j,g,ids)):np.testing.assert_array_equal(before,after)
        permutation=np.array([2,0,1])
        other=filter_visible_support(x,i[permutation],j[permutation],g[permutation],self.pose,self.p,ids)
        self.assertEqual(result['fingerprints']['raw_topology_sha256'],other['fingerprints']['raw_topology_sha256'])

    def test_persistent_ids_observe_and_never_change_mask(self):
        x=np.array([[0,0,1.],[4,0,1.],[.5,0,.8]])
        i=np.array([0,0]);j=np.array([1,2]);g=np.ones((2,3))
        a=filter_visible_support(x,i,j,g,self.pose,self.p,np.array([10,20,30]))
        b=filter_visible_support(x,i,j,g,self.pose,self.p,np.array([71,31,91]))
        np.testing.assert_array_equal(a['kept_raw_edge_indices'],b['kept_raw_edge_indices'])
        self.assertNotEqual(a['fingerprints']['raw_topology_sha256'],b['fingerprints']['raw_topology_sha256'])

    def test_empty_and_invalid_inputs_fail_closed(self):
        x=np.array([[0,0,.5],[.5,0,1.]])
        r=filter_visible_support(x,np.array([],int),np.array([],int),np.empty((0,3)),self.pose,self.p,np.array([1,2]))
        self.assertEqual(r['diagnostics']['visible_edge_count'],0)
        with self.assertRaises(ValueError):
            filter_visible_support(x,np.array([1]),np.array([0]),np.ones((1,3)),self.pose,self.p,np.array([1,2]))
        with self.assertRaises(ValueError):
            filter_visible_support(x,np.array([0,0]),np.array([1,1]),np.ones((2,3)),self.pose,self.p,np.array([1,2]))
        with self.assertRaises(ValueError):self.edge([[np.nan,0,1],[0,0,1]])
        with self.assertRaises(ValueError):self.edge([[0,0,1],[0,0,1]])
        with self.assertRaises(ArithmeticError):
            self.edge(np.array([[0,0,.5],[.5,0,1.]])+1e16,pose=dict(origin=np.full(3,1e16),rotation=np.eye(3)))


if __name__=='__main__':
    unittest.main(verbosity=2)
