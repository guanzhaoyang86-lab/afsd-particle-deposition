# STATIC FIXTURE ADAPTATION FROM test_visible_support_path.py; no monkeypatch.
"""Adversarial components for the fixed-graph continuous visibility audit."""
import unittest
import numpy as np
from visible_support import filter_visible_support
from visible_support_path_tight import audit_visible_support_path,CERTIFIED,FOUND,UNRESOLVED


def zrotation(angle):
    c,s=np.cos(angle),np.sin(angle)
    return np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])


class FrozenPathFixtureAdaptationTests(unittest.TestCase):
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




# STATIC FIXTURE ADAPTATION FROM test_visible_support_numerics.py; no monkeypatch.
"""Limited-domain and numerical ambiguity rejection components."""
import unittest
import numpy as np
from visible_support_numerics import (check_partition_numerical_domain,
    require_visible_support_numerical_domain,VisibleSupportNumericalDomainError,PASS,UNRESOLVED as NUMERICAL_UNRESOLVED)
from visible_support_path_tight import audit_visible_support_path,UNRESOLVED as PATH_UNRESOLVED


class FrozenNumericalFixtureAdaptationTests(unittest.TestCase):
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
        self.assertEqual(r['status'],NUMERICAL_UNRESOLVED)
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
        self.assertEqual(r['status'],NUMERICAL_UNRESOLVED)
        r=check_partition_numerical_domain(np.array([[2.9,0.,.5]]),np.array([[0.,0.,1.]]),2.9,1.)
        self.assertEqual(r['status'],NUMERICAL_UNRESOLVED)

    def test_enclosure_query_itself_is_domain_checked(self):
        x0=np.array([[7.,0.,-1.],[7.,0.,-2.]])
        x1=np.array([[-7.,0.,-1.],[-7.,0.,-2.]])
        r=audit_visible_support_path(x0,x1,self.i,self.j,self.pose,self.pose,self.p,self.ids)
        self.assertEqual(r['status'],PATH_UNRESOLVED)
        self.assertEqual(r['numerical_domain_failures'][0]['kind'],'enclosure')


from visible_support_path_tight import _tight_body_speed, _rotation_path, TIGHT_ROUNDOFF_FACTOR


def arbitrary_rotation(axis,angle):
    axis=np.asarray(axis,float);axis=axis/np.linalg.norm(axis)
    x,y,z=axis;K=np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])
    return np.eye(3)+np.sin(angle)*K+(1-np.cos(angle))*(K@K)


class CancellationAwareBoundTests(unittest.TestCase):
    def bounds(self,x0,x1,R0,R1,o0=None,o1=None):
        o0=np.zeros(3) if o0 is None else o0;o1=np.zeros(3) if o1 is None else o1
        angle,at=_rotation_path(R0,R1)
        rho=np.maximum(np.linalg.norm(x0-o0,axis=1),np.linalg.norm(x1-o1,axis=1))
        old=(np.linalg.norm((x1-x0)-(o1-o0),axis=1)+angle*rho)*(1+1e-11)+1e-12
        scale=max(np.max(np.abs(x0)),np.max(np.abs(x1)),np.max(np.abs(o0)),np.max(np.abs(o1)),rho.max(),1.)
        chosen,margin,info=_tight_body_speed(x0,x1,o0,o1,R0,R1,angle,at,old,1e-12,scale)
        return chosen,margin,info,angle,at

    def assert_derivative_and_chord_enclosed(self,x0,x1,R0,R1,o0=None,o1=None):
        o0=np.zeros(3) if o0 is None else o0;o1=np.zeros(3) if o1 is None else o1
        chosen,margin,info,angle,at=self.bounds(x0,x1,R0,R1,o0,o1)
        relative=R0.T@R1
        skew=.5*np.array([relative[2,1]-relative[1,2],relative[0,2]-relative[2,0],relative[1,0]-relative[0,1]])
        axis=skew/np.linalg.norm(skew);a,b,c=axis
        K=np.array([[0.,-c,b],[c,0.,-a],[-b,a,0.]]);K2=K@K
        v=(x1-x0)-(o1-o0);midr=x0+.5*(x1-x0)-o0-.5*(o1-o0);midq=midr@at(.5)
        # Analytic polynomial derivative oracle; sampled checks supplement,
        # and do not replace, the contract's derivative/curvature derivation.
        for s in np.linspace(0,1,257):
            r=x0+s*(x1-x0)-o0-s*(o1-o0)
            prime=R0@(angle*(np.cos(angle*s)*K+np.sin(angle*s)*K2))
            derivative=v@at(s)+r@prime
            self.assertTrue(np.all(np.linalg.norm(derivative,axis=1)<=chosen))
            self.assertTrue(np.all(np.linalg.norm(r@at(s)-midq,axis=1)<=abs(s-.5)*chosen+margin))
        return chosen,margin,info

    def test_corotating_world_chord_cancellation_and_curvature(self):
        x0=np.array([[.2,.1,.5],[.8,-.1,1.5]])
        R0=np.eye(3);R1=zrotation(.03);x1=x0@R1.T
        chosen,margin,info=self.assert_derivative_and_chord_enclosed(x0,x1,R0,R1)
        self.assertEqual(info['tightened_particle_count'],2)
        self.assertGreater(np.min(np.array(info['per_particle_old_speed_m'])/chosen),10.)
        self.assertGreater(margin,1e-12)
        # World chords are not rigid arcs: a nonzero curvature allowance is required.
        self.assertTrue(np.all(np.array(info['per_particle_curvature_upper_m'])>0.))

    def test_general_axis_with_world_translation(self):
        R0=arbitrary_rotation([1,2,3],.4);R1=R0@arbitrary_rotation([2,-1,4],.13)
        x0=np.array([[.3,.2,.5],[-.5,.8,1.2],[1.2,-.2,.3]])
        x1=x0+np.array([[.02,-.01,.003],[-.03,.015,-.001],[.01,.01,-.004]])
        self.assert_derivative_and_chord_enclosed(x0,x1,R0,R1,np.array([.2,-.1,.3]),np.array([.201,-.102,.304]))

    def test_world_coordinate_rotation_covariance(self):
        x0=np.array([[.3,.2,.5],[-.5,.8,1.2]]);R0=arbitrary_rotation([1,2,3],.4)
        R1=R0@arbitrary_rotation([2,-1,4],.08);x1=x0@arbitrary_rotation([1,2,3],.06).T
        Q=arbitrary_rotation([-3,1,2],.63);shift=np.array([.2,-.3,.1])
        original=self.bounds(x0,x1,R0,R1)[0]
        transformed=self.bounds(x0@Q.T+shift,x1@Q.T+shift,Q@R0,Q@R1,shift,shift)[0]
        np.testing.assert_allclose(original,transformed,rtol=2e-8,atol=1e-11)

    def test_small_nonorthogonal_R0_defect_is_in_frobenius_bound(self):
        R0=np.diag([1+1e-13,1-2e-13,1.]);R1=R0@zrotation(.02)
        x0=np.array([[.2,.1,.5],[.8,-.1,1.5]]);x1=x0@zrotation(.02).T
        _,_,info=self.assert_derivative_and_chord_enclosed(x0,x1,R0,R1)
        self.assertGreater(info['rodrigues_derivative_coefficient_upper'],4.)
        p=dict(tool_radius=3.,bore_width=2.,layer=0.,tool_top=2.)
        r=audit_visible_support_path(x0,x1,np.array([0]),np.array([1]),
            dict(origin=np.zeros(3),rotation=R0),dict(origin=np.zeros(3),rotation=R1),p,np.array([1,2]))
        self.assertEqual(r['status'],CERTIFIED)

    def test_near_cancelled_derivative_gets_absolute_error_allowance(self):
        x0=np.array([[.2,.1,.5],[.8,-.1,1.5]]);R0=np.eye(3);R1=zrotation(1e-4)
        _,_,info,_,_=self.bounds(x0,x0@R1.T,R0,R1)
        err=np.array(info['per_particle_absolute_derivative_roundoff_m'])
        residual=np.array(info['per_particle_midpoint_derivative_m'])
        self.assertTrue(np.all(err>TIGHT_ROUNDOFF_FACTOR*residual*1e6))

    def test_fallback_small_large_and_zero_angles(self):
        x=np.array([[.2,.1,.5],[.8,-.1,1.5]])
        for angle in (0.,1e-8,1.7):
            _,margin,info,_,_=self.bounds(x,x,np.eye(3),zrotation(angle))
            self.assertEqual(info['method'],'INHERITED_TRIANGLE_BOUND_FALLBACK')
            self.assertEqual(margin,1e-12)

    def test_shift_cancellation_inflates_absolute_margin(self):
        x0=np.array([[.2,.1,.5],[.8,-.1,1.5]]);R1=zrotation(.02);x1=x0@R1.T
        base=self.bounds(x0,x1,np.eye(3),R1)
        offset=np.array([1e5,-1e5,1e5])
        shifted=self.bounds(x0+offset,x1+offset,np.eye(3),R1,offset,offset)
        self.assertGreater(shifted[1],base[1]*1e4)
        self.assertTrue(np.all(np.array(shifted[2]['per_particle_absolute_derivative_roundoff_m'])>
            np.array(base[2]['per_particle_absolute_derivative_roundoff_m'])*1e4))

    def test_input_arrays_are_unchanged_in_tight_branch(self):
        x=np.array([[.2,.1,.5],[.8,-.1,1.5]]);R0=np.eye(3);R1=zrotation(.03);x1=x@R1.T
        saved=[a.copy() for a in (x,x1,R0,R1)]
        self.bounds(x,x1,R0,R1)
        for before,after in zip(saved,(x,x1,R0,R1)):np.testing.assert_array_equal(before,after)


if __name__=='__main__':unittest.main(verbosity=2)



