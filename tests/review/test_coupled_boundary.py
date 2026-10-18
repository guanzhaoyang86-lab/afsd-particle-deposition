"""Independent frozen-kick/geometry/heat-interface references, NOT trajectories."""
from pathlib import Path
import inspect
import sys
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph import coupled_boundary as boundary


def parameters(**overrides):
    p=dict(dx=.2,layer=1.,tool_top=2.,tool_radius=2.,bore_width=1.,substrate_depth=1.,
        traverse=0.,omega=0.,base_temperature=300.,tool_temperature=450.,normal_stiffness=20.,
        friction=0.,contact_conductance=10.,bar_length=2.,actuator_depth=.25,feed=.5,release_time=10.)
    p.update(overrides);return p


def one(x=(0.,0.,0.),v=(1.,-2.,-.5),force=(4.,-8.,6.),z0=1.9,time=0.,dt=.25,mass=2.,p=None):
    return boundary.evaluate_boundary(np.array([x]),np.array([v]),np.array([mass]),np.array([400.]),
        time,dt,parameters() if p is None else p,np.array([force]),np.array([z0]),np.array([101]))


class UnifiedBoundary(unittest.TestCase):
    def test_01_uncontacted_axial_reaction_literal(self):
        r=one();np.testing.assert_array_equal(r['velocity'],[[1.5,-3.,-.5]])
        np.testing.assert_array_equal(r['ledger']['actuator_impulse_Ns'],[0.,0.,-1.5])
        self.assertEqual(r['ledger']['actuator_work_J'],.75)
        self.assertEqual(r['ledger']['internal_work_J'],5.5)
        self.assertEqual(r['ledger']['delta_kinetic_J'],6.25)
        self.assertEqual(r['ledger']['contacts'],0);self.assertEqual(len(r['traces']),1)
        self.assertEqual(r['ledger']['kinetic_work_residual_J'],0.)

    def test_02_base_load_is_coupled_with_actuator_not_after_overwrite(self):
        r=one(x=(0.,0.,-.91),p=parameters(normal_stiffness=200.))
        np.testing.assert_allclose(r['velocity'],[[1.5,-3.,-.5]],rtol=0.,atol=1e-12)
        np.testing.assert_allclose(r['ledger']['base_impulse_Ns'],[0.,0.,.5],rtol=0.,atol=1e-12)
        np.testing.assert_allclose(r['ledger']['actuator_impulse_Ns'],[0.,0.,-2.],rtol=0.,atol=1e-12)
        self.assertAlmostEqual(r['ledger']['actuator_work_J'],1.,delta=1e-12)
        self.assertAlmostEqual(r['ledger']['normal_material_work_J'],-.25,delta=1e-12)
        self.assertEqual(r['ledger']['normal_wall_work_J'],0.)
        self.assertAlmostEqual(r['ledger']['delta_kinetic_J'],6.25,delta=1e-12)

    def test_03_release_keeps_initial_actual_velocity_and_no_reaction(self):
        r=one(v=(1.5,-3.,-.5),time=11.)
        np.testing.assert_array_equal(r['velocity'],[[2.,-4.,.25]])
        np.testing.assert_array_equal(r['ledger']['actuator_impulse_Ns'],[0.,0.,0.])
        self.assertEqual(r['ledger']['actuator_work_J'],0.)
        self.assertEqual(r['ledger']['internal_work_J'],8.5625)
        self.assertEqual(r['ledger']['delta_kinetic_J'],8.5625)

    def test_04_two_bore_faces_share_mobility_and_reverse_slip(self):
        # m=1, two muFn=1 contacts; z is their common free tangent.
        # dz acceleration -6 to stop1/6, then-2 for5/6: v=-5/3.
        r=one(x=(.45,.45,1.5),v=(0.,0.,1.),force=(1.,1.,-4.),z0=0.,dt=1.,mass=1.,p=parameters(friction=1.))
        np.testing.assert_allclose(r['velocity'],[[0.,0.,-5/3]],rtol=0.,atol=1e-10)
        self.assertAlmostEqual(r['ledger']['friction_raw_J'],14/9,delta=1e-10)
        self.assertAlmostEqual(r['ledger']['internal_work_J'],22/9,delta=1e-10)
        self.assertAlmostEqual(r['ledger']['delta_kinetic_J'],8/9,delta=1e-10)
        trace=r['traces'][0];np.testing.assert_array_equal([f['feature_id'] for f in trace['contact_rows']],[3,5])
        self.assertAlmostEqual(trace['result']['events'][0]['time_s'],1/6,delta=1e-10)
        self.assertEqual(trace['result']['contact_mobility'][0,1,2,2],1.)

    def test_05_two_faces_and_z_actuator_keep_xy_free_and_account_heat(self):
        r=one(x=(.45,.45,1.5),v=(0.,0.,-.5),force=(1.,1.,0.),z0=1.9,dt=.1,mass=1.,p=parameters(friction=1.))
        np.testing.assert_allclose(r['velocity'],[[0.,0.,-.5]],rtol=0.,atol=1e-12)
        np.testing.assert_allclose(r['ledger']['actuator_impulse_Ns'],[0.,0.,-.2],rtol=0.,atol=1e-12)
        self.assertAlmostEqual(r['ledger']['friction_raw_J'],.1,delta=1e-12)
        self.assertAlmostEqual(r['ledger']['actuator_work_J'],.1,delta=1e-12)
        np.testing.assert_array_equal(r['traces'][0]['free_axes'],[1.,1.,0.])

    def test_06_rotating_bore_same_centre_zero_slip_and_work_torque(self):
        # x(.45,.2), omega2: wallcentre(-.4,.9); wallface(-.4,1).
        # NormalJn=(-.1,0,0), Wwall=.04; opposite internalW=-.04.
        r=one(x=(.45,.2,1.5),v=(-.4,.9,0.),force=(1.,0.,0.),z0=0.,dt=.1,mass=1.,p=parameters(omega=2.,friction=.3))
        np.testing.assert_allclose(r['velocity'],[[-.4,.9,0.]],rtol=0.,atol=1e-12)
        self.assertAlmostEqual(r['ledger']['friction_raw_J'],0.,delta=1e-12)
        self.assertAlmostEqual(r['ledger']['wall_work_J'],.04,delta=1e-12)
        self.assertAlmostEqual(r['ledger']['internal_work_J'],-.04,delta=1e-12)
        self.assertAlmostEqual(2*r['ledger']['tool_angular_impulse_Nms'][2],.04,delta=1e-12)
        row=r['traces'][0]['contact_rows'][0]
        np.testing.assert_allclose(row['wall_velocity'],[-.4,.9,0.],rtol=0.,atol=1e-12)
        np.testing.assert_allclose(row['wall_velocity_at_face'],[-.4,1.,0.],rtol=0.,atol=1e-12)

    def test_07_open_rotating_square_is_not_solid_cylinder(self):
        x=np.array([[0.,0.,1.5],[.1,.1,.5]])
        r=boundary.evaluate_boundary(x,np.zeros((2,3)),np.ones(2),np.full(2,400.),.3,.01,
            parameters(omega=2.),np.zeros((2,3)),np.zeros(2),np.array([7,9]))
        self.assertEqual(len(r['rows']['particle_id']),0)
        np.testing.assert_array_equal(r['velocity'],np.zeros((2,3)))

    def test_08_contact_power_is_only_signed_conduction_not_friction(self):
        x=np.array([[0.,0.,-.91],[.45,.2,1.5]])
        r=boundary.evaluate_boundary(x,np.zeros((2,3)),np.ones(2),np.array([350.,400.]),0.,.01,
            parameters(),np.zeros((2,3)),np.zeros(2),np.array([7,9]))
        np.testing.assert_allclose(r['contact_power_W'],[-20.,20.],rtol=0.,atol=1e-12)
        np.testing.assert_allclose(r['wall_row_W_K'],[.4,.4],rtol=0.,atol=1e-12)
        np.testing.assert_array_equal(r['contact_power_W'],r['contact_conduction_W'])
        np.testing.assert_array_equal(r['friction_raw_J'],[0.,0.])
        self.assertAlmostEqual(r['ledger']['contact_conduction_into_material_J'],0.,delta=1e-12)
        # p has no heat_partition: its only owner is coupled_thermal.
        self.assertNotIn('heat_partition',parameters())

    def test_09_actual_feature_potential_keeps_corner_faces_separate(self):
        r=boundary.normal_snapshot(np.array([[.45,.45,1.5],[0.,0.,-.91]]),np.array([8,2]),0.,parameters())
        self.assertEqual(set(r['active_contact_ids']),{'8:3','8:5','2:0'})
        self.assertAlmostEqual(r['potential_by_contact_id_J']['8:3'],.025,delta=1e-12)
        self.assertAlmostEqual(r['potential_by_contact_id_J']['8:5'],.025,delta=1e-12)
        self.assertAlmostEqual(r['potential_by_contact_id_J']['2:0'],.001,delta=1e-12)
        self.assertAlmostEqual(r['normal_potential_J'],.051,delta=1e-12)

    def test_10_finite_edge_disappearance_is_exposed_not_heat_fixed(self):
        # Known unsupported old finite-side edge: the caller must reject/report
        # disappearing nonzero potential, not call this a valid energy closure.
        before=boundary.normal_snapshot(np.array([[.45,.2,1.+1e-8]]),np.array([8]),0.,parameters())
        after=boundary.normal_snapshot(np.array([[.45,.2,1.-1e-8]]),np.array([8]),0.,parameters())
        self.assertEqual(before['active_contact_ids'],['8:3']);self.assertEqual(after['active_contact_ids'],[])
        self.assertAlmostEqual(before['normal_potential_J'],.025,delta=1e-12)
        self.assertEqual(after['normal_potential_J'],0.)

    def test_11_label_free_api_and_original_material_selection(self):
        signature=inspect.signature(boundary.evaluate_boundary)
        self.assertFalse(any('tracer' in k or 'as_rs' in k for k in signature.parameters))
        driven=one(z0=1.9);free=one(z0=0.)
        np.testing.assert_array_equal(driven['traces'][0]['free_axes'],[1.,1.,0.])
        np.testing.assert_array_equal(free['traces'][0]['free_axes'],[1.,1.,1.])
        self.assertEqual(free['velocity'][0,2],.25)

    def test_12_permutation_preserves_id_contact_and_impulses(self):
        x=np.array([[.45,.2,1.5],[0.,0.,-.91]]);ids=np.array([91,4]);v=np.zeros((2,3))
        args=dict(time=0.,dt=.01,p=parameters(),force=np.zeros((2,3)),z0=np.zeros(2))
        a=boundary.evaluate_boundary(x,v,np.ones(2),np.full(2,400.),particle_ids=ids,**args)
        b=boundary.evaluate_boundary(x[::-1],v,np.ones(2),np.full(2,400.),particle_ids=ids[::-1],**args)
        np.testing.assert_array_equal(a['velocity'],b['velocity'][::-1])
        self.assertEqual({f['contact_id'] for t in a['traces'] for f in t['contact_rows']},
                         {f['contact_id'] for t in b['traces'] for f in t['contact_rows']})

    def test_13_failure_keeps_particle_feature_and_original_exception_context(self):
        failure=RuntimeError('INJECTED_CONE_FAILURE');failure.cone_context={'independent_failure_marker':31}
        with patch.object(boundary,'kick',side_effect=failure):
            with self.assertRaises(boundary.BoundaryEvaluationError) as raised:
                one(x=(.45,.45,1.5))
        context=raised.exception.contact_context
        self.assertEqual(context['particle_id'],101);self.assertEqual(context['contact_ids'],['101:3','101:5'])
        self.assertEqual(context['failed_cone'],{'independent_failure_marker':31})
        self.assertFalse(context['restart_allowed'])
        np.testing.assert_array_equal(context['free_axes'],[1.,1.,0.])

    def test_14_bad_physical_inputs_and_inconsistent_constraint_fail(self):
        for p in (parameters(normal_stiffness=-1.),parameters(contact_conductance=-1.),parameters(friction=-1.)):
            with self.assertRaises(ValueError):one(p=p)
        with self.assertRaises(ValueError):one(mass=0.)
        with self.assertRaisesRegex(boundary.BoundaryEvaluationError,'inconsistent prescribed initial velocity'):
            one(v=(0.,0.,0.))


if __name__=='__main__':unittest.main(verbosity=2)
