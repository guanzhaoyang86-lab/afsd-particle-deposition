"""Small manufactured observation fixtures, never AFSD validation."""
import copy
import unittest
import numpy as np
from deposition_stage_observer import observe_deposition_stages


def fixture():
    x=np.array([[[.75,0,.25],[1.5,0,-.25]],
                [[1.125,0,-.125],[1.5,0,-.25]],
                [[1.25,0,-.25],[1.5,0,-.25]],
                [[1.5,0,-.375],[1.5,0,-.25]]])
    return dict(time_s=np.array([0.,.25,.5,.75]),particle_id=np.array([11,19]),
        component_id=np.array([1,1]),mass_kg=np.array([2.,3.]),position_m=x,
        plastic_strain=np.array([[0.,0.],[.1,.9],[.2,1.],[.3,1.1]]),
        tool_origin_m=np.zeros((4,3)),tool_rotation=np.tile(np.eye(3),(4,1,1)),
        eligible_ids=[11],crossing_events=[dict(particle_id=11,time_s=.2,
            tool_local_position_m=[.9,0,0],mass_kg=2.)],dx_m=.5,bore_width_m=2.,tool_radius_m=3.,
        shoulder_contact=np.array([[False,True],[True,True],[False,True],[False,True]]))


class ObserverFixtures(unittest.TestCase):
    def test_initial_lower_material_never_promoted(self):
        r=observe_deposition_stages(**fixture())
        self.assertEqual(r['initial_ineligible_bar_mass_kg'],3.)
        self.assertEqual(r['ineligible_bar_actual_shoulder_contact_count'],1)
        self.assertEqual(r['final']['crossed_centre_under_annular_shoulder_mass_kg'],2.)
        self.assertFalse(r['complete_deposition_established'])
        bad=fixture();bad['crossing_events'].append(dict(particle_id=19,time_s=.2,
            tool_local_position_m=[.9,0,0],mass_kg=3.))
        with self.assertRaises(ValueError):observe_deposition_stages(**bad)

    def test_half_spacing_and_sampled_duration(self):
        r=observe_deposition_stages(**fixture());h=r['history'];p=r['crossed_particle_rows'][0]
        self.assertEqual(h[0]['eligible_centre_crossed_count'],0)
        self.assertEqual(h[1]['crossed_centre_under_annular_shoulder_count'],1)
        self.assertEqual(h[1]['crossed_laterally_half_spacing_separated_count'],0)
        self.assertEqual(h[2]['crossed_laterally_half_spacing_separated_count'],1)
        self.assertEqual(p['maximum_sampled_consecutive_separated_duration_s'],.25)
        self.assertAlmostEqual(p['post_crossing_plastic_increment'],.2)
        self.assertEqual(p['post_crossing_baseline_sample_time_s'],.25)

    def test_pure_arrays_permutation_and_pose_covariance(self):
        a=fixture();before=copy.deepcopy(a)
        for value in a.values():
            if isinstance(value,np.ndarray):value.flags.writeable=False
        r=observe_deposition_stages(**a)
        for name,value in a.items():
            if isinstance(value,np.ndarray):np.testing.assert_array_equal(value,before[name])
        b=fixture();order=np.array([1,0])
        for name in ('particle_id','component_id','mass_kg'):b[name]=b[name][order]
        for name in ('position_m','plastic_strain','shoulder_contact'):b[name]=b[name][:,order]
        self.assertEqual(r,observe_deposition_stages(**b))
        rotation=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
        b=fixture();shift=np.array([4.,5.,6.])
        b['position_m']=b['position_m']@rotation.T+shift
        b['tool_origin_m'][:]=shift;b['tool_rotation'][:]=rotation
        s=observe_deposition_stages(**b)
        self.assertEqual(r['final'],s['final'])
        self.assertEqual(r['crossed_particle_rows'][0]['tool_local_position_m'],s['crossed_particle_rows'][0]['tool_local_position_m'])

    def test_missing_contact_is_unavailable_not_measured_zero(self):
        a=fixture();a.pop('shoulder_contact');r=observe_deposition_stages(**a)
        self.assertFalse(r['final']['contact_observation_available'])
        a=fixture();a['shoulder_contact_available']=np.array([False,False,False,False])
        r=observe_deposition_stages(**a)
        self.assertEqual(r['crossed_particle_rows'][0]['actual_shoulder_contact_sample_count'],0)

    def test_invalid_ids_aperture_and_negative_ep_not_clipped(self):
        for kind in ('duplicate','bad_time','on_edge'):
            a=fixture()
            if kind=='duplicate':a['particle_id']=np.array([11,11])
            elif kind=='bad_time':a['time_s'][2]=a['time_s'][1]
            else:a['crossing_events'][0]['tool_local_position_m']=[1.,0.,0.]
            with self.assertRaises(ValueError):observe_deposition_stages(**a)
        a=fixture();a['plastic_strain'][-1,0]=.05
        r=observe_deposition_stages(**a)
        self.assertAlmostEqual(r['crossed_particle_rows'][0]['post_crossing_plastic_increment'],-.05)
        self.assertEqual(r['final']['crossed_separated_with_post_crossing_plastic_increment_count'],0)


if __name__=='__main__':unittest.main(verbosity=2)
