"""Independent force-driven stop references; no 3D scout or prescribed motion.

Literal targets and SI bounds were registered in PROCESS_STOP_BENCHMARK_CONTRACT
before execution. The scalar driver below discovers its event from a private
incoming-branch kinematic trial, not from a supplied reference stop time.
"""
import copy
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph import hht_zero_event as tested
from afsd_backend.pysph.hht_path import endpoint


def inputs(other=-2.,slip=(0.,0.,0.)):
    return dict(mass_kg=1.,normal=np.array([0.,0.,1.]),normal_force_N=1.,friction=1.,
        other_force_N=np.array([other,0.,-1.]),wall_acceleration_m_s2=np.zeros(3),
        event_slip_m_s=np.array(slip,float),velocity_scale_m_s=3.,incoming_direction=np.array([1.,0.,0.]))


def emit(name,result):
    print('HHT_ZERO_EVENT_EVIDENCE '+json.dumps(dict(name=name,
        classification='FORCE_DRIVEN_COMPONENT_BENCHMARK_NOT_3D_PROCESS',**result),sort_keys=True,allow_nan=False),flush=True)


def scalar_driver(other,outer):
    """Independent scalar force program composing the existing HHT endpoint map.

Only x/v results returned by that map are used. No state is overwritten to
match the event reference or to force a small velocity to zero.
"""
    x=np.zeros(3);v=np.array([1.,0.,0.]);t=0.;ft=-1.;acc=np.array([other+ft,0.,0.])
    totals=dict(raw_heat_J=0.,other_work_J=0.,contact_work_J=0.,contact_impulse_Ns=0.,other_impulse_Ns=0.)
    events=[];mode='SLIDING';sign=1.
    def commit(h,force_t,canonical):
        nonlocal x,v,t
        xa,va=endpoint(x,v,canonical,canonical,h,.9)
        dx=float(xa[0]-x[0])
        totals['raw_heat_J']+=-force_t*dx
        totals['contact_work_J']+=force_t*dx
        totals['other_work_J']+=other*dx
        totals['contact_impulse_Ns']+=h*force_t
        totals['other_impulse_Ns']+=h*other
        x,v,t=xa,va,t+h
    for step in range(int(round(1./outer))):
        remaining=min(outer,1.-t)
        xa,va=endpoint(x,v,acc,acc,remaining,.9)
        if mode=='SLIDING' and sign*v[0]>0 and sign*va[0]<=0:
            # Dense velocity is linear for a constant physical-force branch.
            slope=(va[0]-v[0])/remaining
            hroot=-v[0]/slope
            if not 0<hroot<=remaining:raise AssertionError('private root outside requested interval')
            old_acc=acc.copy();commit(hroot,ft,acc);remaining-=hroot
            before=(x.copy(),v.copy())
            data=inputs(other,v.copy());data['incoming_direction']=np.array([sign,0.,0.])
            result=tested.zero_event_limits(**data)
            if not np.array_equal(x,before[0]) or not np.array_equal(v,before[1]):raise AssertionError('event projected physical state')
            jump=tested.acceleration_storage_jump(old_acc,result['right_acceleration_m_s2'],1.,hroot,.9)
            events.append(dict(time_s=t,position_m=float(x[0]),velocity_m_s=float(v[0]),
                mode=result['classification'],left_acceleration_m_s2=float(old_acc[0]),
                right_acceleration_m_s2=float(result['right_acceleration_m_s2'][0]),
                storage_reference_dt_s=hroot,numerical_storage_jump_J=jump['signed_jump_J']))
            acc=result['right_acceleration_m_s2'].copy();ft=float(result['right_tangential_force_N'][0])
            mode='STICK' if result['stick_feasible'] else 'SLIDING'
            if not result['stick_feasible']:sign=float(result['outgoing_direction'][0])
        if remaining>0:commit(remaining,ft,acc)
    return dict(time_s=t,position_m=float(x[0]),velocity_m_s=float(v[0]),
        delta_kinetic_J=.5*float(v@v)-.5,events=events,**totals)


class ZeroEventReferences(unittest.TestCase):
    def close(self,actual,expected):
        np.testing.assert_allclose(actual,expected,rtol=0.,atol=1e-12)

    def test_01_force_driven_reversal_locates_event_and_balances_raw_heat(self):
        for outer in (1.,.5,.25):
            with self.subTest(outer=outer):
                r=scalar_driver(-2.,outer)
                self.assertEqual(len(r['events']),1);e=r['events'][0]
                self.close(e['time_s'],1./3.);self.close(e['position_m'],1./6.);self.close(e['velocity_m_s'],0.)
                self.assertEqual(e['mode'],'RESTART_SLIDING')
                self.close(e['left_acceleration_m_s2'],-3.);self.close(e['right_acceleration_m_s2'],-1.)
                for key,value in dict(time_s=1.,position_m=-1./18.,velocity_m_s=-2./3.,raw_heat_J=7./18.,
                    other_work_J=1./9.,contact_work_J=-7./18.,contact_impulse_Ns=1./3.,
                    other_impulse_Ns=-2.,delta_kinetic_J=-5./18.).items():self.close(r[key],value)
                self.close(r['delta_kinetic_J'],r['other_work_J']-r['raw_heat_J'])
                self.close(r['velocity_m_s']-1.,r['contact_impulse_Ns']+r['other_impulse_Ns'])
                if outer==1.:self.close(e['numerical_storage_jump_J'],-1./900.)
                emit('analytic_reversal',dict(outer_dt_s=outer,**r))

    def test_02_feasible_stick_has_static_impulse_but_no_extra_heat(self):
        for outer in (1.,.5,.25):
            with self.subTest(outer=outer):
                r=scalar_driver(-.5,outer);self.assertEqual(len(r['events']),1);e=r['events'][0]
                self.assertEqual(e['mode'],'STICK_FEASIBLE');self.close(e['time_s'],2./3.)
                for key,value in dict(position_m=1./3.,velocity_m_s=0.,raw_heat_J=1./3.,other_work_J=-1./6.,
                    contact_impulse_Ns=-.5,other_impulse_Ns=-.5,delta_kinetic_J=-.5).items():self.close(r[key],value)
                if outer==1.:self.close(e['numerical_storage_jump_J'],-1./800.)
                emit('analytic_feasible_stick',dict(outer_dt_s=outer,**r))

    def test_03_static_force_limits_wall_acceleration_and_marginal_capacity(self):
        for other,classification,right_ft,right_a in ((-2.,'RESTART_SLIDING',1.,-1.),
                (-.5,'STICK_FEASIBLE',.5,0.),(-1.,'STICK_FEASIBLE',1.,0.),
                (-1.-1e-8,'RESTART_SLIDING',1.,-1e-8)):
            with self.subTest(other=other):
                r=tested.zero_event_limits(**inputs(other))
                self.assertEqual(r['classification'],classification)
                self.close(r['right_tangential_force_N'],[right_ft,0.,0.])
                self.close(r['right_acceleration_m_s2'],[right_a,0.,0.])
                self.close(r['event_physical_impulse_Ns'],np.zeros(3))
                self.close(r['event_physical_work_J'],0.);self.close(r['event_friction_heat_J'],0.)
                self.assertFalse(r['accepted']);self.assertFalse(r['restart_allowed']);self.assertFalse(r['velocity_was_modified'])
                if other==-1.:self.assertTrue(r['feasibility_margin_is_roundoff_marginal'])
                if other==-2.:self.assertTrue(r['reverse_relative_to_incoming'])
        data=inputs(-.5);data['wall_acceleration_m_s2']=np.array([.5,0.,0.])
        r=tested.zero_event_limits(**data);self.assertTrue(r['stick_feasible'])
        self.close(r['right_tangential_force_N'],[1.,0.,0.]);self.close(r['right_acceleration_m_s2'],[.5,0.,0.])
        self.close(r['right_relative_tangential_acceleration_m_s2'],np.zeros(3))

    def test_04_complete_vector_zero_and_unsupported_scope_guards(self):
        # Scalar x=0 is not full tangential zero while y remains nonzero.
        for slip in ((0.,1.,0.),(1e-9,0.,0.),(0.,0.,1e-9)):
            with self.subTest(slip=slip):
                with self.assertRaises(ValueError):tested.zero_event_limits(**inputs(slip=slip))
        for extra in (dict(other_contact_count=1),dict(free_axes=(1,0,1)),dict(other_contact_count=True)):
            with self.subTest(extra=extra):
                with self.assertRaises(ValueError):tested.zero_event_limits(**inputs(),**extra)
        raw=np.array([1e-15,0.,0.]);r=tested.zero_event_limits(**inputs(slip=raw))
        np.testing.assert_array_equal(r['event_slip_m_s'],raw)

    def test_05_one_sided_linear_and_exact_second_order_direction(self):
        d=np.array([1.,2.,0.]);unit=d/np.sqrt(5.)
        cases=[(np.array([np.zeros(3),d,np.zeros(3)]),0,1,unit),
               (np.array([d,-d,np.zeros(3)]),1,1,unit),
               (np.array([np.zeros(3),np.zeros(3),d]),0,2,unit),
               (np.array([d,-2*d,d]),1,2,unit),
               (np.array([np.zeros(3),-d,np.zeros(3)]),0,1,-unit)]
        for coeff,end,order,direction in cases:
            with self.subTest(endpoint=end,order=order):
                before=coeff.copy();r=tested.one_sided_slip_direction(coeff,end)
                self.assertEqual(r['leading_order'],order);self.close(r['direction'],direction)
                np.testing.assert_array_equal(coeff,before);self.assertFalse(r['velocity_was_modified'])
        # Resolved tiny derivatives are not erased by an arbitrary speed floor.
        coeff=np.zeros((3,3));coeff[1]=1e-15*d
        r=tested.one_sided_slip_direction(coeff,0);self.assertEqual(r['leading_order'],1);self.close(r['direction'],unit)

    def test_06_one_sided_nonzero_full_endpoint_and_unresolved_cancel_reject(self):
        coeff=np.array([[-1.,1.,0.],[1.,0.,0.],[0.,0.,0.]])
        with self.assertRaises(ValueError):tested.one_sided_slip_direction(coeff,1)
        with self.assertRaises(ValueError):tested.one_sided_slip_direction(np.zeros((3,3)),0)
        coeff=np.zeros((3,3));coeff[:,0]=[1.-1e-14,-2.+1e-14,1.]
        with self.assertRaises(ValueError):tested.one_sided_slip_direction(coeff,1)
        with self.assertRaises(ValueError):tested.one_sided_slip_direction(np.zeros((4,3)),0)
        with self.assertRaises(ValueError):tested.one_sided_slip_direction(np.zeros((3,3)),2)
        with self.assertRaises(ValueError):tested.one_sided_slip_direction(np.zeros((3,3)),0,tolerance=1e-6)

    def test_07_signed_storage_is_not_physical_event_work_or_heat(self):
        minus=np.array([-3.,0.,0.]);plus=np.array([-1.,0.,0.]);before=(minus.copy(),plus.copy())
        r=tested.acceleration_storage_jump(minus,plus,1.,1./3.,.9)
        self.close(r['signed_jump_J'],-1./900.);self.close(r['d'],.0025)
        self.assertFalse(r['full_algorithmic_energy_jump_certified']);self.assertFalse(r['numerical_terms_are_material_heat'])
        np.testing.assert_array_equal(minus,before[0]);np.testing.assert_array_equal(plus,before[1])
        self.close(tested.acceleration_storage_jump(plus,minus,1.,1./3.,.9)['signed_jump_J'],1./900.)
        self.close(tested.acceleration_storage_jump(minus,plus,1.,1./3.,1.)['signed_jump_J'],0.)
        multi=tested.acceleration_storage_jump(np.array([minus,plus]),np.zeros((2,3)),np.array([1.,2.]),.5,.9)
        self.close(multi['signed_jump_J'],-.5*.5**2*.0025*11.)

    def test_08_input_purity_and_invalid_inputs(self):
        data=inputs();before=copy.deepcopy(data);tested.zero_event_limits(**data)
        for key in data:np.testing.assert_array_equal(data[key],before[key])
        bad=[('mass_kg',0.),('mass_kg',float('nan')),('normal_force_N',-1.),('friction',-.1),
             ('normal',np.array([0.,0.,2.])),('other_force_N',np.array([float('inf'),0.,0.])),
             ('other_force_N',np.array([1e308,1e308,-1.])),
             ('incoming_direction',np.zeros(3)),('velocity_scale_m_s',-1.)]
        for name,value in bad:
            with self.subTest(field=name):
                case=inputs();case[name]=value;original=copy.deepcopy(case)
                with self.assertRaises(ValueError):tested.zero_event_limits(**case)
                for key in case:np.testing.assert_array_equal(case[key],original[key])
        for args in (([0.,0.,0.],[0.,0.,0.],1.,0.),([0.,0.,0.],[0.,0.,0.],-1.,1.),
                     ([float('nan'),0.,0.],[0.,0.,0.],1.,1.)):
            with self.assertRaises(ValueError):tested.acceleration_storage_jump(*args)


if __name__=='__main__':unittest.main(verbosity=2)
