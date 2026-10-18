"""Independent literal journals; never invoke mechanics/contact/material code."""
from pathlib import Path
from types import SimpleNamespace
import copy
import json
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from afsd_bridge_audit import audit


def fixture(kind='free'):
    """Hand-balanced unit-mass cases; this is not an AFSD load programme.

free: F=2ex, v0=0 -> v1=2ex, W=2J.
actuator: vz=1 constant, F=3ez, Jc=-3ez, Wi=3J, Wc=-3J.
friction: Fn=1ez balanced by F=-ez, Jt=-.5ex, v:1->.5,
          wall=.25ex, Q=.25J, Wwall=-.125J, deltaK=-.375J.
plastic: q=3Pa, dep=.2, V=1m3 => P=.6J, return=.06J,
         chi=.5, C=1J/K => T1=300.3K.
pair: two C=1J/K nodes, G=1W/K, dt=1s, T0=(300,400),
      exact BE T1=(1000/3,1100/3), q_into0=100/3J.
"""
    count=2 if kind=='pair' else 1
    p=dict(shear=1.,rate_floor=1.,contact_conductance=0.,dx=1.,
        base_temperature=400.,tool_temperature=500.,conductivity=.5,
        plastic_heat_fraction=.5,heat_partition=.5,cp=1.,friction=.5)
    pos=np.zeros((count,3));pos[:,0]=np.arange(count)
    v0=np.zeros((count,3));v1=v0.copy();forces=v0.copy();free=np.ones((count,3));reaction=v0.copy()
    T0=np.full(count,300.);T1=T0.copy();ep0=np.zeros(count);ep1=ep0.copy()
    stress=np.zeros((count,3,3));plastic=np.zeros(count);returned=np.zeros(count);rate=np.zeros(count)
    if kind=='free':v1[0,0]=2.;forces[0,0]=2.
    if kind=='actuator':
        v0[0,2]=1.;v1[0,2]=1.;forces[0,2]=3.;free[0,2]=0.;reaction[0,2]=-3.
    if kind in ('friction','mixed'):v0[0,0]=1.;v1[0,0]=.5;forces[0,2]=-1.;T1[0]=300.125
    if kind in ('plastic','mixed'):
        ep0[0]=.1;ep1[0]=.3;stress[0]=np.diag([2.,-1.,-1.]);plastic[0]=.6;returned[0]=.06;rate[0]=2.;T1[0]=300.3
    if kind=='mixed':T1[0]=300.4;p['heat_partition']=.4
    if kind=='pair':T0[:]=[300.,400.];T1[:]=[1000/3,1100/3]
    vm=.5*(v0+v1);m=np.ones(count);rho=m.copy();ids=np.arange(11,11+count,dtype=np.int64)
    pa=SimpleNamespace(m=m,rho=rho,u=v0[:,0],v=v0[:,1],w=v0[:,2],x=pos[:,0],y=pos[:,1],z=pos[:,2],
        temperature=T0,ep=ep0,particle_id=ids)
    totals={key:0. for key in ('friction_raw_J','normal_material_work_J','normal_wall_work_J','normal_relative_work_J',
        'wall_work_J','actuator_work_J','internal_work_J','delta_kinetic_J')}
    traces=[];qpart=np.zeros(count);roundpart=np.zeros(count)
    for i in range(count):
        contact=kind in ('friction','mixed');nc=int(contact);normals=np.array([[0.,0.,1.]]) if contact else np.zeros((0,3))
        walls=np.array([[.25,0.,0.]]) if contact else np.zeros((0,3));fn=np.ones(nc)
        jn=normals.copy();jt=np.array([[-.5,0.,0.]]) if contact else np.zeros((0,3))
        faces=[dict(contact_id=f'{int(ids[i])}:0',feature_id=0,particle_id=int(ids[i]),particle_index=i,
            normal=np.array([0.,0.,1.]),wall_velocity=np.array([.25,0.,0.]),normal_force_N=1.,wall_temperature_K=400.)] if contact else []
        frame=walls[0].copy() if contact else v0[i].copy()
        relative0=v0[i]-frame;relative1=v1[i]-frame;wr=walls-frame
        speed=max(np.linalg.norm(relative0),float(np.linalg.norm(wr)) if contact else 0.,
                  np.linalg.norm(free[i]*(forces[i]+jn.sum(axis=0))),np.finfo(float).tiny)
        tangent=np.array([[.5,0.,0.]]) if contact else np.zeros((0,3))
        q=np.array([.25]) if contact else np.zeros(0)
        bound=256*np.finfo(float).eps*np.linalg.norm(jt,axis=1)*(speed+np.linalg.norm(relative0)+
            np.linalg.norm(relative1)+2*np.linalg.norm(wr,axis=1))
        wi=2. if kind=='free' else 3. if kind=='actuator' else 0.
        wc=-3. if kind=='actuator' else 0.;ww=-.125 if contact else 0.;dk=2. if kind=='free' else -.375 if contact else 0.
        interval=dict(start_s=0.,dt_s=1.,kind='LITERAL_REFERENCE',v_before=v0[i].copy(),v_after=v1[i].copy(),
            normal_impulse=jn,tangential_impulse=jt,constraint_impulse=reaction[i].copy(),
            friction_heat_raw_J=q,friction_rounding_bound_J=bound,tangential_relative_midpoint=tangent,
            reference_frame_velocity=frame,relative_speed_scale_m_s=speed,other_work_J=wi,normal_work_J=np.zeros(nc),
            constraint_work_J=wc,wall_work_J=ww,delta_kinetic_J=dk)
        result=dict(velocity=v1[i].copy(),normal_impulse=jn.copy(),tangential_impulse=jt.copy(),
            constraint_impulse=reaction[i].copy(),friction_heat_J=q.copy(),other_work_J=wi,
            constraint_work_J=wc,wall_work_J=ww,delta_kinetic_J=dk,intervals=[interval],events=[])
        traces.append(dict(particle_index=i,particle_id=int(ids[i]),mass_kg=1.,material_velocity_before_m_s=v0[i].copy(),
            internal_force_N=forces[i].copy(),free_axes=free[i].copy(),prescribed_velocity=v0[i].copy(),
            actuator_active=kind=='actuator',dt_s=1.,time_mid_s=.5,contact_rows=faces,result=result))
        qpart[i]=q.sum();roundpart[i]=bound.sum()
        for key,val in [('friction_raw_J',q.sum()),('wall_work_J',ww),('actuator_work_J',wc),('internal_work_J',wi),('delta_kinetic_J',dk)]:totals[key]+=val
    first=np.array([0],int) if kind=='pair' else np.array([],int)
    second=np.array([1],int) if kind=='pair' else np.array([],int)
    gradient=np.array([[1.,0.,0.]]) if kind=='pair' else np.zeros((0,3))
    conductance=np.ones(1) if kind=='pair' else np.zeros(0)
    pair=T1[second]-T1[first];node=np.array([pair[0],-pair[0]]) if kind=='pair' else np.zeros(count)
    partition=p['heat_partition']
    particle=dict(plastic_work_J=plastic.copy(),plastic_heat_J=.5*plastic,plastic_unthermalized_J=.5*plastic,
        friction_raw_J=qpart.copy(),friction_material_heat_J=partition*qpart,friction_wall_heat_J=(1-partition)*qpart,
        contact_in_J=np.zeros(count),internal_conduction_J=node)
    total=.5*plastic+partition*qpart+node;measured=T1-T0
    particle.update(total_heat_input_J=total,measured_thermal_increment_J=measured,mc_delta_T_residual_J=measured-total)
    ledger={key:float(value.sum()) for key,value in particle.items()};ledger['wall_heat_out_J']=float(((1-partition)*qpart).sum())
    e=dict(classification='PROCESS_GLOBAL_EVENT_AVERAGE_BE_THERMAL_APPROXIMATE',accepted=False,restart_allowed=False,
        production_authorized=False,temperature_layers='FKS_TMID; PAIR_AND_WALL_CONDUCTION_T1_BE; INCREMENT_FROM_T0',
        particle_id=ids.copy(),time_before_s=0.,time_after_s=1.,time_mid_s=.5,dt_s=1.,
        map=dict(position_m=pos+vm,velocity_m_s=v1,mean_velocity_m_s=vm,temperature_K=T1),
        midpoint=dict(density_kg_m3=rho.copy(),position_m=pos.copy(),i=first,j=second,g=gradient),
        material=dict(ep=ep1,deviatoric_stress_Pa=stress,plastic_work_J=plastic,return_dissipation_J=returned,strain_rate_per_s=rate),
        internal_force_N=forces,boundary=dict(velocity=v1.copy(),work_velocity=vm.copy(),traces=traces,
            friction_raw_J=qpart.copy(),wall_row_W_K=np.zeros(count),contact_power_W=np.zeros(count),ledger=totals),
        thermal_boundary_certificate=dict(wall_row_W_K=np.zeros(count),wall_weighted_temperature_W=np.zeros(count),
            friction_raw_J=qpart.copy(),friction_rounding_bound_J=roundpart.copy()),
        thermal=dict(temperature=T1.copy(),pair_conductance_W_K=conductance,pair_heat_into_i_J=pair,
            particle_energy_J=particle,ledger=ledger))
    return pa,p,e


class BridgeAuditTests(unittest.TestCase):
    def check_pass(self,kind):
        pa,p,e=fixture(kind);before=copy.deepcopy((pa,p,e));r=audit(pa,p,0.,1.,e)
        self.assertTrue(r['passed'],r);self.assertFalse(r['accepted']);self.assertFalse(r['can_publish'])
        json.dumps(r,allow_nan=False)
        self.assertEqual(repr((pa,p,e)),repr(before))
        return r

    def test_free_no_contact_impulse_and_work(self):
        r=self.check_pass('free');self.assertEqual(r['summary']['contact_count'],0)
        self.assertEqual(r['summary']['contact_work_totals_J']['internal_work_J'],2.)

    def test_actuator_reaction_and_signed_work(self):
        r=self.check_pass('actuator');self.assertEqual(r['summary']['constraint_work_J'],-3.)
        self.assertEqual(r['summary']['nonzero_constraint_impulse_particle_count'],1)

    def test_positive_same_history_plastic_source(self):
        r=self.check_pass('plastic');self.assertEqual(r['summary']['plastic_particle_ids'],[11])
        self.assertAlmostEqual(r['summary']['plastic_work_J'],.6,places=14)
        self.assertFalse(r['summary']['return_dissipation_is_material_heat'])

    def test_moving_wall_friction_full_work_identity(self):
        r=self.check_pass('friction');self.assertEqual(r['summary']['raw_friction_J'],.25)
        self.assertEqual(r['summary']['contact_work_totals_J']['wall_work_J'],-.125)

    def test_pair_conduction_independent_two_node_solution(self):
        r=self.check_pass('pair');self.assertAlmostEqual(r['summary']['thermal_pair_net_J'],0.)

    def test_wrong_impulse_or_constraint_work_rejected(self):
        for field in ['constraint_work_J','constraint_impulse']:
            pa,p,e=fixture('actuator');e['boundary']['traces'][0]['result']['intervals'][0][field]*=-1
            self.assertFalse(audit(pa,p,0.,1.,e)['passed'])

    def test_negative_heat_or_inflated_certificate_rejected(self):
        for field in ['friction_heat_raw_J','friction_rounding_bound_J']:
            pa,p,e=fixture('friction');e['boundary']['traces'][0]['result']['intervals'][0][field]*=(-1 if field=='friction_heat_raw_J' else 100.)
            self.assertFalse(audit(pa,p,0.,1.,e)['passed'])

    def test_double_plastic_heat_or_return_heat_rejected(self):
        for field in ['plastic_heat_J','total_heat_input_J']:
            pa,p,e=fixture('plastic');e['thermal']['particle_energy_J'][field]+=.06
            self.assertFalse(audit(pa,p,0.,1.,e)['passed'])

    def test_wrong_pair_direction_or_material_work_rejected(self):
        pa,p,e=fixture('pair');e['thermal']['pair_heat_into_i_J']*=-1
        self.assertFalse(audit(pa,p,0.,1.,e)['passed'])
        pa,p,e=fixture('plastic');e['material']['plastic_work_J']*=2
        self.assertFalse(audit(pa,p,0.,1.,e)['passed'])

    def test_missing_nonfinite_identity_and_false_authority_rejected(self):
        for change in ['missing','nonfinite','id','authority']:
            pa,p,e=fixture('free')
            if change=='missing':del e['material']['ep']
            if change=='nonfinite':e['map']['temperature_K'][0]=np.nan
            if change=='id':e['particle_id'][0]=99
            if change=='authority':e['accepted']=True
            r=audit(pa,p,0.,1.,e);self.assertFalse(r['passed']);json.dumps(r,allow_nan=False)

    def test_nonbinary_near_zero_balance_different_legal_reduction_order(self):
        pa,p,e=fixture('mixed');particle=e['thermal']['particle_energy_J']
        before=particle['mc_delta_T_residual_J'].copy()
        # Independent reassociation of the same operands, not a heat source.
        reordered=((particle['measured_thermal_increment_J']-particle['plastic_heat_J'])-
            particle['friction_material_heat_J'])-particle['contact_in_J']-particle['internal_conduction_J']
        self.assertNotEqual(float(reordered[0]),float(before[0]))
        particle['mc_delta_T_residual_J']=reordered
        e['thermal']['ledger']['mc_delta_T_residual_J']=float(reordered.sum())
        r=audit(pa,p,0.,1.,e);self.assertTrue(r['passed'],r)
        # Neither the actual T nor any heat input was modified. A macroscopic
        # forged residual still fails the narrowly scoped arithmetic bound.
        particle['mc_delta_T_residual_J']+=1e-6
        self.assertFalse(audit(pa,p,0.,1.,e)['passed'])


if __name__=='__main__':unittest.main(verbosity=2)
