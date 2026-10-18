"""Independent slow-boundary wiring and mocked publication rollback checks.

No trajectory/integrator step is run. Literal SI wall/contact references do
not call the slow independent reference or contact routines under test.
Roundoff bounds: boundary force2e-11N, geometry2e-18m, wall velocity2e-20m/s;
current-SPH permutation force5e-12N, acceleration2e-8m/s2. Synthetic states
in publication tests are mock evidence only, never physical output.
"""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))

import probe_hht_slow as driver
from afsd_backend.pysph.mechanics import positions, velocities
from afsd_backend.pysph.hht_coupled import initialize
from afsd_backend.pysph.updated_trial import particle_hash


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def record(name,**values):
    print('HHT_SLOW_RUNNER_EVIDENCE '+json.dumps(dict(name=name,
          trajectory_advanced=False,**values),sort_keys=True,allow_nan=False),flush=True)


class SlowRunner(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frozen_paths=[ROOT/name for name in (
            'scripts/probe_hht_coupled.py','scripts/reference_hht_block.py',
            'src/afsd_backend/pysph/hht_coupled.py','docs/PROCESS_HHT_FREEZE.json',
            'docs/PROCESS_HHT_REPORT.md',
            'outputs/process_hht_20261003/sweep_repair1/grid_032/accepted_history.json',
            'outputs/process_hht_20261003/sweep_repair1/grid_032/stage_032/accepted.npz')]
        cls.before={str(path):digest(path) for path in cls.frozen_paths}

    @classmethod
    def tearDownClass(cls):
        after={str(path):digest(path) for path in cls.frozen_paths}
        if after!=cls.before: raise AssertionError('old source/history/checkpoint changed during wiring tests')
        record('old_evidence_untouched',sha256=after)

    def test_01_actual_pa_literal_fields_and_all_free_actuator_inventory(self):
        pa,p=driver.setup()
        ijk=np.array([(i,j,k) for i in range(3) for j in range(3) for k in range(3)])
        xyz=.005*ijk;vol=.005**3*np.prod(np.where(ijk==1,1.,.5),axis=1)
        vel=np.zeros((27,3));vel[:,2]=.02
        np.testing.assert_array_equal(positions(pa),xyz)
        np.testing.assert_array_equal(velocities(pa),vel)
        np.testing.assert_array_equal(pa.m,2700.*vol)
        np.testing.assert_array_equal(pa.rho,np.full(27,2700.))
        np.testing.assert_allclose(pa.h,np.full(27,.0065),atol=1e-18,rtol=0.)
        np.testing.assert_array_equal(pa.S.reshape(27,3,3),np.tile(np.diag([2e4,-1e4,-1e4]),(27,1,1)))
        np.testing.assert_array_equal(pa.temperature,np.repeat([400.,401.,402.],9))
        for field in ('actuator','component_id','ep','plastic_work','return_dissipation'):
            np.testing.assert_array_equal(getattr(pa,field),np.zeros(27),err_msg=field)
        np.testing.assert_array_equal(pa.particle_id,np.arange(1,28))
        np.testing.assert_array_equal(np.column_stack((pa.x0,pa.y0,pa.z0)),xyz)
        self.assertEqual(p['young'],70e9);self.assertEqual(p['hht_weight'],.9)
        before=particle_hash(pa);state=initialize(pa,p,driver.boundary_for(p))
        self.assertEqual(particle_hash(pa),before)
        self.assertGreater(float(np.linalg.norm(state['acceleration_m_s2'])),0.)
        np.testing.assert_allclose(state['acceleration_m_s2'],state['force_N']/pa.m[:,None],atol=0.,rtol=0.)
        np.testing.assert_array_equal(state['boundary_snapshot']['noncontact_force_N'],np.zeros((27,3)))
        record('actual_pa_startup',physical_dofs=81,actuators=0,nonzero_initial_force=True,
               physical_state_unchanged_by_initialization=True)

    def test_02_literal_wall_program_force_and_no_microsecond_jump(self):
        pa,p=driver.setup();boundary=driver.boundary_for(p);before=particle_hash(pa)
        layers=[0.,1e-6-1e-12,1e-6,1e-6+1e-12,20e-6,.001]
        for t in layers:
            rows=boundary(pa,t)['contacts'];self.assertEqual(len(rows),18)
            omega=10.*np.pi
            stroke=2e-7*np.sin(.5*omega*t)**2
            wall_dx=1e-7*omega*np.sin(omega*t)
            for row in rows:
                face=int(row['contact_id'].split(':')[0]);sign=1. if face==0 else -1.
                np.testing.assert_allclose(row['wall_velocity_m_s'],[sign*wall_dx,0.,.2],atol=2e-20,rtol=0.)
                self.assertAlmostEqual(row['penetration_m'],2e-6+stroke,delta=2e-18)
                self.assertAlmostEqual(row['normal_force_N'],7e6*(2e-6+stroke),delta=2e-11)
                np.testing.assert_allclose(row['tangential_force_N'],[0.,0.,.2*7e6*(2e-6+stroke)],atol=2e-11,rtol=0.)
                self.assertEqual(row['conductance_W_K'],.025)
                self.assertEqual(row['wall_temperature_K'],[410.,420.][face])
        self.assertEqual(particle_hash(pa),before)
        left=boundary(pa,1e-6-1e-12)['contacts'];right=boundary(pa,1e-6+1e-12)['contacts']
        for a,b in zip(left,right):
            self.assertGreater(a['tangential_force_N'][2],0.)
            self.assertGreater(b['tangential_force_N'][2],0.)
            self.assertEqual(a['wall_velocity_m_s'][2],b['wall_velocity_m_s'][2])
        record('wall_program',normal_frequency_Hz=5.,constant_tangent_m_s=.2,
               microsecond_control_jump=False)

    def test_03_full_vector_coulomb_accepts_curved_reversal_and_low_slip(self):
        pa,p=driver.setup();boundary=driver.boundary_for(p)
        # Each case is a direct wall-law probe, not an evolved particle state.
        for vy,vz in ((.03,.21),(.03,.2),(0.,.20000001),(0.,.19999999)):
            pa.v[:]=vy;pa.w[:]=vz;before=particle_hash(pa)
            rows=boundary(pa,0.)['contacts'];u=np.array([0.,vy,vz-.2])
            expected=-2.8*u/np.sqrt(vy*vy+(vz-.2)**2)
            for row in rows:
                np.testing.assert_allclose(row['tangential_force_N'],expected,atol=2e-11,rtol=0.)
                self.assertLess(float(row['tangential_force_N']@u),0.)
                self.assertAlmostEqual(float(np.linalg.norm(row['tangential_force_N'])),2.8,delta=2e-11)
            self.assertEqual(particle_hash(pa),before)
        record('full_vector_sliding',positive_z_with_transverse_slip='ACCEPTED',
               zero_z_with_transverse_slip='ACCEPTED',slip_below_point05='ACCEPTED_RAW_COULOMB')

    def test_04_exact_zero_and_normal_opening_reject_without_particle_mutation(self):
        pa,p=driver.setup();journal=[];boundary=driver.boundary_for(p,journal)
        pa.v[:]=0.;pa.w[:]=.2;before=particle_hash(pa)
        with self.assertRaisesRegex(ValueError,'EXACT_ZERO_SLIP_BOUNDARY'):boundary(pa,0.)
        self.assertEqual(particle_hash(pa),before)
        self.assertEqual(len(journal[-1]['contacts']),18)
        self.assertFalse(journal[-1]['accepted'])
        pa.w[:]=.02;pa.x[0]=3e-6;before=particle_hash(pa)
        with self.assertRaisesRegex(ValueError,'OPENING_BOUNDARY'):boundary(pa,0.)
        self.assertEqual(particle_hash(pa),before)
        failed=[row for row in journal[-1]['contacts'] if row['penetration_m']<0]
        self.assertEqual([row['particle_id'] for row in failed],[1])
        record('unsupported_domain',exact_zero_slip='REJECTED',opening='REJECTED',
               gap_not_clipped=True,physical_fields_not_overwritten=True)

    def test_05_labels_and_particle_reorder_preserve_contact_and_total_initial_force(self):
        pa,p=driver.setup();boundary=driver.boundary_for(p)
        order=(7*np.arange(27)+4)%27;changed=pa.extract_particles(order.astype(np.uint32))
        changed.tracer_type[:]=(3*np.arange(27)+5)%7;changed.wholebar_tag[:]=np.arange(27)[::-1]%11
        state0=initialize(pa,p,boundary);state1=initialize(changed,p,boundary)
        np.testing.assert_allclose(state1['force_N'],state0['force_N'][order],atol=5e-12,rtol=0.)
        np.testing.assert_allclose(state1['acceleration_m_s2'],state0['acceleration_m_s2'][order],atol=2e-8,rtol=0.)
        rows0={row['contact_id']:row for row in state0['boundary_snapshot']['contacts']}
        rows1={row['contact_id']:row for row in state1['boundary_snapshot']['contacts']}
        self.assertEqual(set(rows0),set(rows1))
        for cid,row in rows0.items():
            altered=rows1[cid]
            self.assertEqual(int(changed.particle_id[altered['particle_index']]),row['particle_id'])
            for key in ('normal','normal_force_N','tangential_force_N','penetration_m',
                        'wall_velocity_m_s','conductance_W_K','wall_temperature_K'):
                np.testing.assert_array_equal(altered[key],row[key],err_msg=key)
        record('labels_and_reorder',identity_bound_contacts=18,mechanical_invariance=True,
               current_position_not_used_to_relabel_cohort=True)

    def test_06_second_publication_failure_retains_only_first_mock_accepted_stage(self):
        logs=ROOT/'logs/process_slow_20261003';logs.mkdir(exist_ok=True)
        out=Path(tempfile.mkdtemp(prefix='runner_mock_publication_',dir=logs))/'trajectory'
        source={'MOCK_SOURCE':'NOT_A_PHYSICAL_EXECUTION'};real_write=driver.write
        calls={'advance':0,'publication':0}
        def mock_advance(pa,state,p,boundary,dt):
            calls['advance']+=1
            candidate=pa.extract_particles(np.arange(27,dtype=np.uint32))
            new_state=copy.deepcopy(state);new_state['time_s']=state['time_s']+dt
            new_state['state_hash']='MOCK_CANDIDATE_'+str(calls['advance'])
            diagnostic=dict(evaluations=1,elapsed_s=.001,residuals={},
                maximum_sampled_coulomb_relative_defect=0.,normal_temporal_residual_J=0.,publication_checks={})
            return candidate,new_state,diagnostic
        def mock_integrals(cumulative,*args): cumulative['friction_raw']+=1.
        def mock_frame(pa,state,integrals):
            return dict(classification='MOCK_NOT_PHYSICAL',time_s=state['time_s'],
                        raw=integrals['friction_raw'].copy())
        def mock_save(directory,pa,state,p,diag,integrals):
            return dict(time_s=state['time_s'],checkpoint='MOCK_NOT_WRITTEN',
                        checkpoint_sha256='MOCK',state_hash=state['state_hash'])
        def fail_second_history(path,value):
            if Path(path).name=='accepted_history.json':
                calls['publication']+=1
                if calls['publication']==2:raise OSError('INJECTED_PUBLICATION_FAILURE')
            return real_write(path,value)
        with mock.patch.multiple(driver,advance=mock_advance,checks=lambda *args:None,
                add_integrals=mock_integrals,save_stage=mock_save,frame=mock_frame,
                inventory=lambda:source,write=fail_second_history):
            result=driver.trajectory(out,0,20e-6,source,maximum_steps=2)
        self.assertEqual(result['status'],'STOPPED_LOCAL_SCOUT_NOT_PROCESS_FEASIBLE')
        self.assertEqual(result['published_steps'],1)
        self.assertEqual(result['time_s'],250e-9)
        rejection=driver.read(out/'REJECTED.json')
        self.assertEqual(rejection['last_accepted_time_s'],250e-9)
        self.assertEqual(rejection['last_accepted_state_hash'],'MOCK_CANDIDATE_1')
        self.assertFalse(rejection['accepted']);self.assertFalse(rejection['restart_allowed'])
        self.assertEqual(rejection['candidate_interval_s'],[250e-9,500e-9])
        history=driver.read(out/'accepted_history.json')
        self.assertEqual(len(history),1);self.assertEqual(history[0]['time_s'],250e-9)
        frames=driver.read(out/'frames.json');self.assertEqual(len(frames),2)
        self.assertEqual(frames[-1]['raw'],[1.,1.])
        self.assertEqual(calls,dict(advance=2,publication=2))
        real_write(out.parent/'MOCK_SCOPE.json',dict(classification='MOCK_PUBLICATION_TEST_NOT_PHYSICS',
            trajectory_advanced=False,independent_mechanical_accuracy='NOT_RUN'))
        record('publication_rollback',artifacts=str(out.parent.relative_to(ROOT)),
               synthetic_candidates=2,published_mock_stages=1,last_accepted_time_s=250e-9,
               physical_trajectory_advanced=False)


if __name__=='__main__':unittest.main(verbosity=2)
