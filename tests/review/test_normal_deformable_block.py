"""Short independent gates for free-normal SPH; no full release trajectory.

Rest-state references are literal impulse formulas, not normal/contact calls.
Each instance advances at most five steps; full 120/240-step runs are separate.
"""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import probe_normal_deformable_block as c
from afsd_backend.pysph.mechanics import positions, velocities, corrected_operators
from afsd_backend.pysph.output import restore, checkpoint
from afsd_scaffold.vtk_legacy import read_point_cloud

MASS = 3.375e-7
KN = 657500.
DT = 5e-9


def rest_reference(x):
    top = np.isclose(x[:,2],.00175,rtol=0.,atol=1e-17)
    wall = np.column_stack((-20*x[top,1],20*x[top,0],np.zeros(16)))
    jt = .2*DT*wall/np.linalg.norm(wall,axis=1)[:,None]
    impulse = jt.copy();impulse[:,2] = -DT
    v = np.zeros((64,3));v[top] = impulse/MASS
    return dict(top=top,wall=wall,jt=jt,impulse=impulse,velocity=v,position=x+.5*DT*v,
        tangential_kinetic=float(np.sum(jt*jt)/(2*MASS)),
        kinetic=float(np.sum(impulse*impulse)/(2*MASS)),wall_work=float(np.sum(jt*wall)))


class NormalDeformableShort(unittest.TestCase):
    def test_01_unscaled_cube_and_fortyeight_free_particles(self):
        b=c.Block(DT);x=positions(b.pa)
        self.assertEqual(len(x),64)
        self.assertEqual(np.linalg.matrix_rank(x-x.mean(axis=0)),3)
        for axis in range(3):self.assertEqual(len(np.unique(x[:,axis])),4)
        np.testing.assert_allclose(b.pa.m,MASS,rtol=0.,atol=1e-21)
        self.assertAlmostEqual(b.pa.m.sum(),2.16e-5,delta=1e-18)
        self.assertEqual(np.count_nonzero(b.pa.boundary_role!=1),48)
        self.assertEqual(np.count_nonzero(b.pa.boundary_role==1),16)
        b.advance()
        for trace in b.last_traces:
            role=b.pa.boundary_role[trace['particle_index']]
            self.assertEqual(trace['free_axes'],[0.,0.,0.] if role==1 else [1.,1.,1.])

    def test_02_literal_initial_overlap_force_potential_and_source(self):
        b=c.Block(DT);normal=b.record()['normal_contact']
        np.testing.assert_allclose(normal['signed_penetration_m'],1/KN,rtol=0.,atol=1e-18)
        np.testing.assert_allclose(normal['normal_force_N'],1.,rtol=0.,atol=1e-12)
        np.testing.assert_allclose(normal['potential_J'],1/(2*KN),rtol=0.,atol=1e-18)
        self.assertAlmostEqual(b.initial_normal_potential_J,16/(2*KN),delta=1e-17)
        self.assertEqual(np.count_nonzero(normal['active']),16)
        cfg=json.loads((ROOT/'configs/benchmarks/normal_deformable_block.json').read_text(encoding='utf-8'))
        self.assertFalse(cfg['production_authorized'])
        for item in cfg['parameters'].values():
            self.assertTrue(all(item.get(k) is not None for k in ('value','unit','source','status')))

    def test_03_first_kick_independent_free_normal_solution(self):
        b=c.Block(DT);expected=rest_reference(positions(b.pa).copy());b.advance()
        np.testing.assert_allclose(velocities(b.pa),expected['velocity'],rtol=0.,atol=1e-12)
        np.testing.assert_allclose(positions(b.pa),expected['position'],rtol=0.,atol=2e-17)
        np.testing.assert_array_equal(b.pa.S,np.zeros(64*9))
        np.testing.assert_array_equal(b.pa.ep,np.zeros(64))
        np.testing.assert_array_equal(b.pa.rho,np.full(64,2700.))
        self.assertAlmostEqual(b.record()['kinetic_J'],expected['kinetic'],delta=1e-18)
        self.assertTrue(np.all(b.pa.w[expected['top']]<0.))

    def test_04_top_has_no_clamp_impulse_and_friction_does_not_absorb_normal_work(self):
        b=c.Block(DT);reference=rest_reference(positions(b.pa).copy());b.advance()
        topids=set(b.pa.particle_id[reference['top']]);q=0.;wall=0.;normalwork=0.
        for trace in b.last_traces:
            if trace['particle_id'] not in topids:continue
            result=trace['result']
            np.testing.assert_allclose(result['normal_impulse'],[[0.,0.,-DT]],rtol=0.,atol=1e-20)
            np.testing.assert_allclose(result['constraint_impulse'],[0.,0.,0.],rtol=0.,atol=1e-17)
            q+=float(np.sum(result['friction_heat_J']));wall+=result['wall_work_J']
            normalwork+=sum(float(np.sum(a['normal_work_J'])) for a in result['intervals'])
        self.assertAlmostEqual(q,reference['wall_work']-reference['tangential_kinetic'],delta=1e-18)
        self.assertAlmostEqual(wall,reference['wall_work'],delta=1e-18)
        self.assertAlmostEqual(normalwork,16*DT*DT/(2*MASS),delta=1e-18)
        self.assertGreater(normalwork,0.)

    def test_05_gap_and_normal_force_evolve_not_prescribed_one_newton(self):
        b=c.Block(DT);before=b.record()['normal_contact'];b.advance();after=b.record()['normal_contact']
        expected_delta=1/KN-DT*DT/(2*MASS)
        np.testing.assert_allclose(after['signed_penetration_m'],expected_delta,rtol=0.,atol=2e-18)
        np.testing.assert_allclose(after['normal_force_N'],KN*expected_delta,rtol=0.,atol=2e-12)
        np.testing.assert_allclose(after['potential_J'],.5*KN*expected_delta**2,rtol=0.,atol=3e-18)
        self.assertTrue(np.all(after['normal_force_N']<before['normal_force_N']))
        self.assertTrue(np.all(after['gap_m']>before['gap_m']))

    def test_06_normal_potential_work_residual_is_raw_not_added_to_heat(self):
        b=c.Block(DT);b.advance();r=b.record();j=b.last_journal;cumulative=b.cumulative
        du=float(np.sum(j['normal_after']['potential_J'])-np.sum(j['normal_before']['potential_J']))
        self.assertAlmostEqual(cumulative['normal_temporal_residual_J'],
            du+cumulative['normal_work_J']-cumulative['normal_wall_work_J'],delta=1e-18)
        self.assertEqual(cumulative['normal_wall_work_J'],0.)
        identity=(r['kinetic_J']+cumulative['friction_raw_J']-cumulative['wall_work_J']-
            cumulative['constraint_work_J']-cumulative['internal_work_J']-cumulative['normal_relative_work_J'])
        self.assertAlmostEqual(identity,0.,delta=1e-18)
        expected=rest_reference(np.column_stack((b.pa.x0,b.pa.y0,b.pa.z0)))
        self.assertAlmostEqual(cumulative['friction_raw_J'],expected['wall_work']-expected['tangential_kinetic'],delta=1e-18)

    def test_07_three_steps_transmit_compression_and_shear_to_free_interior(self):
        b=c.Block(DT)
        for _ in range(3):b.advance()
        inner=b.pa.boundary_role==0;top=b.pa.boundary_role==2;bottom=b.pa.boundary_role==1
        self.assertGreater(float(np.linalg.norm(velocities(b.pa)[inner])),1e-10)
        self.assertGreater(float(np.max(np.abs(b.pa.w[inner]))),1e-10)
        self.assertTrue(np.all(b.pa.w[top]<0.))
        np.testing.assert_array_equal(velocities(b.pa)[bottom],np.zeros((16,3)))
        s=b.pa.S.reshape(-1,3,3)
        self.assertGreater(float(np.max(np.abs(s[:,2,2]))),1.)
        self.assertGreater(float(np.max(np.abs(s[:,0,2]))),1.)
        self.assertGreater(float(np.max(np.abs(s[:,1,2]))),1.)
        self.assertGreater(float(np.max(np.abs(b.pa.rho-2700.))),1e-8)

    def test_08_labels_swap_and_remove_cannot_change_dynamics_or_ledger(self):
        a,b,d=[c.Block(DT) for _ in range(3)]
        b.pa.tracer_type[:]=99-b.pa.tracer_type;b.pa.wholebar_tag[:]=17;b.pa.as_rs_tag[:]=-5
        d.pa.tracer_type[:]=0;d.pa.wholebar_tag[:]=0;d.pa.as_rs_tag[:]=0
        # The top-layer marker is no longer a constraint or force selector;
        # only bottom role 1 remains a physical boundary definition.
        d.pa.boundary_role[d.pa.boundary_role==2]=0
        for _ in range(3):a.advance();b.advance();d.advance()
        for other in (b,d):
            for key in ('x','y','z','u','v','w','rho','S','ep','plastic_work'):
                np.testing.assert_array_equal(getattr(a.pa,key),getattr(other.pa,key))
            self.assertEqual(a.cumulative,other.cumulative)
            self.assertEqual(a.record()['normal_potential_J'],other.record()['normal_potential_J'])

    def test_09_actual_nonzero_checkpoint_and_two_step_continuation(self):
        a=c.Block(DT);a.advance();a.advance()
        self.assertGreater(np.linalg.norm(a.pa.S),1.)
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'accepted.npz';a.checkpoint(path);b=c.Block.from_checkpoint(path)
            for _ in range(2):a.advance();b.advance()
        for key in a.pa.properties:np.testing.assert_array_equal(getattr(a.pa,key),getattr(b.pa,key))
        self.assertEqual(a.cumulative,b.cumulative);self.assertEqual(a.time,b.time)
        self.assertEqual(a.record()['normal_potential_J'],b.record()['normal_potential_J'])
        self.assertEqual(a.initial_normal_potential_J,b.initial_normal_potential_J)

    def test_10_nonzero_particle_permutation_preserves_id_and_normal_state(self):
        a=c.Block(DT);a.advance();a.advance()
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'accepted.npz';a.checkpoint(path);pa,state=restore(path)
            permutation=np.random.default_rng(73).permutation(64).astype(np.int32)
            b=c.Block(DT,pa=pa.extract_particles(permutation),state=state)
            for _ in range(2):a.advance();b.advance()
        order=np.argsort(b.pa.particle_id)
        for key in ('particle_id','boundary_role','tracer_type','m'):
            np.testing.assert_array_equal(getattr(a.pa,key),getattr(b.pa,key)[order])
        for key in ('x','y','z','u','v','w','rho','ep'):
            np.testing.assert_allclose(getattr(a.pa,key),getattr(b.pa,key)[order],rtol=0.,atol=1e-12)
        np.testing.assert_allclose(a.pa.S.reshape(64,9),b.pa.S.reshape(64,9)[order],rtol=0.,atol=1e-6)
        for key in a.cumulative:self.assertAlmostEqual(a.cumulative[key],b.cumulative[key],delta=1e-12)

    def test_11_acceptance_failure_rolls_back_material_and_normal_ledger(self):
        b=c.Block(DT);b.advance();b.advance()
        before={key:getattr(b.pa,key).copy() for key in b.pa.properties}
        ledger=copy.deepcopy(b.cumulative);old_time=b.time;old_step=b.step
        journal,traces=b.last_journal,b.last_traces;normal=b.record()['normal_potential_J']
        validate=b._validate_candidate
        def reject_after_normal_accounting(increment):
            validate(increment)
            raise ValueError('INJECTED_POST_NORMAL_LEDGER_REJECTION')
        with patch.object(b,'_validate_candidate',side_effect=reject_after_normal_accounting):
            with self.assertRaisesRegex(ValueError,'INJECTED_POST_NORMAL_LEDGER_REJECTION'):b.advance()
        self.assertFalse(np.array_equal(b.rejected_candidate['properties']['S'],before['S']))
        for key,value in before.items():np.testing.assert_array_equal(getattr(b.pa,key),value)
        self.assertEqual(b.cumulative,ledger);self.assertEqual(b.time,old_time);self.assertEqual(b.step,old_step)
        self.assertIs(b.last_journal,journal);self.assertIs(b.last_traces,traces)
        self.assertEqual(b.record()['normal_potential_J'],normal)
        self.assertEqual(b.rejected_candidate['classification'],'REJECTED')
        self.assertFalse(b.rejected_candidate['restart_allowed'])
        self.assertIn('normal_energy',b.rejected_candidate['journal'])

    def test_12_rejected_checkpoint_cannot_restart(self):
        b=c.Block(DT);b.advance()
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'accepted.npz';b.checkpoint(path);pa,state=restore(path)
            self.assertFalse(state['production_authorized'])
            state['classification']='REJECTED';state['restart_allowed']=False
            checkpoint(Path(td)/'rejected.npz',pa,state)
            with self.assertRaises(ValueError):c.Block.from_checkpoint(Path(td)/'rejected.npz')

    def test_13_acoustic_and_normal_stiffness_bounds_both_enforced(self):
        b=c.Block(DT)
        bulk=2*26.3e9*1.33/(3*.34)
        acoustic=.1*.0005*1.3/np.sqrt((bulk+4*26.3e9/3)/2700.)
        normal=.1*np.sqrt(MASS/KN)
        self.assertAlmostEqual(b.stable_dt(),min(acoustic,normal),delta=1e-20)
        self.assertLessEqual(DT,b.stable_dt())
        with self.assertRaises(ValueError):c.Block(1e-6)

    def test_14_vtk_reports_normal_contact_not_fabricated_heat(self):
        b=c.Block(DT);b.advance();b.advance()
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'frame.vtk';b.vtk(path);vtk=read_point_cloud(path)
        np.testing.assert_array_equal(vtk['points'],positions(b.pa))
        np.testing.assert_array_equal(vtk['scalars']['particle_id'],b.pa.particle_id)
        self.assertIn('density_kg_per_m3',vtk['scalars']);self.assertIn('stress_yz_Pa',vtk['scalars'])
        self.assertFalse(any('temperature' in key.lower() for key in vtk['scalars']))
        self.assertNotIn('internal_energy_J_per_kg',vtk['scalars'])
        for key in ('normal_force_N','normal_potential_J','contact_active','signed_gap_m'):
            self.assertIn(key,vtk['scalars'])
        self.assertTrue(np.all(np.asarray(vtk['scalars']['normal_force_N'])[b.pa.boundary_role!=2]==0.))

    def test_15_corrected_neighborhoods_remain_true_three_dimensional(self):
        b=c.Block(DT);b.advance();b.advance();x=positions(b.pa)
        i,j,gradient,L,drho,B=corrected_operators(b.pa)
        for k in range(64):
            neighbors=np.r_[j[i==k],i[j==k]]
            self.assertEqual(np.linalg.matrix_rank(x[neighbors]-x[k]),3)
        self.assertLess(float(np.max(np.linalg.cond(np.linalg.inv(B)))),1e6)
        self.assertGreater(float(np.max(np.abs(L[:,2,2]))),0.)


if __name__=='__main__':unittest.main(verbosity=2)
