"""Short full-chain regression after local contact gates; no production run."""
import json
import tempfile
import time
import sys
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.config_adapter import load_settings
from afsd_backend.pysph.particles import create_geometry
from afsd_backend.pysph.solver import SolidSolver
from afsd_backend.pysph.output import checkpoint,restore
from test_mechanics_review import keyed_property
from review_artifacts import artifact_path

class EventSolverSmoke(unittest.TestCase):
    def setUp(self):
        _,self.p,_=load_settings(ROOT/'configs/debug/starved.json')

    def test_legacy_first_failure_and_accepted_snapshot_unchanged(self):
        solver=SolidSolver(create_geometry(self.p),self.p)
        for _ in range(31): solver.advance()
        reference,state=restore(ROOT/'outputs/coupling_20260930/phase1_verified/last_accepted.npz')
        for field in reference.properties:
            np.testing.assert_array_equal(getattr(solver.pa,field),getattr(reference,field),err_msg=field)
        self.assertEqual(solver.time,state['t'])
        self.assertEqual(solver.ledger,state['ledger'])
        with self.assertRaisesRegex(ValueError,'negative midpoint friction'):
            solver.advance()

    def test_event_variant_label_sort_restart_full_chain(self):
        pa,pb,pc=(create_geometry(self.p) for _ in range(3))
        pb.tracer_type[:]=np.where(pb.tracer_type==1,2,1);pb.wholebar_tag[:]=91
        pc.remove_property('tracer_type');pc.remove_property('wholebar_tag')
        order=np.random.default_rng(828).permutation(pa.get_number_of_particles()).astype(np.uint32)
        pd=pa.extract_particles(order)
        a,b,c,d=(SolidSolver(p,self.p,contact_variant='event_aware_v1') for p in (pa,pb,pc,pd))
        start=time.perf_counter()
        a.advance(); seconds=time.perf_counter()-start
        for solver in (b,c,d):solver.advance()
        for field in ('x','y','z','u','v','w','rho','S','temperature','ep'):
            np.testing.assert_array_equal(getattr(pa,field),getattr(pb,field))
            np.testing.assert_array_equal(getattr(pa,field),getattr(pc,field))
            np.testing.assert_allclose(keyed_property(pa,field),keyed_property(pd,field),
                                       rtol=1e-11,atol=1e-12)
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'accepted.npz';checkpoint(path,pa,a.state())
            pe,state=restore(path)
            e=SolidSolver(pe,self.p,state,contact_variant='event_aware_v1')
            a.advance();e.advance()
            for field in pa.properties:
                np.testing.assert_array_equal(getattr(pa,field),getattr(pe,field))
            self.assertEqual(a.ledger,e.ledger)
        self.assertEqual(a.diagnostics()['alarms'],[])
        report=dict(status='PASS',N=pa.get_number_of_particles(),first_step_wall_s=seconds,
            initial_dt_s=a.last['dt_s'],estimated_steps_to_2us=int(np.ceil(2e-6/a.last['dt_s'])),
            linear_estimated_2us_wall_s=seconds*np.ceil(2e-6/a.last['dt_s']),
            warning='timing initial state only; later stick/coupled events may increase cost')
        artifact_path('preflight_timing.json').write_text(json.dumps(report,indent=2)+'\n')
        print('PREFLIGHT',json.dumps(report))

    def test_checkpoint_cannot_silently_switch_contact_model(self):
        solver=SolidSolver(create_geometry(self.p),self.p)
        with self.assertRaisesRegex(ValueError,'variant mismatch'):
            SolidSolver(solver.pa,self.p,solver.state(),contact_variant='event_aware_v1')

if __name__=='__main__':unittest.main(verbosity=2)

