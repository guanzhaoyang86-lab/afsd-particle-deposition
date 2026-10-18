"""New two-face blocker: independent closed KKT solution, retained as FAIL."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.contact_kick import cone_solve,ContactConvergenceError
from afsd_backend.pysph.kick_diagnostics import clean
from review_artifacts import artifact_path

class Contact439OpenRisk(unittest.TestCase):
    def test_two_face_frozen_cone_must_reach_known_kkt_solution(self):
        path=ROOT/'outputs/contact454_20260930/fixed_dx05_cfl01/REJECTED_contact_debug.json'
        data=json.loads(path.read_text())['failed_local_context']['failed_cone']
        base=np.array(data['base_tangent']);T=np.array(data['tangent_bases'])
        M=np.array(data['mobility']);limit=np.array(data['limits_N'])
        coefficient=data['coefficient'];lam=coefficient*M[0,0]
        np.testing.assert_allclose(M,np.eye(3)*M[0,0],rtol=1e-14,atol=0.)
        self.assertLess(abs(float(T[0,:,0]@T[1,:,0])),1e-14)
        np.testing.assert_allclose(T[0,:,1],T[1,:,1],atol=1e-14)
        self.assertAlmostEqual(base[0,1],base[1,1],delta=1e-24)
        # Independent solution: first disk interior, second disk at (-limit,0).
        # Shared z force is carried entirely by first disk; no UUT solve used.
        reference=np.array([[-base[0,0]/lam,-base[0,1]/lam],[-limit[1],0.]])
        self.assertLess(np.linalg.norm(reference[0]),limit[0])
        A=np.column_stack([T[0,:,0],T[0,:,1],T[1,:,0],T[1,:,1]])
        gradient=base.ravel()+coefficient*(A.T@M@A)@reference.ravel()
        np.testing.assert_allclose(gradient[:2],0.,atol=1e-23)
        self.assertAlmostEqual(gradient[3],0.,delta=1e-23)
        self.assertGreater(gradient[2],0.) # outward normal multiplier positive
        expected=np.einsum('nki,ni->nk',T,reference)
        report=dict(reference='closed four-variable convex KKT solution',
            particle_id=439,reference_tangent_force_N=reference,
            reference_cartesian_force_N=expected,kkt_gradient=gradient,
            first_disk_margin_N=float(limit[0]-np.linalg.norm(reference[0])),
            second_disk_multiplier=float(gradient[2]/limit[1]))
        try:
            actual=cone_solve(base,T,M,coefficient,limit)
        except ContactConvergenceError as error:
            report.update(status='FAIL',error=str(error),
                          failed_certificate=getattr(error,'cone_context',None))
            artifact_path('contact439_open_risk.json').write_text(json.dumps(clean(report),indent=2)+'\n')
            self.fail(str(error)+'; feasible closed KKT reference exists')
        np.testing.assert_allclose(actual,expected,atol=1e-11,rtol=1e-10)
        report.update(status='PASS',actual_cartesian_force_N=actual,
                      force_error_N=float(np.linalg.norm(actual-expected)))
        artifact_path('contact439_open_risk.json').write_text(json.dumps(clean(report),indent=2)+'\n')

if __name__=='__main__':unittest.main(verbosity=2)

