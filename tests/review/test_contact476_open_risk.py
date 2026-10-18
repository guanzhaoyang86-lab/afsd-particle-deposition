"""Retain the actual E3 near-stop cone failure; no new physical tolerance."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.contact_cone import solve, ConeCertificateError


class Contact476OpenRisk(unittest.TestCase):
    def test_frozen_near_stop_problem_must_return_a_certified_force(self):
        path = ROOT/'outputs/engineering_e3_20261001/dx05_cfl01/REJECTED_contact_debug.json'
        data = json.loads(path.read_text(encoding='utf-8'))['failed_local_context']['failed_cone']
        base = np.asarray(data['base_tangent'])
        basis = np.asarray(data['tangent_bases'])
        mobility = np.asarray(data['mobility'])
        caps = np.asarray(data['limits_N'])
        coefficient = data['coefficient']
        # This is the existing caller's 2e-12 cone tolerance, NOT a relaxed
        # experiment/force-accuracy requirement chosen to match the candidate.
        try:
            forces = solve(base, basis, mobility, coefficient, caps, tolerance=2e-12)
        except ConeCertificateError as exc:
            self.fail('frozen real candidate476 cannot advance; '+str(exc))
        tangent = np.einsum('nki,nk->ni', basis, forces)
        scale_f = float(caps.max())
        scale_v = max(float(np.linalg.norm(base)), coefficient*float(np.linalg.norm(mobility, 2))*scale_f)
        gradient = base+coefficient*np.einsum('nki,k->ni', basis, mobility@forces.sum(axis=0))
        gap = float(np.sum(gradient*tangent)+np.sum(caps*np.linalg.norm(gradient, axis=1)))
        self.assertLessEqual(np.max(np.linalg.norm(tangent, axis=1)-caps)/scale_f, 2e-12)
        self.assertLessEqual(abs(gap)/(scale_f*scale_v), 2e-12)
        # No exact individual force split is forced by this caller-level gate.
        # Independent high-precision comparison is a separate diagnostic.


if __name__ == '__main__':
    unittest.main(verbosity=2)
