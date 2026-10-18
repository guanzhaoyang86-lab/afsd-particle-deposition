"""Acceptance guards: exact rejection semantics, not extra trajectory accuracy."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
import probe_hht_coupled as driver


class PublicationGuards(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.old,cls.p=driver.setup()
        cls.state=driver.initialize(cls.old,cls.p,driver.boundary_for(cls.p,False))
        cls.candidate,cls.next_state,cls.diag=driver.advance(cls.old,cls.state,cls.p,driver.boundary_for(cls.p,False),1e-9)

    def test_saved_temperature_not_only_mapped_temperature(self):
        changed=self.candidate.extract_particles(np.arange(27,dtype=np.uint32))
        changed.temperature[13]+=1e-8
        # Mapped BE diagnostics are unchanged and fine. The actual saved node
        # has an extra m*cp*1e-8 J, independently well above1e-12 J.
        self.assertGreater(changed.m[13]*self.p['cp']*1e-8,1e-12)
        with self.assertRaisesRegex(ValueError,'saved temperature thermal identity'):
            driver.physical_checks(changed,self.old,self.old,self.next_state,copy.deepcopy(self.diag),self.p)

    def test_signed_slip_margin_includes_quadratic_interior_vertex(self):
        diagnostic=copy.deepcopy(self.diag)
        # Expected signed slip .1-.36*s+.36*s^2 is .1 at both endpoints,
        # but .01 at s=.5; the registered .05 margin must reject it.
        diagnostic['contact']['contacts'][0]['slip_coefficients_m_s']=np.array([[0.,0.,-.1],[0.,0.,.36],[0.,0.,-.36]])
        with self.assertRaisesRegex(ValueError,'interior registered sliding branch margin'):
            driver.physical_checks(self.candidate,self.old,self.old,self.next_state,diagnostic,self.p)

    def test_unpublished_outer_guard_failure_retains_previous_time(self):
        with tempfile.TemporaryDirectory(prefix='hht-publication-') as temporary:
            destination=Path(temporary)/'trial'
            with patch.object(driver,'advance',return_value=(self.candidate,self.next_state,self.diag)), \
                 patch.object(driver,'physical_checks',side_effect=ValueError('INJECTED_OUTER_GUARD')):
                with self.assertRaisesRegex(ValueError,'INJECTED_OUTER_GUARD'):
                    driver.trajectory(destination,8,1)
            rejected=json.loads((destination/'REJECTED.json').read_text())
            self.assertEqual(rejected['last_accepted_time_s'],0.)
            self.assertEqual(rejected['published_steps'],0)
            self.assertFalse(rejected['restart_allowed'])
            self.assertFalse((destination/'accepted_history.json').exists())
            self.assertEqual(list(destination.glob('stage_*')),[])


if __name__=='__main__':unittest.main(verbosity=2)
