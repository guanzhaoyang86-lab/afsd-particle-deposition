"""Retain the new observed coupled-kick failure as RED; no production rerun."""
import json
import sys
import unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from afsd_backend.pysph.contact_kick import kick, ContactConvergenceError
from afsd_backend.pysph.kick_diagnostics import clean
from review_artifacts import artifact_path

class NewCandidateRed(unittest.TestCase):
    def test_observed_454_single_face_stop_must_converge(self):
        source=ROOT/'outputs/coupling_20260930/event_dx05_cfl01/REJECTED_contact_debug.json'
        c=json.loads(source.read_text())['failed_local_context']
        face=c['faces'][0]
        normal=np.array(face['normal']);wall=np.array(face['wall_velocity'])
        tangent=np.eye(3)-np.outer(normal,normal)
        v=np.array(c['velocity_before']);other=np.array(c['internal_force'])
        fn=2.63e12*.0005**2*face['penetration']
        capacity=.3*fn;mass=c['mass_kg']
        u=tangent@(v-wall);drive=tangent@other
        speed=np.linalg.norm(u);f=np.linalg.norm(drive)
        b=(1+float(u@drive)/(speed*f))/2;d=1-b
        tstop=mass*speed*(b/(capacity-f)+d/(capacity+f))
        heat=mass*capacity*speed**2*(b*b/(2*(capacity-f))+b*d/capacity+d*d/(2*(capacity+f)))
        self.assertLess(f,capacity);self.assertLess(tstop,c['dt_s'])
        report=dict(reference='independent constant-force vector Coulomb polar integral',
            source=str(source.relative_to(ROOT)),particle_id=454,t_stop_s=tstop,
            friction_heat_J=heat,static_friction_capacity_N=capacity,tangential_other_force_N=f)
        try:
            result=kick(v,mass,other,[normal],[fn],[wall],.3,c['dt_s'],
                free_axes=c['free_axes'],prescribed_velocity=c['prescribed_velocity'])
        except ContactConvergenceError as error:
            report.update(status='FAIL',error=str(error),failed_cone=getattr(error,'cone_context',None))
            artifact_path('candidate454.json').write_text(json.dumps(clean(report),indent=2)+'\n')
            self.fail(str(error)+'; independent stopping reference is feasible')
        self.assertLess(np.linalg.norm(tangent@(result['velocity']-wall)),1e-10)
        self.assertLess(abs(float(result['friction_heat_J'].sum())/heat-1),1e-3)
        self.assertTrue(result['events'],'stop event must be identified, not only nonnegative heat')
        actual_stop=result['events'][0]['time_s']
        self.assertLess(abs(actual_stop/tstop-1),1e-3)
        report.update(status='PASS',actual_stop_s=actual_stop,
            actual_friction_heat_J=float(result['friction_heat_J'].sum()),
            final_tangent_speed_m_s=float(np.linalg.norm(tangent@(result['velocity']-wall))),
            intervals=len(result['intervals']),events=result['events'],
            kinetic_work_residual_J=result['kinetic_work_residual_J'])
        artifact_path('candidate454.json').write_text(json.dumps(clean(report),indent=2)+'\n')

if __name__=='__main__':unittest.main(verbosity=2)

