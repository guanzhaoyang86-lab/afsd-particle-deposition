"""Independent scalar reversal and actual-updater closed elastic loop.
--legacy-red intentionally FAILS; --cycle emits convergence evidence.
"""
import argparse
import json
import math
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.contact import WallContact
from afsd_backend.pysph.material import update
from afsd_backend.pysph.mechanics import bulk_modulus, pressure

# Hand-integrated piecewise Newton/Coulomb equations, not values from the UUT.
REVERSAL_REFERENCE = dict(t_stop=1/3, u_final=-2/3, Q_friction=7/18,
                          W_other=1/9, delta_K=-5/18)

class ElasticOnly:
    """Disable plasticity only in this prescribed material-point experiment."""
    def yield_stress(self, ep, rate, temperature):
        return np.full_like(np.asarray(ep, float), np.inf)

def legacy_reversal():
    dt, mass = 1., 1.
    force, _, _, _ = WallContact(10., 1.).force(0., [0.,0.,1.], [1.,0.,0.], .1, mass, dt)
    u_final = 1.+dt*(-2.+force[0])/mass
    midpoint = (1.+u_final)/2
    return dict(u_final=u_final, Q_friction=-dt*force[0]*midpoint,
                W_other=-2.*dt*midpoint, delta_K=.5*mass*(u_final*u_final-1.))

def elastic_cycle(n):
    G, rho0, mass, amplitude, Jc = 26.3e9, 2700., 1., .001, .99
    K = bulk_modulus(G, .33); V0 = mass/rho0
    A = np.diag([1., -.5, -.5])
    endpoints = [np.zeros(3), amplitude*np.diag(A),
                 amplitude*np.diag(A)+math.log(Jc)/3,
                 np.full(3, math.log(Jc)/3), np.zeros(3)]
    s = np.zeros((1,3,3)); ep = np.zeros(1); F = np.eye(3)
    work = 0.; deviatoric_work = 0.; hydrostatic_work = 0.; es = []
    max_stress = 0.; records = []
    dt = 1./n
    for segment in range(4):
        lo, hi = endpoints[segment:segment+2]
        L = np.diag(hi-lo)[None]
        for k in range(n):
            log0 = lo+(hi-lo)*k/n
            log1 = lo+(hi-lo)*(k+1)/n
            logmid = (log0+log1)/2
            Jmid = math.exp(float(logmid.sum()))
            rho_mid = rho0/Jmid
            old = s.copy()
            s, ep, heat, numerical, _ = update(s, ep, L, np.array([293.15]), dt, G, ElasticOnly())
            if np.any(ep != 0) or np.any(heat != 0) or np.any(numerical != 0):
                raise AssertionError('elastic-only material experiment produced plasticity')
            smid = .5*(old+s)
            p = float(pressure(np.array([rho_mid]), rho0, K)[0])
            dw_s = float(np.sum(smid*L))*V0*Jmid*dt
            dw_p = -p*float(np.trace(L[0]))*V0*Jmid*dt
            deviatoric_work += dw_s; hydrostatic_work += dw_p; work += dw_s+dw_p
            F = np.diag(np.exp((hi-lo)*dt))@F
            Jend = math.exp(float(log1.sum()))
            es_now = V0*Jend*float(np.sum(s*s))/(4*G)
            es.append(es_now)
            max_stress = max(max_stress, float(np.linalg.norm(s)))
            records.append(dict(segment=segment+1, substep=k+1, J=Jend,
                integral_sigma_D_V_dt_J=work, Es_diagnostic_J=es_now,
                stress_Pa=s[0].tolist(), F=F.tolist()))
    # Independent analytic integral: loading at V0 minus unloading at Jc V0.
    reference = (1-Jc)*G*1.5*amplitude**2*V0
    return dict(substeps_per_segment=n, dt_s=dt, duration_s=4., mass_kg=mass,
        Jc=Jc, G_Pa=G, shear_amplitude=amplitude, work_J=work,
        deviatoric_work_J=deviatoric_work, hydrostatic_work_J=hydrostatic_work,
        independent_work_reference_J=reference, Es_final_J=es[-1],
        Es_peak_J=max(es), final_stress_Pa=s[0].tolist(), final_F=F.tolist(),
        final_stress_relative=float(np.linalg.norm(s))/max_stress,
        final_F_error=float(np.linalg.norm(F-np.eye(3))),
        internal_energy_variable=None,
        internal_energy_variable_status='NOT_SOLVED: pa.e is upstream scratch; updater has no internal energy variable',
        plasticity='DISABLED_BY_TEST_ONLY_INFINITE_YIELD_LAW',
        thermal='OFF', contact='OFF', records=records)

class LegacyRed(unittest.TestCase):
    def test_original_friction_must_resolve_true_reversal(self):
        measured = legacy_reversal()
        print('LEGACY_RED', json.dumps(measured), 'INDEPENDENT_REFERENCE', REVERSAL_REFERENCE)
        for key, value in measured.items():
            self.assertAlmostEqual(value, REVERSAL_REFERENCE[key], delta=1e-10, msg=key)

class ElasticLoop(unittest.TestCase):
    def test_actual_updater_nonzero_loop_is_not_an_Es_only_error(self):
        for n in (16, 32, 64, 128):
            result = elastic_cycle(n)
            self.assertLess(result['final_F_error'], 1e-11)
            self.assertLess(result['final_stress_relative'], 1e-11)
            self.assertAlmostEqual(result['work_J'], (1-.99)*26.3e9*1.5*.001**2/2700, delta=1e-9)
            self.assertGreater(result['work_J'], .1)
            print('CLOSED_CYCLE', json.dumps({k:v for k,v in result.items() if k!='records'}))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--legacy-red', action='store_true')
    parser.add_argument('--cycle', action='store_true')
    parser.add_argument('--output')
    args = parser.parse_args()
    if args.output:
        path = Path(args.output)
        if path.exists(): raise FileExistsError(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps([elastic_cycle(n) for n in (16,32,64,128)], indent=2)+'\n')
    selected = LegacyRed if args.legacy_red else ElasticLoop
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(selected))
    raise SystemExit(0 if result.wasSuccessful() else 1)

