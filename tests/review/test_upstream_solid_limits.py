"""Independent RED probes of the pinned upstream solid-mechanics equations.

These are adoption guards, not expected-failure decorators and not AFSD tests.
The upstream implementation is called directly without copying or patching it.
With the currently pinned PySPH source, both physical assertions must FAIL.

Run with the isolated interpreter:
    python -B tests/review/test_upstream_solid_limits.py
    python -B tests/review/test_upstream_solid_limits.py --plastic-red
    python -B tests/review/test_upstream_solid_limits.py --energy-red

Tolerance preregistration belongs to the engineering-approximation contract:
relative J2 yield error <= 1e-10; stationary power <= 1e-12 J/(kg s).
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
UPSTREAM = ROOT / "third_party" / "PySPH"
sys.path.insert(0, str(UPSTREAM))

from pysph.sph.solid_mech import basic, hvi


def source_evidence(module):
    path = Path(module.__file__).resolve()
    expected = (UPSTREAM / "pysph" / "sph" / "solid_mech" / path.name).resolve()
    if path != expected:
        raise AssertionError(f"Not the pinned upstream source: {path}")
    return {
        "source_path": path.relative_to(ROOT).as_posix(),
        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


class UpstreamSolidLimits(unittest.TestCase):
    def test_pure_xz_shear_must_return_to_three_dimensional_yield(self):
        # Independent 3D J2 invariant: q = sqrt(3) * abs(tau) for pure shear.
        # Input q = 2Y, so perfect plastic projection must return q = Y.
        yield_pa = 90.0e6
        tau_pa = 2.0 * yield_pa / math.sqrt(3.0)
        stress = {name: np.zeros(1) for name in
                  ("s00", "s01", "s02", "s11", "s12", "s22")}
        stress["s02"][0] = tau_pa
        equation = hvi.VonMisesPlasticity2D(
            dest="material", sources=None, flow_stress=yield_pa
        )
        equation.loop(
            d_idx=0, d_s00=stress["s00"], d_s01=stress["s01"],
            d_s02=stress["s02"], d_s11=stress["s11"],
            d_s12=stress["s12"], d_s22=stress["s22"],
        )
        squared_norm = sum(float(stress[k][0]) ** 2
                           for k in ("s00", "s11", "s22"))
        squared_norm += 2.0 * sum(float(stress[k][0]) ** 2
                                 for k in ("s01", "s02", "s12"))
        actual_q = math.sqrt(1.5 * squared_norm)
        relative_error = abs(actual_q / yield_pa - 1.0)
        print("UPSTREAM_PLASTIC_RED " + json.dumps({
            **source_evidence(hvi),
            "classification": "UPSTREAM_ADOPTION_GUARD_NOT_AFSD",
            "q_trial_Pa": 2.0 * yield_pa,
            "expected_q_Pa": yield_pa,
            "actual_q_Pa": actual_q,
            "relative_error": relative_error,
            "tolerance_relative": 1.0e-10,
            "reference": "three-dimensional J2 perfect-plastic yield surface",
        }), flush=True)
        self.assertLessEqual(
            relative_error, 1.0e-10,
            "Pure xz shear remains overyield: this upstream limiter is not "
            "an admissible three-dimensional plasticity replacement.",
        )

    def test_stationary_pressure_pair_must_have_zero_power(self):
        # Pressure work has a velocity factor. With both particles stationary,
        # pair pressure work and artificial-viscosity work are exactly zero.
        # Values are a dimensional manufactured fixture, not AFSD parameters.
        one = np.array([1.0])
        pressure = np.array([2.0])
        energy_rate = np.zeros(1)
        equation = basic.EnergyEquationWithStress(
            dest="material", sources=["material"],
            alpha=0.0, beta=0.0, eta=0.01,
        )
        equation.initialize(d_idx=0, d_ae=energy_rate)
        equation.loop(
            d_idx=0, s_idx=0, s_m=one, d_rho=one, s_rho=one,
            d_p=pressure, s_p=pressure, d_cs=one, s_cs=one,
            d_ae=energy_rate,
            XIJ=np.array([1.0, 0.0, 0.0]),
            VIJ=np.zeros(3),
            DWIJ=np.array([-1.0, 0.0, 0.0]),
            HIJ=1.0, R2IJ=1.0, RHOIJ1=1.0,
        )
        # Zero stress and zero velocity gradient make the separate deviatoric
        # post_loop term zero; no reference value is obtained from that method.
        actual = float(energy_rate[0])
        print("UPSTREAM_ENERGY_RED " + json.dumps({
            **source_evidence(basic),
            "classification": "UPSTREAM_ADOPTION_GUARD_NOT_AFSD",
            "mass_kg": 1.0,
            "density_kg_m3": 1.0,
            "pressure_Pa": 2.0,
            "relative_velocity_m_s": [0.0, 0.0, 0.0],
            "expected_specific_power_J_kg_s": 0.0,
            "actual_d_ae_raw": actual,
            "actual_units_warning": "upstream output lacks its velocity-gradient "
                                    "factor; do not interpret as valid power",
            "tolerance_absolute_J_kg_s": 1.0e-12,
            "reference": "zero mechanical power at zero velocity",
        }), flush=True)
        self.assertAlmostEqual(
            actual, 0.0, delta=1.0e-12,
            msg="A stationary positive-pressure pair must not create energy.",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--plastic-red", action="store_true")
    selection.add_argument("--energy-red", action="store_true")
    args = parser.parse_args()
    if args.plastic_red:
        suite = unittest.TestSuite([UpstreamSolidLimits(
            "test_pure_xz_shear_must_return_to_three_dimensional_yield")])
    elif args.energy_red:
        suite = unittest.TestSuite([UpstreamSolidLimits(
            "test_stationary_pressure_pair_must_have_zero_power")])
    else:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(
            UpstreamSolidLimits)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(0 if result.wasSuccessful() else 1)

