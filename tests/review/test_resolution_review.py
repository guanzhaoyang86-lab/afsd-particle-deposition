import copy
from pathlib import Path
import unittest
import numpy as np
from afsd_backend.pysph.config_adapter import load_settings
from afsd_backend.pysph.particles import create_geometry
from afsd_backend.pysph.solver import SolidSolver

ROOT = Path(__file__).resolve().parents[2]


class ResolutionBoundaryReview(unittest.TestCase):
    def test_wall_heat_reservoir_limits_actual_dt(self):
        from afsd_backend.pysph.particles import make_array
        _, p, _ = load_settings(ROOT/'configs/debug/starved.json')
        p['contact_conductance'] = 1e20  # numerical stiffness probe, NOT a calibrated case
        x = [[0., 0., -p['substrate_depth'] + .4*p['dx']]]
        pa = make_array(x, np.zeros(1, int), p)
        solver = SolidSolver(pa, p)
        oracle = .25*(p['rho0']*p['dx']**3)*p['cp']/(p['contact_conductance']*p['dx']**2)
        self.assertLessEqual(solver.stable_dt(), oracle*(1+1e-12))

    def test_same_physical_actuator_and_area_stiffness(self):
        _, p, _ = load_settings(ROOT/'configs/debug/starved.json')
        expected_volume = .0025**2*.0005
        expected_stiffness = .05*26.3e9/.0005
        for dx in (.0005, .00025):
            q = copy.deepcopy(p); q['dx'] = dx
            pa = create_geometry(q)
            solver = SolidSolver(pa, q)
            with self.subTest(dx=dx, quantity='actuator volume'):
                actual = np.sum(pa.m[pa.actuator == 1]/pa.rho[pa.actuator == 1])
                self.assertAlmostEqual(actual/expected_volume, 1, places=12)
            with self.subTest(dx=dx, quantity='normal stiffness per area'):
                self.assertAlmostEqual((solver.kn/dx**2)/expected_stiffness, 1, places=12)
            self.assertAlmostEqual(np.sum(pa.m), .0002217375, places=14)
