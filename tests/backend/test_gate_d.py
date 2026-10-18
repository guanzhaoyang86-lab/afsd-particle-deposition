import tempfile
import unittest
from pathlib import Path
import numpy as np
from afsd_backend.pysph.config_adapter import load_settings
from afsd_backend.pysph.particles import create_geometry
from afsd_backend.pysph.output import checkpoint, restore
from afsd_backend.pysph.integrator import HostLeapfrog
from afsd_backend.pysph.solver import SolidSolver

ROOT = Path(__file__).resolve().parents[2]


class GateD(unittest.TestCase):
    def setUp(self):
        _, self.p, _ = load_settings(ROOT/'configs/debug/starved.json')

    def test_leapfrog_constant_acceleration_and_order(self):
        def oscillator(dt):
            pa = create_geometry(self.p).extract_particles(np.array([0], dtype=np.uint32))
            pa.x[:] = 1; pa.u[:] = 0
            def acceleration(t, step): pa.au[:] = -pa.x
            stepper = HostLeapfrog(pa, acceleration)
            for n in range(round(1/dt)): stepper.advance(n*dt, dt)
            return abs(pa.x[0]-np.cos(1))
        self.assertGreater(oscillator(.04)/oscillator(.02), 3.9)
        pa = create_geometry(self.p)
        initial = pa.x.copy(); initial_v = pa.u.copy()
        def constant(t, dt): pa.au[:] = 2
        HostLeapfrog(pa, constant).advance(0, .001)
        np.testing.assert_allclose(pa.x, initial+.001*initial_v+.001**2, atol=1e-15)

    def test_coupled_restart_and_passive_tracers(self):
        pa = create_geometry(self.p); pb = create_geometry(self.p)
        pb.tracer_type[:] = 99; pb.wholebar_tag[:] = 77
        a = SolidSolver(pa, self.p); b = SolidSolver(pb, self.p)
        for _ in range(3): a.advance(); b.advance()
        for field in ('x', 'y', 'z', 'u', 'v', 'w', 'rho', 'S', 'temperature', 'ep'):
            np.testing.assert_array_equal(getattr(pa, field), getattr(pb, field))
        with tempfile.TemporaryDirectory() as d:
            checkpoint(Path(d)/'checkpoint.npz', pa, a.state())
            pc, state = restore(Path(d)/'checkpoint.npz'); c = SolidSolver(pc, self.p, state)
            a.advance(); c.advance()
            for field in ('x', 'y', 'z', 'u', 'v', 'w', 'rho', 'S', 'temperature', 'ep'):
                np.testing.assert_array_equal(getattr(pa, field), getattr(pc, field))
            self.assertEqual(a.ledger, c.ledger)
        self.assertEqual(a.diagnostics()['alarms'], [])
        self.assertGreater(a.last['pair_count'], 1000)

    def test_nan_failure(self):
        solver = SolidSolver(create_geometry(self.p), self.p)
        solver.pa.ep[0] = np.nan
        with self.assertRaises(RuntimeError): solver.diagnostics()
