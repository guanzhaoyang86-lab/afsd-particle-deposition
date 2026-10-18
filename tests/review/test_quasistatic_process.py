"""Independent initialization and pure coupled-map invariants; not AFSD validation."""
from pathlib import Path
import json
import sys
import unittest
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.quasistatic_process import (
    independent_initial_state, QuasistaticTrial, _geometry, pose)
from afsd_backend.pysph.mechanics import positions, operators


def parameters():
    p = json.loads((ROOT/'outputs/contact_material_search_20261004/root/resolved_config.json').read_text())['parameters']
    p['bore_width'] = .002504
    return p


class QuasistaticProcessTests(unittest.TestCase):
    def test_independent_initial_equilibrium_and_physical_constraints(self):
        p = parameters(); pa, state = independent_initial_state(p)
        trial = QuasistaticTrial(pa, state, p)
        x = positions(pa)
        c = trial.contact.evaluate(x, x, pose(p, 0), pose(p, 0), 1)
        self.assertLess(np.max(np.abs(c['force_N'])), 1e-10)
        self.assertLess(c['normal_potential1_J'].sum(), 1e-20)
        self.assertTrue(np.all(pa.S == 0)); self.assertTrue(np.all(pa.ep == 0))
        self.assertTrue(np.all(pa.u == 0)); self.assertTrue(np.all(pa.w == 0))
        self.assertEqual(int(state['free'].sum()), 809)
        self.assertTrue(np.all(state['free'][pa.actuator.astype(bool), :2]))
        self.assertTrue(np.all(~state['free'][pa.actuator.astype(bool), 2]))

    def test_density_gradient_is_affine_complete_in_current_geometry(self):
        p = parameters(); pa, _ = independent_initial_state(p)
        i, j, g, _, _ = operators(pa); x = positions(pa)
        A = np.array([[.001, .002, -.001], [.003, -.002, .001], [0, -.002, .001]])
        rho = p['rho0']*(1+.001*np.sin(np.arange(len(x))))
        H, B = _geometry(x, x@A.T, pa.m, rho, i, j, g)
        np.testing.assert_allclose(H, np.broadcast_to(A, H.shape), atol=1e-16)

    def test_rest_thermal_coupling_is_pure_and_passive_labels_do_not_enter(self):
        p = parameters()
        p.update(feed=0., traverse=0., omega=0., base_temperature=723.15, tool_temperature=723.15)
        pa, state = independent_initial_state(p)
        before = {k: v.get_npy_array().copy() for k, v in pa.properties.items()}
        trial = QuasistaticTrial(pa, state, p)
        r = trial.evaluate(np.zeros((len(pa.x), 3)), .0005)
        for key, value in before.items():
            np.testing.assert_array_equal(pa.properties[key].get_npy_array(), value)
        self.assertLess(r['force_residual_N'], 1e-10)
        np.testing.assert_allclose(r['temperature'], pa.temperature, rtol=0., atol=1e-12)
        pa.tracer_type[:] = 91; pa.wholebar_tag[:] = -4; pa.material_id[:] = 777
        other = QuasistaticTrial(pa, state, p).evaluate(np.zeros((len(pa.x), 3)), .0005)
        for key in ('force', 'temperature', 'stress', 'ep', 'density'):
            np.testing.assert_array_equal(r[key], other[key])
        np.testing.assert_array_equal(state['contact_history'], np.zeros((len(pa.x), 8, 3)))

    def test_prescribed_displacement_violation_rejects_without_mutation(self):
        p = parameters(); pa, state = independent_initial_state(p)
        trial = QuasistaticTrial(pa, state, p)
        before = positions(pa)
        with self.assertRaisesRegex(ValueError, 'prescribed physical axes'):
            trial.evaluate(np.zeros_like(before), .0005)
        np.testing.assert_array_equal(positions(pa), before)


if __name__ == '__main__': unittest.main()
