import unittest
import numpy as np
from afsd_backend.pysph.kernel import HyperbolicSpline
from afsd_backend.pysph.thermal import conduction, temperature_step, effusivity_partition
from afsd_backend.pysph.contact import WallContact
from afsd_backend.pysph.geometry import surfaces, local_xy


class GateB(unittest.TestCase):
    def test_conduction_exact_two_body_decay(self):
        x = np.array([[0., 0, 0], [1., 0, 0]])
        i, j = np.array([0]), np.array([1])
        g = HyperbolicSpline().gradient(x[i]-x[j], 1)
        T = np.array([400., 300.]); m = np.array([2., 2.]); rho = np.ones(2)
        power, row = conduction(m, rho, T, 1., x, i, j, g)
        coefficient = row[0]; dt = .0001
        initial = np.dot(m, T)
        for _ in range(1000):
            power, _ = conduction(m, rho, T, 1., x, i, j, g)
            T += dt*power/m
        expected = 100*np.exp(-2*coefficient/m[0]*dt*1000)
        self.assertAlmostEqual(np.dot(m, T), initial, places=9)
        self.assertLess(abs((T[0]-T[1])/expected-1), 2e-4)
        uniform, _ = conduction(m, rho, np.ones(2)*350, 1., x, i, j, g)
        np.testing.assert_equal(uniform, 0)

    def test_thermal_domain_and_plastic_heat(self):
        result = temperature_step(np.array([400.]), np.array([900.]), np.array([1.]), 900., 293.15, 855.15)
        self.assertEqual(result[0], 401.)
        with self.assertRaises(ValueError):
            temperature_step(result, np.array([1e9]), np.ones(1), 900., 293., 855.)
        self.assertEqual(effusivity_partition(1, 1, 1, 1, 1, 1), .5)

    def test_contact_work_heat_three_dimensions(self):
        wall = WallContact(1e5, .3)
        for n in (np.array([1., 0, 0]), np.array([0., 1, 0]), np.array([0., 0, 1]), np.ones(3)/np.sqrt(3)):
            v = np.cross(n, [.5, .7, .9]); vw = np.array([.3, -.2, .1])
            f, heat, potential, overlap = wall.force(.009, n, v, .01, 1., 1e-5)
            self.assertAlmostEqual(np.dot(f, n), 100.)
            self.assertGreater(heat, 0)
            self.assertAlmostEqual(np.dot(f, v+vw)+heat, np.dot(f, vw), places=10)
            self.assertAlmostEqual(potential, .05)
        f, heat, potential, overlap = wall.force(.02, [0, 0, 1], [1, 0, 0], .01, 1., 1e-5)
        np.testing.assert_equal(f, 0); self.assertEqual(heat, 0)

    def test_open_bore_rotates_flat_shoulder(self):
        p = dict(dx=.001, traverse=0., omega=np.pi/2, layer=.0015, substrate_depth=.002, bore_width=.004, tool_radius=.004, tool_top=.01, base_temperature=293., tool_temperature=700.)
        x = np.array([[0., 0, .0014], [.003, 0, .0014], [.0017, 0, .004]])
        initial = list(surfaces(x, 0, p))
        self.assertNotIn(0, initial[1][0]); self.assertIn(1, initial[1][0])
        np.testing.assert_equal(initial[1][2][0], [0, 0, -1])
        local, R = local_xy(x, .5, 0., np.pi/2)
        self.assertAlmostEqual(local[2, 0], .0017/np.sqrt(2))
        self.assertGreater(sum(len(s[0]) for s in initial[2:]), 0)

    def test_three_dimensional_fourier_conduction_refinement(self):
        # Interior periodic Fourier mode; explicit infinite-lattice stencil.
        errors = []
        for dx in (.04, .02, .01):
            offsets = np.array(np.meshgrid(*([np.arange(-3, 4)]*3), indexing='ij')).reshape(3, -1).T*dx
            radius = np.linalg.norm(offsets, axis=1)
            neighbors = offsets[(radius > 0) & (radius < 2.6*dx)]
            x = np.vstack((np.zeros((1, 3)), neighbors))
            m = np.ones(len(x))*dx**3; rho = np.ones(len(x))
            T = 350+np.cos(2*np.pi*x[:, 0]); i = np.zeros(len(neighbors), int); j = np.arange(1, len(x))
            g = HyperbolicSpline().gradient(x[i]-x[j], 1.3*dx)
            power, _ = conduction(m, rho, T, 1., x, i, j, g)
            numerical = power[0]/m[0]
            errors.append(abs(numerical/(-4*np.pi**2)-1))
        self.assertLess(errors[-1], .025)
        self.assertLess(errors[-1], errors[0])


if __name__ == '__main__':
    unittest.main(verbosity=2)
