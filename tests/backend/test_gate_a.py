import unittest
import numpy as np
from pysph.base.utils import get_particle_array
from afsd_backend.pysph.kernel import HyperbolicSpline
from afsd_backend.pysph.material import FKS, update, equivalent, rotation
from afsd_backend.pysph.mechanics import pairs, operators, internal_force, pressure, volumetric_energy


class GateA(unittest.TestCase):
    def test_kernel_integral_derivative(self):
        k = HyperbolicSpline()
        r = np.linspace(0, 2, 20001)
        self.assertAlmostEqual(np.trapezoid(4*np.pi*r*r*k.value(r, 1), r), 1, places=8)
        for q in (0.2, 0.99, 1.01, 1.8):
            numerical = (k.value(q+1e-6, 1)-k.value(q-1e-6, 1))/2e-6
            self.assertAlmostEqual(k.gradient([q, 0, 0], 1)[0], numerical, places=8)
        with self.assertRaises(ValueError):
            k.gradient([0, 0, 0], 1)

    def test_nnps_is_three_dimensional(self):
        rng = np.random.default_rng(71)
        x = rng.random((60, 3))
        pa = get_particle_array(x=x[:, 0], y=x[:, 1], z=x[:, 2], h=np.full(60, .14))
        i, j = pairs(pa)
        oracle = {(a, b) for a in range(60) for b in range(a+1, 60) if np.linalg.norm(x[a]-x[b]) < .28}
        self.assertEqual(set(zip(i, j)), oracle)
        self.assertTrue(any(abs(x[a, 2]-x[b, 2]) > .05 for a, b in oracle))

    def test_gradient_and_pair_momentum(self):
        x = np.array(np.meshgrid(*([np.arange(-5, 6)]*3), indexing='ij')).reshape(3, -1).T.astype(float)
        L = np.array([[.1, .2, .3], [.4, -.2, .1], [.2, .3, .1]])
        v = x@L.T
        pa = get_particle_array(x=x[:, 0], y=x[:, 1], z=x[:, 2], u=v[:, 0], v=v[:, 1], w=v[:, 2], h=np.full(len(x), 1.3), m=np.ones(len(x)), rho=np.ones(len(x)))
        i, j, g, grad, drho = operators(pa)
        center = np.flatnonzero(np.linalg.norm(x, axis=1) == 0)[0]
        np.testing.assert_allclose(grad[center], L, rtol=.025, atol=1e-12)
        stress = np.broadcast_to(np.diag([2., 3., 4.]), grad.shape)
        np.testing.assert_allclose(internal_force(pa, stress, i, j, g).sum(axis=0), 0, atol=1e-11)

    def test_fks_independent_and_invalid(self):
        law = FKS()
        ts = (623.15-293.15)/562
        expected = (317.6e6+1e6*np.arctan(70*.2))*(1+12.44*ts**5.65*np.log(10))*(1-(1+np.exp(-13*ts))**(-100))
        self.assertAlmostEqual(float(law.yield_stress(.2, 10, 623.15))/expected, 1, places=12)
        self.assertLess(law.yield_stress(.2, 10, 723.15), law.yield_stress(.2, 10, 623.15))
        with self.assertRaises(ValueError):
            law.yield_stress(0, 0, 900)
        with self.assertRaises(ValueError):
            FKS(b1=12.44e6, rate_floor=.001).yield_stress(0, .001, 700)

    def test_prestressed_rigid_rotation(self):
        s = np.array([[[2e6, 1e6, .2e6], [1e6, -1e6, .3e6], [.2e6, .3e6, -1e6]]])
        initial = s.copy()
        spin = np.array([[[0., -2, 0], [2, 0, -1], [0, 1, 0]]])
        ep = np.zeros(1)
        for _ in range(100):
            s, ep, work, diss, rate = update(s, ep, spin, np.array([293.15]), .001, 26.3e9, FKS())
        R = rotation(spin*.1)
        np.testing.assert_allclose(s, R@initial@R.swapaxes(-1, -2), rtol=2e-12)
        np.testing.assert_allclose(equivalent(s), equivalent(initial), rtol=2e-12)
        self.assertEqual(float(ep[0]), 0)

    def test_plastic_all_three_shear_planes(self):
        for a, b in ((0, 1), (0, 2), (1, 2)):
            L = np.zeros((1, 3, 3)); L[0, a, b] = L[0, b, a] = 100
            s, ep, work, diss, rate = update(np.zeros_like(L), np.zeros(1), L, np.array([723.15]), 1e-5, 26.3e9, FKS())
            y = FKS().yield_stress(ep, rate, np.array([723.15]))
            np.testing.assert_allclose(equivalent(s), y, rtol=1e-10)
            self.assertGreater(ep[0], 0); self.assertGreater(work[0], 0)
            self.assertAlmostEqual(np.trace(s[0]), 0)
            qtrial = np.sqrt(3)*2*26.3e9*100*1e-5
            self.assertAlmostEqual((qtrial-equivalent(s)[0])/(3*26.3e9), ep[0], places=13)

    def test_eos_energy_derivative(self):
        r = np.array([2600., 2700., 2800.]); eps = .01
        de = (volumetric_energy(r+eps, 2700, 68e9)-volumetric_energy(r-eps, 2700, 68e9))/(2*eps)
        np.testing.assert_allclose(de, pressure(r, 2700, 68e9)/r**2, atol=2e-6, rtol=1e-7)

    def test_elastic_unload_and_volumetric_split(self):
        law = FKS(); G = 26.3e9
        L = np.array([np.diag([10., -5., -5.])])
        s, ep, work, diss, rate = update(np.zeros((1, 3, 3)), np.zeros(1), L, np.array([293.15]), 1e-6, G, law)
        unloaded, ep2, *_ = update(s, ep, -L, np.array([293.15]), 1e-6, G, law)
        np.testing.assert_allclose(unloaded, 0, atol=1e-9); np.testing.assert_equal(ep2, 0)
        hydro, *_ = update(np.zeros_like(s), ep, np.array([np.eye(3)*100]), np.array([293.15]), 1e-6, G, law)
        np.testing.assert_equal(hydro, 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
