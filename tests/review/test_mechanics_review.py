"""Independent review oracles; never construct expected values with the UUT.

Normal suite: python -B tests/review/test_mechanics_review.py
Open-risk evidence (deliberately fails until the method is repaired): append
--open-risks. Tolerances are frozen in docs/VALIDATION_CONTRACT.md before running.
"""
import math
from decimal import Decimal, localcontext
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))

from pysph.base.utils import get_particle_array
from afsd_backend.pysph.config_adapter import load_settings
from afsd_backend.pysph.kernel import HyperbolicSpline
from afsd_backend.pysph.material import FKS, update
from afsd_backend.pysph.mechanics import (
    pairs, operators, corrected_operators, internal_force, pressure, volumetric_energy,
)
from afsd_backend.pysph.particles import create_geometry, identity_map
from afsd_backend.pysph.output import checkpoint, restore
from afsd_backend.pysph.solver import SolidSolver
from afsd_backend.pysph.geometry import as_rs


def scalar_paper_fks(ep, rate, temperature):
    """Direct P0 Eqs 4-8/Table 2 oracle, using Python math not material.py."""
    theta = (temperature - 293.15) / 562.0
    hardening = 317600000.0 + 1000000.0 * math.atan(70.0 * ep)
    rate_factor = 1.0 + 12.44 * theta ** 5.65 * math.log(max(rate, 1.0))
    # Direct high-precision bracket evaluation avoids oracle cancellation at
    # 800 K; production uses a different double-precision log/expm1 algorithm.
    with localcontext() as context:
        context.prec = 40
        tstar = (Decimal(str(temperature))-Decimal('293.15'))/Decimal('562')
        softening = float(1-(1+(-13*tstar).exp())**(-100))
    return hardening * rate_factor * softening


def axis_angle_matrix(axis, angle):
    """Independent closed-form axis/angle entries, no UUT rotation helper."""
    x, y, z = np.asarray(axis) / np.linalg.norm(axis)
    c, s, d = math.cos(angle), math.sin(angle), 1 - math.cos(angle)
    return np.array([[c+x*x*d, x*y*d-z*s, x*z*d+y*s],
                     [y*x*d+z*s, c+y*y*d, y*z*d-x*s],
                     [z*x*d-y*s, z*y*d+x*s, c+z*z*d]])


def seeded_array(p):
    """Nonzero, admissible state keyed by immutable ID, not by storage order."""
    pa = create_geometry(p)
    a = pa.particle_id.astype(float)
    pa.u[:] += 0.0002 * np.sin(a)
    pa.v[:] += 0.0003 * np.cos(a)
    pa.w[:] += 0.0001 * np.sin(0.4*a)
    pa.rho[:] *= 1 + 1e-7*np.sin(a)
    pa.temperature[:] += 0.2*np.cos(a)
    pa.ep[:] = 1e-6*(1+np.sin(a))
    s = pa.S.reshape(-1, 3, 3)
    s[:, 0, 0] = 1e4*np.sin(a)
    s[:, 1, 1] = -3e3*np.cos(a)
    s[:, 2, 2] = -s[:, 0, 0]-s[:, 1, 1]
    for j, k in ((0, 1), (0, 2), (1, 2)):
        s[:, j, k] = s[:, k, j] = 500*(j+k+1)*np.sin(a)
    pa.D[:] = np.arange(pa.D.size)*1e-6
    pa.plastic_work[:] = a*1e-15
    pa.return_dissipation[:] = a*2e-16
    pa.feed_crossed[:] = (a % 7 == 0).astype(float)
    return pa


def keyed_property(pa, name):
    values = pa.properties[name].get_npy_array()
    stride = len(values)//pa.get_number_of_particles()
    return values.reshape(-1, stride)[np.argsort(pa.particle_id)]


class ReviewMechanics(unittest.TestCase):
    def test_kernel_piecewise_3d_normalization_and_independent_derivative(self):
        kernel = HyperbolicSpline()
        # Four Gauss points integrate r^2 times either cubic branch exactly.
        node, weight = np.polynomial.legendre.leggauss(4)
        for h in (0.0007, 0.13, 1.9):
            integral = 0.
            for left, right in ((0., h), (h, 2*h)):
                r = (node+1)*(right-left)/2+left
                integral += np.sum(weight*4*math.pi*r*r*kernel.value(r, h))*(right-left)/2
            self.assertLess(abs(integral-1), 1e-8)
        h = 0.43
        direction = np.array([2., -3., 6.])/7
        for q in (0.2, 0.9, 1., 1.5, 2.1):
            df = 3*q*q-6 if q < 1 else -3*(2-q)**2 if q < 2 else 0.
            expected = (15/(62*math.pi))*df/(h**4)*direction
            np.testing.assert_allclose(kernel.gradient(direction*q*h, h), expected, rtol=1e-12, atol=1e-14)
        self.assertEqual(float(kernel.value(2*h, h)), 0.)

    def test_non_coplanar_nnps_and_genuine_out_of_plane_response(self):
        x = np.array([[0, 0, 0], [.1, .2, .3], [.1, .2, 1.2], [.2, -.2, -.3], [-.3, .1, .2]])
        n = len(x); h = .4
        pa = get_particle_array(x=x[:, 0], y=x[:, 1], z=x[:, 2], h=np.full(n, h), m=np.ones(n), rho=np.ones(n))
        expected = {(i, j) for i in range(n) for j in range(i+1, n) if sum((x[i]-x[j])**2) < (2*h)**2}
        i, j = pairs(pa)
        self.assertEqual(set(zip(i, j)), expected)
        # Pure x-z shear stress must exert z force: a planar copy cannot pass.
        stress = np.zeros((n, 3, 3)); stress[:, 0, 2] = stress[:, 2, 0] = 2.
        i, j, g, _, _ = operators(pa)
        f = internal_force(pa, stress, i, j, g)
        self.assertGreater(float(np.linalg.norm(f[:, 2])), 0.)
        self.assertLess(float(np.linalg.norm(f.sum(axis=0))), 1e-12*float(np.linalg.norm(f, axis=1).sum()))

    def test_fks_against_independent_scalar_formula_and_explicit_domain(self):
        law = FKS()
        for ep, rate, temperature in ((0., 0., 293.15), (.2, .4, 600.), (.4, 10., 723.15), (.9, 100., 800.)):
            expected = scalar_paper_fks(ep, rate, temperature)
            actual = float(law.yield_stress(ep, rate, temperature))
            self.assertLess(abs(actual/expected-1), 5e-12)
        for temperature in (293.14, 855.15, 900.):
            with self.assertRaises(ValueError): law.yield_stress(0., 1., temperature)
        self.assertEqual(float(law.yield_stress(.2, 0., 700.)), float(law.yield_stress(.2, 1., 700.)))

    def test_prestressed_arbitrary_axis_objectivity_independent_rotation(self):
        axis = np.array([2., -1., 3.]); axis /= np.linalg.norm(axis)
        wx, wy, wz = 2.3*axis
        spin = np.array([[[0, -wz, wy], [wz, 0, -wx], [-wy, wx, 0]]])
        s0 = np.array([[[2e5, 4e4, -7e4], [4e4, -1.1e5, 3e4], [-7e4, 3e4, -9e4]]])
        s = s0.copy(); ep = np.zeros(1); heat = 0.
        dt = .002
        for _ in range(250):
            s, ep, work, _, _ = update(s, ep, spin, np.array([293.15]), dt, 26.3e9, FKS())
            heat += work[0]
        R = axis_angle_matrix(axis, 2.3*.5)
        expected = R@s0[0]@R.T
        self.assertLess(np.linalg.norm(s[0]-expected)/np.linalg.norm(s0), 1e-11)
        self.assertLess(abs(np.sum(s*s)/np.sum(s0*s0)-1), 1e-11)
        self.assertLess(abs(np.trace(s[0]))/np.linalg.norm(s0), 1e-11)
        self.assertEqual(float(ep[0]), 0.); self.assertEqual(float(heat), 0.)

    def test_three_dimensional_j2_perfect_plastic_analytic_oracle(self):
        G, dt, yield_value = 26.3e9, 1e-4, 1e6
        law = FKS(a1=yield_value, a2=0., b1=0., c1=0.)
        D = np.array([[[9., 4., 7.], [4., -3., -2.], [7., -2., -6.]]])
        trial = 2*G*dt*D
        q_trial = math.sqrt(1.5*float(np.sum(trial*trial)))
        dep = (q_trial-yield_value)/(3*G)
        expected = trial*(yield_value/q_trial)
        s, ep, work, numerical, _ = update(np.zeros_like(D), np.zeros(1), D, np.array([293.15]), dt, G, law)
        np.testing.assert_allclose(s, expected, rtol=2e-12)
        self.assertLess(abs(ep[0]/dep-1), 2e-12)
        self.assertLess(abs(work[0]/(yield_value*dep)-1), 2e-12)
        self.assertLess(abs(numerical[0]/(1.5*G*dep*dep)-1), 2e-12)

    def test_fks_j2_return_independent_scalar_bisection(self):
        G, dt, temperature = 26.3e9, 1e-5, 723.15
        D = np.array([[[10., 20., -30.], [20., 40., 50.], [-30., 50., -50.]]])
        rate = math.sqrt((2/3)*float(np.sum(D*D)))
        trial = 2*G*dt*D
        q_trial = math.sqrt(1.5*float(np.sum(trial*trial)))
        lo, hi, ep0 = 0., q_trial/(3*G), .12
        for _ in range(100):
            increment = (lo+hi)/2
            f = q_trial-3*G*increment-scalar_paper_fks(ep0+increment, rate, temperature)
            if f > 0: lo = increment
            else: hi = increment
        dep = (lo+hi)/2
        q_expected = scalar_paper_fks(ep0+dep, rate, temperature)
        s, ep, _, _, actual_rate = update(np.zeros_like(D), np.array([ep0]), D, np.array([temperature]), dt, G, FKS())
        np.testing.assert_allclose(s, trial*(q_expected/q_trial), rtol=1e-10)
        self.assertLess(abs(ep[0]-ep0-dep), 1e-12)
        self.assertLess(abs(actual_rate[0]/rate-1), 1e-12)

    def test_bulk_pressure_and_energy_against_scalar_primitive(self):
        for rho in (2590., 2700., 2810.):
            K, rho0 = 68e9, 2700.
            expected_p = K*(rho-rho0)/rho0
            expected_e = (K/rho0)*math.log(rho/rho0)+K/rho-K/rho0
            self.assertLess(abs(float(pressure(np.array([rho]), rho0, K)[0])-expected_p), 2e-5)
            self.assertLess(abs(float(volumetric_energy(np.array([rho]), rho0, K)[0])-expected_e), 1e-8)

    def test_full_nonzero_state_sort_and_checkpoint_identity(self):
        _, p, _ = load_settings(ROOT/'configs/debug/starved.json')
        pa = seeded_array(p)
        before = {k: keyed_property(pa, k).copy() for k in pa.properties}
        identity = identity_map(pa)
        order = np.random.default_rng(411).permutation(pa.get_number_of_particles()).astype(np.uint32)
        pb = pa.extract_particles(order)
        for key, value in before.items(): np.testing.assert_array_equal(keyed_property(pb, key), value)
        with tempfile.TemporaryDirectory(dir=ROOT/'logs/review') as tmp:
            path = Path(tmp)/'nonzero_state.npz'
            checkpoint(path, pb, {'t': 0., 'next_id': int(pa.particle_id.max())+1})
            pc, state = restore(path)
            for key, value in before.items(): np.testing.assert_array_equal(keyed_property(pc, key), value)
            self.assertEqual(identity_map(pc), identity)
            self.assertEqual(state['next_id'], int(pa.particle_id.max())+1)

    def test_label_swap_and_removal_do_not_change_physics(self):
        _, p, _ = load_settings(ROOT/'configs/debug/starved.json')
        pa, pb, pc = (seeded_array(p) for _ in range(3))
        pb.tracer_type[:] = np.where(pb.tracer_type == 1, 2, np.where(pb.tracer_type == 2, 1, 0))
        pb.wholebar_tag[:] = np.where(pb.wholebar_tag == 1, 2, np.where(pb.wholebar_tag == 2, 1, 0))
        pc.remove_property('tracer_type'); pc.remove_property('wholebar_tag')
        solvers = [SolidSolver(a, p) for a in (pa, pb, pc)]
        for _ in range(4):
            for solver in solvers: solver.advance()
        for field in ('x', 'y', 'z', 'u', 'v', 'w', 'rho', 'S', 'temperature', 'ep', 'D', 'plastic_work'):
            np.testing.assert_array_equal(getattr(pa, field), getattr(pb, field))
            np.testing.assert_array_equal(getattr(pa, field), getattr(pc, field))
        self.assertEqual(solvers[0].ledger, solvers[1].ledger)
        self.assertEqual(solvers[0].ledger, solvers[2].ledger)

    def test_reordered_nonzero_state_dynamics_and_restart(self):
        _, p, _ = load_settings(ROOT/'configs/debug/starved.json')
        pa = seeded_array(p)
        order = np.random.default_rng(733).permutation(pa.get_number_of_particles()).astype(np.uint32)
        pb = pa.extract_particles(order)
        a, b = SolidSolver(pa, p), SolidSolver(pb, p)
        for _ in range(3): a.advance(); b.advance()
        for field in ('x', 'y', 'z', 'u', 'v', 'w', 'rho', 'S', 'temperature', 'ep', 'D', 'plastic_work'):
            va, vb = keyed_property(pa, field), keyed_property(pb, field)
            scale = max(float(np.linalg.norm(va)), 1e-30)
            self.assertLess(float(np.linalg.norm(va-vb))/scale, 1e-12, field)
        with tempfile.TemporaryDirectory(dir=ROOT/'logs/review') as tmp:
            path = Path(tmp)/'restart.npz'; checkpoint(path, pb, b.state())
            pc, state = restore(path); c = SolidSolver(pc, p, state)
            b.advance(); c.advance()
            for field in pb.properties:
                np.testing.assert_array_equal(pb.properties[field].get_npy_array(), pc.properties[field].get_npy_array())
            self.assertEqual(b.ledger, c.ledger)

    def test_as_rs_follows_kinematics_not_labels(self):
        y = np.array([-2., -1., 0., 1., 2.])
        for omega, traverse in ((31.4, .00212), (-31.4, .00212), (31.4, -.00212)):
            r = np.column_stack((np.ones(5), y, np.zeros(5)))
            vrot = np.cross(np.array([0., 0., omega]), r)
            expected = np.sign(vrot@np.array([traverse, 0., 0.])).astype(int)
            np.testing.assert_array_equal(as_rs(y, omega, traverse), expected)

    def test_corrected_irregular_prestressed_rigid_rotation_chain(self):
        x = np.array([[0, 0, 0], [.2, .1, .3], [-.1, .25, .15], [.3, -.1, -.2], [-.2, -.15, .05]])
        omega = np.array([1., 2., -3.])
        v = np.cross(np.broadcast_to(omega, x.shape), x)
        pa = get_particle_array(x=x[:, 0], y=x[:, 1], z=x[:, 2], u=v[:, 0], v=v[:, 1], w=v[:, 2], h=np.full(len(x), .35), m=np.full(len(x), .02), rho=np.ones(len(x)))
        _, _, _, L, drho, _ = corrected_operators(pa)
        spin = np.array([[0., 3., 2.], [-3., 0., -1.], [-2., 1., 0.]])
        np.testing.assert_allclose(L, np.broadcast_to(spin, L.shape), rtol=1e-11, atol=1e-12)
        self.assertLess(float(np.max(np.abs(drho))), 1e-12)
        s0 = np.broadcast_to(np.array([[2e5, 4e4, -7e4], [4e4, -1.1e5, 3e4], [-7e4, 3e4, -9e4]]), L.shape).copy()
        dt = 1e-4
        s, ep, heat, _, _ = update(s0, np.zeros(len(x)), L, np.full(len(x), 293.15), dt, 26.3e9, FKS())
        R = axis_angle_matrix(omega, np.linalg.norm(omega)*dt)
        reference = np.broadcast_to(R@s0[0]@R.T, s.shape)
        self.assertLess(np.linalg.norm(s-reference)/np.linalg.norm(s0), 1e-10)
        np.testing.assert_array_equal(ep, 0.)
        np.testing.assert_array_equal(heat, 0.)

    def test_corrected_unequal_density_work_and_angular_momentum(self):
        x = np.array([[0, 0, 0], [.2, .1, .3], [-.1, .25, .15], [.3, -.1, -.2], [-.2, -.15, .05]], float)
        rng = np.random.default_rng(872)
        v = rng.normal(size=x.shape)
        pa = get_particle_array(x=x[:, 0], y=x[:, 1], z=x[:, 2], u=v[:, 0], v=v[:, 1], w=v[:, 2], h=np.full(len(x), .35), m=np.array([2., 3., 4., 2., 3.]), rho=np.array([1., 1.2, .8, .9, 1.1]))
        stress = rng.normal(size=(len(x), 3, 3)); stress += stress.swapaxes(1, 2)
        i, j, g, L, _, B = corrected_operators(pa)
        force = internal_force(pa, stress, i, j, g, correction=B)
        kinetic_power = float(np.sum(force*v))
        internal_power = float(np.sum((pa.m/pa.rho)*np.sum(stress*L, axis=(1, 2))))
        defect = abs(kinetic_power+internal_power)/max(abs(kinetic_power), abs(internal_power))
        self.assertLess(defect, 1e-12)
        force_scale = float(np.linalg.norm(force, axis=1).sum())
        self.assertLess(float(np.linalg.norm(force.sum(axis=0)))/force_scale, 1e-12)
        torque = np.cross(x, force).sum(axis=0)
        self.assertLess(float(np.linalg.norm(torque))/(force_scale*np.linalg.norm(x, axis=1).max()), 1e-12)

    def test_correction_fails_closed_for_planar_or_isolated_neighborhoods(self):
        for x in (np.zeros((1, 3)), np.array([[0., 0., 0.], [.1, .2, 0.], [-.1, .2, 0.], [.2, .3, 0.]])):
            pa = get_particle_array(x=x[:, 0], y=x[:, 1], z=x[:, 2], h=np.full(len(x), .3), m=np.ones(len(x)), rho=np.ones(len(x)))
            with self.assertRaisesRegex(ValueError, 'degenerate/ill-conditioned'):
                corrected_operators(pa)


class OpenMethodRisks(unittest.TestCase):
    """Unsatisfied physical invariants are retained as FAIL evidence, not xfail."""
    def test_rigid_rotation_sph_boundary_must_not_create_strain(self):
        x = np.array([[0, 0, 0], [.2, .1, .3], [-.1, .25, .15], [.3, -.1, -.2], [-.2, -.15, .05]])
        omega = np.array([1., 2., -3.])
        v = np.cross(np.broadcast_to(omega, x.shape), x)
        pa = get_particle_array(x=x[:, 0], y=x[:, 1], z=x[:, 2], u=v[:, 0], v=v[:, 1], w=v[:, 2], h=np.full(len(x), .35), m=np.full(len(x), .02), rho=np.ones(len(x)))
        _, _, _, L, _ = operators(pa)
        D = (L+L.swapaxes(1, 2))/2
        defect = float(np.linalg.norm(D)/np.linalg.norm(L))
        print('OPEN_RISK rigid_rotation_spurious_symmetric_gradient_relative=', defect)
        self.assertLess(defect, 1e-11, 'uncorrected SPH boundary gradient is not first-order complete')

    def test_unequal_density_table1_operators_must_be_power_conjugate(self):
        x = np.array([[0., 0., 0.], [.3, .1, -.2]])
        v = np.array([[1., 2., -.5], [-.7, .2, .8]])
        pa = get_particle_array(x=x[:, 0], y=x[:, 1], z=x[:, 2], u=v[:, 0], v=v[:, 1], w=v[:, 2], h=np.full(2, .4), m=np.array([2., 3.]), rho=np.array([1., 1.2]))
        sigma = np.array([[[3., .4, -.2], [.4, 2., .5], [-.2, .5, 1.]], [[2., -.3, .6], [-.3, 4., .1], [.6, .1, -1.]]])
        i, j, g, L, _ = operators(pa)
        force = internal_force(pa, sigma, i, j, g)
        kinetic_power = float(np.sum(force*v))
        internal_power = float(np.sum((pa.m/pa.rho)*np.sum(sigma*L, axis=(1, 2))))
        relative = abs(kinetic_power+internal_power)/max(abs(kinetic_power), abs(internal_power))
        print('OPEN_RISK unequal_density_kinetic_power=', kinetic_power, 'internal_power=', internal_power, 'relative_defect=', relative)
        self.assertLess(relative, 1e-12, 'Table1 volume-gradient and rho^-2 stress force are not conjugate at unequal density')


class UnresolvedThermodynamics(unittest.TestCase):
    """Known uncorrected energy/free-energy incompatibility; not a safety gate."""
    def test_prestressed_compression_needs_consistent_free_energy(self):
        # At rho=rho0, pressure=0 and D=-I/s: sigma:D=0 because trace(S)=0.
        # Jaumann hypoelasticity gives Sdot=2G dev(D)=0. But the currently
        # declared Es=m/rho*S:S/(4G) has Esdot=-3 Es/s under drho/dt=3rho/s.
        # Decreasing dt does not remove this instantaneous incompatibility.
        # This is not an invitation to inject compensating heat into a ledger.
        G, rho0, mass = 26.3e9, 2700., 1.
        s0 = np.array([np.diag([2e6, -1e6, -1e6])])
        D = np.array([-np.eye(3)])
        dt = 1e-6
        s1, ep, heat, _, _ = update(s0, np.zeros(1), D, np.array([293.15]), dt, G, FKS())
        rho1 = rho0*math.exp(3*dt)
        es0 = mass/rho0*float(np.sum(s0*s0))/(4*G)
        es1 = mass/rho1*float(np.sum(s1*s1))/(4*G)
        declared_rate = (es1-es0)/dt
        independent_mechanical_power = mass/rho0*float(np.sum(s0*D))
        defect = abs(declared_rate-independent_mechanical_power)/(3*es0)
        print('OPEN_RISK hypoelastic_compression Es0_J=', es0,
              'declared_Es_rate_J_s=', declared_rate,
              'continuum_deviatoric_power_J_s=', independent_mechanical_power,
              'expected_exact_Es_rate_J_s=', -3*es0,
              'relative_nonconjugacy=', defect)
        np.testing.assert_array_equal(s1, s0)
        np.testing.assert_array_equal(ep, 0.)
        np.testing.assert_array_equal(heat, 0.)
        self.assertLess(defect, 1e-10, 'finite-volume stored shear energy is not consistent with constant-G Jaumann hypoelasticity under compression; requires constitutive/free-energy work, not ledger compensation')


if __name__ == '__main__':
    open_risks = '--open-risks' in sys.argv
    if '--hypoelastic-open-risk' in sys.argv:
        selected = UnresolvedThermodynamics
    else:
        selected = OpenMethodRisks if open_risks else ReviewMechanics
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(selected)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
