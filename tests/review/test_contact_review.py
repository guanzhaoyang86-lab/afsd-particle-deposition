"""Independent N05/N06 benchmarks; references are hand-derived, not solver output.

Tolerances frozen before execution in docs/VALIDATION_CONTRACT.md v1.
"""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.contact import WallContact
from afsd_backend.pysph.thermal import conduction


def wall_parameters(**overrides):
    p = dict(dx=.2, traverse=0., omega=0., layer=1., substrate_depth=1.,
             bore_width=.1, tool_radius=1., tool_top=2.,
             base_temperature=350., tool_temperature=350.,
             contact_conductance=0., heat_partition=.4)
    p.update(overrides)
    return p


class ContactReview(unittest.TestCase):
    def test_instantaneous_force_api_three_dimensional(self):
        wall = WallContact(1e5, .3)
        normal = np.array([1., 2., 3.])/np.sqrt(14.)
        tangent = np.array([2., -1., 0.])/np.sqrt(5.)
        force, power, potential, penetration = wall.force(.009, normal, tangent, .01, 2., .001)
        np.testing.assert_allclose(force, 100.*normal-30.*tangent, rtol=1e-12, atol=1e-12)
        self.assertAlmostEqual(power, 30., places=10)
        self.assertAlmostEqual(potential, .05, places=12)
        self.assertAlmostEqual(penetration, .001, places=12)

    def _kick(self, moving=False, dt=.1):
        p = wall_parameters(traverse=.5 if moving else 0.)
        # One compressed plane; tangent x, normal +/- z. No thermal gradient.
        x = np.array([[.15, 0., .91]]) if moving else np.array([[0., 0., -.91]])
        v = np.array([[1.5 if moving else 1., 0., 0.]])
        mass = np.array([2.]); temperature = np.array([350.])
        wall = WallContact(1e5, .3)
        force, _, _ = wall.evaluate(x, v, mass, temperature, 0., dt, p)
        after = v+dt*force/mass[:, None]
        _, heat, ledger = wall.evaluate(x, v, mass, temperature, 0., dt, p, heat_velocity=.5*(v+after))
        # Independent kinetic-energy difference in tangential DOF and motor work.
        delta_kinetic = .5*2.*(after[0, 0]**2-v[0, 0]**2)
        wall_work = dt*force[0, 0]*(.5 if moving else 0.)
        dissipated = wall_work-delta_kinetic
        deposited = dt*(heat.sum()+ledger['wall_heat_rate'])
        np.testing.assert_allclose(deposited, dissipated, rtol=1e-10, atol=1e-12)
        return dissipated, deposited, after[0, 0]

    def test_stationary_wall_impulse_cap_discrete_heat(self):
        dissipation, heat, speed = self._kick()
        self.assertAlmostEqual(dissipation, 1., places=12)
        self.assertAlmostEqual(speed, 0., places=12)

    def test_moving_wall_impulse_cap_discrete_heat(self):
        dissipation, heat, speed = self._kick(moving=True)
        self.assertAlmostEqual(dissipation, 1., places=12)
        self.assertAlmostEqual(speed, .5, places=12)

    def test_sliding_not_impulse_capped_discrete_heat(self):
        self._kick(dt=.0001)

    def test_negative_midpoint_dissipation_rejected(self):
        wall = WallContact(1e5, .3); p = wall_parameters()
        with self.assertRaisesRegex(ValueError, 'negative.*dissipation'):
            wall.evaluate(np.array([[0., 0., -.91]]), np.array([[1., 0., 0.]]),
                          np.array([2.]), np.array([350.]), 0., .1, p,
                          heat_velocity=np.array([[-1., 0., 0.]]))

    def test_deep_tool_penetration_fails_closed(self):
        p = wall_parameters(dx=.0005, layer=.0015, substrate_depth=.0015,
                            bore_width=.002496, tool_radius=.0027, tool_top=.004)
        with self.assertRaisesRegex(ValueError, 'rigid tool'):
            WallContact(1e5, .3).evaluate(np.array([[.002, 0., .003]]), np.zeros((1, 3)),
                                        np.ones(1), np.array([350.]), 0., .001, p)

    def test_open_square_bore_not_rejected_or_filled(self):
        p = wall_parameters(dx=.02)
        force, heat, ledger = WallContact(1e5, .3).evaluate(np.array([[0., 0., 1.5]]),
            np.zeros((1, 3)), np.ones(1), np.array([350.]), 0., .001, p)
        np.testing.assert_array_equal(force, 0.)
        self.assertEqual(ledger['contacts'], 0)

    def test_contact_heat_and_row_have_si_units(self):
        p = wall_parameters(contact_conductance=10., base_temperature=300.)
        wall = WallContact(1e5, .3)
        x = np.array([[0., 0., -.91]])
        _, heat, ledger = wall.evaluate(x, np.zeros((1, 3)), np.array([2.]), np.array([350.]), 0., .001, p)
        # h*A = 10 W/(m2 K)*(.2m)^2 = .4 W/K; flux=-20W.
        np.testing.assert_allclose(heat, [-20.], rtol=1e-12, atol=1e-12)
        self.assertAlmostEqual(ledger['wall_heat_rate'], 20., places=12)
        np.testing.assert_allclose(wall.conductance_row(x, 0., p), [.4], rtol=1e-12)

    def test_corotating_particle_has_no_spurious_bore_friction(self):
        # Both velocities are evaluated at the same particle centre. Projecting
        # only the wall velocity would create spurious omega cross gap slip.
        p = wall_parameters(dx=.02, omega=2., traverse=.3)
        t = .4
        x = np.array([[.3*t+.045*np.cos(.8), .045*np.sin(.8), 1.5]])
        velocity = np.array([[.3-2*x[0, 1], 2*(x[0, 0]-.3*t), 0.]])
        force, heat, ledger = WallContact(1e5, .3).evaluate(x, velocity, np.array([2.]),
            np.array([350.]), t, .001, p)
        normal = np.array([-np.cos(.8), -np.sin(.8), 0.])
        np.testing.assert_allclose(force[0], 500.*normal, rtol=1e-12, atol=1e-12)
        self.assertAlmostEqual(ledger['friction_power'], 0., places=12)
        self.assertAlmostEqual(float(heat[0]), 0., places=12)

    def test_conduction_independent_coefficient_sign_and_conservation(self):
        # Independent hand calculation: C=(2*3)/(4*5)*(4*6*8)/(6+8)*(.25*2)/.25^2=1152/35 W/K.
        x = np.array([[0., 0., 0.], [.25, 0., 0.]])
        power, row = conduction(np.array([2., 3.]), np.array([4., 5.]), np.array([400., 300.]),
                               np.array([6., 8.]), x, np.array([0]), np.array([1]), np.array([[2., 0., 0.]]))
        np.testing.assert_allclose(row, [1152/35, 1152/35], rtol=1e-12)
        np.testing.assert_allclose(power, [-115200/35, 115200/35], rtol=1e-12)
        self.assertEqual(float(power.sum()), 0.)


if __name__ == '__main__':
    unittest.main(verbosity=2)
