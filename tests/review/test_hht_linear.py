"""Preregistered smooth HHT algebra/accuracy; never an AFSD acceptance test."""
import copy
import json
import time
import unittest
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from afsd_backend.pysph import hht_linear as tested


def direct_matrices(omega, mass, h, weight):
    """Independent simultaneous x,v,a equations, not the tested reduced solve."""
    gamma = 1.5-weight
    beta = (2.-weight)**2/4.
    left = np.array([[1., 0., -beta*h*h],
                     [0., 1., -gamma*h], [weight*omega**2, 0., 1.]])
    right = np.array([[1., h, (.5-beta)*h*h],
                      [0., 1., (1.-gamma)*h],
                      [-(1.-weight)*omega**2, 0., 0.]])
    forcing = np.zeros((3, 2))
    forcing[2] = [(1.-weight)/mass, weight/mass]
    return np.linalg.solve(left, right), np.linalg.solve(left, forcing)


def independent_terms(m, k, x0, v0, a0, x1, v1, a1, f0, f1, h, w, old_dx, old_df):
    """Independent physical/algorithmic quadratic forms; no tested ledger use."""
    delta = 1.-w
    coefficient = delta*(.5-delta)
    dx, da = x1-x0, a1-a0
    E0 = sum(.5*m*v0*v0)+.5*np.dot(x0, k@x0)
    E1 = sum(.5*m*v1*v1)+.5*np.dot(x1, k@x1)
    extra0 = h*h*delta*delta/8.*sum(m*a0*a0)+coefficient/2.*(old_dx@k@old_dx)
    extra1 = h*h*delta*delta/8.*sum(m*a1*a1)+coefficient/2.*(dx@k@dx)
    Da = h*h*delta**3/4.*sum(m*da*da)
    Dx = coefficient/2.*((dx-old_dx)@k@(dx-old_dx))
    W = .5*sum(dx*(f0+f1))
    Wf = coefficient*sum(dx*(f1-f0-old_df))
    defect = m*a0-f0+k@x0+delta*(old_df-k@old_dx)
    Wh = (.5-delta)*sum(dx*defect)
    return dict(physical_energy0_J=E0, physical_energy1_J=E1,
                algorithmic_energy0_J=E0+extra0, algorithmic_energy1_J=E1+extra1,
                algorithmic_acceleration_loss_J=Da, algorithmic_increment_loss_J=Dx,
                external_trapezoidal_work_J=W, algorithmic_load_work_J=Wf,
                acceleration_history_work_J=Wh,
                identity_residual_J=E1+extra1-E0-extra0+Da+Dx-W-Wf-Wh)


def emit(name, values):
    print('HHT_LINEAR_EVIDENCE '+json.dumps(dict(name=name, **values),
          sort_keys=True, allow_nan=False), flush=True)


class SmoothLinearHHT(unittest.TestCase):
    def check_ledger(self, m, k, old, following, f0, f1, h, w, old_dx, old_df):
        independent = independent_terms(m, k, *old, following.x, following.v,
                                        following.acc, f0, f1, h, w, old_dx, old_df)
        for key, value in independent.items():
            self.assertAlmostEqual(following.ledger[key], value, delta=2e-12, msg=key)
        self.assertLessEqual(abs(independent['identity_residual_J']), 5e-11)
        self.assertGreaterEqual(independent['algorithmic_acceleration_loss_J'], -1e-14)
        self.assertGreaterEqual(independent['algorithmic_increment_loss_J'], -1e-14)
        self.assertFalse(following.ledger['algorithmic_terms_are_material_heat'])
        return independent

    def integrate(self, m, k, x, v, force, h, w, duration=10.):
        x, v = np.array(x, float), np.array(v, float)
        a = (force(0.)-k@x)/m  # independent smooth acceleration initialization
        matrices = [direct_matrices(np.sqrt(k[j, j]/m[j]), m[j], h, w)
                    for j in range(len(m))]
        direct = np.array([x, v, a]).T
        old_dx = np.zeros_like(x)
        old_df = np.zeros_like(x)
        max_errors = np.zeros(3)
        identity_sum, losses, work_sum, path_work = 0., 0., 0., 0.
        max_history = 0.
        first_energy = None
        physical_increases = 0
        for index in range(int(round(duration/h))):
            f0, f1 = force(index*h), force((index+1)*h)
            old = (x.copy(), v.copy(), a.copy())
            result = tested.step(m, k, x, v, a, f0, f1, h, w, old_dx, old_df)
            values = self.check_ledger(m, k, old, result, f0, f1, h, w, old_dx, old_df)
            if first_energy is None:
                first_energy = values['algorithmic_energy0_J']
            identity_sum += values['identity_residual_J']
            losses += values['algorithmic_acceleration_loss_J']+values['algorithmic_increment_loss_J']
            work_sum += values['external_trapezoidal_work_J']+values['algorithmic_load_work_J']
            path_work += result.ledger['exact_linear_load_work_J']
            max_history = max(max_history, result.ledger['acceleration_history_defect_max_N'])
            physical_increases += int(values['physical_energy1_J'] > values['physical_energy0_J']+1e-14)
            if np.all(f0 == 0.) and np.all(f1 == 0.):
                self.assertLessEqual(values['algorithmic_energy1_J']-values['algorithmic_energy0_J'], 5e-11)
            for j, (matrix, load) in enumerate(matrices):
                direct[j] = matrix@direct[j]+load@np.array([f0[j], f1[j]])
            calculated = np.array([result.x, result.v, result.acc]).T
            max_errors = np.maximum(max_errors, np.max(np.abs(calculated-direct), axis=0))
            x, v, a = result.x, result.v, result.acc
            old_dx, old_df = result.displacement_increment, result.load_increment
        np.testing.assert_array_less(max_errors, [2e-12, 2e-11, 2e-9])
        self.assertLessEqual(max_history, 1e-10)
        self.assertLessEqual(abs(identity_sum), 2e-9)
        self.assertAlmostEqual(result.ledger['algorithmic_energy1_J']-first_energy+losses,
                               work_sum, delta=2e-9)
        return result, dict(direct_errors_x_v_acc=max_errors.tolist(),
            maximum_history_defect_N=max_history, cumulative_identity_residual_J=identity_sum,
            algorithmic_loss_J=losses, algorithmic_external_work_J=work_sum,
            conditional_cubic_path_external_work_J=path_work,
            physical_energy_increase_steps=physical_increases)

    def test_01_two_mode_common_final_time(self):
        start = time.perf_counter()
        m, omega = np.array([2., 3.]), np.array([1., 100.])
        k = np.diag(m*omega**2)
        x0, v0 = np.array([.3, .002]), np.array([-.2, .12])
        analytic_x = x0*np.cos(omega*10.)+v0/omega*np.sin(omega*10.)
        analytic_v = -omega*x0*np.sin(omega*10.)+v0*np.cos(omega*10.)
        initial_amp = np.hypot(x0, v0/omega)
        attenuation = {}
        for w in (.9, 1.):
            errors = []
            for h in (.1, .05, .025):
                result, evidence = self.integrate(m, k, x0, v0, lambda t: np.zeros(2), h, w)
                x_error = abs(result.x[0]-analytic_x[0])
                v_error = abs(result.v[0]-analytic_v[0])
                self.assertLessEqual(x_error, .6*h*h)
                self.assertLessEqual(v_error, .6*h*h)
                errors.append(float(np.hypot(x_error, v_error)))
                ratios = np.hypot(result.x, result.v/omega)/initial_amp
                attenuation[(w, h)] = ratios[1]
                phase = np.angle(np.exp(1j*(np.arctan2(-result.v/omega, result.x)-
                                      np.arctan2(-analytic_v/omega, analytic_x))))
                radii = [float(max(abs(np.linalg.eigvals(direct_matrices(o, mass, h, w)[0]))))
                         for o, mass in zip(omega, m)]
                self.assertLessEqual(max(radii), 1.+1e-12)
                if w == 1.:
                    initial_energy = .5*sum(m*v0*v0)+.5*x0@k@x0
                    self.assertAlmostEqual(result.ledger['physical_energy1_J'], initial_energy, delta=2e-9)
                emit('two_modes', dict(weight=w, dt_s=h, final_time_s=10.,
                    low_x_error_m=float(x_error), low_v_error_m_s=float(v_error),
                    modal_attenuation=ratios.tolist(), wrapped_modal_phase_error_rad=phase.tolist(),
                    phase_diagnostic_only_below_amplitude_ratio=1e-12,
                    amplification_spectral_radius=radii,
                    final_physical_energy_J=result.ledger['physical_energy1_J'],
                    final_algorithmic_energy_J=result.ledger['algorithmic_energy1_J'], **evidence))
            observed = [errors[i]/errors[i+1] for i in range(2)]
            self.assertTrue(all(3. <= ratio <= 5. for ratio in observed), observed)
            emit('low_mode_refinement', dict(weight=w, joint_errors=errors, refinement_ratios=observed))
        for h in (.1, .05, .025):
            self.assertLess(attenuation[(.9, h)], attenuation[(1., h)])
        emit('two_mode_cost', dict(wall_seconds=time.perf_counter()-start,
                                   complete_steps=1400, physical_run=False))

    def test_02_free_and_ramp_force(self):
        m, k = np.array([2.]), np.zeros((1, 1))
        for w in (.9, 1.):
            for h in (.1, .05, .025):
                free, _ = self.integrate(m, k, [.3], [-.2], lambda t: np.zeros(1), h, w)
                np.testing.assert_allclose(free.x, [-1.7], atol=2e-12, rtol=0)
                np.testing.assert_allclose(free.v, [-.2], atol=2e-12, rtol=0)
                ramp, evidence = self.integrate(m, k, [.3], [-.2],
                                               lambda t: np.array([.7+.2*t]), h, w)
                exact_x = .3-.2*10.+.7/4.*100.+.2/12.*1000.
                exact_v = -.2+.7/2.*10.+.2/4.*100.
                errors = [float(abs(ramp.x[0]-exact_x)), float(abs(ramp.v[0]-exact_v))]
                self.assertLessEqual(max(errors), .6*h*h)
                emit('free_and_ramp_load', dict(weight=w, dt_s=h,
                                               ramp_errors_x_m_v_m_s=errors, **evidence))

    def test_03_moving_spring_load_and_exact_path_work(self):
        m, k = np.array([2.]), np.array([[2.]])
        frequency, target_amplitude = .4, .2
        particular = target_amplitude/(1.-frequency**2)
        def analytic(t):
            u = .3*np.cos(t)+(-.2-particular*frequency)*np.sin(t)+particular*np.sin(frequency*t)
            v = -.3*np.sin(t)+(-.2-particular*frequency)*np.cos(t)+particular*frequency*np.cos(frequency*t)
            return u, v
        exact_x, exact_v = analytic(10.)
        exact_energy = exact_x**2+exact_v**2
        for w in (.9, 1.):
            for h in (.1, .05, .025):
                result, evidence = self.integrate(m, k, [.3], [-.2],
                    lambda t: np.array([2.*target_amplitude*np.sin(frequency*t)]), h, w)
                errors = [float(abs(result.x[0]-exact_x)), float(abs(result.v[0]-exact_v))]
                self.assertLessEqual(max(errors), .6*h*h)
                energy_error = abs(result.ledger['physical_energy1_J']-exact_energy)
                work_error = abs(evidence['conditional_cubic_path_external_work_J']-(exact_energy-.13))
                self.assertLessEqual(energy_error, .8*h*h)
                self.assertLessEqual(work_error, .8*h*h)
                emit('moving_spring_load', dict(weight=w, dt_s=h,
                    errors_x_m_v_m_s=errors, physical_energy_error_J=energy_error,
                    conditional_path_work_error_J=work_error, **evidence))
        # Four-point Gaussian integration is independent of the closed form;
        # it exactly integrates the registered linear-force/quadratic-velocity product.
        initial_x, initial_v, initial_a = np.array([.3]), np.array([-.2]), np.array([-.3])
        result = tested.step(m, k, initial_x, initial_v, initial_a, [0.], [.4], .1)
        nodes, weights = np.polynomial.legendre.leggauss(4)
        integral = 0.
        for s, quadrature_weight in zip((nodes+1.)/2., weights/2.):
            velocity = (initial_v+s*(result.v-initial_v)+6.*s*(1.-s)*
                        (result.displacement_increment/.1-(initial_v+result.v)/2.))
            integral += .1*quadrature_weight*float((s*np.array([.4]))@velocity)
        self.assertAlmostEqual(integral, result.ledger['exact_linear_load_work_J'], delta=2e-12)

    def test_04_restart_detachment_and_explicit_history_defect(self):
        m, k = np.array([2., 3.]), np.array([[7., -2.], [-2., 5.]])
        x, v = np.array([.3, -.2]), np.array([-.4, .1])
        f0, f1 = np.array([.1, .4]), np.array([.2, .3])
        a = tested.initialize_acceleration(m, k, x, f0)
        args = [m, k, x, v, a, f0, f1]
        pristine = copy.deepcopy(args)
        first = tested.step(*args, .05)
        for before, after in zip(pristine, args):
            np.testing.assert_array_equal(before, after)
        checkpoint = copy.deepcopy(first)
        live = tested.step(m, k, first.x, first.v, first.acc, f1, f0, .05, .9,
                           first.displacement_increment, first.load_increment)
        resumed = tested.step(m, k, checkpoint.x, checkpoint.v, checkpoint.acc, f1, f0, .05, .9,
                              checkpoint.displacement_increment, checkpoint.load_increment)
        for name in ('x', 'v', 'acc', 'displacement_increment', 'load_increment'):
            np.testing.assert_array_equal(getattr(live, name), getattr(resumed, name))
        first.x[:] = 8.
        np.testing.assert_array_equal(x, pristine[2])
        inconsistent = tested.step(m, k, x, v, a+1., f0, f1, .05)
        self.assertGreater(inconsistent.ledger['acceleration_history_defect_max_N'], 1.)
        self.assertGreater(abs(inconsistent.ledger['acceleration_history_work_J']), 1e-6)
        self.check_ledger(m, k, (x, v, a+1.), inconsistent, f0, f1, .05, .9,
                          np.zeros(2), np.zeros(2))
        emit('restart_and_bad_history', dict(restart='BITWISE_EQUAL',
            bad_history_defect_N=inconsistent.ledger['acceleration_history_defect_N'].tolist(),
            bad_history_signed_work_J=inconsistent.ledger['acceleration_history_work_J'],
            bad_history_identity_residual_J=inconsistent.ledger['identity_residual_J']))

    def test_05_reject_invalid_inputs(self):
        arguments = dict(mass=[2.], stiffness=[[2.]], x0=[.3], v0=[-.2],
                         acc0=[-.3], f0=[0.], f1=[0.], dt=.1)
        for key, value in [('mass', [0.]), ('stiffness', [[-1.]]), ('dt', 0.),
                           ('weight', .5), ('x0', [float('nan')]), ('acc0', [1., 2.])]:
            bad = dict(arguments, **{key: value})
            with self.assertRaises(ValueError, msg=key):
                tested.step(**bad)


if __name__ == '__main__':
    unittest.main(verbosity=2)
