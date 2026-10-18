"""Independent frozen elastic-density block tests; not AFSD process steps."""
from pathlib import Path
import hashlib
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from afsd_backend.pysph.stage_continuity_preconditioner import ContinuityStagePreconditioner
from afsd_backend.pysph.structural_mobility import ConstrainedResponse
from afsd_backend.pysph.stage_iteration import FIELDS
from afsd_backend.pysph.particles import make_array
from afsd_backend.pysph.mechanics import corrected_operators, internal_force


def compression(constrained=False):
    # Literal m=4,V=2,Kb=3,E=(1,0,0), K=E^T V Kb E=diag(6,0,0).
    response = ConstrainedResponse([4., 4., 4.], np.diag([6., 0., 0.]), 2.,
                                   [not constrained, True, True])
    return ContinuityStagePreconditioner.from_response(response, np.array([17]),
        [[1., 0., 0.]], [2.], [5.], 3.)


def fields():
    return dict(position_m=np.zeros((1, 3)), velocity_m_s=np.zeros((1, 3)),
        density_kg_m3=np.array([5.]), temperature_K=np.array([400.]),
        mean_velocity_m_s=np.zeros((1, 3)))


def copied(y):
    return {key: value.copy() for key, value in y.items()}


def particle_hash(pa):
    h = hashlib.sha256()
    for name in sorted(pa.properties):
        h.update(name.encode()); h.update(pa.properties[name].get_npy_array().tobytes())
    return h.hexdigest()


class CompressionBlockTests(unittest.TestCase):
    def close(self, a, b):
        np.testing.assert_allclose(a, b, rtol=2e-12, atol=2e-13)

    def test_01_literal_two_by_two_pressure_continuity_schur(self):
        c = fields(); mapping = copied(c)
        mapping['mean_velocity_m_s'][0, 0] = .5; mapping['density_kg_m3'][0] += 2.
        result = compression().apply(c, mapping); delta = result['correction']; info = result['diagnostics']
        # Equations db - .15 dr=.5; 10 db+dr=2 => db=.32, dr=-1.2.
        self.close(delta['mean_velocity_m_s'][0], [.32, 0., 0.])
        self.close(delta['density_kg_m3'], [-1.2])
        self.close(info['inertial_rhs_Ns'], [2., 0., 0.])
        self.close(info['pressure_density_rhs_Ns'], [1.2, 0., 0.])
        self.close(info['rhs_Ns'], [3.2, 0., 0.])
        self.close(delta['position_m'], [[-.36, 0., 0.]])
        self.close(delta['velocity_m_s'], [[-.36, 0., 0.]])
        self.close(info['density_backsubstitution_kg_m3'], [-3.2])
        self.assertLessEqual(info['backward_error_inf'], 1e-12)

    def test_02_nonzero_constraint_density_backsubstitution_and_reaction(self):
        c = fields(); mapping = copied(c)
        mapping['mean_velocity_m_s'][0, 0] = .5; mapping['density_kg_m3'][0] += 2.
        got = compression(True).apply(c, mapping)
        self.close(got['correction']['mean_velocity_m_s'], [[.5, 0., 0.]])
        self.close(got['correction']['density_kg_m3'], [-3.])
        # A11=10; constrained .5 => reaction=5-3.2=1.8 Ns.
        self.close(got['diagnostics']['algebraic_reaction_impulse_Ns'], [1.8, 0., 0.])
        self.close(got['correction']['position_m'], np.zeros((1, 3)))

    def test_03_literal_coupled_root_in_one_correction(self):
        # The independent real map db=a+.15*r; r=q-10*b has root b=.32,r=-1.2.
        c = fields(); c['mean_velocity_m_s'][0, 0] = .7; c['density_kg_m3'][0] = 4.
        def mapping(y):
            result = copied(y)
            density_increment = y['density_kg_m3'][0]-5.
            vbar = .5+.15*density_increment
            result['mean_velocity_m_s'][0, 0] = vbar
            result['density_kg_m3'][0] = 7.-10*y['mean_velocity_m_s'][0, 0]
            result['position_m'][0, 0] = 2*vbar
            result['velocity_m_s'][0, 0] = 2*vbar
            return result
        root = compression().apply(c, mapping(c))['proposal']
        self.close(root['mean_velocity_m_s'], [[.32, 0., 0.]])
        self.close(root['density_kg_m3'], [3.8])
        for key in FIELDS:
            self.close(root[key], mapping(root)[key])

    def test_04_event_average_and_zero_residual_preserved(self):
        c = fields(); c['mean_velocity_m_s'][:, 0] = -1/18
        c['velocity_m_s'][:, 0] = -2/3; c['position_m'][:, 0] = -1/18
        result = compression().apply(c, copied(c))
        for key in FIELDS:
            np.testing.assert_array_equal(result['proposal'][key], c[key])
            np.testing.assert_array_equal(result['correction'][key], np.zeros_like(c[key]))
        self.assertNotEqual(result['proposal']['mean_velocity_m_s'][0, 0], 1/6)

    def test_05_no_false_zero_and_unmodified_kinematic_defect(self):
        c = fields(); mapping = copied(c)
        mapping['density_kg_m3'][0] += 3.; mapping['position_m'][0, 0] = .9
        mapping['velocity_m_s'][0, 0] = -.3; mapping['temperature_K'][0] += 2.
        result = compression().apply(c, mapping); delta = result['correction']
        b, r = delta['mean_velocity_m_s'][0, 0], delta['density_kg_m3'][0]
        self.close(b-.15*r, 0.); self.close(10*b+r, 3.)
        self.close(delta['position_m'][0, 0]-2*b, .9)
        self.close(delta['velocity_m_s'][0, 0]-2*b, -.3)
        np.testing.assert_array_equal(delta['temperature_K'], [2.])

    def test_06_no_density_clipping_or_metadata_alias(self):
        c = fields(); mapping = copied(c); mapping['mean_velocity_m_s'][0, 0] = 10.
        p = compression(); got = p.apply(c, mapping)
        self.assertLess(got['proposal']['density_kg_m3'][0], 0.)
        got['diagnostics']['density_reference_kg_m3'][:] = -1.
        np.testing.assert_array_equal(p.rho_reference, [5.])
        np.testing.assert_array_equal(c['density_kg_m3'], [5.])
        self.assertFalse(got['diagnostics']['accepted'])

    def test_07_dimensionally_transformed_compression_is_covariant(self):
        # Pure manufactured units conversion: L'=5L, t'=3t, m'=7m, T'=2T.
        c = fields(); mapping = copied(c)
        mapping['mean_velocity_m_s'][0, 0] = .5; mapping['density_kg_m3'][0] += 2.
        scale = dict(position_m=5., velocity_m_s=5/3, mean_velocity_m_s=5/3,
                     density_kg_m3=7/125, temperature_K=2.)
        response = ConstrainedResponse(np.full(3, 4.*7), np.diag([6., 0., 0.])*7/9,
                                       2.*3, np.ones(3, bool))
        p = ContinuityStagePreconditioner.from_response(response, np.array([17]),
            [[1/5, 0., 0.]], [2.*125], [5.*7/125], 3.*7/(5*9))
        base = compression().apply(c, mapping)
        transformed = p.apply({k: v*scale[k] for k, v in c.items()},
                              {k: v*scale[k] for k, v in mapping.items()})
        for key in FIELDS:
            self.close(transformed['correction'][key], base['correction'][key]*scale[key])


class ActualFrozenSPHBlockTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        xyz = np.array(np.meshgrid([-.0005, 0., .0005], [-.0005, 0., .0005],
                                   [-.0005, 0., .0005], indexing='ij')).reshape(3, -1).T
        cls.p = dict(dx=.0005, h_ratio=1.3, rho0=2700., width=.0015,
            marker_low=-.001, marker_high=.001, bar_length=.01, actuator_depth=.0005,
            temperature=400., traverse=.00212, feed=.00212, shear=26.3e9,
            poisson=.33, release_time=1.)
        cls.pa = make_array(xyz, np.zeros(27, int), cls.p)
        cls.before = particle_hash(cls.pa); cls.dt = 2e-7
        cls.pc = ContinuityStagePreconditioner(cls.pa, cls.p, cls.dt, 0.)
        cls.c = dict(position_m=xyz, velocity_m_s=np.zeros_like(xyz),
            mean_velocity_m_s=np.zeros_like(xyz), density_kg_m3=cls.pa.rho.copy(),
            temperature_K=cls.pa.temperature.copy())

    def independent_force_divergence(self, y):
        # Existing force path plus literal constitutive formula, NEVER new E/K.
        q = self.pa.extract_particles(np.arange(27, dtype=np.uint32))
        q.u[:], q.v[:], q.w[:] = y['mean_velocity_m_s'].T
        i, j, g, l, _, correction = corrected_operators(q)
        divergence = np.trace(l, axis1=1, axis2=2)
        gshear, nu, rho0 = self.p['shear'], self.p['poisson'], self.p['rho0']
        bulk = 2*gshear*(1+nu)/(3*(1-2*nu))
        midpoint_s = gshear*self.dt*(.5*(l+l.swapaxes(1, 2))-divergence[:, None, None]*np.eye(3)/3)
        midpoint_p = bulk*(y['density_kg_m3']-rho0)/(2*rho0)
        sigma = midpoint_s-midpoint_p[:, None, None]*np.eye(3)
        return internal_force(q, sigma, i, j, g, correction=correction), divergence

    def test_08_actual_three_dimensional_coupled_map_root(self):
        # Affine velocity root is specified independently; body forcing is
        # manufactured with the old operator, not the preconditioner matrix.
        a = np.array([[1.1, -.3, .2], [.5, -.7, .1], [-.2, .4, .6]])
        target = copied(self.c)
        target['mean_velocity_m_s'] = self.c['position_m']@a.T+np.array([.001, -.002, .003])
        target['density_kg_m3'] = np.full(27, self.p['rho0']*(1-self.dt*np.trace(a)))
        target['position_m'] = self.c['position_m']+self.dt*target['mean_velocity_m_s']
        target['velocity_m_s'] = 2*target['mean_velocity_m_s']
        fstar, _ = self.independent_force_divergence(target)
        forcing = target['mean_velocity_m_s']-.5*self.dt*fstar/self.pa.m[:, None]
        def real_map(y):
            force, divergence = self.independent_force_divergence(y)
            bar = forcing+.5*self.dt*force/self.pa.m[:, None]
            return dict(position_m=self.c['position_m']+self.dt*bar,
                mean_velocity_m_s=bar, velocity_m_s=2*bar,
                density_kg_m3=self.p['rho0']*(1-self.dt*divergence),
                temperature_K=self.c['temperature_K'].copy())
        first = self.pc.apply(self.c, real_map(self.c)); root = first['proposal']
        mapped = real_map(root)
        for key in FIELDS:
            scale = max(np.max(np.abs(target[key])), 1e-20)
            self.assertLessEqual(np.max(np.abs(root[key]-target[key]))/scale, 2e-11, key)
            self.assertLessEqual(np.max(np.abs(mapped[key]-root[key]))/scale, 2e-11, key)
        self.assertLessEqual(first['diagnostics']['backward_error_inf'], 1e-12)

    def test_09_six_rigid_modes_do_not_create_density_increment(self):
        for axis in np.eye(3):
            for mode in (np.broadcast_to(axis, (27, 3)), np.cross(axis, self.c['position_m'])):
                mapping = copied(self.c); mapping['mean_velocity_m_s'] = mode.copy()
                got = self.pc.apply(self.c, mapping)
                self.assertLessEqual(np.linalg.norm(got['correction']['mean_velocity_m_s']-mode),
                                     2e-11*np.linalg.norm(mode))
                # Dimensionless density difference, scaled by dt*rho*|v|/dx.
                scale = self.dt*self.p['rho0']*np.linalg.norm(mode)/self.p['dx']
                self.assertLessEqual(np.linalg.norm(got['correction']['density_kg_m3']), 2e-11*scale)

    def test_10_passive_labels_preserve_schur_and_source_history(self):
        q = self.pa.extract_particles(np.arange(27, dtype=np.uint32))
        q.tracer_type[:] = np.arange(27) % 3; q.wholebar_tag[:] = np.arange(27) % 4
        other = ContinuityStagePreconditioner(q, self.p, self.dt, 0.)
        mapping = copied(self.c); mapping['density_kg_m3'] += np.sin(np.arange(27))
        first, second = self.pc.apply(self.c, mapping), other.apply(self.c, mapping)
        for key in FIELDS:
            np.testing.assert_array_equal(first['correction'][key], second['correction'][key])
        self.assertEqual(self.before, particle_hash(self.pa))
        np.testing.assert_array_equal(other.response.mass_per_dof, np.repeat(q.m, 3))


if __name__ == '__main__':
    unittest.main(verbosity=2)
