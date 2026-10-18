"""Frozen 3D elastic structural preconditioner, NOT a process integrator.

The physical particle masses and elastic moduli are unchanged.  This tangent
omits prestress, evolving geometry, FKS/plasticity, thermal and contact terms.
No accepted history, particle position, velocity or force is committed here.
"""
import numpy as np

CLASSIFICATION = 'FROZEN_CONFIGURATION_ELASTIC_PRECONDITIONER_NOT_PROCESS_SOLVER'


def _real(value, name, shape=None):
    raw = np.asarray(value)
    if raw.dtype.kind not in 'iuf' or (shape is not None and raw.shape != shape):
        raise ValueError('invalid real shape/type: '+name)
    out = np.array(raw, dtype=float, copy=True)
    if not np.isfinite(out).all():
        raise ValueError('nonfinite '+name)
    return out


def _positive(value, name):
    result = _real(value, name, ())
    if result <= 0:
        raise ValueError('positive '+name+' required')
    return float(result)


def _readonly(value):
    result = np.array(value, copy=True)
    result.flags.writeable = False
    return result


class ConstrainedResponse:
    """A=M+h^2 K/4; impulse -> endpoint velocity increment [m/s].

``mass_per_dof`` repeats each actual particle mass on its three axes.
``free`` selects prescribed versus free axes; it never removes an unrequested
mode. Reactions are impulses [N s] in the ORIGINAL full DOF inventory.
This small dense component solves its exact A, not a pseudoinverse.
"""
    def __init__(self, mass_per_dof, stiffness, dt, free):
        m = _real(mass_per_dof, 'mass_per_dof')
        if m.ndim != 1 or not len(m) or np.any(m <= 0):
            raise ValueError('positive nonempty mass_per_dof required')
        n = len(m)
        k = _real(stiffness, 'stiffness', (n, n))
        h = _positive(dt, 'dt')
        mask = np.asarray(free)
        if (mask.shape != (n,) or mask.dtype.kind not in 'biuf' or
                not np.isfinite(mask).all() or not np.all((mask == 0) | (mask == 1))):
            raise ValueError('free must be an exact boolean/zero-one DOF vector')
        scale = max(float(np.linalg.norm(k, 2)), np.finfo(float).tiny)
        asymmetry = float(np.linalg.norm(k-k.T, 2))
        if asymmetry > 1e-12*scale:
            raise ValueError('nonsymmetric stiffness')
        eigenvalues = np.linalg.eigvalsh(k)
        if eigenvalues[0] < -1e-12*scale:
            raise ValueError('negative elastic stiffness; not clipped')
        matrix = np.diag(m)+h*h*k/4.
        if not np.isfinite(matrix).all():
            raise ValueError('nonfinite inertial structural matrix')
        selected = mask.astype(bool)
        ff = matrix[np.ix_(selected, selected)]
        if len(ff):
            try:
                np.linalg.cholesky(ff)
            except np.linalg.LinAlgError as error:
                raise ValueError('free inertial matrix not positive definite; no pinv') from error
        self.mass_per_dof = _readonly(m)
        self.stiffness = _readonly(k)
        self.matrix = _readonly(matrix)
        self.free = _readonly(selected)
        self.dt = h
        self.diagnostics = dict(classification=CLASSIFICATION,
            production_authorized=False, full_tangent=False,
            stiffness_asymmetry_norm_N_m=asymmetry,
            stiffness_min_eigenvalue_N_m=float(eigenvalues[0]),
            stiffness_near_zero_eigenvalues=int(np.count_nonzero(np.abs(eigenvalues) <= 1e-10*scale)),
            stiffness_zero_mode_relative_cutoff=1e-10,
            free_matrix_condition=(float(np.linalg.cond(ff)) if len(ff) else None),
            free_dofs=int(selected.sum()), prescribed_dofs=int((~selected).sum()),
            inertia_does_not_certify_stiffness_coercivity=True)

    def response(self, impulse, prescribed_delta_velocity=None):
        """Solve A dv = impulse + reaction, with caller-specified constrained dv.

Prescribed input contains only constrained-axis increments; nonzero values on
free axes reject rather than silently overriding an independently solved DOF.
"""
        n = len(self.mass_per_dof)
        rhs = _real(impulse, 'impulse', (n,))
        dv = (np.zeros(n) if prescribed_delta_velocity is None else
              _real(prescribed_delta_velocity, 'prescribed_delta_velocity', (n,)))
        if np.any(dv[self.free] != 0):
            raise ValueError('prescribed increment on a free axis')
        free, constrained = self.free, ~self.free
        if np.any(free):
            aff = self.matrix[np.ix_(free, free)]
            shifted = rhs[free]-self.matrix[np.ix_(free, constrained)]@dv[constrained]
            dv[free] = np.linalg.solve(aff, shifted)
        reaction = self.matrix@dv-rhs
        if not np.isfinite(dv).all() or not np.isfinite(reaction).all():
            raise ValueError('nonfinite structural response; no clipping')
        return dict(delta_velocity_m_s=dv, reaction_impulse_Ns=reaction,
            free_equation_residual_Ns=reaction[free].copy(),
            classification=CLASSIFICATION, accepted=False, full_tangent=False)

    def contact_response(self, contact_map, at='endpoint'):
        """G Z Aff^-1 Z^T G^T [1/kg], or half this for midpoint velocity.

Columns of G^T distribute contact impulse to the same power-conjugate DOFs.
Geometry/normal derivatives and nonzero prescribed velocities are NOT encoded
in W. The latter's affine contribution is obtained with response separately.
"""
        n = len(self.mass_per_dof)
        g = _real(contact_map, 'contact_map')
        if g.ndim != 2 or g.shape[1] != n:
            raise ValueError('contact_map requires rows by full DOFs')
        if at not in ('endpoint', 'midpoint'):
            raise ValueError('explicit endpoint or midpoint convention required')
        if np.any(self.free):
            gf = g[:, self.free]
            w = gf@np.linalg.solve(self.matrix[np.ix_(self.free, self.free)], gf.T)
        else:
            w = np.zeros((len(g), len(g)))
        result = w if at == 'endpoint' else .5*w
        if not np.isfinite(result).all():
            raise ValueError('nonfinite contact structural response')
        return result


class FrozenElastic:
    """Detached current-geometry 3D gradient and paired elastic stiffness."""
    def __init__(self, gradient, stiffness, volume, mass, x, particle_ids, shear, bulk, pairs):
        self.gradient_map = _readonly(gradient)
        self.stiffness = _readonly(stiffness)
        self.volume = _readonly(volume)
        self.mass = _readonly(mass)
        self.position = _readonly(x)
        self.particle_ids = _readonly(particle_ids)
        self.shear, self.bulk = shear, bulk
        self.pair_count = pairs
        self.classification = CLASSIFICATION

    def gradient(self, displacement):
        u = _real(displacement, 'displacement', self.position.shape)
        return (self.gradient_map@u.ravel()).reshape(-1, 3, 3)

    def stress(self, displacement):
        gradient = self.gradient(displacement)
        trace = np.trace(gradient, axis1=1, axis2=2)
        return (self.shear*(gradient+gradient.swapaxes(1, 2)) +
                (self.bulk-2*self.shear/3)*trace[:, None, None]*np.eye(3))

    def internal_force(self, displacement):
        u = _real(displacement, 'displacement', self.position.shape)
        return (-self.stiffness@u.ravel()).reshape(self.position.shape)

    def energy(self, displacement):
        u = _real(displacement, 'displacement', self.position.shape).ravel()
        return .5*float(u@self.stiffness@u)

    def response_operator(self, dt, free_axes):
        mask = np.asarray(free_axes)
        if mask.shape != self.position.shape:
            raise ValueError('free_axes requires particle by three axes')
        return ConstrainedResponse(np.repeat(self.mass, 3), self.stiffness, dt, mask.ravel())


def build_frozen_elastic(pa, shear, bulk):
    """Reuse actual PySPH NNPS, corrected B, m/rho; never mutate pa.

K is the derivative ONLY of frozen-gradient, unprestressed linear elastic
stress.  It is not the derivative of the full updated-configuration force.
"""
    from .mechanics import corrected_operators, positions
    shear = _positive(shear, 'shear'); bulk = _positive(bulk, 'bulk')
    n = pa.get_number_of_particles()
    if not n:
        raise ValueError('nonempty particle array required')
    private = pa.extract_particles(np.arange(n, dtype=np.uint32))
    if (not np.isfinite(private.m).all() or not np.isfinite(private.rho).all() or
            np.any(private.m <= 0) or np.any(private.rho <= 0)):
        raise ValueError('positive finite particle mass/density required')
    i, j, g, _, _, correction = corrected_operators(private)
    volume = private.m/private.rho
    gi = np.einsum('nji,nj->ni', correction[i], g)
    gj = np.einsum('nji,nj->ni', correction[j], g)
    gradient = np.zeros((n, 3, 3, 3*n))
    for component in range(3):
        for axis in range(3):
            np.add.at(gradient, (i, component, axis, 3*j+component), volume[j]*gi[:, axis])
            np.add.at(gradient, (i, component, axis, 3*i+component), -volume[j]*gi[:, axis])
            np.add.at(gradient, (j, component, axis, 3*j+component), volume[i]*gj[:, axis])
            np.add.at(gradient, (j, component, axis, 3*i+component), -volume[i]*gj[:, axis])
    d = gradient.reshape(n, 9, 3*n)
    elasticity = np.zeros((9, 9))
    for a in range(3):
        for b in range(3):
            elasticity[3*a+b, 3*a+b] += shear
            elasticity[3*a+b, 3*b+a] += shear
            if a == b:
                for c in range(3):
                    elasticity[3*a+b, 3*c+c] += bulk-2*shear/3
    stiffness = np.einsum('nik,nij,njl->kl', d, volume[:, None, None]*elasticity, d,
                          optimize=True)
    if not np.isfinite(d).all() or not np.isfinite(stiffness).all():
        raise ValueError('nonfinite frozen elastic operator')
    return FrozenElastic(d.reshape(9*n, 3*n), stiffness, volume, private.m,
                         positions(private), private.particle_id, shear, bulk, len(i))
