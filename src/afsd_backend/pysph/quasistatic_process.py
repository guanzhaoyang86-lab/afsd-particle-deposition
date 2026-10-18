"""Independent hot-start, end-point quasi-static research model.

No inertial checkpoint may be passed off as this model's initial state.  The
frozen tangent is a search direction only; all gates use the complete trial.
"""
import time
import numpy as np
from .particles import create_geometry
from .mechanics import (positions, operators, internal_force, bulk_modulus,
                        pressure, volumetric_energy)
from .material import FKS, update
from .implicit_thermal import advance_implicit_thermal

CLASSIFICATION = 'HOT_START_QUASISTATIC_ELASTOPLASTIC_AFSD_APPROXIMATE'


def pose(p, t):
    a = p['omega']*t
    c, s = np.cos(a), np.sin(a)
    return dict(origin=np.array([p['traverse']*t, 0., p['layer']]),
                rotation=np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]]))


def independent_initial_state(p):
    """Create an unloaded geometry; no previous dynamic state is consumed."""
    if p['bore_width'] != 0.002504:
        raise ValueError('this registered initial geometry requires 2.504 mm bore')
    pa = create_geometry(p)
    # create_geometry also serves old dynamics and presets actuator velocity.
    # This model initializes at rest before any prescribed increment is applied.
    pa.u[:] = 0.; pa.v[:] = 0.; pa.w[:] = 0.
    free = np.ones((len(pa.x), 3), dtype=bool)
    bottom = (pa.component_id == 0) & (pa.z < -p['substrate_depth']+p['dx'])
    free[bottom] = False
    free[pa.actuator.astype(bool), 2] = False
    return pa, dict(time_s=0., step=0, classification=CLASSIFICATION,
                    contact_history=np.zeros((len(pa.x), 8, 3)), free=free,
                    crossing_ids=[], independent_initialization=True)


def _geometry(x, displacement, mass, rho, i, j, g):
    """Rebuild the volume-dependent correction on the same current neighbors."""
    volume = mass/rho
    moment = np.zeros((len(x), 3, 3)); raw = np.zeros_like(moment)
    pair_m = (x[j]-x[i])[:, :, None]*g[:, None, :]
    pair_u = (displacement[j]-displacement[i])[:, :, None]*g[:, None, :]
    for target, neighbor in ((i, j), (j, i)):
        np.add.at(moment, target, volume[neighbor, None, None]*pair_m)
        np.add.at(raw, target, volume[neighbor, None, None]*pair_u)
    ev = np.linalg.eigvalsh(moment)
    if np.any(ev[:, 0] <= 0) or np.any(ev[:, -1] > 1e6*ev[:, 0]):
        raise ValueError('singular/ill-conditioned current 3D neighborhood')
    correction = np.linalg.inv(moment)
    return raw@correction, correction


class QuasistaticTrial:
    def __init__(self, pa, state, p):
        from .quasistatic_contact import QuasistaticContact
        if state.get('classification') != CLASSIFICATION:
            raise ValueError('not a quasi-static history')
        self.pa, self.state, self.p = pa, state, dict(p)
        self.x0 = positions(pa)
        self.s0 = pa.S.reshape(-1, 3, 3).copy()
        self.bulk = bulk_modulus(p['shear'], p['poisson'])
        self.law = FKS(**{k: p[k] for k in FKS.__dataclass_fields__})
        self.contact = QuasistaticContact(tool_radius=p['tool_radius'],
            bore_width=p['bore_width'], tool_height=p['tool_top']-p['layer'],
            particle_radius=p['dx']/2, normal_stiffness=p['normal_stiffness'],
            tangential_stiffness=p['normal_stiffness'], friction=p['friction'],
            base_z=-p['substrate_depth'])

    def prescribed(self, dt):
        u = np.zeros_like(self.x0)
        u[self.pa.actuator.astype(bool), 2] = -self.p['feed']*dt
        return u

    def evaluate(self, displacement, dt):
        start = time.perf_counter()
        p, pa = self.p, self.pa
        u = np.array(displacement, dtype=float, copy=True)
        if u.shape != self.x0.shape or not np.isfinite(u).all() or dt <= 0:
            raise ValueError('finite displacement and positive interval required')
        if np.max(np.linalg.norm(u, axis=1)) > .1*p['dx']:
            raise ValueError('increment exceeds registered 0.1 dx domain')
        fixed = ~self.state['free']
        if not np.allclose(u[fixed], self.prescribed(dt)[fixed], atol=1e-16, rtol=0):
            raise ValueError('prescribed physical axes violated')
        x1 = self.x0+u
        cp = pa.extract_particles(np.arange(len(u), dtype=np.uint32))
        cp.x[:], cp.y[:], cp.z[:] = x1.T
        cp.u[:], cp.v[:], cp.w[:] = (u/dt).T
        i, j, g, _, _ = operators(cp)  # exactly one real 3D NNPS per trial
        rho = pa.rho.copy()
        for density_iteration in range(40):
            H, B = _geometry(x1, u, pa.m, rho, i, j, g)
            denominator = 1+np.trace(H, axis1=1, axis2=2)
            if np.any(denominator <= 0):
                raise ValueError('nonpositive BE density denominator')
            target = pa.rho/denominator
            density_residual = float(np.max(np.abs(target-rho)/p['rho0']))
            if density_residual <= 1e-10:
                break
            rho = target
        else:
            raise ValueError('density/correction fixed point did not converge')
        if np.max(np.abs(rho/p['rho0']-1)) > p['density_limit']:
            raise ValueError('registered density domain exceeded')
        cp.rho[:] = rho
        L, volume = H/dt, pa.m/rho
        c = self.contact.evaluate(self.x0, x1, pose(p, self.state['time_s']),
            pose(p, self.state['time_s']+dt), dt, self.state['contact_history'])
        # Finite patches use the same reference area each; overlap multiplicity
        # is an explicit approximation, consistently shared by mechanics/heat.
        active = np.asarray(c['active'], bool)
        conductance = active.astype(float)*p['contact_conductance']*p['dx']**2
        bath_T = np.full(8, p['tool_temperature']); bath_T[-1] = p['base_temperature']
        wall = conductance.sum(axis=1)
        bath = conductance@bath_T
        T = pa.temperature.copy()
        for thermal_iteration in range(30):
            material = update(self.s0, pa.ep, L, .5*(pa.temperature+T), dt, p['shear'], self.law)
            S, ep, wp, return_loss, rate = material
            thermal = advance_implicit_thermal(pa.m, rho, pa.temperature,
                p['conductivity'], x1, i, j, g, dt, p['cp'], volume*wp,
                np.asarray(c['friction_dissipation_J']).sum(axis=1),
                np.zeros(len(u)), wall, bath, p['plastic_heat_fraction'],
                p['heat_partition'], p['room'], p['melt'])
            thermal_residual = float(np.max(np.abs(thermal['temperature']-T)))
            if thermal_residual <= 1e-8:
                break
            T = thermal['temperature']
        else:
            raise ValueError('thermal/material fixed point did not converge')
        # Keep the actual feedback temperature that produced S; discrepancy
        # from the final heat solution is bounded above, not silently hidden.
        sigma = S-pressure(rho, p['rho0'], self.bulk)[:, None, None]*np.eye(3)
        fint = internal_force(cp, sigma, i, j, g, correction=B)
        force = fint+np.asarray(c['force_N'])
        scale = max(1., float(np.max(np.abs(fint))), float(np.max(np.abs(c['force_N']))))
        limit = 1e-5+1e-7*scale
        free_residual = float(np.max(np.abs(force[self.state['free']])))
        virtual = float(np.sum(fint*u)+np.sum(volume[:, None, None]*sigma*H))
        work = float(np.sum(volume[:, None, None]*sigma*H))
        if abs(virtual) > 1e-12+1e-10*abs(work):
            raise ValueError('paired mechanical virtual work failure')
        result = dict(displacement=u, position=x1, density=rho, gradient=L,
            correction=B, volume=volume, i=i, j=j, g=g, contact=c,
            material=material, temperature=thermal['temperature'],
            material_temperature=.5*(pa.temperature+T), thermal=thermal,
            stress=S, ep=ep, rate=rate, force=force, internal_force=fint,
            force_residual_N=free_residual, force_limit_N=limit,
            scaled_residual=free_residual/limit, density_residual=density_residual,
            thermal_residual_K=thermal_residual, virtual_work_residual_J=virtual,
            material_endpoint_work_J=work, density_iterations=density_iteration+1,
            thermal_iterations=thermal_iteration+1, dt_s=dt)
        result['wall_seconds'] = time.perf_counter()-start
        return result

    def tangent(self, result):
        """Neighbor-local assembly + local contact secant; NOT full Jacobian."""
        from .quasistatic_tangent import (gradient_map, material_tangent,
                                          assemble_frozen_force_tangent)
        p, r = self.p, result
        dt = r['dt_s']
        D = gradient_map(r['volume'], r['i'], r['j'], r['g'], r['correction'])
        material = material_tangent(self.s0, self.pa.ep, r['gradient'],
            r['material_temperature'], dt, p['shear'], self.law, material=r['material'])
        J = assemble_frozen_force_tangent(D, r['volume'], material['tangent_Pa_s'],
            self.pa.rho, p['rho0'], self.bulk, dt, r['displacement'])
        blocks = np.zeros((len(self.x0), 3, 3))
        # Batch perturb independent particle contacts, six contact calls total.
        eps = 1e-9
        for axis in range(3):
            plus = r['position'].copy(); minus = plus.copy()
            plus[:, axis] += eps; minus[:, axis] -= eps
            args = (pose(p, self.state['time_s']), pose(p, self.state['time_s']+dt),
                    dt, self.state['contact_history'])
            a = self.contact.evaluate(self.x0, plus, *args, check_path=False)
            b = self.contact.evaluate(self.x0, minus, *args, check_path=False)
            blocks[:, :, axis] = (a['force_N']-b['force_N'])/(2*eps)
        for particle, block in enumerate(blocks):
            rows = slice(3*particle, 3*particle+3)
            J[rows, rows] += block
        return J


def newton_direction(trial, result):
    """No mass shift, pseudoinverse, or artificial constraints on free modes."""
    free = trial.state['free'].ravel()
    J = trial.tangent(result)
    A = J[np.ix_(free, free)]
    rhs = -result['force'].ravel()[free]
    delta = np.zeros(free.shape)
    delta[free] = np.linalg.solve(A, rhs)
    backward_error = float(np.max(np.abs(A@delta[free]-rhs)))
    if backward_error > 1e-9*max(1., float(np.max(np.abs(rhs)))):
        raise ValueError('quasi-static linear solve backward error')
    return delta.reshape(trial.x0.shape), dict(
        linear_backward_error_N=backward_error, nnz=int(np.count_nonzero(A)),
        coercivity='NOT_CERTIFIED_BY_LU; INITIAL_ELASTIC_RANK_AUDIT_REQUIRED', inertia=False)


def energy_diagnostics(trial, r):
    """Raw mechanical storage and work ledger, no residual compensation."""
    from .quasistatic_energy import material_energy_breakdown
    p, pa, c = trial.p, trial.pa, r['contact']
    V0, V1 = pa.m/pa.rho, r['volume']
    U0 = np.sum(V0*np.sum(trial.s0**2, axis=(1, 2))/(4*p['shear']))
    U1 = np.sum(V1*np.sum(r['stress']**2, axis=(1, 2))/(4*p['shear']))
    bulk0 = np.sum(pa.m*volumetric_energy(pa.rho, p['rho0'], trial.bulk))
    bulk1 = np.sum(pa.m*volumetric_energy(r['density'], p['rho0'], trial.bulk))
    plastic = float(np.sum(V1*r['material'][2]))
    returned = float(np.sum(V1*r['material'][3]))
    # This raw diagnostic retains the known hypoelastic/BE work-storage gap.
    material_residual = r['material_endpoint_work_J']-(U1-U0+bulk1-bulk0)-plastic
    contact_normal_error = float(np.sum(c['normal_endpoint_work_error_J']))
    reaction = -r['force']*(~trial.state['free'])
    prescribed_work = float(np.sum(reaction*r['displacement']))
    tool_work = float(np.sum(c['wall_work_J']))
    normal_change = float(np.sum(c['normal_potential1_J']-c['normal_potential0_J']))
    tangent_change = float(np.sum(c['tangential_potential1_J']-c['tangential_potential0_J']))
    slide = float(np.sum(c['friction_dissipation_J']))
    friction_algorithm = float(np.sum(c['algorithmic_dissipation_J']))
    release = float(np.sum(c['release_energy_J']))
    global_raw = prescribed_work+tool_work-(U1-U0+bulk1-bulk0+normal_change+tangent_change+plastic+slide)
    scale = sum(abs(v) for v in (prescribed_work, tool_work, U1-U0,
                bulk1-bulk0, normal_change, tangent_change, plastic, slide))
    breakdown = material_energy_breakdown(trial.s0, r['stress'], pa.rho,
        r['density'], pa.m, r['gradient'], r['dt_s'], p['shear'], trial.bulk,
        p['rho0'], r['material'][2], r['material'][3])
    return dict(material_storage_increment_J=float(U1-U0+bulk1-bulk0),
        plastic_work_J=plastic, material_return_algorithmic_J=returned,
        material_work_storage_raw_residual_J=float(material_residual),
        contact_normal_endpoint_error_J=contact_normal_error,
        contact_tangential_algorithmic_J=friction_algorithm,
        contact_release_energy_J=release, normal_storage_increment_J=normal_change,
        tangential_storage_increment_J=tangent_change, friction_sliding_J=slide,
        prescribed_boundary_work_J=prescribed_work, tool_work_J=tool_work,
        global_raw_residual_J=float(global_raw), absolute_work_scale_J=scale,
        registered_energy_limit_J=1e-7+.05*scale,
        material_breakdown=breakdown,
        contact_normal_error_absolute_sum_J=float(np.abs(c['normal_endpoint_work_error_J']).sum()),
        prescribed_reaction_N=reaction)


def accepted_copy(trial, r):
    """Only explicit acceptance produces a detached committable state."""
    if not r['contact']['sweep'].get('certified', False):
        raise ValueError('uncertified swept geometry cannot be accepted')
    if r['scaled_residual'] > 1:
        raise ValueError('unconverged free forces')
    energy = energy_diagnostics(trial, r)
    if abs(energy['global_raw_residual_J']) > energy['registered_energy_limit_J']:
        raise ValueError('registered raw energy diagnostic failed')
    penetration = np.asarray(r['contact']['penetration_m'])
    if np.max(penetration) > .1*trial.p['dx']:
        raise ValueError('registered penalty penetration exceeded')
    pa = trial.pa.extract_particles(np.arange(len(trial.pa.x), dtype=np.uint32))
    pa.x[:], pa.y[:], pa.z[:] = r['position'].T
    pa.u[:], pa.v[:], pa.w[:] = (r['displacement']/r['dt_s']).T
    pa.rho[:] = r['density']; pa.S[:] = r['stress'].ravel(); pa.ep[:] = r['ep']
    pa.D[:] = (.5*(r['gradient']+r['gradient'].swapaxes(1, 2))).ravel()
    pa.temperature[:] = r['temperature']; pa.strain_rate[:] = r['rate']
    pa.ep_rate[:] = (r['ep']-trial.pa.ep)/r['dt_s']
    pa.plastic_work[:] += r['volume']*r['material'][2]
    pa.return_dissipation[:] += r['volume']*r['material'][3]
    if abs(pa.m.sum()/trial.pa.m.sum()-1) > 1e-12:
        raise ValueError('mass changed')
    state = dict(trial.state, time_s=trial.state['time_s']+r['dt_s'],
        step=trial.state['step']+1, contact_history=r['contact']['history'].copy())
    # A crossing requires a bar particle previously ABOVE the shoulder;
    # initial material below it cannot create a spurious deposition claim.
    above = trial.x0[:, 2] >= trial.p['layer']
    below = r['position'][:, 2] < trial.p['layer']
    cross = (pa.component_id == 1) & above & below
    ids = set(state['crossing_ids']); ids.update(int(v) for v in pa.particle_id[cross])
    state['crossing_ids'] = sorted(ids)
    state['energy'] = energy
    return pa, state
