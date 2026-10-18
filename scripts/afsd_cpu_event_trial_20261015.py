"""Explicit CPU one-sided event trial; original CPU operators are inherited.

Only pre-registered start-blocked incoming pairs may be held on their old side
while solving the bounded event-time bracket. The same selected i/j/g enter
all density, material, force, tangent and heat operators. This is an explicit
bounded activation-delay discretization, not an accepted endpoint filter fix.
"""
import ast
import inspect
import textwrap
import time
import numpy as np
import afsd_deposition_remote_setup_20261007 as parent
from afsd_deposition_remote_setup_20261007 import (operators, _geometry, filter_visible_support,
    update, pressure, internal_force, edge_fingerprint, SPATIAL_CLASSIFICATION,
    require_visible_support_numerical_domain, require_remote_filter_consistency,
    advance_implicit_thermal)
from afsd_event_graph_20261010 import restrict_activation_graph

ROOT, CLASSIFICATION = parent.ROOT, parent.CLASSIFICATION
DEFAULT_CONFIG, NUMERICAL_PROTOCOL = parent.DEFAULT_CONFIG, parent.NUMERICAL_PROTOCOL
COMPUTE_PROTOCOL = 'ORIGINAL_CPU_THERMOMECHANICAL_PROVIDERS_20261015_V1'
require, identity_digest = parent.require, parent.identity_digest
validate_state, load_configuration = parent.validate_state, parent.load_configuration
advance_candidate = parent.advance_candidate
EVENT_PROTOCOL = 'COUPLED_BOUNDED_ACTIVATION_EVENT_20261010_V1'


class ProcessTrial(parent.ProcessTrial):
    def __init__(self, pa, state, dt, incoming_pid_pairs=()):
        super().__init__(pa, state, dt)
        self.incoming_pid_pairs = tuple(tuple(v for v in pair) for pair in incoming_pid_pairs)

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
        raw_i, raw_j, raw_g = i, j, g
        numerical_domain=require_visible_support_numerical_domain(x1,i,j,
            self.motion.pose(self.state['time_s']+dt),p,pa.particle_id)
        graph = filter_visible_support(x1, i, j, g,
            self.motion.pose(self.state['time_s']+dt), p, pa.particle_id)
        require_remote_filter_consistency(numerical_domain, graph)
        if self.incoming_pid_pairs:
            graph = restrict_activation_graph(self.x0, x1, raw_i, raw_j, raw_g,
                self.motion.pose(self.state['time_s']), self.motion.pose(self.state['time_s']+dt),
                p, pa.particle_id, self.incoming_pid_pairs)
        i, j, g = graph['i'], graph['j'], graph['g']
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
        c = self.contact.evaluate(self.x0, x1, self.motion.pose(self.state['time_s']),
            self.motion.pose(self.state['time_s']+dt), dt, self.state['contact_history'])
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
        result['graph_fingerprint']=edge_fingerprint(i,j,pa.particle_id)
        result.update(visibility_graph=graph, raw_i=raw_i, raw_j=raw_j, raw_g=raw_g,
                      spatial_classification=SPATIAL_CLASSIFICATION,visibility_numerical_domain=numerical_domain)
        result['wall_seconds'] = time.perf_counter()-start
        return result


def operator_reuse_audit():
    expected = ast.parse(textwrap.dedent(inspect.getsource(parent.ProcessTrial.evaluate)))
    actual = ast.parse(textwrap.dedent(inspect.getsource(ProcessTrial.evaluate)))
    extra = [n for n in actual.body[0].body if isinstance(n, ast.If)
             and isinstance(n.test, ast.Attribute) and n.test.attr == 'incoming_pid_pairs']
    require(len(extra) == 1, 'exactly one explicit common-graph selection extension required')
    actual.body[0].body.remove(extra[0])
    require(ast.dump(actual) == ast.dump(expected), 'physical equations changed beyond common event graph')
    require(ProcessTrial.tangent is parent.ProcessTrial.tangent and advance_candidate is parent.advance_candidate,
            'event trial must preserve tangent and fresh physical advance')
    for name in ('operators', '_geometry', 'filter_visible_support', 'update',
                 'pressure', 'internal_force', 'edge_fingerprint', 'advance_implicit_thermal',
                 'require_visible_support_numerical_domain', 'require_remote_filter_consistency'):
        require(globals()[name] is getattr(parent, name), 'original CPU provider changed: '+name)
    return dict(status='PASS_CPU_EXPLICIT_ONE_SIDED_EVENT_GRAPH_ADAPTATION', event_protocol=EVENT_PROTOCOL,
        physical_equations_unchanged=True, common_graph_selection_extensions=1,
        runtime_original_global_rebinding=False, inherited=parent.operator_reuse_audit())
