"""Versioned visible-support variant of the fixed PySPH quasi-static model.

The evaluate body is a static adaptation of quasistatic_process.py: the only
physical change is the common visible graph before every volume/correction,
material-gradient, paired-force and thermal assembly. No frozen source is
patched or rebound at runtime. This is a new ENGINEERING_APPROXIMATE spatial
model, not retrospective acceptance of any old failure.
"""
import copy
import hashlib
import time
import numpy as np
import run_discrete_contact_process as core
from afsd_backend.pysph.discrete_process import DiscreteContactTrial, CLASSIFICATION as OLD_CLASSIFICATION
from afsd_backend.pysph.quasistatic_process import _geometry, pose
from afsd_backend.pysph.mechanics import operators, internal_force, pressure
from afsd_backend.pysph.material import update
from afsd_backend.pysph.implicit_thermal import advance_implicit_thermal
from batched_contact import BatchedChainRuleContact
from visible_support import filter_visible_support
from visible_support_numerics import require_visible_support_numerical_domain

CLASSIFICATION = 'HOT_QS_VISIBLE_SUPPORT_CHAIN_RULE_CONTACT_APPROXIMATE_V1'


def edge_fingerprint(i, j, particle_ids):
    ids=np.asarray(particle_ids)
    edges=np.sort(np.column_stack((ids[i],ids[j])).astype('<i8'),axis=1)
    if len(edges): edges=edges[np.lexsort((edges[:,1],edges[:,0]))]
    return hashlib.sha256(edges.tobytes()).hexdigest()


def explicit_transition(state, parent_sha256, *, diagnostic_only=False):
    """Copy all history; only an auditable spatial-model label is changed."""
    if state.get('classification') != OLD_CLASSIFICATION or len(parent_sha256)!=64:
        raise ValueError('explicit transition requires the qualified old spatial model and parent hash')
    if state.get('restart_allowed') is False and not diagnostic_only:
        raise ValueError('unaccepted reference may only be used for a named map diagnostic')
    out=copy.deepcopy(state)
    out.update(classification=CLASSIFICATION, spatial_model_transition=dict(
        from_classification=OLD_CLASSIFICATION, to_classification=CLASSIFICATION,
        parent_checkpoint_sha256=parent_sha256, time_s=state['time_s'], step=state['step'],
        fields_reset=False, interval_energy_reset=False, observer_reset=False,
        change='One common finite-tool-visible support graph before all spatial operators',
        diagnostic_only=bool(diagnostic_only),
        qualification='UNACCEPTED_FIXED_PATH_MAP_DIAGNOSTIC_ONLY' if diagnostic_only else
            'REQUIRES_PARENT_EQUIVALENCE_AND_INDEPENDENT_BENCHMARKS'))
    return out


def accepted_copy(trial, result):
    if trial.state['spatial_model_transition'].get('diagnostic_only'):
        raise ValueError('a diagnostic-only model state cannot advance an accepted process')
    return core.accepted_copy(trial, result)


class VisibleSupportTrial(DiscreteContactTrial):
    def __init__(self, pa, state, p):
        if state.get('classification')!=CLASSIFICATION or not state.get('spatial_model_transition'):
            raise ValueError('explicit visible-support model history required')
        adapter=dict(state, classification=OLD_CLASSIFICATION)
        super().__init__(pa, adapter, p)
        self.state=state
        self.contact=BatchedChainRuleContact(tool_radius=p['tool_radius'],bore_width=p['bore_width'],
            tool_height=p['tool_top']-p['layer'],particle_radius=p['dx']/2,
            normal_stiffness=p['normal_stiffness'],tangential_stiffness=p['normal_stiffness'],
            friction=p['friction'],base_z=-p['substrate_depth'])

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
            pose(p,self.state['time_s']+dt),p,pa.particle_id)
        graph = filter_visible_support(x1, i, j, g,
            pose(p, self.state['time_s']+dt), p, pa.particle_id)
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
        result['graph_fingerprint']=edge_fingerprint(i,j,pa.particle_id)
        result.update(visibility_graph=graph, raw_i=raw_i, raw_j=raw_j, raw_g=raw_g,
                      spatial_classification=CLASSIFICATION,visibility_numerical_domain=numerical_domain)
        result['wall_seconds'] = time.perf_counter()-start
        return result

    def tangent(self, result):
        graph=result['visibility_graph']
        for key in ('i','j','g'):
            if not np.array_equal(result[key],graph[key]):
                raise ValueError('tangent must use the exact residual visible graph')
        return super().tangent(result)
