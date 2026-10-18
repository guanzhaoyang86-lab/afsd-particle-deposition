"""Independent 918-particle ramped hot-QS diagnostic; no old restart adapter.

Only the boundary history changes. evaluate is statically copied from the
frozen VisibleSupportTrial with four pose calls replaced by self.motion.pose.
Its spatial CLASSIFICATION name is an alias retaining the original value.
The original QS tangent is statically copied with two pose calls replaced;
the pre-existing visible-graph guard and an explicit interval check precede it.
Its package-relative tangent import becomes an absolute import in scripts/.
Everything else, including 1e-9 contact FD, material floors, thermal and force
thresholds, is unchanged. static_adaptation_audit checks those AST boundaries.

advance_candidate deliberately does NOT call discrete_process.accepted_copy:
that helper hardcodes the old constant-speed pose in its strict observer.
Here fresh force -> original QS accepted_copy -> complete interval energy ->
strict observer all share one IntervalMotion. Detached states always remain
accepted=False/restart_allowed=False. A supervising diagnostic still owns
scalar/geometry/path/time-screen qualification; this helper is not publication.
"""
from pathlib import Path
import ast
import copy
import time
import numpy as np
import expanded_hot_domain as expanded
from ramped_hot_motion import RampedHotMotion,IntervalMotion,RAMP_DURATION_S,MOTION_CLASSIFICATION
from visible_quasistatic import (operators,_geometry,filter_visible_support,
    require_visible_support_numerical_domain,update,advance_implicit_thermal,
    pressure,internal_force,edge_fingerprint,CLASSIFICATION as SPATIAL_CLASSIFICATION)
from afsd_backend.pysph.mechanics import positions,bulk_modulus
from afsd_backend.pysph.material import FKS
from afsd_backend.pysph import quasistatic_process as original_qs
from afsd_backend.pysph.interval_acceptance import (initial_bore_crossing_state,
    validate_bore_crossing_state,validate_interval_energy_state,require_interval_energy,
    observe_square_bore_crossings)

ROOT=expanded.ROOT
CLASSIFICATION='INDEPENDENT_RAMPED_20MS_HOT_QS_VISIBLE_DIAGNOSTIC_V1'
SCOPE='NEW_918_PARTICLE_HOT_ZERO_TIME_SMOOTH_VELOCITY_STARTUP_ASSUMPTION_V1'
require=expanded.require


def source_hashes():
    result=expanded.source_hashes()
    for name in ('ramped_hot_motion.py','ramped_hot_trial.py'):
        path=ROOT/'scripts'/name;result[path.relative_to(ROOT).as_posix()]=expanded.sha(path)
    return result


def independent_initial_state():
    """Allocate a new t=0 918-particle factory state; no checkpoint is read."""
    p,nominal_provenance=expanded.parameters(.008)
    pa,state=expanded.independent_initial_state(p,nominal_provenance,variant='visible')
    require(len(pa.x)==918 and int(state['free'].sum())==1961,'registered independent 8 mm layout differs')
    history=RampedHotMotion(p,RAMP_DURATION_S)
    nominal=copy.deepcopy(nominal_provenance)
    provenance=dict(scope=SCOPE,nominal_parameter_reference=nominal,
        boundary_condition_history=history.to_dict(),
        changed_boundary_history_fields=['omega','traverse','feed'],
        nominal_values_are_not_instantaneous_startup_velocities=True,
        physical_qualification='NOT_RUN_NEW_STARTUP_BOUNDARY_ASSUMPTION',
        old_state_read=False,parent_checkpoint=None,mass_scale=1.,time_scale=1.,modulus_scale=1.)
    initialization=copy.deepcopy(state['independent_domain_initialization'])
    initialization.update(kind=SCOPE,id_namespace='independent_ramped_hot_8mm',
        physical_qualification='NOT_RUN_NEW_STARTUP_BOUNDARY_ASSUMPTION',
        factory_scope=expanded.SCOPE,motion_protocol=history.to_dict())
    state.update(classification=CLASSIFICATION,model_scope=SCOPE,diagnostic_only=True,
        accepted=False,restart_allowed=False,parameters=copy.deepcopy(p),
        nominal_parameters=copy.deepcopy(p),nominal_parameter_provenance=nominal,
        parameter_provenance=provenance,independent_domain_initialization=initialization,
        source_sha256=source_hashes(),ramp_component_source_sha256=source_hashes(),
        ramp_motion=history.to_dict(),ramped_boundary_history=[],
        ramped_physical_qualification='NOT_RUN',accepted_time_advanced_s=0.)
    # At t=0 both poses coincide, but rebuild this observer explicitly under the
    # new motion object rather than inheriting the old observer's provenance.
    state['bore_observation']=initial_bore_crossing_state(positions(pa),pa.particle_id,
        pa.component_id,pa.m,0.,history.global_pose(0.),p['bore_width'])
    for key in ('spatial_model_transition','branch_origin_sha256','expanded_initial_qualification'):
        require(key not in state,'independent startup must not fabricate an old qualification/origin')
    validate_interval_energy_state(state['interval_energy']);validate_bore_crossing_state(state['bore_observation'])
    return pa,state


def validate_ramped_state(state):
    require(state.get('classification')==CLASSIFICATION and state.get('model_scope')==SCOPE and
        state.get('diagnostic_only') is True and state.get('independent_initialization') is True and
        state.get('accepted') is False and state.get('restart_allowed') is False,
        'independent unaccepted ramped diagnostic state required')
    require('spatial_model_transition' not in state and 'expanded_initial_qualification' not in state and
        state['independent_domain_initialization']['kind']==SCOPE and
        state['independent_domain_initialization']['parent_checkpoint'] is None and
        state['independent_domain_initialization']['old_state_read'] is False,
        'old process state or old qualification cannot become ramped startup')
    p=state['nominal_parameters'];expected,nominal_provenance=expanded.parameters(.008)
    require(p==expected and state['parameters']==p and state['nominal_parameter_provenance']==nominal_provenance,
        'registered nominal material/geometric inputs changed')
    motion=RampedHotMotion(p,RAMP_DURATION_S)
    require(state['ramp_motion']==motion.to_dict() and
        state['parameter_provenance']['boundary_condition_history']==motion.to_dict() and
        state['parameter_provenance']['scope']==SCOPE,'ramped motion provenance differs')
    validate_interval_energy_state(state['interval_energy']);validate_bore_crossing_state(state['bore_observation'])
    rows=state['ramped_boundary_history']
    require(len(rows)==state['step']==state['interval_energy']['steps'] and
        len(state['interval_energy']['records'])==state['step'],'complete independent ramp/energy prefix required')
    t=0.
    for row in rows:
        require(row['t0_s']==t and row==motion.step(t,row['dt_s']).to_dict(),'ramp interval history is broken or changed')
        t=row['t1_s']
    require(t==state['time_s']==state['bore_observation']['last_time_s'],'independent ramp time/prefix differs')
    return motion


class RampedHotTrial:
    def __init__(self,pa,state,dt):
        history=validate_ramped_state(state)
        require(len(pa.x)==918 and int(np.asarray(state['free']).sum())==1961,'918-particle fixed-N diagnostic required')
        self.pa,self.state,self.p=pa,state,dict(state['nominal_parameters'])
        self.x0=positions(pa);self.s0=pa.S.reshape(-1,3,3).copy()
        self.bulk=bulk_modulus(self.p['shear'],self.p['poisson'])
        self.law=FKS(**{k:self.p[k] for k in FKS.__dataclass_fields__})
        self.motion=history.step(state['time_s'],dt)
        self.contact=expanded.BatchedChainRuleContact(tool_radius=self.p['tool_radius'],bore_width=self.p['bore_width'],
            tool_height=self.p['tool_top']-self.p['layer'],particle_radius=self.p['dx']/2,
            normal_stiffness=self.p['normal_stiffness'],tangential_stiffness=self.p['normal_stiffness'],
            friction=self.p['friction'],base_z=-self.p['substrate_depth'])

    def _require_interval(self,dt):
        require(isinstance(self.motion,IntervalMotion) and float(dt)==self.motion.dt_s and
            self.state['time_s']==self.motion.t0 and self.state['ramp_motion']==self.motion.history.to_dict(),
            'one explicit matching motion interval required for every trial call')

    def prescribed(self,dt):
        self._require_interval(dt)
        u=np.zeros_like(self.x0)
        u[self.pa.actuator.astype(bool),2]=self.motion.feed_displacement_m
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
        raw_i, raw_j, raw_g = i, j, g
        numerical_domain=require_visible_support_numerical_domain(x1,i,j,
            self.motion.pose(self.state['time_s']+dt),p,pa.particle_id)
        graph = filter_visible_support(x1, i, j, g,
            self.motion.pose(self.state['time_s']+dt), p, pa.particle_id)
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


    def tangent(self, result):
        """Neighbor-local assembly + local contact secant; NOT full Jacobian."""
        from afsd_backend.pysph.quasistatic_tangent import (gradient_map, material_tangent,
                                          assemble_frozen_force_tangent)
        self._require_interval(result['dt_s'])
        graph=result['visibility_graph']
        for key in ('i','j','g'):
            if not np.array_equal(result[key],graph[key]):
                raise ValueError('tangent must use the exact residual visible graph')
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
            args = (self.motion.pose(self.state['time_s']), self.motion.pose(self.state['time_s']+dt),
                    dt, self.state['contact_history'])
            a = self.contact.evaluate(self.x0, plus, *args, check_path=False)
            b = self.contact.evaluate(self.x0, minus, *args, check_path=False)
            blocks[:, :, axis] = (a['force_N']-b['force_N'])/(2*eps)
        for particle, block in enumerate(blocks):
            rows = slice(3*particle, 3*particle+3)
            J[rows, rows] += block
        return J



def make_trial(pa,state,dt):return RampedHotTrial(pa,state,dt)


def advance_candidate(trial,result):
    """Detached fresh-qualified advance only; no restart or publication."""
    require(isinstance(trial,RampedHotTrial),'explicit ramped trial required')
    trial._require_interval(result['dt_s']);validate_ramped_state(trial.state)
    fresh=trial.evaluate(result['displacement'],result['dt_s'])
    require(np.array_equal(result['i'],fresh['i']) and np.array_equal(result['j'],fresh['j']) and
        result['graph_fingerprint']==fresh['graph_fingerprint'],'fresh graph changed at identical displacement')
    energy=original_qs.energy_diagnostics(trial,fresh)
    interval=require_interval_energy(trial.state['interval_energy'],energy)
    pa,state=original_qs.accepted_copy(trial,fresh)
    observer=observe_square_bore_crossings(trial.state['bore_observation'],trial.x0,positions(pa),
        pa.particle_id,pa.component_id,pa.m,trial.motion.t0,trial.motion.t1,trial.motion.pose)
    state.update(interval_energy=interval,bore_observation=observer['next_state'],
        crossing_ids=observer['crossing_ids'],crossing_mass_kg=observer['crossing_mass_kg'],
        ramped_boundary_history=copy.deepcopy(trial.state['ramped_boundary_history'])+[trial.motion.to_dict()],
        last_ramped_interval=trial.motion.to_dict(),accepted=False,restart_allowed=False,diagnostic_only=True,
        ramped_physical_qualification='FRESH_AND_FULL_PREFIX_ENERGY_ONLY_REQUIRES_EXTERNAL_GEOMETRY_AND_TIME_SCREEN',
        accepted_time_advanced_s=0.)
    require(state['interval_energy']['records'][:-1]==trial.state['interval_energy']['records'] and
        state['step']==trial.state['step']+1 and state['time_s']==trial.motion.t1,
        'ramped candidate reset history or advanced a different interval')
    validate_ramped_state(state)
    return pa,state,fresh


def static_adaptation_audit():
    """AST check; no code generation, monkeypatching or physical evaluation."""
    def method(path,cls,name):
        tree=ast.parse(path.read_text(encoding='utf-8'))
        owner=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
        return copy.deepcopy(next(n for n in owner.body if isinstance(n,ast.FunctionDef) and n.name==name))
    class IntendedPoseDelta(ast.NodeTransformer):
        def visit_ImportFrom(self,node):
            if node.level==1 and node.module=='quasistatic_tangent':
                node.level=0;node.module='afsd_backend.pysph.quasistatic_tangent'
            return node
        def visit_Call(self,node):
            self.generic_visit(node)
            if isinstance(node.func,ast.Name) and node.func.id=='pose':
                require(len(node.args)==2 and isinstance(node.args[0],ast.Name) and node.args[0].id=='p',
                    'unexpected original pose call')
                node.func=ast.Attribute(value=ast.Attribute(value=ast.Name(id='self',ctx=ast.Load()),
                    attr='motion',ctx=ast.Load()),attr='pose',ctx=ast.Load());node.args=node.args[1:]
            return node
        def visit_Name(self,node):
            if node.id=='CLASSIFICATION':node.id='SPATIAL_CLASSIFICATION'
            return node
    here=Path(__file__).resolve();frozen=ROOT/'src/afsd_backend/pysph/quasistatic_process.py'
    old=method(ROOT/'scripts/visible_quasistatic.py','VisibleSupportTrial','evaluate')
    new=method(here,'RampedHotTrial','evaluate');expected=IntendedPoseDelta().visit(old)
    require(ast.dump(expected,include_attributes=False)==ast.dump(new,include_attributes=False),
        'evaluate differs beyond four interval poses and unchanged spatial-classification alias')
    old_t=IntendedPoseDelta().visit(method(frozen,'QuasistaticTrial','tangent'))
    new_t=method(here,'RampedHotTrial','tangent')
    # Three added statements: interval check, actual graph assignment, graph loop.
    inserts=[n for n in new_t.body if isinstance(n,(ast.Assign,ast.Expr,ast.For)) and (
        (isinstance(n,ast.Assign) and any(isinstance(k,ast.Name) and k.id=='graph' for k in n.targets)) or
        (isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Attribute) and n.value.func.attr=='_require_interval') or
        (isinstance(n,ast.For) and isinstance(n.target,ast.Name) and n.target.id=='key'))]
    require(len(inserts)==3,'tangent guard additions differ')
    for node in inserts:new_t.body.remove(node)
    require(ast.dump(old_t,include_attributes=False)==ast.dump(new_t,include_attributes=False),
        'tangent differs beyond interval poses/guards and absolute import relocation')
    require(SPATIAL_CLASSIFICATION=='HOT_QS_VISIBLE_SUPPORT_CHAIN_RULE_CONTACT_APPROXIMATE_V1',
        'original spatial classification alias changed')
    return dict(status='PASS_STATIC_RAMPED_BOUNDARY_ONLY_EVALUATE_AND_TANGENT_ADAPTATION',
        evaluate_pose_replacements=4,tangent_pose_replacements=2,contact_fd_epsilon_m=1e-9,
        runtime_monkeypatch=False,physics_evaluations=0,nnps_evaluations=0,
        non_pose_deltas=['unchanged original spatial CLASSIFICATION alias',
            'tangent adds same actual visible graph check plus explicit interval check',
            'package-relative tangent import becomes absolute in the scripts module'])
