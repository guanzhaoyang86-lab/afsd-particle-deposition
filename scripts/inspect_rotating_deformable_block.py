"""Independent saved-array audit of the 64-particle rotating deformable block.

Only JSON, NumPy and the generic VTK reader are imported.  Constitutive and
contact implementations are deliberately not reference-answer generators.
Passing algebraic identities is NOT a finite-strain energy certification.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_scaffold.vtk_legacy import read_point_cloud

CLASSIFICATION = 'ENGINEERING_APPROXIMATE_ROTATING_DEFORMABLE_BLOCK_NOT_AFSD'
MASS = 2700*(.0005)**3
G = 26.3e9
BULK = 2*G*1.33/(3*.34)
END = 2e-7


def require(condition, message):
    if not condition: raise ValueError(message)


def array(value, shape, name):
    a = np.asarray(value, float)
    require(a.shape == shape and np.isfinite(a).all(), 'shape/nonfinite '+name)
    return a


def close(actual, expected, atol, name):
    a, b = np.asarray(actual, float), np.asarray(expected, float)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(), 'shape/nonfinite '+name)
    defect = float(np.max(np.abs(a-b))) if a.size else 0.
    require(defect <= atol, name+': defect '+str(defect)+' > '+str(atol))
    return defect


def load(path):
    def reject(value): raise ValueError('nonfinite JSON constant '+value)
    return json.loads(Path(path).read_text(encoding='utf-8'), parse_constant=reject)


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def file_hash(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reconstruct_gradient(velocity, i, j, g, volume, B):
    """Literal pair sum, independent of the mechanics implementation."""
    velocity, volume = np.asarray(velocity), np.asarray(volume)
    raw = np.zeros((len(velocity), 3, 3))
    for a, b, grad in zip(i, j, g):
        difference = velocity[b]-velocity[a]
        raw[a] += volume[b]*np.outer(difference, grad)
        raw[b] += volume[a]*np.outer(difference, grad)
    return raw@B


def audit_operator(x, velocity, mass, density, i, j, g, B, L, force, sigma):
    """Audit 3D support, radial hyperbolic gradient and adjoint pair arithmetic."""
    x, velocity, mass, density = map(np.asarray, (x, velocity, mass, density))
    count = len(x); i, j = np.asarray(i, int), np.asarray(j, int)
    require(i.ndim == j.ndim == 1 and i.shape == j.shape and len(i) > 0, 'pair shape')
    require(np.all((i >= 0)&(j < count)&(i < j)), 'canonical pair indices')
    require(len(set(zip(i.tolist(), j.tolist()))) == len(i), 'duplicate pair')
    g, B = array(g, (len(i), 3), 'radial gradient'), array(B, (count, 3, 3), 'correction B')
    volume = mass/density
    require(np.all(volume > 0), 'positive particle volumes')
    displacement = x[i]-x[j]
    distance = np.linalg.norm(displacement, axis=1)
    require(np.all(distance > 0), 'distinct support particles')
    h = .00065; q = distance/h
    derivative = np.where(q < 1, 3*q*q-6, np.where(q < 2, -3*(2-q)**2, 0.))
    radial_reference = (15/(62*np.pi*h**4)*derivative/distance)[:, None]*displacement
    close(g, radial_reference, 1e-10*max(1., float(np.max(np.abs(radial_reference)))), 'independent radial kernel derivative')
    expected_pairs = {(a,b) for a in range(count) for b in range(a+1,count) if np.linalg.norm(x[a]-x[b]) < 2*h}
    require(set(zip(i.tolist(),j.tolist())) == expected_pairs, 'complete 3D radial neighbor inventory')
    moment = np.zeros_like(B); reconstructed_force = np.zeros((count, 3))
    for a, b, grad in zip(i, j, g):
        dyad = np.outer(x[b]-x[a], grad)
        moment[a] += volume[b]*dyad; moment[b] += volume[a]*dyad
        pair = volume[a]*volume[b]*(sigma[a]@B[a].T@grad+sigma[b]@B[b].T@grad)
        reconstructed_force[a] += pair; reconstructed_force[b] -= pair
    condition = np.linalg.cond(moment)
    require(np.isfinite(condition).all() and np.max(condition) <= 1e6, 'ill-conditioned correction moment')
    close(moment@B, np.repeat(np.eye(3)[None], count, axis=0), 1e-10, 'correction inverse moment')
    ranks = []
    for a in range(count):
        neighbors = np.r_[j[i == a], i[j == a]]
        ranks.append(int(np.linalg.matrix_rank(x[neighbors]-x[a])))
    require(min(ranks) == 3, 'neighborhood not spatially three-dimensional')
    computed_L = reconstruct_gradient(velocity, i, j, g, volume, B)
    close(L, computed_L, 1e-10, 'independent corrected velocity gradient')
    force_scale = float(np.linalg.norm(force, axis=1).sum())
    close(force, reconstructed_force, 1e-12+1e-10*force_scale, 'independent adjoint force')
    netforce = np.sum(force, axis=0)
    torque_each = np.cross(x, force); nettorque = torque_each.sum(axis=0)
    close(netforce, np.zeros(3), 1e-12+1e-10*force_scale, 'total internal force')
    close(nettorque, np.zeros(3), 1e-14+1e-10*np.linalg.norm(torque_each, axis=1).sum(), 'total internal torque')
    power = float(np.sum(force*velocity)); stress_power = float(np.einsum('n,nij,nij', volume, sigma, computed_L))
    power_scale = float(np.sum(np.linalg.norm(force, axis=1)*np.linalg.norm(velocity, axis=1)))
    defect = close(power+stress_power, 0., 1e-12+1e-10*power_scale, 'discrete power conjugacy')
    return dict(maximum_condition=float(np.max(condition)), minimum_neighbor_rank=min(ranks),
        power_residual_W=defect, internal_power_W=power, stress_power_W=stress_power,
        net_internal_force_N=netforce.tolist(), net_internal_torque_N_m=nettorque.tolist())


def audit_interval(interval, before, force, wall, free_axes, mass, dt, start, initial):
    """Signed impulse work using independently recomputed actual interval velocities."""
    nc = len(wall); wall = np.asarray(wall, float).reshape(nc,3)
    h = float(interval['dt_s'])
    require(np.isfinite(h) and 0 < h <= dt, 'invalid interval duration')
    close(interval['start_s'], start, 1e-21, 'interval chronology')
    vb, va = array(interval['v_before'], (3,), 'interval before'), array(interval['v_after'], (3,), 'interval after')
    close(vb, before, 1e-12, 'interval velocity continuity')
    jn = np.asarray(interval['normal_impulse'], float).reshape(nc,3)
    jt = np.asarray(interval['tangential_impulse'], float).reshape(nc,3)
    reaction = array(interval['constraint_impulse'], (3,), 'constraint impulse')
    close(jn, np.tile([0.,0.,-h], (nc,1)), 1e-18, 'prescribed normal impulse')
    close(jt[:,2], np.zeros(nc), 1e-18, 'tangential normal component')
    close(mass*(va-vb), h*force+jn.sum(axis=0)+jt.sum(axis=0)+reaction, 1e-18, 'actual kick momentum')
    free = np.asarray(free_axes)
    close((1-free)*va, np.zeros(3), 1e-12, 'fixed velocity')
    close(free*reaction, np.zeros(3), 1e-18, 'constraint impulse on a free axis')
    midpoint = .5*(va+vb); relative = midpoint[None,:]-wall; relative[:,2] = 0.
    close(np.asarray(interval['tangential_relative_midpoint']).reshape(nc,3), relative, 1e-12, 'relative midpoint')
    q_each = -np.einsum('ij,ij->i', jt, relative)
    saved_q = array(interval['friction_heat_raw_J'], (nc,), 'raw friction')
    close(saved_q, q_each, 1e-18, 'signed friction work')
    require(np.all(saved_q[q_each < 0] < 0), 'negative raw friction clipped or made absolute')
    frame = wall[0] if nc else initial
    close(interval['reference_frame_velocity'], frame, 1e-12, 'relative reference frame')
    normal_force = np.array([0.,0.,-float(nc)])
    scale = max(float(np.linalg.norm(initial-frame)), float(np.linalg.norm(wall-frame)) if nc else 0.,
                dt*float(np.linalg.norm(free*(force+normal_force)/mass)), np.finfo(float).tiny)
    close(interval['relative_speed_scale_m_s'], scale, 1e-12*max(1.,scale), 'unmodified relative speed scale')
    bound = 256*np.finfo(float).eps*np.linalg.norm(jt,axis=1)*(scale+np.linalg.norm(vb-frame)+np.linalg.norm(va-frame)+2*np.linalg.norm(wall-frame,axis=1))
    close(interval['friction_rounding_bound_J'], bound, max(1e-35, float(np.max(bound,initial=0))*1e-12), 'original friction rounding bound')
    require(np.all(saved_q >= -bound), 'negative accepted friction')
    # This manufactured benchmark is sliding; neither positive Q alone nor a
    # velocity overwritten to the wall is sufficient to satisfy Coulomb's law.
    if nc:
        speed = np.linalg.norm(relative,axis=1)
        require(np.all(speed > .01), 'left manufactured sliding domain')
        close(jt, -.2*h*relative/speed[:,None], 1e-18, 'Coulomb direction and cap')
    other = float(h*force@midpoint); wall_work = float(np.sum((jn+jt)*wall))
    constraint = float(reaction@midpoint); normal = np.einsum('ij,j->i',jn,midpoint)
    dk = .5*mass*float(va@va-vb@vb); q = float(q_each.sum())
    scale_energy = abs(other)+abs(wall_work)+abs(constraint)+abs(q)
    tolerance = 1e-15+1e-9*scale_energy
    for key, expected in [('other_work_J',other),('wall_work_J',wall_work),('constraint_work_J',constraint),
                          ('normal_work_J',normal),('delta_kinetic_J',dk)]:
        close(interval[key], expected, tolerance, 'interval '+key)
    residual = dk-other-constraint-float(np.sum((jn+jt)*midpoint))
    close(interval['kinetic_work_residual_J'], residual, tolerance, 'saved kinetic work residual')
    close(residual, 0., tolerance, 'contact work identity')
    return dict(v_after=va, h=h, work_velocity_integral=h*midpoint, internal=other, wall=wall_work,
        constraint=constraint, normal=float(normal.sum()), q=q, dk=dk, jn=jn.sum(axis=0), jt=jt.sum(axis=0),
        reaction=reaction, negative=int(np.sum(q_each<0)), minimum_raw=float(np.min(q_each)) if nc else None)


def independent_first_step(x, dt):
    """Closed-form rest-state impulse; no dynamics implementation imported."""
    x = np.asarray(x); v = np.zeros_like(x); top = np.isclose(x[:,2], .00175, rtol=0., atol=1e-16)
    wall = np.column_stack((-20*x[top,1],20*x[top,0],np.zeros(np.sum(top))))
    impulse = .2*dt*wall/np.linalg.norm(wall,axis=1)[:,None]
    v[top] = impulse/MASS
    work = float(np.sum(impulse*wall)); kinetic = float(np.sum(impulse*impulse)/(2*MASS))
    return dict(position=x+.5*dt*v, velocity=v, wall_work=work, kinetic=kinetic, q=work-kinetic)


SOURCE_PATHS = {'scripts/probe_rotating_deformable_block.py',
    'configs/benchmarks/rotating_deformable_block.json','docs/ROTATING_DEFORMABLE_CONTRACT.md',
    *('src/afsd_backend/pysph/'+name+'.py' for name in ('particles','kernel','mechanics','material',
       'geometry','rigid_contact_kinematics','contact_kick','contact_cone','integrator','output')),
    'src/afsd_scaffold/vtk_legacy.py','third_party/PySPH/pysph/sph/integrator.py',
    'third_party/PySPH/pysph/sph/integrator_step.py','third_party/PySPH/pysph/base/nnps.py',
    'third_party/PySPH/pysph/base/nnps_base.pyx',
    'third_party/PySPH/pysph/base/linked_list_nnps.pyx'}
PARAMETERS = dict(dx=.0005,h_ratio=1.3,rho0=2700,shear=26.3e9,poisson=.33,temperature=293.15,
    omega=20,traverse=0,layer=.00195,tool_top=.02,tool_radius=.02,bore_width=.008,substrate_depth=.01,
    base_temperature=293.15,tool_temperature=293.15,normal_force=1,friction=.2,cfl=.1,width=.002,
    bar_length=.002,actuator_depth=0,marker_low=0,marker_high=.002,feed=0,end_time=2e-7,
    a1=317600000,a2=1000000,a3=70,b1=12.44,b2=5.65,c1=13,c2=.01,room=293.15,melt=855.15,
    reference_rate=1,rate_floor=1)


def audit_configuration(config):
    require(config['classification'] == CLASSIFICATION and config['production_authorized'] is False, 'configuration classification')
    require(config['parameters'] == PARAMETERS, 'frozen physical parameters changed')
    require(config['dt_s'] in (5e-9,2.5e-9), 'unregistered timestep')
    require(config['physical_scaling'] == 'NONE' and config['tags_are_observation_only'] is True, 'scaling/passive tags policy')
    require(config['spatial_operator'] == 'corrected_work_conjugate', 'spatial operator policy')
    require(config['normal_law'] == 'PRESCRIBED 1 N PER TOP PARTICLE; NORMAL DOF CLAMPED; NOT PENALTY', 'normal law policy')
    status = config['field_status']
    require(status['temperature'] == 'PRESCRIBED_INPUT_NOT_SOLVED' and
        status['Es'] == 'DIAGNOSTIC_ONLY_NOT_CERTIFIED_FREE_ENERGY' and
        status['friction_dissipation'] == 'SOLVED_RAW_NOT_DEPOSITED_AS_HEAT', 'field honesty policy')
    source = config['source_configuration']
    require(source['classification'] == CLASSIFICATION and source['production_authorized'] is False, 'source policy')
    require({k:v['value'] for k,v in source['parameters'].items()} == PARAMETERS, 'source parameter values')
    require(all(v.get('unit') and v.get('source') and v.get('status') for v in source['parameters'].values()), 'source/unit/status missing')
    require(source['dt_cases_s'] == [5e-9,2.5e-9] and source['run_timeout_s'] == 110 and source['external_timeout_s'] == 120, 'registered budget changed')
    return float(config['dt_s'])


def audit_row(row, pid, x, t):
    require(row['particle_id'] == pid and row['feature_id'] == 1 and row['contact_id'] == str(pid)+':1', 'persistent shoulder identity')
    close(row['application_point'], x, 2e-17, 'force application point')
    close(row['contact_point'], [x[0],x[1],.00195], 2e-17, 'shoulder projection')
    close(row['normal'], [0.,0.,-1.], 1e-15, 'shoulder normal')
    close(row['body_origin'], [0.,0.,.00195], 1e-15, 'tool origin')
    close(row['translation_velocity'], [0.,0.,0.], 0., 'tool translation')
    close(row['omega'], [0.,0.,20.], 0., 'tool angular velocity')
    close(row['gap'], .00195-x[2], 2e-17, 'shoulder gap')
    close(row['penetration'], .00025-(.00195-x[2]), 2e-17, 'shoulder overlap')
    angle = 20*t; c,s = np.cos(angle),np.sin(angle)
    local = [c*x[0]+s*x[1], -s*x[0]+c*x[1]]
    require(max(map(abs,local)) >= .004 and np.linalg.norm(local) <= .02, 'finite open-bore shoulder domain')
    wall = np.array([-20*x[1],20*x[0],0.])
    close(row['wall_velocity'], wall, 1e-12, 'independent rigid wall velocity')
    close(row['wall_velocity_at_face'], wall, 1e-12, 'shoulder special-case face velocity')
    return wall


def audit_trace(trace, pid, role, x, vb, va, force, dt, step):
    require(trace['particle_id'] == pid and trace['step'] == step, 'trace particle/step identity')
    require(trace['stage'] == 'MIDPOINT_FROZEN_GEOMETRY_KICK', 'trace stage')
    close(trace['time_mid_s'], (step-.5)*dt, 1e-21, 'contact time layer')
    close(trace['position_m'], x, 2e-17, 'trace midpoint position')
    close(trace['material_velocity_before_m_s'], vb, 1e-12, 'trace precontact velocity')
    close(trace['internal_force_N'], force, 1e-12, 'trace internal force')
    close(trace['mass_kg'], MASS, 0., 'unscaled mass')
    free = [0,0,0] if role == 1 else [1,1,0] if role == 2 else [1,1,1]
    require(trace['free_axes'] == free, 'material boundary role/free axes')
    rows = trace['contact_rows']
    require(len(rows) == (1 if role == 2 else 0), 'contact occurs only on top layer')
    wall = np.asarray([audit_row(row,pid,x,(step-.5)*dt) for row in rows]).reshape(-1,3)
    result = trace['result']; close(result['velocity'], va, 1e-12, 'kick final velocity')
    close(result['constraint_projector'], np.diag(free), 0., 'constraint projection')
    require(bool(result['intervals']), 'accepted kick interval missing')
    sums = dict(internal=0.,wall=0.,constraint=0.,normal=0.,q=0.,dk=0.,negative=0)
    impulses = dict(jn=np.zeros(3),jt=np.zeros(3),reaction=np.zeros(3))
    work_velocity = np.zeros(3); elapsed = 0.; before = vb; minq = None
    for interval in result['intervals']:
        a = audit_interval(interval,before,force,wall,free,MASS,dt,elapsed,vb)
        before = a['v_after']; elapsed += a['h']; work_velocity += a['work_velocity_integral']
        for key in sums: sums[key] += a[key]
        for key in impulses: impulses[key] += a[key]
        if a['minimum_raw'] is not None: minq = a['minimum_raw'] if minq is None else min(minq,a['minimum_raw'])
    close(elapsed, dt, 1e-21, 'whole kick duration'); close(before, va, 1e-12, 'whole kick endpoint')
    close(result['constraint_impulse'], impulses['reaction'], 1e-18, 'aggregate constraint impulse')
    for key, expected in [('other_work_J',sums['internal']),('wall_work_J',sums['wall']),
            ('constraint_work_J',sums['constraint']),('delta_kinetic_J',sums['dk'])]:
        close(result[key], expected, 1e-15, 'aggregate '+key)
    for key, expected in [('normal_impulse',impulses['jn']),('tangential_impulse',impulses['jt'])]:
        saved = np.asarray(result[key]).reshape(len(rows),3)
        close(saved.sum(axis=0), expected, 1e-18, 'aggregate '+key)
    close(np.sum(result['friction_heat_J']), sums['q'], 1e-18, 'aggregate raw friction')
    close(result['kinetic_work_residual_J'], 0., 1e-15, 'aggregate kinetic work')
    return dict(**sums,**impulses,work_velocity=work_velocity/dt,minimum_raw=minq)


def audit_vtk(path, record, original):
    vtk = read_point_cloud(path)
    require(CLASSIFICATION in vtk['title'] and 'ISOTHERMAL_INPUT_ONLY_NO_SOLVED_TEMPERATURE' in vtk['title'], 'VTK classification')
    expected = {'particle_id','material_id','component_id','tracer_type','wholebar_tag','actuator','boundary_role',
        'as_rs_tag','mass_kg','density_kg_per_m3','pressure_Pa','equivalent_stress_Pa','equivalent_plastic_strain',
        'smoothing_length_m',*('stress_'+a+'_Pa' for a in ('xx','yy','zz','xy','xz','yz'))}
    if record['step']: expected |= {'equivalent_deviatoric_strain_rate_per_s','equivalent_plastic_strain_rate_per_s'}
    require(set(vtk['scalars']) == expected and set(vtk['vectors']) == {'velocity_m_per_s','initial_position_m'}, 'fabricated/omitted VTK fields')
    close(vtk['points'], record['position_m'], 0., 'VTK position')
    close(vtk['vectors']['velocity_m_per_s'], record['velocity_m_s'], 0., 'VTK velocity')
    close(vtk['vectors']['initial_position_m'], original, 0., 'VTK persistent initial position')
    for key in ('particle_id','tracer_type','boundary_role','as_rs_tag'): require(vtk['scalars'][key] == record[key], 'VTK '+key)
    fields = dict(mass_kg='mass_kg',density_kg_per_m3='density_kg_m3',pressure_Pa='pressure_Pa',
        equivalent_stress_Pa='equivalent_stress_Pa',equivalent_plastic_strain='ep')
    if record['step']: fields.update(equivalent_deviatoric_strain_rate_per_s='strain_rate_per_s',equivalent_plastic_strain_rate_per_s='ep_rate_per_s')
    for target, source in fields.items(): close(vtk['scalars'][target],record[source],0.,'VTK '+target)
    for a,b,key in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
        close(vtk['scalars']['stress_'+key+'_Pa'],np.asarray(record['stress_Pa'])[:,a,b],0.,'VTK stress '+key)
    close(vtk['scalars']['smoothing_length_m'],np.full(64,.00065),1e-18,'VTK fixed smoothing length')
    close(vtk['scalars']['material_id'],np.ones(64),0.,'same material for all tracers')
    for key in ('component_id','wholebar_tag','actuator'): close(vtk['scalars'][key],np.zeros(64),0.,'VTK '+key)


def audit_history_record(r, original, ids, roles, step, dt):
    require(r['step'] == step and r['particle_id'] == ids and r['boundary_role'] == roles, 'record identity/role/step')
    require(r['tracer_type'] == [p%3 for p in ids] and r['as_rs_tag'] == [0]*64, 'passive label history')
    close(r['time_s'], step*dt, 1e-21, 'record time')
    x = array(r['position_m'], (64,3), 'position'); v = array(r['velocity_m_s'], (64,3), 'velocity')
    mass = array(r['mass_kg'], (64,), 'mass'); rho = array(r['density_kg_m3'], (64,), 'density')
    close(mass, np.full(64,MASS), 1e-21, 'unscaled mass')
    require(abs(mass.sum()-2.16e-5) <= 2.16e-17, 'mass conservation')
    require(np.max(np.abs(rho/2700-1)) <= .01, 'density bound')
    require(np.max(np.linalg.norm(v,axis=1)) <= 1., 'speed bound')
    require(np.max(np.linalg.norm(x-original,axis=1)) <= 1e-5, 'displacement bound')
    fixed = np.zeros((64,3),bool); roles_a = np.asarray(roles)
    fixed[roles_a==1] = True; fixed[roles_a==2,2] = True
    close((x-original)[fixed],np.zeros(np.sum(fixed)),2e-17,'fixed position')
    close(v[fixed],np.zeros(np.sum(fixed)),1e-12,'fixed velocity')
    s = array(r['deviatoric_stress_Pa'],(64,3,3),'deviatoric stress')
    close(s,s.swapaxes(1,2),1e-6,'symmetric stress')
    close(np.trace(s,axis1=1,axis2=2),np.zeros(64),1e-6,'deviatoric stress trace')
    pressure = BULK*(rho/2700-1)
    close(r['pressure_Pa'],pressure,1e-6,'EOS pressure')
    close(r['stress_Pa'],s-pressure[:,None,None]*np.eye(3),1e-6,'total stress')
    close(r['equivalent_stress_Pa'],np.sqrt(1.5*np.sum(s*s,axis=(1,2))),1e-6,'von Mises stress')
    ep = array(r['ep'],(64,),'plastic strain'); require(np.min(ep) >= 0.,'negative plastic strain')
    close(r['tool_origin_m'],[0.,0.,.00195],0.,'tool origin history')
    close(r['tool_angle_rad'],20*step*dt,1e-20,'tool angle history')
    kinetic = float(.5*np.sum(mass[:,None]*v*v))
    shear = float(np.sum(mass/rho*np.sum(s*s,axis=(1,2))/(4*G)))
    ratio = rho/2700
    volume = float(np.sum(mass*BULK/2700*(np.log(ratio)+1/ratio-1)))
    close(r['kinetic_J'],kinetic,1e-18,'kinetic energy')
    close(r['Es_diagnostic_J'],shear+volume,1e-16,'Es diagnostic (not free energy)')
    close(r['total_mass_kg'],mass.sum(),1e-18,'reported mass')
    c = r['cumulative']
    residual = kinetic+shear+volume+c['plastic_work_J']+c['friction_raw_J']-c['wall_work_J']-c['constraint_work_J']
    close(r['energy_diagnostic_residual_J'],residual,1e-16,'reported energy defect')
    warning = abs(residual)>1e-12+.01*(abs(c['wall_work_J'])+abs(c['constraint_work_J'])+c['friction_raw_J'])
    require(r['energy_status'] == ('ENERGY_WARNING' if warning else 'DIAGNOSTIC_WITHIN_SCREEN_NOT_CERTIFIED_FREE_ENERGY'),'energy warning label')
    return dict(x=x,v=v,mass=mass,rho=rho,s=s,ep=ep,K=kinetic,Es=shear+volume,residual=residual,warning=warning)


def audit_case(directory):
    directory = Path(directory); manifest = load(directory/'manifest.json')
    require(manifest['status']=='COMPLETE_DEFORMABLE_COMPONENT_NOT_AFSD','full registered run not complete')
    require(manifest['production_authorized'] is False and manifest['mode']=='FULL_REGISTERED_COMPONENT','component scope')
    cfg = manifest['configuration']; dt = audit_configuration(cfg); count = round(END/dt)
    require(canonical_hash(cfg)==manifest['configuration_sha256'],'configuration hash')
    require(cfg['source_configuration']==load(ROOT/'configs/benchmarks/rotating_deformable_block.json'),'current configuration content')
    require(set(manifest['source_sha256'])==SOURCE_PATHS,'complete source inventory')
    for name,digest in manifest['source_sha256'].items(): require(file_hash(ROOT/name)==digest,'source hash '+name)
    expected_files = {'history.json','steps.jsonl','last_accepted.npz',*(f'frame_{i:04d}.vtk' for i in range(count+1))}
    require(set(manifest['artifact_sha256'])==expected_files,'artifact inventory')
    for name,digest in manifest['artifact_sha256'].items(): require(file_hash(directory/name)==digest,'artifact hash '+name)
    require(manifest['accepted_steps']==count and manifest['target_steps']==count,'accepted step count')
    close(manifest['accepted_time_s'],END,1e-21,'complete common endpoint')
    require(manifest['wall_seconds']<=110,'internal wall budget')
    history = load(directory/'history.json')
    journal = [json.loads(line) for line in (directory/'steps.jsonl').read_text(encoding='utf-8').splitlines()]
    require(len(history)==count+1 and len(journal)==count,'history/journal completeness')
    original = np.array(np.meshgrid(.01+(np.arange(4)-1.5)*.0005,(np.arange(4)-1.5)*.0005,
        (np.arange(4)+.5)*.0005,indexing='ij')).reshape(3,-1).T
    ids = list(range(1,65)); roles = np.where(original[:,2]==.00025,1,np.where(original[:,2]==.00175,2,0)).tolist()
    state = audit_history_record(history[0],original,ids,roles,0,dt)
    close(state['x'],original,0.,'initial geometry'); close(state['v'],np.zeros((64,3)),0.,'initial velocity')
    close(state['s'],np.zeros((64,3,3)),0.,'initial stress'); close(state['rho'],np.full(64,2700.),0.,'initial density')
    close(state['ep'],np.zeros(64),0.,'initial plastic strain')
    audit_vtk(directory/'frame_0000.vtk',history[0],original)
    sums = {key:0. for key in history[0]['cumulative']}
    require(all(v==0 for v in history[0]['cumulative'].values()),'zero initial ledger')
    maximum_condition=0.; maximum_power=0.; maximum_force=0.; maximum_torque=0.
    maximum_energy=0.; maximum_free_w=0.; negative=0; minimum_q=None; warnings=0
    total_external_impulse=np.zeros(3); total_external_angular=np.zeros(3)
    for step,(row,record) in enumerate(zip(journal,history[1:]),1):
        require(row['step']==step and row['particle_id']==ids and row['boundary_role']==roles,'journal identity')
        close(row['dt_s'],dt,0.,'journal dt'); close(row['time_mid_s'],(step-.5)*dt,1e-21,'journal time')
        current = audit_history_record(record,original,ids,roles,step,dt)
        xmid = array(row['x_mid_m'],(64,3),'midpoint position')
        vb = array(row['v_before_m_s'],(64,3),'before velocity'); va = array(row['v_after_m_s'],(64,3),'after velocity')
        close(vb,state['v'],1e-12,'velocity history continuity'); close(va,current['v'],1e-12,'accepted velocity')
        close(xmid,state['x']+.5*dt*vb,2e-17,'actual first drift')
        close(current['x'],xmid+.5*dt*va,2e-17,'actual second drift')
        close(row['x_before_m'],state['x'],0.,'saved pre-step geometry'); close(row['x_after_m'],current['x'],0.,'saved post-step geometry')
        rho = array(row['rho_before_kg_m3'],(64,),'operator density')
        close(rho,state['rho'],1e-12,'density history continuity')
        close(row['force_density_kg_m3'],rho,0.,'force uses same density as moment')
        close(row['mass_kg'],current['mass'],0.,'journal mass')
        i,j=np.asarray(row['i'],int),np.asarray(row['j'],int); g=np.asarray(row['g']); B=np.asarray(row['B'])
        sigma=array(row['sigma_mid_Pa'],(64,3,3),'midpoint stress'); force=array(row['internal_force_N'],(64,3),'internal force')
        Lbefore=reconstruct_gradient(vb,i,j,g,current['mass']/rho,B)
        close(row['L_before_per_s'],Lbefore,1e-10,'precontact gradient')
        drho=-rho*np.trace(Lbefore,axis1=1,axis2=2)
        close(row['drho_kg_m3_s'],drho,1e-6,'continuity rate')
        close(current['rho'],rho+dt*drho,1e-10,'actual continuity kick')
        close(row['rho_after_kg_m3'],current['rho'],0.,'saved final density')
        close(row['rho_mid_kg_m3'],rho+.5*dt*drho,1e-10,'midpoint pressure density')
        close(row['S_before_Pa'],state['s'],1e-6,'stress history continuity')
        close(row['S_after_Pa'],current['s'],1e-6,'stress accepted candidate')
        close(sigma,.5*(state['s']+current['s'])-BULK*((rho+.5*dt*drho)/2700-1)[:,None,None]*np.eye(3),1e-6,'midpoint stress used by force')
        acoustic=.1*.00065/(np.sqrt((BULK+4*G/3)/2700)+np.max(np.linalg.norm(vb,axis=1)))
        require(dt<=acoustic,'acoustic timestep condition')
        traces=row['traces']; require(len(traces)==64,'all-particle kick trace')
        workv=np.zeros((64,3)); increments=dict.fromkeys(sums,0.); external=np.zeros((64,3)); contact=np.zeros((64,3))
        for a,trace in enumerate(traces):
            t=audit_trace(trace,ids[a],roles[a],xmid[a],vb[a],va[a],force[a],dt,step)
            workv[a]=t['work_velocity']; external[a]=t['jn']+t['jt']+t['reaction']; contact[a]=t['jn']+t['jt']
            for dest,source in [('friction_raw_J','q'),('wall_work_J','wall'),('constraint_work_J','constraint'),
                ('internal_work_J','internal'),('normal_work_J','normal'),('delta_kinetic_J','dk')]: increments[dest]+=t[source]
            negative+=t['negative']
            if t['minimum_raw'] is not None: minimum_q=t['minimum_raw'] if minimum_q is None else min(minimum_q,t['minimum_raw'])
        close(row['work_velocity_m_s'],workv,1e-12,'actual subinterval work-average velocity')
        close(row['endpoint_midpoint_velocity_m_s'],.5*(vb+va),1e-12,'endpoint midpoint distinct diagnostic')
        operator=audit_operator(xmid,workv,current['mass'],rho,i,j,g,B,np.asarray(row['L_midpoint_per_s']),force,sigma)
        maximum_condition=max(maximum_condition,operator['maximum_condition'])
        maximum_power=max(maximum_power,operator['power_residual_W'])
        maximum_force=max(maximum_force,float(np.linalg.norm(force,axis=1).max()))
        maximum_torque=max(maximum_torque,float(np.linalg.norm(operator['net_internal_torque_N_m'])))
        increments['material_old_work_J']=float(dt*np.einsum('n,nij,nij',current['mass']/rho,sigma,Lbefore))
        increments['midpoint_stress_work_J']=dt*operator['stress_power_W']
        # This registered short trajectory stays elastic. Do not infer plastic
        # constitutive validation from its zero plastic work.
        close(current['ep'],np.zeros(64),0.,'registered elastic trajectory domain')
        close(row['plastic_work_J'],np.zeros(64),0.,'elastic plastic work')
        close(row['return_dissipation_J'],np.zeros(64),0.,'elastic return diagnostic')
        old_rate=np.sqrt(2/3*np.sum((.5*(Lbefore+Lbefore.swapaxes(1,2))-
            np.trace(Lbefore,axis1=1,axis2=2)[:,None,None]*np.eye(3)/3)**2,axis=(1,2)))
        close(record['strain_rate_per_s'],old_rate,1e-10,'precontact strain-rate field')
        close(record['ep_rate_per_s'],np.zeros(64),0.,'elastic plastic strain-rate field')
        close(increments['internal_work_J']+increments['midpoint_stress_work_J'],0.,1e-15,'actual spatial work identity')
        step_dp=np.sum(current['mass'][:,None]*(va-vb),axis=0)
        step_dl=np.sum(current['mass'][:,None]*(np.cross(current['x'],va)-np.cross(state['x'],vb)),axis=0)
        angular=np.cross(xmid,external).sum(axis=0)
        close(step_dp,external.sum(axis=0),1e-17,'whole-body external impulse')
        close(step_dl,angular,1e-19,'whole-body external angular impulse')
        close(increments['wall_work_J'],20*np.cross(xmid,contact).sum(axis=0)[2],1e-15,'rotating wall power torque identity')
        total_external_impulse+=external.sum(axis=0); total_external_angular+=angular
        for key,value in increments.items():
            close(row['increments'][key],value,1e-15,'saved increment '+key)
            sums[key]+=value; close(record['cumulative'][key],sums[key],1e-14,'cumulative '+key)
        close(current['K']+sums['friction_raw_J']-sums['wall_work_J']-sums['constraint_work_J']-sums['internal_work_J'],0.,1e-15,'global kinetic/contact balance')
        if step==1:
            ref=independent_first_step(original,dt)
            close(current['x'],ref['position'],2e-17,'analytic first-step position')
            close(va,ref['velocity'],1e-12,'analytic first-step velocity')
            for key,value in [('wall_work_J',ref['wall_work']),('friction_raw_J',ref['q'])]: close(sums[key],value,1e-18,'analytic first-step '+key)
        maximum_energy=max(maximum_energy,abs(current['residual'])); warnings+=int(current['warning'])
        maximum_free_w=max(maximum_free_w,float(np.max(np.abs(va[np.asarray(roles)==0,2]))))
        audit_vtk(directory/f'frame_{step:04d}.vtk',record,original)
        state=current
    require(maximum_force>1e-6 and maximum_free_w>1e-8,'no meaningful internal 3D force transmission')
    require(np.max(np.abs(state['s'][:,0,2]))>1 and np.max(np.abs(state['s'][:,1,2]))>1,'missing xz/yz stress response')
    require(manifest['final_record']==history[-1],'manifest final state versus history')
    return dict(status='PASS_DEFORMABLE_READBACK_NOT_AFSD',classification=CLASSIFICATION,
        production_authorized=False,case=str(directory),steps=count,end_time_s=END,dt_s=dt,frames_checked=count+1,
        source_sha256=manifest['source_sha256'],configuration_sha256=manifest['configuration_sha256'],
        reader_sha256=file_hash(Path(__file__)),manifest_sha256=file_hash(directory/'manifest.json'),
        minimum_neighbor_rank=3,maximum_condition=maximum_condition,maximum_power_residual_W=maximum_power,
        maximum_internal_force_N=maximum_force,maximum_net_internal_torque_Nm=maximum_torque,
        maximum_free_z_velocity_m_s=maximum_free_w,minimum_raw_friction_J=minimum_q,negative_raw_count=negative,
        cumulative=sums,maximum_energy_diagnostic_residual_J=maximum_energy,energy_warning_frames=warnings,
        energy_certified=False,final_energy_diagnostic_residual_J=state['residual'],
        final_xz_stress_max_Pa=float(np.max(np.abs(state['s'][:,0,2]))),final_yz_stress_max_Pa=float(np.max(np.abs(state['s'][:,1,2]))),
        final_ep_max=float(state['ep'].max()),total_external_impulse_Ns=total_external_impulse.tolist(),
        total_external_angular_impulse_Nms=total_external_angular.tolist(),
        checkpoint_validation='HASH_CHECKED_HERE; real nonzero-state restart is separate unit test',
        continuum_accuracy='NOT_VALIDATED; no independent full continuum trajectory reference')


def compare_cases(coarse,fine):
    a,b=audit_case(coarse),audit_case(fine)
    require(a['dt_s']==5e-9 and b['dt_s']==2.5e-9 and a['end_time_s']==b['end_time_s']==END,'same registered endpoints')
    require(a['source_sha256']==b['source_sha256'],'same physics source')
    ca=load(Path(coarse)/'manifest.json')['configuration']; cb=load(Path(fine)/'manifest.json')['configuration']
    require({k:v for k,v in ca.items() if k!='dt_s'}=={k:v for k,v in cb.items() if k!='dt_s'},'same physical problem')
    ha,hb=load(Path(coarse)/'history.json'),load(Path(fine)/'history.json')
    metrics={}; limits=dict(position_m=1e-7,velocity_m_s=1e-2,stress_Pa=1e5,density_kg_m3=1.)
    for field,limit in limits.items():
        defects=[float(np.max(np.abs(np.asarray(x[field])-np.asarray(hb[2*k][field])))) for k,x in enumerate(ha)]
        metrics[field]=dict(common_history_max_absolute=max(defects),endpoint_max_absolute=defects[-1],screen_limit=limit,within_screen=max(defects)<=limit)
    for field in ('friction_raw_J','wall_work_J'):
        defects=[abs(x['cumulative'][field]-hb[2*k]['cumulative'][field]) for k,x in enumerate(ha)]
        metrics[field]=dict(common_history_max_absolute=max(defects),endpoint_max_absolute=defects[-1],screen_limit=1e-8,within_screen=max(defects)<=1e-8)
    return dict(status='COMPARISON_COMPLETE_NOT_CONVERGENCE_CERTIFICATE',production_authorized=False,
        screening_status='WITHIN_SCREEN' if all(v['within_screen'] for v in metrics.values()) else 'SENSITIVE_WARNING',
        common_end_s=END,common_frames=len(ha),metrics=metrics,cases=[a,b],
        note='Same physical geometry/parameters. Two timesteps are not a measured order of convergence or paper validation.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path); parser.add_argument('--coarse',type=Path); parser.add_argument('--fine',type=Path)
    parser.add_argument('--output',type=Path,required=True); args=parser.parse_args()
    if args.output.exists(): raise FileExistsError('fresh readback output required')
    if args.case and (args.coarse or args.fine): parser.error('use --case OR --coarse and --fine')
    if not args.case and not (args.coarse and args.fine): parser.error('case or comparison pair required')
    try:
        result=audit_case(args.case) if args.case else compare_cases(args.coarse,args.fine); code=0
    except Exception as error:
        result=dict(status='FAIL_READBACK',error=str(error),production_authorized=False,reader_sha256=file_hash(Path(__file__))); code=2
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,allow_nan=False)); return code


if __name__=='__main__': raise SystemExit(main())
