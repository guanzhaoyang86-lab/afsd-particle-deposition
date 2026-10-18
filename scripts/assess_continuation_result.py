"""Independent saved-state continuation audit; no physical advancement.

Only published accepted_history entries count. References use saved arrays,
explicit work/impulse/heat algebra, and fixed-plane crossing definitions.
PySPH is used only to deserialize checkpoints, never to solve an equation.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import sys
import time
import numpy as np
from pysph.solver.utils import load

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from afsd_scaffold.vtk_legacy import read_point_cloud

CLASSIFICATION = 'UPDATED_STAGE_ANDERSON_CONTINUATION_APPROXIMATE'
DEFAULT = ROOT/'outputs/nonlinear_continuation_20261003/continuation'
TOL = dict(position_m=1e-12, velocity_m_s=1e-9, density_kg_m3=1e-10,
           temperature_K=5e-12, mean_velocity_m_s=1e-9)


def require(value, message):
    if not value: raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    with (gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz')
          else Path(path).open('r', encoding='utf-8')) as handle:
        return json.load(handle)


def array(value):
    a = np.asarray(value, dtype=float)
    require(np.isfinite(a).all(), 'nonfinite evidence')
    return a


def close(a, b, name, atol=1e-18, rtol=1e-11):
    a, b = array(a), array(b)
    require(a.shape == b.shape and np.all(np.abs(a-b) <= atol+rtol*np.abs(b)), name)


def exact(a, b, name):
    close(a, b, name, 0., 0.)


def clean(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [clean(v) for v in value]
    return value


def xyz(pa): return np.column_stack((pa.x, pa.y, pa.z))
def velocity(pa): return np.column_stack((pa.u, pa.v, pa.w))


def state_hash(pa):
    h = hashlib.sha256()
    for key in sorted(pa.properties):
        a = pa.properties[key].get_npy_array()
        h.update(key.encode()); h.update(a.dtype.str.encode())
        h.update(str(a.shape).encode()); h.update(a.tobytes())
    return h.hexdigest()


def stage(case, entry, before, before_time, p, manifest):
    directory = (case/entry['snapshot']).parent
    require(directory.resolve().is_relative_to(case), 'unsafe published stage path')
    require(sha(directory/'accepted.npz') == entry['snapshot_sha256'], 'published checkpoint hash')
    data = load(str(directory/'accepted.npz')); after, s = data['arrays']['solid'], data['solver_data']
    require(s['classification'] == CLASSIFICATION and s['accepted'] is True and
            s['restart_allowed'] is True and s['production_authorized'] is False, 'checkpoint acceptance scope')
    require(s['source_sha256'] == manifest['source_sha256'] and s['parameters'] == p and
            s['configuration_sha256'] == manifest['configuration_sha256'], 'checkpoint provenance')
    require(s['particle_state_sha256'] == entry['particle_state_sha256'] == state_hash(after), 'full property history hash')
    exact(s['particle_ids'], after.particle_id, 'checkpoint ID metadata')
    n = len(before.x); m = before.m; ids = before.particle_id.astype(np.int64)
    require(len(np.unique(ids)) == n, 'unique persistent IDs')
    for key in ('particle_id', 'component_id', 'material_id', 'tracer_type', 'wholebar_tag',
                'as_rs_tag', 'actuator', 'm', 'h', 'x0', 'y0', 'z0'):
        exact(getattr(after, key), getattr(before, key), 'immutable '+key)
    payload = read(directory/'accepted_trial.json.gz'); c, e = payload['candidate'], payload['evaluation']
    require(e['accepted'] is False and e['restart_allowed'] is False and
            e['material_scope'] == 'UNCHANGED_FKS_APPROXIMATE_PARAMETERS', 'wrong trial provenance')
    require(e['accepted_state_hash'] == state_hash(before), 'same accepted history input')
    exact(ids, e['particle_id'], 'trial ID correspondence')
    dt = e['dt_s']; require(dt == manifest['settings']['dt_s'] == 5e-8, 'unscaled physical dt')
    exact(e['time_before_s'], before_time, 'stage start'); exact(e['time_after_s'], before_time+dt, 'stage end')
    exact(s['t'], e['time_after_s'], 'checkpoint time'); exact(entry['time_s'], s['t'], 'published time')
    require(entry['stage_index'] == s['stage_index'], 'stage metadata index')
    c = {k: array(v) for k, v in c.items()}; mapping = {k: array(v) for k, v in e['map'].items()}
    norms = {k: float(np.max(np.abs(c[k]-mapping[k])/(before.rho if k == 'density_kg_m3' else 1.))) for k in TOL}
    scaled = max(norms[k]/TOL[k] for k in TOL); require(scaled <= 1., 'accepted global root residual')
    for actual, value, name in ((xyz(after), c['position_m'], 'position'),
        (velocity(after), c['velocity_m_s'], 'velocity'), (after.rho, c['density_kg_m3'], 'rho'),
        (after.temperature, c['temperature_K'], 'temperature'),
        (after.S, np.asarray(e['material']['deviatoric_stress_Pa']).ravel(), 'stress'),
        (after.ep, e['material']['ep'], 'plastic strain')): exact(actual, value, 'saved candidate '+name)
    mid = e['midpoint']; xm = array(mid['position_m']); rho = array(mid['density_kg_m3']); V = m/rho
    exact(xm, .5*(xyz(before)+xyz(after)), 'midpoint position')
    exact(rho, .5*(before.rho+after.rho), 'midpoint density')
    exact(mid['temperature_K'], .5*(before.temperature+after.temperature), 'midpoint temperature')
    exact(mid['mass_kg'], m, 'unscaled midpoint mass')
    L = array(mid['velocity_gradient_per_s']); sigma = array(mid['stress_Pa']); F = array(e['internal_force_N'])
    mean = c['mean_velocity_m_s']; power = float(np.sum(F*mean)+np.sum(V[:, None, None]*sigma*L))
    scale = float(np.sum(np.linalg.norm(F, axis=1)*np.linalg.norm(mean, axis=1)))
    require(abs(power) <= 1e-12+1e-10*scale, 'independent spatial power')
    netF, netTorque = F.sum(axis=0), np.cross(xm, F).sum(axis=0)
    require(np.linalg.norm(netF) <= 1e-12+1e-10*np.linalg.norm(F, axis=1).sum(), 'net internal force')
    require(np.linalg.norm(netTorque) <= 1e-14+1e-10*np.linalg.norm(np.cross(xm,F),axis=1).sum(), 'net internal torque')
    div = np.trace(L, axis1=1, axis2=2)
    close(mapping['density_kg_m3'], before.rho*(1-.5*dt*div)/(1+.5*dt*div), 'continuity', atol=1e-12)
    dep = after.ep-before.ep; require(np.all(dep >= 0), 'nondecreasing plastic history')
    S = after.S.reshape(n,3,3); P = np.sqrt(1.5*np.sum(S*S, axis=(1,2)))*dep*V
    R = 1.5*p['shear']*dep*dep*V
    close(e['material']['plastic_work_J'], P, 'plastic work increment')
    close(e['material']['return_dissipation_J'], R, 'return diagnostic increment')
    close(after.plastic_work, before.plastic_work+P, 'single cumulative plastic increment')
    close(after.return_dissipation, before.return_dissipation+R, 'single cumulative return increment')
    close(after.ep_rate, dep/dt, 'plastic rate', atol=1e-12)
    exact(after.D, (.5*(L+L.swapaxes(1,2))).ravel(), 'saved symmetric strain rate')
    b = e['boundary']; Q = np.zeros(n); vbmean = np.zeros((n,3)); v_end = vbmean.copy()
    flux = np.zeros(n); interval_count = contact_count = 0; seen = set(); momentum_max = 0.
    total_wall = total_other = total_constraint = total_normal_relative = 0.; minimum_Q = None
    impulses = np.zeros((n,3)); v0 = velocity(before)
    cap = (before.z0 > p['bar_length']-p['actuator_depth']) & (before.z0 <= p['bar_length'])
    for trace in b['traces']:
        k = int(trace['particle_index']); require(k not in seen and trace['particle_id'] == ids[k], 'contact trace ID')
        seen.add(k); rows = trace['contact_rows']; nc = len(rows); contact_count += nc
        normals = np.asarray([r['normal'] for r in rows], float).reshape(nc,3)
        walls = np.asarray([r['wall_velocity'] for r in rows], float).reshape(nc,3)
        fn = np.asarray([r['normal_force_N'] for r in rows], float)
        free = np.array([1.,1.,0. if cap[k] and e['time_mid_s'] < p['release_time'] else 1.])
        exact(free, trace['free_axes'], 'only original top axial constraint')
        for row in rows:
            require(row['particle_id'] == ids[k] and row['contact_id'] == f"{ids[k]}:{row['feature_id']}", 'contact persistent feature ID')
            flux[k] += p['contact_conductance']*p['dx']**2*(row['wall_temperature_K']-mid['temperature_K'][k])
        previous = v0[k].copy(); covered = 0.; Jn = np.zeros((nc,3)); Jt = Jn.copy(); reaction = np.zeros(3)
        frame = walls[0] if nc else v0[k]
        speed = max(np.linalg.norm(v0[k]-frame), np.linalg.norm(walls-frame) if nc else 0.,
            dt*np.linalg.norm(free/m[k]*(F[k]+(fn[:,None]*normals).sum(axis=0))), np.finfo(float).tiny)
        for a in trace['result']['intervals']:
            h = float(a['dt_s']); require(h > 0, 'positive contact interval'); interval_count += 1
            close(a['start_s'], covered, 'continuous event time', atol=128*np.finfo(float).eps*dt)
            va, vz = array(a['v_before']), array(a['v_after']); vm = .5*(va+vz)
            close(va, previous, 'continuous event velocity', atol=1e-12)
            jn = array(a['normal_impulse']).reshape(nc,3); jt = array(a['tangential_impulse']).reshape(nc,3)
            jr = array(a['constraint_impulse']); close(jn,h*fn[:,None]*normals,'normal impulse')
            ut = vm-walls; ut -= np.sum(ut*normals,axis=1)[:,None]*normals
            q = -np.sum(jt*ut,axis=1); close(a['friction_heat_raw_J'],q,'raw signed friction',atol=1e-24,rtol=1e-10)
            rounding = 256*np.finfo(float).eps*np.linalg.norm(jt,axis=1)*(speed+np.linalg.norm(va-frame)+np.linalg.norm(vz-frame)+2*np.linalg.norm(walls-frame,axis=1))
            require(np.all(q >= -rounding-1e-29), 'negative friction beyond independent roundoff')
            defect = m[k]*(vz-va)-h*F[k]-jn.sum(axis=0)-jt.sum(axis=0)-jr
            momentum_max = max(momentum_max,float(np.abs(defect).max()))
            require(np.abs(defect).max() <= 1e-17 and np.abs(free*jr).max() <= 1e-17, 'interval momentum/free reaction')
            vbmean[k] += h*vm/dt; Q[k] += q.sum(); Jn += jn; Jt += jt; reaction += jr
            total_wall += float(np.sum((jn+jt)*walls)); total_other += float(h*F[k]@vm)
            total_constraint += float(jr@vm); total_normal_relative += float(np.sum(jn*(vm-walls)))
            if nc: minimum_Q = float(q.min()) if minimum_Q is None else min(minimum_Q,float(q.min()))
            previous = vz; covered += h
        close(covered,dt,'event interval coverage',atol=128*np.finfo(float).eps*dt)
        close(np.asarray(trace['result']['normal_impulse']).reshape(nc,3),Jn,'accumulated normal impulse')
        close(np.asarray(trace['result']['tangential_impulse']).reshape(nc,3),Jt,'accumulated tangent impulse')
        close(trace['result']['constraint_impulse'],reaction,'accumulated constraint impulse')
        v_end[k] = previous; impulses[k] = Jn.sum(axis=0)+Jt.sum(axis=0)+reaction
    require(len(seen) == n and interval_count == b['ledger']['contact_subintervals'] and
            contact_count == b['ledger']['contacts'], 'complete contact/particle count')
    close(b['friction_raw_J'],Q,'per particle friction'); close(b['work_velocity'],vbmean,'event average',atol=1e-12)
    close(mapping['velocity_m_s'],v_end,'endpoint map',atol=1e-12)
    close(mapping['mean_velocity_m_s'],vbmean,'mean map',atol=1e-12)
    close(mapping['position_m'],xyz(before)+dt*vbmean,'true event-integrated position')
    for key, expected in [('friction_raw_J',Q.sum()),('wall_work_J',total_wall),('internal_work_J',total_other),
        ('actuator_work_J',total_constraint),('normal_relative_work_J',total_normal_relative)]:
        close(b['ledger'][key],expected,'work ledger '+key,rtol=1e-9)
    dKmap = .5*np.sum(m[:,None]*(v_end*v_end-v0*v0))
    contact_work_defect = dKmap+Q.sum()-total_wall-total_constraint-total_other-total_normal_relative
    workscale = sum(abs(v) for v in (dKmap,Q.sum(),total_wall,total_constraint,total_other,total_normal_relative))
    require(abs(contact_work_defect) <= 1e-15+1e-9*workscale,'contact work identity')
    close(b['contact_power_W'],flux,'signed wall conduction',atol=1e-14)
    i,j=np.asarray(mid['i'],int),np.asarray(mid['j'],int); g=array(mid['g']); rij=xm[i]-xm[j]
    coefficient=2*p['conductivity']*V[i]*V[j]*(-np.sum(rij*g,axis=1))/np.sum(rij*rij,axis=1)
    require(np.all(coefficient>=0),'conduction direction'); pair=coefficient*(before.temperature[j]-before.temperature[i])
    conduction=np.zeros(n);np.add.at(conduction,i,pair);np.add.at(conduction,j,-pair)
    heat=p['plastic_heat_fraction']*P+p['heat_partition']*Q+dt*(flux+conduction)
    close(e['thermal']['pair_power_W'],conduction,'pair conduction',atol=1e-14)
    close(e['thermal']['particle_energy_J']['total_heat_input_J'],heat,'independent heat, no return/residual')
    close(mapping['temperature_K'],before.temperature+heat/(m*p['cp']),'thermal map',atol=1e-12)
    heat_residual=float(np.sum(m*p['cp']*(after.temperature-before.temperature)-heat))
    require(abs(heat_residual)<=1e-12+1e-9*np.abs(heat).sum(),'accepted mcDT balance')
    transport=read(directory/'transport.json'); order=np.argsort(ids); bar=before.component_id==1
    down=bar&(before.z>=p['layer'])&(after.z<p['layer']);up=bar&(before.z<p['layer'])&(after.z>=p['layer'])
    for name,mask in [('downward',down),('upward',up)]:
        row=transport['crossing'][name]; exact(row['particle_ids'],sorted(ids[mask].tolist()),'crossing IDs '+name)
        exact(row['particle_count'],int(mask.sum()),'crossing count '+name);close(row['mass_kg'],m[mask].sum(),'crossing mass '+name)
    exact(transport['crossing']['plane_z_m'],p['layer'],'frozen sampling plane')
    exact(transport['per_particle']['particle_id'],ids[order],'transport ID order')
    for key,ref in [('position_before_m',xyz(before)),('position_after_m',xyz(after)),('lab_displacement_m',xyz(after)-xyz(before)),
        ('equivalent_plastic_strain_increment',dep)]:close(transport['per_particle'][key],ref[order],'transport '+key)
    tool_displacement=xyz(after)-xyz(before)-np.array([p['traverse']*dt,0.,0.])
    close(transport['per_particle']['tool_translation_displacement_m'],tool_displacement[order],'translation-only frame',atol=1e-18)
    pose=read(directory/'tool_pose.json');close(pose['centre_m'],[p['traverse']*s['t'],0.,p['layer']],'tool translation')
    close(pose['angle_rad'],p['omega']*s['t'],'tool angle');require(pose==s['tool_pose'],'checkpoint tool pose')
    vtk=read_point_cloud(directory/'accepted.vtk');exact(vtk['points'],xyz(after),'VTK positions')
    exact(vtk['vectors']['velocity_m_per_s'],velocity(after),'VTK velocity')
    pr=2*p['shear']*(1+p['poisson'])/(3*(1-2*p['poisson']))*(after.rho/p['rho0']-1)
    scalar=dict(particle_id=after.particle_id,material_id=after.material_id,component_id=after.component_id,
        tracer_type=after.tracer_type,wholebar_tag=after.wholebar_tag,actuator=after.actuator,as_rs_tag=after.as_rs_tag,
        mass_kg=m,density_kg_per_m3=after.rho,temperature_K=after.temperature,pressure_Pa=pr,
        equivalent_stress_Pa=np.sqrt(1.5*np.sum(S*S,axis=(1,2))),equivalent_plastic_strain=after.ep,
        cumulative_plastic_work_J=after.plastic_work,cumulative_return_diagnostic_J=after.return_dissipation)
    for a,baxis,label in ((0,0,'xx'),(1,1,'yy'),(2,2,'zz'),(0,1,'xy'),(0,2,'xz'),(1,2,'yz')):
        scalar['stress_'+label+'_Pa']=S[:,a,baxis]-(pr if a==baxis else 0.)
    for key,value in scalar.items():exact(vtk['scalars'][key],value,'VTK '+key)
    checks=read(directory/'acceptance.json');close(checks['maximum_scaled_residual'],scaled,'accepted residual metadata',atol=0.)
    freebar=bar&~cap
    report=dict(stage_index=entry['stage_index'],time_after_s=s['t'],dt_s=dt,root_residuals=norms,
        maximum_scaled_residual=scaled,contact_subintervals=interval_count,active_contacts=contact_count,
        interval_momentum_max_Ns=momentum_max,raw_friction_J=float(Q.sum()),minimum_raw_friction_J=minimum_Q,
        contact_work_residual_J=float(contact_work_defect),thermal_residual_J=heat_residual,
        spatial_power_residual_W=power,plastic_increment_J=float(P.sum()),return_diagnostic_increment_J=float(R.sum()),
        plastic_increment_ids=ids[dep>0].tolist(),continued_plastic_ids=ids[(before.ep>0)&(dep>0)].tolist(),
        free_bar_max_displacement_m=float(np.linalg.norm(xyz(after)[freebar]-xyz(before)[freebar],axis=1).max(initial=0.)),
        downward_mass_kg=float(m[down].sum()),upward_mass_kg=float(m[up].sum()),
        checkpoint_vtk_transport_pose='PASS_SAVED_DATA_CORRESPONDENCE',free_energy_certification='NOT_RUN')
    return after,s,report,e,transport


def audit(case):
    case=Path(case).resolve();require(case.is_relative_to(ROOT),'project-local run required')
    manifest=read(case/'manifest.json'); require(manifest['classification']==CLASSIFICATION and
        manifest['mode']=='BOUNDED_CONTINUATION' and manifest['production_authorized'] is False,'case scope')
    require(manifest['status'] in ('COMPLETE_REGISTERED_FOUR_STAGES_NOT_DEPOSITION','BLOCKED_WITH_ACCEPTED_PREFIX',
        'BLOCKED_NO_ACCEPTED_PROGRESS'),'run not finished; do not read mutable evidence')
    watched={str((case/'manifest.json').relative_to(ROOT)):sha(case/'manifest.json')}
    for section in ('source_sha256','input_sha256'):
        require(bool(manifest[section]),'missing provenance')
        for name,digest in manifest[section].items():
            path=(ROOT/name).resolve();require(path.is_relative_to(ROOT) and sha(path)==digest,'source/input drift '+name)
            watched[str(path.relative_to(ROOT))]=digest
    for name,digest in manifest['artifact_sha256'].items():
        path=(case/name).resolve();require(path.is_relative_to(case) and sha(path)==digest,'artifact drift '+name)
        watched[str(path.relative_to(ROOT))]=digest
    resolved=read(case/'resolved_config.json');p=resolved['parameters']
    config=hashlib.sha256(json.dumps(dict(classification=CLASSIFICATION,parameters=p,settings=manifest['settings']),
        sort_keys=True,allow_nan=False).encode()).hexdigest()
    require(config==manifest['configuration_sha256']==resolved['configuration_sha256'],'configuration hash')
    require(manifest['settings']==resolved['settings'],'settings mismatch')
    ref=read(case/'last_accepted_reference.json');parent=(ROOT/ref['path']).resolve()
    require(parent.is_relative_to(ROOT) and sha(parent)==ref['sha256'],'parent checkpoint')
    source=load(str(parent));pa=source['arrays']['solid'];old=source['solver_data']
    require(old['accepted'] is True and old['count']==62 and old['parameters']==p,'true unchanged source step62')
    time_s=old['t'];initial_xyz=xyz(pa);initial_ep=pa.ep.copy();initial_hash=state_hash(pa)
    history=read(case/'accepted_history.json');require(len(history)==manifest['accepted_steps'],'published accepted count')
    reports=[];totalQ=totalP=0.;total_intervals=0;boundary_sum={};thermal_sum={};down_total=up_total=0.
    for number,entry in enumerate(history,1):
        require(entry['stage_index']==number,'contiguous published history')
        pa,s,row,e,transport=stage(case,entry,pa,time_s,p,manifest);time_s=s['t'];reports.append(row)
        totalQ+=row['raw_friction_J'];totalP+=row['plastic_increment_J'];total_intervals+=row['contact_subintervals']
        down_total+=row['downward_mass_kg'];up_total+=row['upward_mass_kg']
        for kind,collection in [('boundary',boundary_sum),('thermal',thermal_sum)]:
            for key,value in e[kind]['ledger'].items():
                if kind=='boundary' and key in ('normal_potential_J','max_penetration_m','contacts','minimum_raw_friction_J'):continue
                collection[key]=np.asarray(collection.get(key,0.))+np.asarray(value)
            for key,value in collection.items():close(s['ledger'][kind][key],value,'cumulative '+kind+' '+key)
        close(s['transport_cumulative']['downward_event_mass_kg'],down_total,'cumulative down mass')
        close(s['transport_cumulative']['upward_event_mass_kg'],up_total,'cumulative up mass')
        require(s['transport_cumulative']==transport['cumulative'],'transport checkpoint cumulative')
    exact(manifest['accepted_time_s'],time_s,'manifest actual accepted time')
    require(len(history)<=manifest['settings']['maximum_steps'] and
            manifest['total_evaluations']<=manifest['settings']['total_evaluation_budget'],'registered process limits')
    rejected=[]
    published={entry['snapshot'] for entry in history}
    for path in case.rglob('REJECTED_candidate.json.gz'):
        item=read(path);require(item['classification']=='REJECTED' and item['accepted'] is False and
            item['restart_allowed'] is False,'rejected candidate mislabelled')
        rejected.append(dict(path=str(path.relative_to(case)),error=item['error'],scope='EXCLUDED_FROM_ACCEPTED_HISTORY'))
    orphan=[str(path.relative_to(case)) for path in case.rglob('accepted.npz') if str(path.relative_to(case)).replace('\\','/') not in published]
    bar=pa.component_id==1;cap=(pa.z0>p['bar_length']-p['actuator_depth'])&(pa.z0<=p['bar_length'])
    free=bar&~cap; displacement=xyz(pa)-initial_xyz
    report=dict(status='PASS_ACCEPTED_CONTINUATION_READBACK_NOT_DEPOSITION' if history else 'PASS_NO_ACCEPTED_STAGES_DOWNSTREAM_NOT_RUN',
        source_status=manifest['status'],accepted_steps=len(history),source_time_s=old['t'],accepted_time_s=time_s,
        physical_increment_s=time_s-old['t'],actual_driver_wall_seconds=manifest['wall_seconds'],
        nonlinear_evaluations=manifest['total_evaluations'],particle_count=len(pa.x),source_files_checked=len(manifest['source_sha256']),
        artifact_files_checked=len(manifest['artifact_sha256']),total_contact_subintervals=total_intervals,
        original_snapshot_hash=initial_hash,final_snapshot_hash=state_hash(pa),stages=reports,rejected_candidates=rejected,
        orphan_snapshots_not_counted=orphan,total_raw_friction_J=totalQ,total_plastic_increment_J=totalP,
        free_bar_max_net_displacement_m=float(np.linalg.norm(displacement[free],axis=1).max(initial=0.)),
        free_bar_mass_weighted_net_displacement_m=(np.sum(pa.m[free,None]*displacement[free],axis=0)/pa.m[free].sum()).tolist(),
        fixed_plane_z_m=p['layer'],downward_event_mass_kg=down_total,upward_event_mass_kg=up_total,
        total_plastic_increase_ids=pa.particle_id[pa.ep>initial_ep].tolist(),
        initially_plastic_and_increased_ids=pa.particle_id[(initial_ep>0)&(pa.ep>initial_ep)].tolist(),
        deposition_status='NOT_ESTABLISHED',macro_transport_status='NOT_ESTABLISHED',
        restart_recompute_equivalence='NOT_RUN',same_end_time_sensitivity='NOT_RUN',production_authorized=False,
        reference='INDEPENDENT_SAVED_ARRAY_ALGEBRA_NO_PHYSICS_SOLVER_CALLS')
    for name,digest in watched.items():require(sha(ROOT/name)==digest,'input changed during readback '+name)
    report['all_input_hashes_unchanged']=True
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--case',type=Path,default=DEFAULT)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    path=args.output.resolve();require(path.is_relative_to(ROOT) and not path.exists(),'fresh output required')
    start=time.perf_counter();code=0
    try: result=audit(args.case)
    except Exception as error: result=dict(status='FAIL_INDEPENDENT_CONTINUATION_READBACK',error=repr(error),production_authorized=False);code=2
    result.update(reader_sha256=sha(__file__),readback_wall_seconds=time.perf_counter()-start)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8') as handle:json.dump(clean(result),handle,indent=2,allow_nan=False)
    print(json.dumps(clean(result),allow_nan=False),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())
