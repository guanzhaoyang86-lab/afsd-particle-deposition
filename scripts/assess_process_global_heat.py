"""Read-only final-journal audit with independent vector/heat algebra.

Imports no AFSD solver, acceptance, heat, material or contact implementation.
Upstream load reads the immutable parent. Never advances a physical state.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from pysph.solver.utils import load

ROOT = Path(__file__).resolve().parents[1]
TOLS = dict(position_m=1e-12, velocity_m_s=1e-9, density_kg_m3=1e-10,
            temperature_K=5e-12, mean_velocity_m_s=1e-9)


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def require(condition, message):
    if not condition: raise ValueError(message)


def close(a, b, atol, rtol, message):
    a = np.asarray(a); b = np.asarray(b)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
        and np.all(np.abs(a-b) <= atol+rtol*np.abs(b)), message)


def audit(case):
    started = time.perf_counter(); manifest_path = case/'manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    require(manifest['status'] != 'RUNNING_UNCOMMITTED', 'WAIT_FINAL_MANIFEST; no partial journal audit')
    watched = {manifest_path: sha(manifest_path)}
    for rel, digest in manifest['source_sha256'].items():
        path = (ROOT/rel).resolve(); require(path.is_relative_to(ROOT) and sha(path) == digest, 'source drift '+rel)
        watched[path] = digest
    for rel, digest in manifest['artifact_sha256'].items():
        path = (case/rel).resolve(); require(path.is_relative_to(case) and sha(path) == digest, 'artifact drift '+rel)
        watched[path] = digest
    parent = ROOT/manifest['parent_reference']['path']
    require(sha(parent) == manifest['parent_reference']['sha256'], 'parent drift')
    watched[parent] = sha(parent)
    old = load(str(parent))['arrays']['solid']
    accepted = manifest['accepted_steps'] == 1 and manifest['accepted_progress'] is True
    history = json.loads((case/'accepted_history.json').read_text(encoding='utf-8'))
    require(len(history) == (1 if accepted else 0), 'published history/status mismatch')
    journal = case/('accepted_trial.json.gz' if accepted else 'latest_trial.json.gz')
    with gzip.open(journal, 'rt', encoding='utf-8') as stream: saved = json.load(stream)
    require(saved['evaluation'] is not None, 'NO_COMPLETED_FINAL_CANDIDATE; local heat audit unavailable')
    e = saved['evaluation']; c = {k: np.asarray(v) for k, v in saved['candidate'].items()}
    mapped = {k: np.asarray(v) for k, v in e['map'].items()}
    p = json.loads((case/'resolved_config.json').read_text(encoding='utf-8'))['parameters']
    require(e['accepted'] is False and e['restart_allowed'] is False, 'journal is a pure evaluated trial')
    require(e['classification'] == 'PROCESS_GLOBAL_EVENT_AVERAGE_BE_THERMAL_APPROXIMATE', 'wrong heat variant')
    ids = old.particle_id; n = len(ids); dt = e['dt_s']; mid = e['midpoint']; boundary = e['boundary']
    close(e['particle_id'], ids, 0., 0., 'ID mismatch')
    xmid = np.asarray(mid['position_m']); rho = np.asarray(mid['density_kg_m3'])
    i = np.asarray(mid['i'], int); j = np.asarray(mid['j'], int)
    gradient = np.asarray(mid['g']); m = old.m; cp = p['cp']; capacity = m*cp
    require(n == 342 and np.linalg.matrix_rank(xmid-xmid.mean(axis=0)) == 3, 'actual 3D parent inventory required')
    # Independently reconstruct radial gradient from the scalar polynomial.
    rij = xmid[i]-xmid[j]; radius = np.linalg.norm(rij, axis=1); h = .5*(old.h[i]+old.h[j]); q = radius/h
    derivative = np.where(q < 1., 3*q*q-6., np.where(q < 2., -3*(2-q)**2, 0.))
    ref_gradient = (15/(62*np.pi*h**4)*derivative/radius)[:, None]*rij
    close(gradient, ref_gradient, 1e-3, 1e-13, 'radial kernel gradient mismatch')
    coefficient = (m[i]/rho[i])*(m[j]/rho[j])*2*p['conductivity']*(-np.sum(rij*ref_gradient, axis=1))/radius**2
    require(np.all(coefficient >= 0) and len(np.unique(np.sort(np.column_stack((i,j)), axis=1),axis=0)) == len(i),
            'invalid undirected conductive graph')
    close(e['thermal']['pair_conductance_W_K'], coefficient, 1e-18, 1e-12, 'conductance units/coefficient')
    H = np.zeros(n); B = np.zeros(n); Q = np.zeros(n); Qround = np.zeros(n); seen = set()
    interval_count = 0; active_faces = 0; min_raw = None; max_q_error = 0.; max_impulse = 0.
    max_map_momentum = 0.; max_candidate_momentum = 0.; max_candidate_mean = 0.
    max_contact_identity = 0.; max_wall_velocity_error = 0.; candidate_momentum_ok = True
    summed = dict(kinetic=0., friction=0., wall=0., actuator=0., internal=0., normal_relative=0.)
    eps = np.finfo(float).eps; v0all = np.column_stack((old.u, old.v, old.w))
    force = np.asarray(e['internal_force_N'])
    for trace in boundary['traces']:
        a = trace['particle_index']; require(a not in seen and trace['particle_id'] == ids[a], 'trace ID duplicate/mismatch')
        seen.add(a); faces = trace['contact_rows']; nc = len(faces); active_faces += nc
        normals = np.asarray([face['normal'] for face in faces]).reshape(nc,3)
        walls = np.asarray([face['wall_velocity'] for face in faces]).reshape(nc,3)
        fn = np.asarray([face['normal_force_N'] for face in faces])
        for face in faces:
            fid = int(face['feature_id']); Tw = p['base_temperature'] if fid == 0 else p['tool_temperature']
            require(face['wall_temperature_K'] == Tw, 'wall reservoir wrong')
            hc = p['contact_conductance']*p['dx']**2; H[a] += hc; B[a] += hc*Tw
            wall_ref = np.zeros(3) if fid == 0 else np.array([p['traverse']-p['omega']*xmid[a,1],
                p['omega']*(xmid[a,0]-p['traverse']*e['time_mid_s']), 0.])
            close(face['wall_velocity'], wall_ref, 1e-14, 1e-13, 'power-collocated rigid wall velocity')
            max_wall_velocity_error = max(max_wall_velocity_error, float(np.max(np.abs(wall_ref-face['wall_velocity']))))
        result = trace['result']; free = np.asarray(trace['free_axes']); prescribed = np.asarray(trace['prescribed_velocity'])
        selected = old.z0[a] > p['bar_length']-p['actuator_depth'] and old.z0[a] <= p['bar_length'] and e['time_mid_s'] < p['release_time']
        require(np.array_equal(free, [1.,1.,0.] if selected else [1.,1.,1.]), 'actuator original-position selection')
        frame = walls[0] if nc else v0all[a]; total_other = force[a]+(fn[:,None]*normals).sum(axis=0)
        speed = max(np.linalg.norm(v0all[a]-frame), np.linalg.norm(walls-frame) if nc else 0.,
                    dt*np.linalg.norm(free*total_other/m[a]), np.finfo(float).tiny)
        event_mean = np.zeros(3); total_impulse = np.zeros(3); covered = 0.; previous = v0all[a]
        for interval in result['intervals']:
            interval_count += 1; local_h = interval['dt_s']; va = np.asarray(interval['v_before']); vb = np.asarray(interval['v_after'])
            require(local_h > 0 and abs(interval['start_s']-covered) <= 128*eps*dt, 'interval coverage')
            close(va, previous, 1e-15, 1e-13, 'interval velocity continuity')
            jn = np.asarray(interval['normal_impulse']).reshape(nc,3); jt = np.asarray(interval['tangential_impulse']).reshape(nc,3)
            reaction = np.asarray(interval['constraint_impulse']); vmid = .5*(va+vb)
            rel = .5*((va-frame)+(vb-frame))-(walls-frame)
            ut = rel-normals*np.sum(normals*rel,axis=1)[:,None]
            qref = -np.sum(jt*ut,axis=1); qraw = np.asarray(interval['friction_heat_raw_J'])
            bound = 256*eps*np.linalg.norm(jt,axis=1)*(speed+np.linalg.norm(va-frame)+
                np.linalg.norm(vb-frame)+2*np.linalg.norm(walls-frame,axis=1))
            close(interval['friction_rounding_bound_J'], bound, 1e-30, 1e-12, 'roundoff bound changed')
            # Arithmetic budget includes independent projection, no heat clipping.
            require(np.all(np.abs(qraw-qref) <= 1e-25+1e-12*np.abs(qref)+bound), 'independent friction impulse work')
            require(np.all(qraw >= -np.asarray(interval['friction_rounding_bound_J'])), 'individual negative raw dissipation')
            Q[a] += qraw.sum(); Qround[a] += bound.sum()
            close(jn, local_h*fn[:,None]*normals, 1e-25, 1e-12, 'normal impulse')
            impulse = jn.sum(axis=0)+jt.sum(axis=0)+reaction
            residual = m[a]*(vb-va)-local_h*force[a]-impulse
            max_impulse = max(max_impulse,float(np.max(np.abs(residual))))
            require(np.max(np.abs(residual)) <= 1e-17+m[a]*1e-9, 'subinterval momentum')
            require(np.max(np.abs(free*reaction)) <= 1e-17, 'free-axis reaction')
            require(np.max(np.abs((1-free)*(vb-prescribed))) <= 1e-12, 'subinterval actuator velocity')
            internal = local_h*float(force[a]@vmid); actuator = float(reaction@vmid)
            wall_work = float(np.sum((jn+jt)*walls)); normal_relative = float(np.sum(jn*(vmid-walls)))
            kinetic = .5*m[a]*float((vb-va)@(vb+va))
            identity = kinetic+float(qraw.sum())-internal-actuator-wall_work-normal_relative
            max_contact_identity = max(max_contact_identity, abs(identity))
            scale = sum(abs(v) for v in (kinetic,float(qraw.sum()),internal,actuator,wall_work,normal_relative))
            require(abs(identity) <= 1e-15+1e-9*scale, 'independent contact work identity')
            for key,val in [('kinetic',kinetic),('friction',float(qraw.sum())),('internal',internal),
                            ('actuator',actuator),('wall',wall_work),('normal_relative',normal_relative)]: summed[key] += val
            if nc:
                min_raw = float(qraw.min()) if min_raw is None else min(min_raw,float(qraw.min()))
                max_q_error = max(max_q_error,float(np.max(np.abs(qraw-qref))))
            event_mean += local_h*vmid/dt; total_impulse += impulse; covered += local_h; previous = vb
        require(abs(covered-dt) <= 128*eps*dt, 'incomplete event coverage')
        close(event_mean, mapped['mean_velocity_m_s'][a], 1e-15, 1e-12, 'mapped event mean')
        map_defect = m[a]*(mapped['velocity_m_s'][a]-v0all[a])-dt*force[a]-total_impulse
        candidate_defect = m[a]*(c['velocity_m_s'][a]-v0all[a])-dt*force[a]-total_impulse
        max_map_momentum = max(max_map_momentum,float(np.max(np.abs(map_defect))))
        max_candidate_momentum = max(max_candidate_momentum,float(np.max(np.abs(candidate_defect))))
        candidate_momentum_ok &= bool(np.max(np.abs(candidate_defect)) <= m[a]*1e-9+1e-17)
        max_candidate_mean = max(max_candidate_mean,float(np.max(np.abs(event_mean-c['mean_velocity_m_s'][a]))))
    require(len(seen) == n, 'missing particle trace')
    close(boundary['friction_raw_J'], Q, 1e-25, 1e-12, 'particle raw heat accumulation')
    close(boundary['wall_row_W_K'], H, 1e-15, 1e-13, 'wall H dimensions/count')
    stress = np.asarray(e['material']['deviatoric_stress_Pa']); ep = np.asarray(e['material']['ep'])
    plastic = (m/rho)*np.sqrt(1.5*np.sum(stress*stress,axis=(1,2)))*(ep-old.ep)
    require(np.all(ep >= old.ep), 'plastic history decrease')
    close(e['material']['plastic_work_J'], plastic, 1e-20, 1e-10, 'plastic q dep V energy unit')
    th = e['thermal']; readings = {}
    for name, T1 in [('mapped',mapped['temperature_K']),('candidate',c['temperature_K'])]:
        require(np.all(T1 >= p['room']) and np.all(T1 < p['melt']), name+' solid domain')
        pairs = dt*coefficient*(T1[j]-T1[i]); internal = np.zeros(n)
        np.add.at(internal,i,pairs); np.add.at(internal,j,-pairs)
        wall = dt*(B-H*T1); source = p['plastic_heat_fraction']*plastic+p['heat_partition']*Q
        thermal_residual = capacity*(T1-old.temperature)-source-internal-wall
        limit = 1e-12+1e-10*float(np.abs(source+internal+wall).sum())
        total_limit = 1e-12+1e-9*float(np.abs(source+internal+wall).sum())
        readings[name] = dict(maximum_node_residual_J=float(np.max(np.abs(thermal_residual))),
            total_residual_J=float(thermal_residual.sum()),per_node_limit_J=limit,total_limit_J=total_limit,
            balance_PASS=bool(np.max(np.abs(thermal_residual)) <= limit and abs(thermal_residual.sum()) <= total_limit),
            pair_net_J=float(internal.sum()),wall_to_material_J=float(wall.sum()),
            thermal_increment_J=float(np.sum(capacity*(T1-old.temperature))),
            temperature_min_K=float(T1.min()),temperature_max_K=float(T1.max()))
        if name == 'mapped':
            require(readings[name]['balance_PASS'],'mapped independent BE node balance')
            close(th['particle_energy_J']['internal_conduction_J'],internal,1e-25,1e-10,'pair ledger')
            close(th['particle_energy_J']['contact_in_J'],wall,1e-25,1e-10,'wall ledger')
            close(th['particle_energy_J']['plastic_heat_J'],p['plastic_heat_fraction']*plastic,1e-20,1e-10,'plastic heat partition')
            close(th['particle_energy_J']['plastic_unthermalized_J'],(1-p['plastic_heat_fraction'])*plastic,1e-20,1e-10,'unthermalized plastic')
            close(th['particle_energy_J']['friction_material_heat_J'],p['heat_partition']*Q,1e-25,1e-12,'friction material fraction')
            close(th['particle_energy_J']['friction_wall_heat_J'],(1-p['heat_partition'])*Q,1e-25,1e-12,'friction wall fraction')
            close(boundary['contact_power_W'],B-H*T1,1e-18,1e-12,'BE wall flux time layer')
            close(np.asarray(th['ledger']['wall_heat_out_J']),np.asarray(float(((1-p['heat_partition'])*Q-wall).sum())),1e-20,1e-10,'global wall heat out')
    norms = {key: float(np.max(np.abs(c[key]-mapped[key])/(old.rho if key == 'density_kg_m3' else 1.))) for key in TOLS}
    max_scaled = max(norms[k]/TOLS[k] for k in TOLS)
    if accepted:
        require(max_scaled <= 1 and candidate_momentum_ok and max_candidate_mean <= 1e-9 and readings['candidate']['balance_PASS'],
                'published accepted state violates independent candidate identity')
    for path,digest in watched.items(): require(sha(path) == digest, 'watched file changed during independent audit')
    return dict(status='PASS_ACCEPTED_LOCAL_HEAT_REVIEW_NOT_PROCESS_VALIDATION' if accepted else
        'PASS_REJECTED_CANDIDATE_LOCAL_HEAT_REVIEW_NOT_ACCEPTANCE', accepted=accepted,
        physical_acceptance_created_by_this_audit=False, parent_time_s=manifest['source_time_s'],
        physical_last_accepted_time_s=manifest['accepted_time_s'],candidate_requested_end_s=e['time_after_s'],
        final_run_status=manifest['status'],journal=journal.relative_to(ROOT).as_posix(),
        particle_count=n,pair_count=len(i),contact_subintervals=interval_count,active_faces=active_faces,
        friction_raw_sum_J=float(Q.sum()),minimum_raw_subinterval_heat_J=min_raw,
        max_independent_friction_work_difference_J=max_q_error,total_friction_rounding_bound_J=float(Qround.sum()),
        max_subinterval_momentum_residual_Ns=max_impulse,max_contact_work_identity_residual_J=max_contact_identity,
        map_momentum_max_residual_Ns=max_map_momentum,candidate_momentum_max_residual_Ns=max_candidate_momentum,
        candidate_momentum_PASS=candidate_momentum_ok,candidate_mean_velocity_max_residual_m_s=max_candidate_mean,
        candidate_global_raw_residual_norms=norms,candidate_maximum_scaled_raw_residual=max_scaled,
        plastic_work_J=float(plastic.sum()),contact_work_sums=summed,thermal=readings,
        maximum_wall_velocity_difference_m_s=max_wall_velocity_error,
        claimed_scope='INDEPENDENT_FORMULA_AUDIT; AUTHOR_IS_ALSO_VARIANT_IMPLEMENTER; NOT_EXTERNAL_PEER_REVIEW',
        source_sha256={str(Path(__file__).relative_to(ROOT).as_posix()):sha(__file__)},
        input_sha256={path.relative_to(ROOT).as_posix():digest for path,digest in watched.items()},
        wall_seconds=time.perf_counter()-started,production_authorized=False)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();target=a.output.resolve();case=a.case.resolve()
    require(target.is_relative_to(ROOT) and not target.exists() and case.is_relative_to(ROOT),'fresh local audit output required')
    result=audit(case);target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('input_sha256',)},allow_nan=False),flush=True)


if __name__=='__main__':main()
