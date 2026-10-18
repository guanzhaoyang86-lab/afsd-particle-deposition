"""Independent saved-array cross-check; no solver or comparator imports."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def check(ok, message):
    if not ok:
        raise ValueError(message)


def close(a, b, message, atol=1e-14):
    if not np.allclose(a, b, rtol=1e-12, atol=atol):
        raise ValueError(message)


def load(path):
    r = read(path/'report.json')
    check(r['status'] == 'J2_CONE_PATH_COMPLETE' and r['accepted_time_s'] == .2,
          'partial case is not comparable')
    for name, expected in r['artifact_sha256'].items():
        check(Path(name).name == name, 'unsafe artifact name')
        check(hashlib.sha256((path/name).read_bytes()).hexdigest() == expected,
              'artifact changed: '+name)
    with np.load(path/'history.npz', allow_pickle=False) as f:
        h = {k:f[k] for k in f.files}
    check(all(np.isfinite(a).all() for a in h.values()), 'nonfinite saved history')
    times = h['physical_time_s']
    check(np.all(np.diff(times)>0), 'unordered accepted times')
    ix = []
    for k in range(17):
        found = [i for i,t in enumerate(times) if abs(t-k/80.) <= 1e-12]
        check(len(found) == 1, 'ambiguous/missing common time')
        i = found[0]
        check(k == 0 or r['records'][i]['event_kind'] == 'REQUESTED',
              'common endpoint not a requested accepted frame')
        ix.append(i)
    return r, h, np.asarray(ix), read(path/'resolved_config.json')


def signed_work(r, h):
    wall = np.array([row['wall_displacement_m'] for row in r['records']])
    ids = list(h['particle_id'])
    nodes = [ids.index(pid) for pid in h['contact_particle_id']]
    qs, wws, wms, wrs = [], [], [], []
    for i in range(1, len(wall)):
        du = h['displacement_m'][i]-h['displacement_m'][i-1]
        dw = wall[i]-wall[i-1]
        f = (h['contact_force_N'][i]+h['contact_force_N'][i-1])/2
        support = (h['support_force_N'][i]+h['support_force_N'][i-1])/2
        # Deliberately use assembled contact forces rather than stored tau/Q.
        q = np.array([-np.dot(f[j,1:],du[j,1:]-dw[1:]) for j in nodes])
        qs.append(q)
        wws.append(sum(np.dot(f[j],dw) for j in range(len(ids))))
        wms.append(sum(np.dot(f[j],du[j]) for j in range(len(ids))))
        wrs.append(sum(np.dot(support[j],du[j]) for j in range(len(ids))))
    q, ww, wm, wr = map(np.asarray, (qs,wws,wms,wrs))
    q_error = float(np.max(np.abs(q-h['interval_friction_J'][1:])))
    check(q_error <= 1e-14, 'raw signed interval Q mismatch')
    check(np.max(np.abs(np.cumsum(q,axis=0)-h['friction_dissipation_J'][1:])) <= 1e-14,
          'signed cumulative Q mismatch')
    defect = ww-wm-q.sum(axis=1)
    check(np.max(np.abs(defect)) <= 1e-12, 'contact work identity')
    check(q.min() >= -1e-10, 'frozen negative Q gate')
    records = r['records']
    close(np.cumsum(ww), [row['wall_work_J'] for row in records[1:]], 'wall work sum')
    close(np.cumsum(wm), [row['contact_material_work_J'] for row in records[1:]], 'material contact work sum')
    material = np.array([row['stress_work_J']-row['elastic_diagnostic_J']-row['plastic_work_J']
                         for row in records])
    warning_steps = [i for i,value in enumerate(material) if abs(value)>1e-6]
    check(warning_steps == r['material_work_warning_steps'], 'material warning inventory')
    close(material,[row['work_minus_elastic_minus_plastic_J'] for row in records], 'material warning')
    close(h['plastic_work_J'],90e6*h['reference_volume_m3'][None,:]*h['equivalent_plastic_strain'],
          'independent perfect-plastic Y V ep identity')
    return dict(raw_Q_min_J=float(q.min()),raw_negative_components=int(np.count_nonzero(q<0)),
        signed_Q_sum_J=float(q.sum()),Wwall_J=float(ww.sum()),Wmaterial_J=float(wm.sum()),
        Wsupport_J=float(wr.sum()),maximum_Q_reconstruction_error_J=q_error,
        maximum_contact_identity_J=float(np.abs(defect).max()),
        maximum_stress_contact_defect_J=float(np.max(np.abs(np.cumsum(wm+wr)-
            np.array([row['stress_work_J'] for row in records[1:]])))),
        material_warning_J=float(material[-1]),material_warning_steps=warning_steps,
        material_energy_certified=False)


def events(r, h):
    entries, order, groups, phase_groups = {}, [], [], [0,0,0,0]
    for event in r['events']:
        t=event['physical_time_s']
        phase=next(p for p in range(1,5) if t <= .05*p+1e-12)
        keys=[]
        for i,(released,margin) in enumerate(zip(event['released_rows'],event['cone_margin_N'])):
            if not released and abs(margin)<=9e-5:
                prefix=f"{phase}:{int(h['contact_particle_id'][i])}:{int(h['contact_feature_id'][i])}:"
                occurrence=sum(k.startswith(prefix) for k in entries)+1
                key=prefix+str(occurrence)
                entries[key]=t
                keys.append(key)
        check(bool(keys),'no new cone-boundary candidate')
        groups.append(keys)
        order.extend(sorted(keys))
        phase_groups[phase-1]+=1
    return entries,order,groups,phase_groups


def audit(coarse,fine,comparison):
    ra,ha,ia,ca=load(coarse)
    rb,hb,ib,cb=load(fine)
    check(ca.pop('increments_per_phase')==4 and cb.pop('increments_per_phase')==8,'wrong grids')
    check(ca==cb,'physical problem changed')
    for k in ('particle_id','tracer_type','reference_position_m','reference_volume_m3',
              'mass_kg','support_components','contact_particle_id','contact_feature_id'):
        check(np.array_equal(ha[k],hb[k]),'static field changed '+k)
    data=[]
    for r,h,ix in ((ra,ha,ia),(rb,hb,ib)):
        v=h['reference_volume_m3']
        force=h['contact_force_N'][ix].sum(axis=1)
        data.append(dict(displacement_m=h['displacement_m'][ix],contact_resultant_N=force,
            contact_tangent_resultant_N=force[:,1:],stress_Pa=h['stress_Pa'][ix],
            ep=h['equivalent_plastic_strain'][ix],
            volume_mean_ep=np.sum(h['equivalent_plastic_strain'][ix]*v,axis=1)/v.sum(),
            plastic_work_J=h['plastic_work_J'][ix].sum(axis=1),
            friction_work_J=h['friction_dissipation_J'][ix].sum(axis=1),
            wall_work_J=np.array([r['records'][i]['wall_work_J'] for i in ix]),
            stress_work_J=np.array([r['records'][i]['stress_work_J'] for i in ix])))
    # Specification constants, not obtained by importing/calling the comparator.
    specs={'displacement_m':(40e-6,1),'contact_resultant_N':(9000.,1),
        'contact_tangent_resultant_N':(900.,1),'stress_Pa':(90e6,2),
        'ep':(.004,0),'volume_mean_ep':(.004,0),'plastic_work_J':(.36,0),
        'friction_work_J':(.36,0),'wall_work_J':(.36,0),'stress_work_J':(.36,0)}
    result={}
    for name,(scale,rank) in specs.items():
        a,b=data[0][name],data[1][name]
        def norm(x):
            return np.linalg.norm(x,axis=tuple(range(x.ndim-rank,x.ndim))) if rank else np.abs(x)
        denominator=max(norm(a).max(),norm(b).max(),.01*scale)
        per_time=np.array([norm((b-a)[i:i+1]).max() for i in range(17)])
        result[name]=dict(denominator=float(denominator),maximum_absolute_difference=float(per_time.max()),
            maximum_normalized_difference=float(per_time.max()/denominator),
            final_absolute_difference=float(per_time[-1]),final_normalized_difference=float(per_time[-1]/denominator),
            per_time_absolute_difference=per_time.tolist(),per_time_normalized_difference=(per_time/denominator).tolist())
    ea,oa,ga,pa=events(ra,ha)
    eb,ob,gb,pb=events(rb,hb)
    check(set(ea)==set(eb),'event inventory differs')
    shifts={key:eb[key]-ea[key] for key in ea}
    work=[signed_work(ra,ha),signed_work(rb,hb)]
    # Only now load the tested postprocessor's result for cross-comparison.
    published=read(comparison)
    check(published['status']=='COMPARISON_COMPLETE','comparison incomplete')
    check(published['matched_indices']==[ia.tolist(),ib.tolist()],'common frame mapping mismatch')
    for name,row in result.items():
        for key,value in row.items():
            close(value,published['metrics'][name][key],'independent metric mismatch '+name+'/'+key)
        expected='WITHIN_SCREEN' if row['maximum_normalized_difference']<=.05 else 'SENSITIVE_WARNING'
        check(published['metrics'][name]['status']==expected,'screen mismatch '+name)
    close(max(abs(s) for s in shifts.values()),published['events']['maximum_paired_absolute_shift_s'],'event shift')
    for pair in published['events']['paired']:
        close(pair['signed_difference_s'],shifts[pair['key']],'event pair difference')
    check(published['events']['sequence_equal']==(oa==ob),'event order mismatch')
    return dict(status='PASS_INDEPENDENT_ARRAY_CROSSCHECK_NOT_CONTINUUM_VALIDATION',
        imports='stdlib + NumPy only; neither solver nor tested comparator imported/called',
        common_frames=17,accepted_steps=[ra['accepted_steps'],rb['accepted_steps']],
        endpoint_s=[ra['accepted_time_s'],rb['accepted_time_s']],metrics=result,
        events=dict(groups=[len(ga),len(gb)],phase_group_counts=[pa,pb],entries=[len(ea),len(eb)],
            same_sequence=oa==ob,same_grouping=ga==gb,maximum_shift_s=max(abs(s) for s in shifts.values()),
            signed_shifts_s=shifts),signed_work=work,
        final_values=[{k:v[-1].tolist() for k,v in d.items()} for d in data],
        not_validated=['continuum solution','convergence order','material energy closure',
            'normal opening/recontact','thermal','rotating tool','AFSD deposition','paper quantitative agreement'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--coarse',type=Path,required=True)
    p.add_argument('--fine',type=Path,required=True)
    p.add_argument('--comparison',type=Path,required=True)
    args=p.parse_args()
    print(json.dumps(audit(args.coarse,args.fine,args.comparison),indent=2,allow_nan=False))
