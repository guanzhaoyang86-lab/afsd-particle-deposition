"""Q6 audited saved-data comparisons; no solver import, interpolation or re-solve."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from inspect_process_friction_cone import inspect, ROOT


CASES=('oblique_i4','oblique_i8','basis37_i4')


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_verified(path):
    path=Path(path); checked=inspect(path)
    if checked['status']!='PASS_COMPLETE_READBACK' or checked['accepted_time_s']!=.2:
        raise ValueError('comparison requires complete common 0.2 s endpoint')
    report=json.loads((path/'report.json').read_text(encoding='utf-8'))
    config=json.loads((path/'resolved_config.json').read_text(encoding='utf-8'))
    with np.load(path/'history.npz',allow_pickle=False) as source:
        h={key:source[key].copy() for key in source.files}
    if not all(np.isfinite(a).all() for a in h.values()):raise ValueError('nonfinite comparison input')
    return dict(path=path,checked=checked,report=report,config=config,history=h)


def requested_indices(case,times):
    history=case['history']; rows=case['report']['records']; found=[]
    for t in times:
        matches=[i for i,actual in enumerate(history['physical_time_s'])
                 if abs(actual-t)<=1e-14 and (i==0 or rows[i]['event_kind']=='REQUESTED')]
        if len(matches)!=1:raise ValueError('unique saved requested state required; no interpolation')
        found.append(matches[0])
    return np.asarray(found,dtype=int)


def equal_physics(a,b):
    ca,cb=map(lambda c:dict(c['config']),(a,b))
    for c in (ca,cb):
        for field in ('increments_per_phase','basis_angle_deg','tangent_basis_global_yz'):
            c.pop(field)
    if ca!=cb:raise ValueError('changed physical configuration cannot be a sensitivity comparison')
    if a['report']['source_sha256']!=b['report']['source_sha256']:
        raise ValueError('different source versions cannot be a same-physics comparison')
    ah,bh=a['history'],b['history']
    ids=ah['particle_id']; lookup={int(pid):i for i,pid in enumerate(bh['particle_id'])}
    if set(map(int,ids))!=set(lookup):raise ValueError('particle identity mismatch between runs')
    particle=np.array([lookup[int(pid)] for pid in ids])
    pairs=lambda h:list(zip(map(int,h['contact_particle_id']),map(int,h['contact_feature_id'])))
    paira,pairb=pairs(ah),pairs(bh); bypair={pair:i for i,pair in enumerate(pairb)}
    if len(bypair)!=len(pairb) or set(paira)!=set(pairb):raise ValueError('contact identity mismatch between runs')
    contact=np.array([bypair[pair] for pair in paira])
    for key in ('reference_position_m','reference_volume_m3','mass_kg','support_components'):
        if not np.array_equal(ah[key],bh[key][particle]):raise ValueError('changed geometry/mass/boundary: '+key)
    return particle,contact


def compare_pair(a,b):
    particle,contact=equal_physics(a,b)
    times=np.linspace(0.,.2,17); ia,ib=requested_indices(a,times),requested_indices(b,times)
    ah,bh=a['history'],b['history']; ra,rb=a['report']['records'],b['report']['records']
    rows=[]
    for t,i,j in zip(times,ia,ib):
        rows.append(dict(physical_time_s=float(t),
            maximum_displacement_norm_difference_m=float(np.linalg.norm(bh['displacement_m'][j,particle]-ah['displacement_m'][i],axis=1).max()),
            maximum_stress_component_difference_Pa=float(np.max(np.abs(bh['stress_Pa'][j,particle]-ah['stress_Pa'][i]))),
            maximum_tangent_contact_force_norm_difference_N=float(np.linalg.norm(bh['tangent_force_global_yz_N'][j,contact]-ah['tangent_force_global_yz_N'][i],axis=1).max()),
            maximum_contact_force_norm_difference_N=float(np.linalg.norm(bh['contact_force_N'][j,particle]-ah['contact_force_N'][i],axis=1).max()),
            total_friction_difference_J=float(bh['friction_dissipation_J'][j].sum()-ah['friction_dissipation_J'][i].sum()),
            total_wall_work_difference_J=float(rb[j]['wall_work_J']-ra[i]['wall_work_J'])))
    ea,eb=a['report']['events'],b['report']['events']
    def identity(event,permutation=None):
        modes=np.asarray(event['ending_interval_modes']); released=np.asarray(event['released_rows'])
        if permutation is not None:modes,released=modes[permutation],released[permutation]
        return event['kind'],tuple(modes),tuple(released)
    matching=len(ea)==len(eb) and all(identity(x)==identity(y,contact) for x,y in zip(ea,eb))
    events=dict(status='MATCHED_EVENT_INVENTORY' if matching else 'NOT_COMPARABLE_EVENT_INVENTORY',
        count_a=len(ea),count_b=len(eb),maximum_event_time_difference_s=None)
    if matching:
        differences=[y['physical_time_s']-x['physical_time_s'] for x,y in zip(ea,eb)]
        events.update(signed_event_time_differences_s=differences,
            maximum_event_time_difference_s=float(max(map(abs,differences),default=0.)))
    result=dict(case_a=a['path'].name,case_b=b['path'].name,
        difference_convention='case_b minus case_a; nonnegative vector/component norms explicitly named',
        sampling='17 common explicitly saved requested times; no interpolation of event frames',
        common_requested_states=rows,final=rows[-1],event_comparison=events,
        interpretation='Increment sensitivity or coordinate representation check only; no convergence-order or spatial-accuracy claim.')
    return result


def compare(root):
    cases={name:load_verified(Path(root)/name) for name in CASES}
    for name,increments,angle in (('oblique_i4',4,0.),('oblique_i8',8,0.),('basis37_i4',4,37.)):
        r=cases[name]['report']
        if r['increments_per_phase']!=increments or r['basis_angle_deg']!=angle:
            raise ValueError('case name does not match declared program/representation')
    refinement=compare_pair(cases['oblique_i4'],cases['oblique_i8'])
    basis=compare_pair(cases['oblique_i4'],cases['basis37_i4'])
    end=basis['final']; within=(end['maximum_displacement_norm_difference_m']<=2e-9 and abs(end['total_wall_work_difference_J'])<=1e-6)
    basis.update(status='PASS_BASIS_REPRESENTATION' if within else 'FAIL_BASIS_REPRESENTATION',
        preregistered_final_displacement_norm_tolerance_m=2e-9,
        preregistered_total_wall_work_tolerance_J=1e-6,
        caution='Same global geometry/material/wall path; this changes tangent coordinates, NOT physical tool rotation.')
    result=dict(status='COMPARISON_COMPLETE' if within else 'COMPARISON_FAIL',production_authorized=False,
        common_end_time_s=.2,classification='ENGINEERING_APPROXIMATE_FIXED_NORMAL_CONE_COMPARISON',
        increment_sensitivity=refinement,basis_representation=basis,
        source_sha256={str(path.relative_to(ROOT)).replace('\\','/'):digest(path)
                       for path in (Path(__file__),ROOT/'scripts/inspect_process_friction_cone.py')},
        input_sha256={name:{file.name:digest(file) for file in sorted(case['path'].iterdir()) if file.is_file()}
                      for name,case in cases.items()},
        readbacks={name:dict(status=case['checked']['status'],accepted_time_s=case['checked']['accepted_time_s'])
                   for name,case in cases.items()},
        not_run=['spatial_convergence','physical_tool_rotation','dynamic_AFSD','thermal','production'])
    # JSON serialization is the final finite-number guard for scalar diagnostics.
    json.dumps(result,allow_nan=False)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    result=compare(args.root)
    with args.output.open('x',encoding='utf-8') as target:target.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'],increment_final=result['increment_sensitivity']['final'],
        increment_events=result['increment_sensitivity']['event_comparison'],basis_final=result['basis_representation']['final'],
        basis_events=result['basis_representation']['event_comparison'],basis_status=result['basis_representation']['status']),indent=2,allow_nan=False))
    raise SystemExit(0 if result['status']=='COMPARISON_COMPLETE' else 2)
