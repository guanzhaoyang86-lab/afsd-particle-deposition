"""One same-parent original/inverse search comparison; never publish a state."""
from pathlib import Path
import argparse
import ast
import copy
import json
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import run_expanded_hot_campaign_window as campaign
import inverse_direction_ab_candidate as complete
import visibility_equilibrium as original_search
import visibility_inverse_equilibrium as inverse_search
from verified_inverse_direction import verified_inverse_direction
from secant_equilibrium import _full_matrix

ROOT=campaign.ROOT
core,post,old,model=campaign.core,campaign.post,campaign.old,campaign.model
require,rel=campaign.require,campaign.relative
PARENT=ROOT/'outputs/expanded_hot_budget_prefix_20261005/window_001/fine/step_012/accepted.npz'
PARENT_SHA='12ea12433454492d3157d3efd03c6d6178934611c3a484f6a3d2647d76659f3b'
AUDIT=ROOT/'outputs/expanded_hot_budget_prefix_20261005/window_001/readback/summary.json'
AUDIT_SHA='d30f8f1f30e58a5c66a6dd6d8069c482210d2980b2884ba4b58762f23b403373'
FAILED_AUDIT=ROOT/'outputs/expanded_hot_campaign_20261005/window_001/readback/summary.json'
FAILED_AUDIT_SHA='89f999f398a9f3de70535500121bbba3a359ab3b39cbd351a0e8b00960e88658'
ALGEBRA=ROOT/'outputs/verified_inverse_direction_20261005/algebra_v1/receipt.json'
ALGEBRA_SHA='9912b64a09c5ce7432ac07d2c40558ad9033291a4189a33c1b060191a597232f'
CONTRACT=ROOT/'docs/INVERSE_DIRECTION_AB_CONTRACT_20261005.md'
CLASSIFICATION='SAME_PARENT_VERIFIED_INVERSE_NUMERICAL_DIAGNOSTIC_NOT_ACCEPTED_V1'
PASS='PASS_SAME_PARENT_INVERSE_NUMERICAL_SCREEN_NOT_PROCESS_ACCEPTANCE'
DT=7.8125e-6
THREADS=('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')
LIMITS=dict(position_max_m=1e-9,temperature_max_K=1e-4,
    stress_pointwise_relative=1e-3,stress_floor_Pa=1e6,plastic_absolute=1e-7,
    density_relative=1e-6,contact_history_max_m=1e-9,
    energy_absolute_J=1e-12,energy_relative=1e-3)


def adapter_static_check():
    """Compare every original candidate operation, excluding declared metadata."""
    def function(path):
        tree=ast.parse(path.read_text(encoding='utf-8'))
        return next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='candidate')
    before=function(ROOT/'scripts/expanded_hot_continuation_candidate.py')
    after=function(Path(complete.__file__))
    require([a.arg for a in after.args.kwonlyargs]==['solver'],'solver must be the only new explicit argument')
    after.args.kwonlyargs=[];after.args.kw_defaults=[]
    metadata=[n for n in after.body if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call)
        and isinstance(n.value.func,ast.Attribute) and isinstance(n.value.func.value,ast.Name)
        and n.value.func.value.id=='nextstate' and n.value.func.attr=='update']
    require(len(metadata)==1,'exactly one detached metadata update expected')
    allowed={'diagnostic_only','diagnostic_classification','diagnostic_solver','diagnostic_parent_sha256',
        'continuation_parent_sha256','window_origin_sha256','implementation_source_sha256','controller_source_sha256'}
    require({k.arg for k in metadata[0].value.keywords}==allowed and not metadata[0].value.args,
        'undeclared metadata or physical state change')
    after.body.remove(metadata[0])
    require(ast.dump(before,include_attributes=False)==ast.dump(after,include_attributes=False),
        'complete candidate adapter differs outside explicit solver/metadata scope')
    return dict(status='PASS_STATIC_COMPLETE_CANDIDATE_BODY',physics_evaluations=0)


def preflight():
    require(core.sha(PARENT)==PARENT_SHA and core.sha(AUDIT)==AUDIT_SHA,'fixed true parent/readback drift')
    c=campaign.parent_context(PARENT,AUDIT)
    require(c['state']['step']==14 and abs(c['state']['time_s']-62.5e-6)<1e-15 and
        len(c['prefix'])==14,'same accepted 62.5 us / 14-record parent required')
    old_frozen=dict(c['state']['implementation_source_sha256'])
    failed=post.read(FAILED_AUDIT)
    require(core.sha(FAILED_AUDIT)==FAILED_AUDIT_SHA and failed['exit_code']==0 and
        failed['status']=='PASS_EXPANDED_HOT_CAMPAIGN_READBACK_NOT_PHYSICAL_VALIDATION' and
        failed['added_accepted_steps']==0 and failed['accepted_global_steps']==14 and
        failed['stages']['fine']['status']=='FAILED_NO_NEW_ACCEPTED_STATE' and failed['window_completed'] is False,
        'old time-screen failure must remain failed and unaccepted')
    campaign.bind_map({},failed['evidence_sha256'])
    algebra=post.read(ALGEBRA)
    require(core.sha(ALGEBRA)==ALGEBRA_SHA and algebra['status']=='PASS_ALGEBRA_AND_MANUFACTURED_SEARCH_ONLY' and
        algebra['exit_code']==0 and algebra['test_count']==15 and algebra['source_unchanged'] is True and
        algebra['physical_solver_evaluations']==0 and algebra['nnps_evaluations']==0,
        'matching bounded algebra qualification required')
    campaign.bind_map({},algebra['source_sha256'])
    static=adapter_static_check()
    refs=dict(c['parent_evidence']);refs.update(failed['evidence_sha256'])
    refs.update({rel(p):core.sha(p) for p in (AUDIT,FAILED_AUDIT,ALGEBRA,ALGEBRA.parent/'tests.log')})
    costs=[]
    for stage in ('coarse','fine'):
        folder=FAILED_AUDIT.parent.parent/stage
        rows=post.read(folder/'candidate_history.json')
        manifest=post.read(folder/'manifest.json')
        require(manifest['candidate_steps']==len(rows)==2,'fixed completed baseline candidate count differs')
        costs.extend(row['wall_seconds'] for row in rows)
        refs.update({rel(p):core.sha(p) for p in (folder/'candidate_history.json',folder/'manifest.json')})
    measured=max(costs)
    require(np.isfinite(measured) and abs(measured-28.01225140000315)<1e-9,
        'registered measured complete-candidate maximum differs')
    per_root=measured+1.;forecast=30.+1.5*2*per_root
    require(forecast<120. and forecast<=450.,'registered A/B forecast exceeds its scope')
    frozen=dict(c['frozen'])
    extras=[Path(__file__).resolve(),Path(complete.__file__),CONTRACT,
        ROOT/'scripts/verified_inverse_direction.py',ROOT/'scripts/visibility_inverse_equilibrium.py',
        ROOT/'scripts/check_verified_inverse_direction.py',ROOT/'scripts/visibility_equilibrium.py',
        ROOT/'scripts/secant_equilibrium.py',ALGEBRA,ALGEBRA.parent/'tests.log']
    campaign.bind_map(frozen,{rel(p):core.sha(p) for p in extras})
    campaign.bind_map(frozen,algebra['source_sha256'])
    protocol=dict(schema=CLASSIFICATION,status='REGISTERED_NUMERICAL_DIAGNOSTIC',
        parent_checkpoint=rel(PARENT),parent_sha256=PARENT_SHA,parent_readback=rel(AUDIT),parent_readback_sha256=AUDIT_SHA,
        parent_time_s=c['state']['time_s'],parent_global_step=14,inherited_energy_records=14,
        independent_origin_sha256=campaign.BRANCH_ORIGIN_SHA,dt_s=DT,candidate_time_s=c['state']['time_s']+DT,
        root_order=['original','inverse'],root_count=2,initial_direction_probe_count=1,
        source_sha256=frozen,parent_implementation_source_sha256=old_frozen,reference_sha256=refs,
        parameters=c['p'],parameter_provenance=c['provenance'],numerical_comparison_limits=LIMITS,
        exact_candidate_adapter_static_check=static,algebra_receipt=rel(ALGEBRA),algebra_receipt_sha256=ALGEBRA_SHA,
        cost=dict(measured_maximum_complete_candidate_s=measured,per_candidate_allowance_s=1.,
            complete_candidate_cost_s=per_root,overhead_and_initial_probe_allowance_s=30.,
            formula='30 + 1.5 * 2 * (measured_maximum_complete_candidate_s + 1)',
            forecast_s=forecast,forecast_registration_limit_s=120.,global_forecast_limit_s=450.,
            one_tangent_probe_not_reused_by_either_root=True),internal_limit_s=540.,external_limit_s=600.,
        accepted=False,restart_allowed=False,accepted_steps=0,accepted_time_advanced_s=0.,
        physical_solver_evaluations_in_preflight=0,nnps_evaluations_in_preflight=0,
        classification=model.CLASSIFICATION,diagnostic_classification=CLASSIFICATION)
    c.update(frozen=frozen,refs=refs,qualification=campaign.QUALIFICATION,ab_classification=CLASSIFICATION,
        per_root=per_root,protocol=protocol)
    return c


def check_context(c):
    campaign.bind_map({},c['frozen']);campaign.bind_map({},c['refs'])
    require(core.sha(PARENT)==PARENT_SHA and core.sha(AUDIT)==AUDIT_SHA and
        model.parameters(.008)==(c['p'],c['provenance']),'fixed parent or parameter drift')


def initial_direction_probe(c,directory,deadline):
    """One real initial full map/tangent; the two subsequent roots rebuild their own."""
    started=time.perf_counter();pa,state=core.restore(PARENT)
    trial=complete.RecordedTrial(pa,state,c['p']);guess=complete.batch.guess_for(trial,DT)
    before=time.perf_counter();result=trial.evaluate(guess,DT);map_time=time.perf_counter()-before
    original_search.journal(directory,trial,result,'INITIAL_UNCOMMITTED_NOT_EQUILIBRIUM')
    old.save_graph(directory,trial,result)
    free=np.asarray(state['free'],bool);before=time.perf_counter()
    A,H,condition=_full_matrix(trial,result,free.ravel());tangent_time=time.perf_counter()-before
    require(time.perf_counter()<deadline,'initial tangent probe consumed complete diagnostic budget')
    rhs=-np.asarray(result['force'])[free]
    A.setflags(write=False);H.setflags(write=False);rhs.setflags(write=False)
    before=time.perf_counter();dense=np.linalg.solve(A,rhs)
    backward=float(np.max(np.abs(A@dense-rhs)));dense_time=time.perf_counter()-before
    limit=1e-9*max(1.,float(np.max(np.abs(rhs))))
    require(np.isfinite(dense).all() and backward<=limit,'original dense initial direction failed its unchanged gate')
    checked,diagnostics=verified_inverse_direction(A,H,rhs,deadline=deadline)
    np.savez_compressed(directory/'same_initial_linear_system_NOT_RESTARTABLE.npz',
        A=A,H=H,rhs=rhs,free=free,particle_id=pa.particle_id,guess_m=guess,
        original_direction_m=dense,inverse_direction_m=checked,accepted=False,restart_allowed=False)
    info=dict(status='PASS_SAME_INITIAL_LINEAR_BACKWARD_GATES_NOT_ROOT_QUALIFICATION',
        matrix_size=len(rhs),condition_infinity_estimate=condition,graph_fingerprint=result['graph_fingerprint'],
        full_map_wall_s=map_time,tangent_and_initial_inverse_wall_s=tangent_time,
        original_dense_and_backward_wall_s=dense_time,verified_inverse=diagnostics,
        original_linear_backward_error_N=backward,linear_backward_limit_N=limit,
        direction_max_absolute_difference_m=float(np.max(np.abs(dense-checked))),
        direction_norm_relative_difference=float(np.linalg.norm(dense-checked)/max(np.linalg.norm(dense),1e-30)),
        initial_guess_max_m=float(np.max(np.abs(guess))),physical_full_map_evaluations=1,
        wall_seconds=time.perf_counter()-started,reused_by_roots=False,accepted=False,restart_allowed=False,
        caveat='Single initial-system timing includes only direction operations; no whole-process speedup inference.')
    core.write(directory/'comparison.json',info)
    return info,guess


def leaf_differences(a,b,path=''):
    """Report every complete-map leaf; numerical limits are applied separately."""
    rows={}
    if isinstance(a,dict) and isinstance(b,dict):
        for key in sorted(set(a)|set(b)):
            name=path+'.'+str(key) if path else str(key)
            if key not in a or key not in b:rows[name]=dict(same_structure=False)
            else:rows.update(leaf_differences(a[key],b[key],name))
    elif isinstance(a,(list,tuple)) and isinstance(b,(list,tuple)):
        if len(a)!=len(b):rows[path]=dict(same_structure=False)
        else:
            for i,(x,y) in enumerate(zip(a,b)):rows.update(leaf_differences(x,y,f'{path}[{i}]'))
    elif isinstance(a,(np.ndarray,np.number,int,float,bool)) and isinstance(b,(np.ndarray,np.number,int,float,bool)):
        x,y=np.asarray(a),np.asarray(b)
        if x.shape!=y.shape:rows[path]=dict(same_structure=False,original_shape=list(x.shape),inverse_shape=list(y.shape))
        else:
            require(np.isfinite(x).all() and np.isfinite(y).all(),'nonfinite comparison field '+path)
            diff=x.astype(float)-y.astype(float)
            rows[path]=dict(same_structure=True,shape=list(x.shape),bitwise_value_equal=bool(np.array_equal(x,y)),
                maximum_absolute_difference=float(np.max(np.abs(diff))) if diff.size else 0.,
                l2_difference=float(np.linalg.norm(diff.ravel())))
    else:rows[path]=dict(same_structure=type(a)==type(b),equal=core.clean(a)==core.clean(b))
    return rows


def comparison(a,b,c,directory):
    pa,pb=a['pa'],b['pa'];ra,rb=a['result'],b['result'];sa,sb=a['state'],b['state']
    require(np.array_equal(pa.particle_id,pb.particle_id),'A/B particle identity/order differs')
    require(sa['time_s']==sb['time_s']==c['protocol']['candidate_time_s'] and sa['step']==sb['step']==15,
        'A/B endpoints differ')
    require(sa['interval_energy']['records'][:-1]==sb['interval_energy']['records'][:-1]==c['prefix'],
        'A/B full accepted prefix was reset')
    propsa={k:p.get_npy_array() for k,p in pa.properties.items()}
    propsb={k:p.get_npy_array() for k,p in pb.properties.items()}
    core.write(directory/'all_particle_property_differences.json',leaf_differences(propsa,propsb))
    core.write(directory/'all_full_map_leaf_differences.json',leaf_differences(ra,rb))
    core.write(directory/'observer_differences.json',dict(original=sa['bore_observation'],inverse=sb['bore_observation'],
        differences=leaf_differences(sa['bore_observation'],sb['bore_observation'])))
    xerr=np.linalg.norm(ra['position']-rb['position'],axis=1)
    serr=np.linalg.norm((ra['stress']-rb['stress']).reshape(len(pa.x),-1),axis=1)
    snorm=np.linalg.norm(ra['stress'].reshape(len(pa.x),-1),axis=1)
    srel=serr/np.maximum(snorm,LIMITS['stress_floor_Pa'])
    eerr=np.abs(ra['ep']-rb['ep']);terr=np.abs(ra['temperature']-rb['temperature'])
    np.savez_compressed(directory/'pointwise_numerical_differences_NOT_RESTARTABLE.npz',
        particle_id=pa.particle_id,position_difference_m=xerr,stress_difference_Pa=serr,
        original_stress_norm_Pa=snorm,stress_relative=srel,plastic_strain_absolute=eerr,
        temperature_difference_K=terr,accepted=False,restart_allowed=False)
    energy_a=core.energy_diagnostics(a['trial'],ra);energy_b=core.energy_diagnostics(b['trial'],rb)
    energy=[]
    for key in sorted(energy_a):
        if key.endswith('_J') and np.ndim(energy_a[key])==0:
            va,vb=float(energy_a[key]),float(energy_b[key])
            limit=LIMITS['energy_absolute_J']+LIMITS['energy_relative']*max(abs(va),abs(vb))
            energy.append(dict(field=key,original_J=va,inverse_J=vb,absolute_difference_J=abs(va-vb),
                numerical_limit_J=limit,passed=abs(va-vb)<=limit))
    metrics=dict(position_max_m=float(xerr.max()),temperature_max_K=float(terr.max()),
        stress_pointwise_relative=float(srel.max()),plastic_absolute=float(eerr.max()),
        density_relative=float(np.max(np.abs(ra['density']-rb['density'])/np.maximum(np.abs(ra['density']),1.))),
        contact_history_max_m=float(np.max(np.linalg.norm(
            (sa['contact_history']-sb['contact_history']).reshape(-1,3),axis=1))))
    same_graph=ra['graph_fingerprint']==rb['graph_fingerprint'] and all(np.array_equal(ra[k],rb[k]) for k in ('i','j','raw_i','raw_j'))
    branches={key:bool(np.array_equal(value,rb['contact'][key])) for key,value in ra['contact'].items()
        if isinstance(value,np.ndarray) and value.dtype.kind in 'biu'}
    immutable=campaign.immutable_hashes(pa,sa)==campaign.immutable_hashes(pb,sb)==campaign.immutable_hashes(c['pa'],c['state'])
    crossing=sa['bore_observation']['crossing_ids']==sb['bore_observation']['crossing_ids']
    passed=all(metrics[key]<=LIMITS[key] for key in metrics) and same_graph and all(branches.values()) and immutable and crossing and all(v['passed'] for v in energy)
    report=dict(status='PASS_NUMERICAL_COMPARISON' if passed else 'FAIL_NUMERICAL_COMPARISON',passed=passed,
        scope='SAME_DT_SAME_PARENT_ALGORITHM_COMPARISON_NOT_TIME_REFINEMENT_OR_PHYSICAL_VALIDATION',
        numerical_limits=LIMITS,metrics=metrics,maximum_stress_error_particle_id=int(pa.particle_id[np.argmax(srel)]),
        exact_graph_match=same_graph,contact_integer_boolean_branch_match=branches,immutable_fields_match=immutable,
        crossing_ids_match=crossing,energy_comparison=energy,
        original_wall_s=a['wall_seconds'],inverse_wall_s=b['wall_seconds'],
        observed_original_over_inverse_wall_ratio=a['wall_seconds']/b['wall_seconds'],
        timing_caveat='One ordered pair including fresh/scalar/geometry/path/I/O; not a distribution or production speedup.',
        accepted=False,restart_allowed=False,accepted_time_advanced_s=0.)
    core.write(directory/'numerical_screen.json',report)
    return report


def worker(directory):
    started=time.perf_counter();deadline=started+525.;rows=[];code=2
    report=dict(status='RUNNING_NUMERICAL_DIAGNOSTIC',accepted=False,restart_allowed=False,
        accepted_steps=0,accepted_time_advanced_s=0.,parent_sha256=PARENT_SHA)
    def stop():
        try:core.write(directory/'HARD_TIMEOUT.json',dict(status='TIMEOUT_NOT_ACCEPTED',candidate_count=len(rows)))
        finally:os._exit(124)
    timer=threading.Timer(540.,stop);timer.daemon=True;timer.start()
    try:
        require(all(os.environ.get(k)=='1' for k in THREADS),'single BLAS/OpenMP thread required')
        c=preflight();protocol=post.read(directory/'protocol.json')
        require(c['protocol']==protocol,'registered source/parent/numerical limits changed')
        report.update(protocol,status='RUNNING_NUMERICAL_DIAGNOSTIC')
        core.write(directory/'manifest.json',report);core.write(directory/'accepted_history.json',[])
        probe=directory/'initial_direction';probe.mkdir();info,expected_guess=initial_direction_probe(c,probe,deadline)
        report['initial_direction_probe']=info;core.write(directory/'manifest.json',report)
        items={}
        for name,solver in (('original',original_search),('inverse',inverse_search)):
            check_context(c)
            require(time.perf_counter()+1.5*c['per_root']+15<=deadline,'insufficient measured budget for next full candidate')
            pa,state=core.restore(PARENT)
            require(state['implementation_source_sha256']==protocol['parent_implementation_source_sha256'] and
                state['interval_energy']['records']==c['prefix'],'reload changed fixed parent/source/history')
            trial=complete.RecordedTrial(pa,state,c['p'])
            require(np.array_equal(complete.batch.guess_for(trial,DT),expected_guess),'root initial guess differs')
            sd=directory/name;sd.mkdir()
            core.write(sd/'diagnostic_parent_binding.json',dict(parent_checkpoint=rel(PARENT),parent_sha256=PARENT_SHA,
                parent_global_step=14,parent_time_s=state['time_s'],dt_s=DT,controller=solver.CONTROLLER,
                accepted=False,restart_allowed=False))
            item=complete.candidate(pa,state,c,DT,sd,deadline,c['per_root'],solver=solver)
            require(item['state']['accepted'] is False and item['state']['restart_allowed'] is False and
                item['state']['diagnostic_only'] is True and item['state']['diagnostic_parent_sha256']==PARENT_SHA,
                'A/B helper attempted publication or changed direct parent')
            require(np.array_equal(complete.batch.guess_for(item['trial'],DT),expected_guess),'parent arrays changed during root')
            row=old.candidate_row(item,1,PARENT_SHA,directory);row['route']=name;rows.append(row);items[name]=item
            core.write(directory/'candidate_history.json',rows)
            report.update(completed_candidates=len(rows),candidate_history=rows)
            core.write(directory/'manifest.json',report)
            print(json.dumps(dict(route=name,candidate_wall_s=item['wall_seconds'],metrics=item['metrics'],accepted=False)),flush=True)
        screen=comparison(items['original'],items['inverse'],c,directory)
        check_context(c);require(not list(directory.rglob('accepted.npz')),'diagnostic produced restartable publication')
        report.update(numerical_screen=screen,roots_completed=True)
        require(screen['passed'],'pre-registered numerical-only A/B screen failed')
        report['status']=PASS;code=0
    except Exception as error:
        report.update(status='FAIL_INVERSE_AB_NOT_ACCEPTED',error=repr(error),traceback=traceback.format_exc())
    finally:
        timer.cancel();report.update(wall_seconds=time.perf_counter()-started,completed_candidates=len(rows),
            peak_working_set_bytes=core.memory_bytes(),accepted=False,restart_allowed=False,accepted_steps=0,
            accepted_time_advanced_s=0.)
        report['artifacts_sha256']={rel(p):core.sha(p) for p in directory.rglob('*') if p.is_file() and
            p.name not in ('manifest.json','run.log','supervision.json')}
        core.write(directory/'manifest.json',report)
    return code


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    parser.add_argument('--preflight-only',action='store_true');parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args();directory=campaign.project(args.output)
    require(directory.is_relative_to(ROOT/'outputs') and directory!=ROOT/'outputs','output must be a new project outputs directory')
    if args.worker:return worker(directory)
    c=preflight();protocol=c['protocol']
    if args.preflight_only:
        print(json.dumps(dict(status=protocol['status'],parent_sha256=PARENT_SHA,parent_time_s=protocol['parent_time_s'],
            parent_global_step=14,dt_s=DT,forecast_s=protocol['cost']['forecast_s'],
            numerical_comparison_limits=LIMITS,physical_solver_evaluations=0,nnps_evaluations=0),indent=2));return 0
    directory.mkdir(parents=True,exist_ok=False);campaign.write_new(directory/'protocol.json',protocol)
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--output',rel(directory),'--worker']
    start=time.perf_counter();code=124
    with (directory/'run.log').open('x',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=ROOT,env=dict(os.environ,**{k:'1' for k in THREADS}),
            stdout=stream,stderr=subprocess.STDOUT,timeout=600.).returncode
        except subprocess.TimeoutExpired:stream.write('\nEXTERNAL_TIMEOUT_NOT_ACCEPTED\n')
    error=None
    try:check_context(c)
    except Exception as failure:error=repr(failure);code=2
    manifest=directory/'manifest.json'
    campaign.write_new(directory/'supervision.json',dict(command=command,exit_code=code,
        wall_seconds=time.perf_counter()-start,source_inputs_still_match=error is None,source_error=error,
        manifest_sha256=core.sha(manifest) if manifest.exists() else None,log_sha256=core.sha(directory/'run.log'),
        internal_limit_s=540.,external_limit_s=600.,accepted=False,restart_allowed=False))
    print(json.dumps(dict(exit_code=code,output=rel(directory),status=post.read(manifest)['status'] if manifest.exists() else 'NO_MANIFEST')))
    return code


if __name__=='__main__':raise SystemExit(main())
