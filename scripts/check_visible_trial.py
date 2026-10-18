"""Same-input full-map diagnostics; no solve, acceptance, or time advancement."""
from pathlib import Path
import json
import time
import numpy as np
from pysph.solver.utils import load
import run_batched_process as batch
import visible_quasistatic as visible
from afsd_backend.pysph.discrete_contact import ChainRuleDiscreteContact

core,post=batch.core,batch.post
ROOT=core.ROOT
BASE=ROOT/'outputs/visible_support_20261004/qualification_v3'
PARENT=ROOT/'outputs/coupled_cost_20261004/batched/pre_occlusion_response/fine/step_032/accepted.npz'
PARENT_SHA='553f57b1dd5f0555d5a3726f7942705bd485ffd688ef3bf34f16ff342fd9a804'
OLD=ROOT/'outputs/coupled_cost_20261004/batched/post_passage_short_response/coarse'


def require(condition,message):
    if not condition:raise ValueError(message)


def equal_nested(a,b,path):
    if isinstance(a,dict):
        require(isinstance(b,dict) and a.keys()==b.keys(),'nested keys differ '+path)
        for key in a:equal_nested(a[key],b[key],path+'.'+str(key))
    elif isinstance(a,(list,tuple)):
        require(len(a)==len(b),'nested lengths differ '+path)
        for k,(left,right) in enumerate(zip(a,b)):equal_nested(left,right,path+'.'+str(k))
    elif isinstance(a,np.ndarray):
        require(a.shape==np.asarray(b).shape,'array shape differs '+path)
        require(np.array_equal(a,b) if a.dtype.kind in 'biu' else
            np.allclose(a,b,rtol=1e-10,atol=1e-11,equal_nan=False),'nested array differs '+path)
    elif isinstance(a,(float,np.floating)):
        require(np.isclose(a,b,rtol=1e-10,atol=1e-11,equal_nan=False),'nested number differs '+path)
    else:require(a==b,'nested value differs '+path)


def same_map(a,b,old_trial,new_trial):
    a=dict(a,energy_comparison=core.energy_diagnostics(old_trial,a))
    b=dict(b,energy_comparison=core.energy_diagnostics(new_trial,b))
    values=batch.candidate_equivalence(a,b)
    for key in ('i','j','g','correction','volume','rate','internal_force','material_temperature'):
        require(np.array_equal(a[key],b[key]),'unchanged graph changed field '+key)
    for key in ('pair_conductance_W_K','pair_heat_into_i_J','linear_system_residual_J'):
        require(np.array_equal(a['thermal'][key],b['thermal'][key]),'unchanged graph changed heat '+key)
    require(core.clean(a['thermal']['ledger'])==core.clean(b['thermal']['ledger']),'thermal ledger differs')
    for key in ('contact','thermal','material'):equal_nested(a[key],b[key],key)
    return values


def main():
    BASE.mkdir(parents=True,exist_ok=False);started=time.perf_counter()
    files=[Path(__file__),ROOT/'scripts/visible_quasistatic.py',ROOT/'scripts/visible_support.py',
        ROOT/'scripts/visibility_equilibrium.py',ROOT/'scripts/secant_equilibrium.py',
        ROOT/'scripts/audit_solid_occluded_neighbors.py',ROOT/'scripts/batched_contact.py',
        ROOT/'scripts/visible_support_numerics.py',ROOT/'scripts/run_batched_process.py',ROOT/'scripts/run_discrete_post_start.py',
        ROOT/'docs/VISIBLE_SUPPORT_PROCESS_CONTRACT_20261004.md',PARENT,
        OLD/'manifest.json',OLD/'step_019/NOT_RESTARTABLE_candidate.npz',
        OLD/'step_020/NOT_RESTARTABLE_candidate.npz',OLD/'step_020/FULL_UNCOMMITTED_candidate.npz']
    frozen=dict(core.frozen_sources(),**{p.relative_to(ROOT).as_posix():core.sha(p) for p in files})
    report=dict(status='RUNNING_DIAGNOSTIC_ONLY',source_sha256=frozen,accepted_steps=0,
                accepted_time_advanced_s=0.,classification=visible.CLASSIFICATION)
    core.write(BASE/'manifest.json',report)
    require(core.sha(PARENT)==PARENT_SHA,'qualified parent drift')
    pa,state=core.restore(PARENT);p,provenance=core.parameters()
    require(state['step']==431 and abs(state['time_s']-.039375)<1e-15 and
            state['parameters']==p and state['parameter_provenance']==provenance and
            state['source_sha256']==core.frozen_sources(),'parent model differs')
    post.check_sources(state['implementation_source_sha256'])
    translated=visible.explicit_transition(state,PARENT_SHA)
    for key in state:
        if key!='classification':require(core.clean(translated[key])==core.clean(state[key]),'transition reset '+key)
    old=batch.make_trial(pa,state,p);new=visible.VisibleSupportTrial(pa,translated,p)
    dt=.0000625;guess=batch.guess_for(old,dt)
    a=old.evaluate(guess,dt);b=new.evaluate(guess,dt)
    require(b['visibility_graph']['diagnostics']['blocked_edge_count']==0,'registered parent forecast is already occluded')
    differences=same_map(a,b,old,new)
    ta=time.perf_counter();J0=old.tangent(a);old_tangent_wall=time.perf_counter()-ta
    ta=time.perf_counter();J1=new.tangent(b);new_tangent_wall=time.perf_counter()-ta
    require(np.array_equal(J0,J1),'unfiltered parent search tangent changed')
    report['parent_equivalence']=dict(status='PASS_SAME_INPUT_COMPLETE_MAP_AND_TANGENT',
        parent_sha256=PARENT_SHA,dt_s=dt,differences=differences,
        graph=b['visibility_graph']['diagnostics'],old_candidate_wall_s=a['wall_seconds'],
        numerical_domain=b['visibility_numerical_domain'],
        new_candidate_wall_s=b['wall_seconds'],old_tangent_wall_s=old_tangent_wall,new_tangent_wall_s=new_tangent_wall,
        tangent_maximum_difference_N_m=float(np.max(np.abs(J0-J1))),history_reset=False)
    # Deliberate diagnostic read: the public restart API correctly rejects this
    # reference. Keep that guard intact and never pass it to process acceptance.
    rejected=load(str(OLD/'step_019/NOT_RESTARTABLE_candidate.npz'))
    old_pa,old_state=rejected['arrays']['solid'],rejected['solver_data']
    require(old_state.get('restart_allowed') is False,'expected unaccepted diagnostic reference')
    failed=OLD/'step_020/FULL_UNCOMMITTED_candidate.npz'
    with np.load(failed) as evidence:u=evidence['displacement_m'].copy()
    legacy=batch.make_trial(old_pa,old_state,p)
    diagnostic=visible.explicit_transition(old_state,core.sha(OLD/'step_019/NOT_RESTARTABLE_candidate.npz'),diagnostic_only=True)
    changed=visible.VisibleSupportTrial(old_pa,diagnostic,p)
    before=legacy.evaluate(u,.000125);after=changed.evaluate(u,.000125)
    graph=after['visibility_graph'];blocked=graph['blocked_edge_pid_pairs']
    require(blocked.tolist()==[[219,251]],'saved failure no longer has exactly the independently audited obstruction')
    require(np.array_equal(after['i'],before['i'][graph['kept_raw_edge_indices']]) and
            np.array_equal(after['g'],before['g'][graph['kept_raw_edge_indices']]),'candidate operators did not use graph subset')
    require(len(after['thermal']['pair_conductance_W_K'])==len(after['i']),'heat did not use same graph')
    changed.contact=ChainRuleDiscreteContact(tool_radius=p['tool_radius'],bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'],particle_radius=p['dx']/2,normal_stiffness=p['normal_stiffness'],
        tangential_stiffness=p['normal_stiffness'],friction=p['friction'],base_z=-p['substrate_depth'])
    scalar=changed.evaluate(u,.000125);same_map(scalar,after,changed,changed)
    try:visible.accepted_copy(changed,after)
    except ValueError as error:require('diagnostic-only' in str(error),'unexpected diagnostic acceptance error')
    else:raise ValueError('diagnostic-only failed parent was allowed to advance')
    rows={key:float(np.max(np.abs(after[key]-before[key]))) for key in
        ('force','stress','ep','temperature','density','gradient','correction')}
    report['old_failure_map_diagnostic']=dict(status='NEW_SPATIAL_MAP_DIFFERENCE_NOT_SOLVED_OR_ACCEPTED',
        old_four_ms_status=post.read(OLD/'manifest.json')['status'],graph=graph['diagnostics'],differences=rows,
        numerical_domain=after['visibility_numerical_domain'],
        new_scaled_force_residual=after['scaled_residual'],old_scaled_force_residual=before['scaled_residual'],
        new_virtual_work_residual_J=after['virtual_work_residual_J'],new_energy=core.energy_diagnostics(changed,after),
        old_maximum_correction_condition=float(np.linalg.cond(before['correction']).max()),
        new_maximum_correction_condition=float(np.linalg.cond(after['correction']).max()),
        old_energy=core.energy_diagnostics(legacy,before),new_thermal=after['thermal']['ledger'],
        scalar_batch_visible_map='PASS',restart_allowed=False)
    np.savez_compressed(BASE/'failed_path_visible_map.npz',accepted=False,restart_allowed=False,
        particle_id=old_pa.particle_id,raw_i=before['i'],raw_j=before['j'],i=after['i'],j=after['j'],
        g=after['g'],correction=after['correction'],force_N=after['force'],stress_Pa=after['stress'],
        temperature_K=after['temperature'],plastic_strain=after['ep'],density_kg_m3=after['density'])
    post.check_sources(frozen)
    report.update(status='PASS_DIAGNOSTICS_NOT_PROCESS_ACCEPTANCE',full_map_evaluations=5,
                  tangent_evaluations=2,wall_seconds=time.perf_counter()-started)
    core.write(BASE/'manifest.json',report)
    print(json.dumps(dict(status=report['status'],parent_map=differences,blocked_ids=blocked.tolist(),
        new_force_scaled=after['scaled_residual'],wall_seconds=report['wall_seconds'])))


if __name__=='__main__':
    try:main()
    except Exception as error:
        path=BASE/'manifest.json'
        if path.exists():
            failed=post.read(path)
            if failed['status']=='RUNNING_DIAGNOSTIC_ONLY':
                failed.update(status='FAILED_DIAGNOSTIC_NO_ACCEPTANCE',error=repr(error));core.write(path,failed)
        raise
