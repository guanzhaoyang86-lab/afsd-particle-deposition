"""Expanded continuation candidate; unchanged physical and temporal gates.

Static adaptation of visible_tight_path_candidate with explicit independent origin.
"""
import time
import numpy as np
import run_visible_support_process as original

core,post,event,solver,batch=(original.core,original.post,original.event,original.solver,original.batch)
import expanded_hot_domain as model
from expanded_hot_continuation_acceptance import accepted_copy


class RecordedTrial(model.ExpandedVisibleTrial):
    def evaluate(self,displacement,dt):
        result=super().evaluate(displacement,dt)
        self.last_complete_result=result
        return result

ChainRuleDiscreteContact=original.ChainRuleDiscreteContact
audit_contact_shell_finite_enclosure=original.audit_contact_shell_finite_enclosure
from visible_support_path_tight import audit_visible_support_path
require=original.require
ALLOWED_FEATURES=frozenset(('SHOULDER','TOP'))
QUALIFICATION='VISIBLE_SUPPORT_SHOULDER_AND_TOP_EVENT_SCREEN_V1'


def candidate(pa,state,c,dt,directory,deadline,per_root):
    started=time.perf_counter();trial=RecordedTrial(pa,state,c['p']);guess=batch.guess_for(trial,dt)
    attempt=directory/'broyden';attempt.mkdir();fallback=None
    try:r,metrics=solver.solve(trial,dt,guess,attempt,deadline,24)
    except (ValueError,ArithmeticError,np.linalg.LinAlgError) as error:
        require(time.perf_counter()+per_root+15<=deadline,'insufficient budget for one same-parent Newton retry')
        fallback=repr(error);attempt=directory/'newton_retry';attempt.mkdir()
        try:r,metrics=core.solve(trial,dt,guess,attempt,deadline,24)
        except Exception:
            if hasattr(trial,'last_complete_result'):
                evidence=attempt/'LAST_FULL_EVALUATION_NOT_NECESSARILY_BEST_ITERATE';evidence.mkdir()
                solver.journal(evidence,trial,trial.last_complete_result,'last_FULL_UNCOMMITTED_candidate')
                original.save_graph(evidence,trial,trial.last_complete_result)
            raise
        metrics['controller']='ORIGINAL_MODIFIED_NEWTON_SAME_PARENT_VISIBLE_FULL_MAP'
    metrics['fallback_reason']=fallback
    core.journal(directory,trial,r,'SOLVED_UNCOMMITTED_candidate');original.save_graph(directory,trial,r)
    nextpa,nextstate,fresh=accepted_copy(trial,r,c['qualification'])
    require(np.array_equal(r['i'],fresh['i']) and np.array_equal(r['j'],fresh['j']) and
        r['graph_fingerprint']==fresh['graph_fingerprint'],'fresh graph changed at identical displacement')
    require(nextstate['interval_energy']['records'][:len(c['prefix'])]==c['prefix'],'complete old prefix was reset')
    target,evidence=event.save_candidate(directory,trial,nextpa,nextstate,fresh,c['frozen'])
    original.save_graph(directory,trial,fresh)
    scalar=model.ExpandedVisibleTrial(pa,state,c['p']);p=c['p']
    scalar.contact=ChainRuleDiscreteContact(tool_radius=p['tool_radius'],bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'],particle_radius=p['dx']/2,normal_stiffness=p['normal_stiffness'],
        tangential_stiffness=p['normal_stiffness'],friction=p['friction'],base_z=-p['substrate_depth'])
    sr=scalar.evaluate(fresh['displacement'],dt)
    differences=batch.candidate_equivalence(dict(sr,energy_comparison=core.energy_diagnostics(scalar,sr)),
        dict(fresh,energy_comparison=core.energy_diagnostics(trial,fresh)))
    require(sr['graph_fingerprint']==fresh['graph_fingerprint'] and np.array_equal(sr['g'],fresh['g']) and
        np.array_equal(sr['thermal']['pair_conductance_W_K'],fresh['thermal']['pair_conductance_W_K']),
        'scalar map did not use same graph/heat')
    core.write(directory/'scalar_batch_equivalence.json',dict(status='PASS_SAME_FILTERED_COMPLETE_MAP',differences=differences))
    geometry=post.geometry_audit(trial,fresh);core.write(directory/'geometry_audit.json',geometry)
    require(not geometry['transient_contacts_missed_by_endpoints'] and
        set(geometry['new_non_bore_features'])<=ALLOWED_FEATURES,'contact outside registered SHOULDER+TOP domain')
    t=state['time_s'];pose0,pose1=core.pose(p,t),core.pose(p,t+dt)
    shell=audit_contact_shell_finite_enclosure(trial.contact,trial.x0,fresh['position'],pose0,pose1)
    core.write(directory/'continuous_shell_audit.json',shell)
    require(shell['status']==event.CERTIFIED,'continuous contact shell unresolved')
    before=time.perf_counter()
    path=audit_visible_support_path(trial.x0,fresh['position'],fresh['i'],fresh['j'],pose0,pose1,p,pa.particle_id)
    path.update(wall_seconds=time.perf_counter()-before,candidate_sha256=core.sha(target),graph_sha256=core.sha(directory/'visible_graph.npz'))
    core.write(directory/'visible_path_audit.json',path)
    require(path['status']==original.CERTIFIED,'retained support graph not certified over complete step')
    events=[]
    for fraction in (0.,.25,.5,.75,1.):
        endpoint=trial.contact.geometry(trial.x0+fraction*fresh['displacement'],core.pose(p,t+fraction*dt))
        for i,j in np.argwhere(endpoint['active']):
            if j in (0,1,6,7):events.append(dict(fraction=fraction,particle_id=int(pa.particle_id[i]),feature=trial.contact.features[j]))
    require(set(e['feature'] for e in events)<=ALLOWED_FEATURES,'sampled unqualified contact feature')
    core.write(directory/'sampled_nonbore_events.json',events)
    core.write(directory/'contact_event_qualification.json',dict(qualification=QUALIFICATION,
        allowed_nonbore_features=sorted(ALLOWED_FEATURES),actual_nonbore_features=sorted(set(e['feature'] for e in events)),
        status='CANDIDATE_GATES_PASS_REQUIRES_SAME_PARENT_PAIR_SCREEN',physics_changed=False,
        candidate_sha256=core.sha(target),time_s=nextstate['time_s']))
    return dict(pa=nextpa,state=nextstate,trial=trial,result=fresh,metrics=metrics,directory=directory,
        checkpoint=target,evidence=evidence,events=events,wall_seconds=time.perf_counter()-started)
