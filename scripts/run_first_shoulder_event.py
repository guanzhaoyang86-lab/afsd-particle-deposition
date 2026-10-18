"""Bounded first-shoulder approach and same-parent transactional time screen."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import threading
import time
import traceback
import numpy as np
import run_batched_window as window
from contact_shell_audit import audit_contact_shell, CERTIFIED

batch, core, post = window.batch, window.core, window.post
BASE = batch.BASE/'first_shoulder_event'
PARENT_DIR = window.BASE/'refinement/eighth'
CONTRACT = core.ROOT/'docs/FIRST_SHOULDER_EVENT_CONTRACT_20261004.md'
SHELL_RECEIPT = core.ROOT/'logs/coupled_cost_20261004/shell_audit/receipt.json'
SHELL_READBACK = core.ROOT/'outputs/coupled_cost_20261004/shell_readback_20ms.json'
DRIVER_TEST = core.ROOT/'tests/review/test_first_shoulder_event.py'
DRIVER_RECEIPT = core.ROOT/'logs/coupled_cost_20261004/first_shoulder_driver_tests/receipt.json'
PARENT_SHA = '356ef74b474f52cef07455b64fb025b984c0939516b882865e57f24023e33832'
DT, MAX_APPROACH_STEPS = .000125, 16


def require(condition, message):
    if not condition:
        raise ValueError(message)


def checked_file(relative, expected_hash):
    path = (core.ROOT/relative).resolve()
    require(path.is_relative_to(core.ROOT.resolve()), 'evidence path outside project')
    require(core.sha(path) == expected_hash, 'evidence hash drift: '+relative)
    return path


def preflight(stage):
    """Read metadata and construct a trial only; no physical evaluation."""
    p, provenance, baseline = batch.check_inputs()
    meta = post.read(PARENT_DIR/'manifest.json')
    rows = post.read(PARENT_DIR/'accepted_history.json')
    require(meta['status'] == 'PASS_WINDOW_REFINEMENT_NOT_TRANSPORT' and
            meta['temporal_comparison']['pass'] and meta['accepted_steps'] == len(rows) == 112 and
            abs(meta['accepted_time_s']-.020) < 1e-15, 'qualified 20ms eighth parent absent')
    require(meta['provenance'] == provenance, 'parent manifest parameter provenance differs')
    post.check_sources(meta['source_sha256'])
    parent = PARENT_DIR/rows[-1]['checkpoint']
    require(core.sha(parent) == rows[-1]['sha256'] == PARENT_SHA, 'registered 20ms parent hash differs')
    pa, state = core.restore(parent)
    require(state['step'] == 136 and state['time_s'] == rows[-1]['time_s'] and
            abs(state['time_s']-.020) < 1e-15 and state['source_sha256'] == baseline and
            state['implementation_source_sha256'] == meta['source_sha256'] and
            state['parameters'] == p and state['parameter_provenance'] == provenance,
            '20ms checkpoint time/history/source/parameters differ')
    batch.make_trial(pa, state, p)
    receipt, readback = post.read(SHELL_RECEIPT), post.read(SHELL_READBACK)
    require(receipt['exit_code'] == 0 and receipt['physical_time_advanced_s'] == 0,
            'shell fixture gate absent')
    post.check_sources(receipt['source_sha256'])
    require(readback['status'] == CERTIFIED and readback['physical_evaluations'] == 0,
            'saved 20ms shell readback gate absent')
    post.check_sources(readback['source_sha256'])
    for name, digest in readback['input_sha256'].items():
        checked_file(name, digest)
    driver_tests = post.read(DRIVER_RECEIPT)
    require(driver_tests['exit_code'] == 0 and driver_tests['physical_evaluations'] == 0,
            'event controller fixture gate absent')
    post.check_sources(driver_tests['source_sha256'])
    extra_paths = [Path(__file__).resolve(), CONTRACT, core.ROOT/'scripts/contact_shell_audit.py',
                   core.ROOT/'tests/review/test_contact_shell_audit.py', SHELL_RECEIPT, SHELL_READBACK,
                   DRIVER_TEST, DRIVER_RECEIPT]
    frozen = dict(meta['source_sha256'], **{x.relative_to(core.ROOT).as_posix():core.sha(x) for x in extra_paths})
    per_root = max(row['wall_seconds'] for row in rows)
    origin = parent; origin_hash = PARENT_SHA
    event = None; references = {}
    if stage != 'approach':
        approach = post.read(BASE/'approach/manifest.json')
        require(approach['status'] == 'EVENT_CANDIDATE_NOT_ACCEPTED' and
                approach['source_sha256'] == frozen and approach['provenance'] == provenance,
                'matching uncommitted event approach absent')
        event_file = BASE/'approach/event.json'; event = post.read(event_file)
        require(core.sha(event_file) == approach['event_record_sha256'], 'event record drift')
        require(event['new_features'] == ['SHOULDER'] and event['event_particle_ids'],
                'only a first shoulder event may enter this local screen')
        require(event['dt_s'] == DT and event['source_sha256'] == frozen,
                'event interval/source differs')
        parent = checked_file(event['parent_checkpoint'], event['parent_sha256'])
        references[event_file] = core.sha(event_file)
        references[checked_file(event['candidate_checkpoint'], event['candidate_sha256'])] = event['candidate_sha256']
        pa, state = core.restore(parent)
        require(state['time_s'] == event['parent_time_s'] and state['step'] == event['parent_global_step'] and
                .020-1e-15 <= state['time_s'] < .022 and state['parameters'] == p and
                state['parameter_provenance'] == provenance and state['source_sha256'] == baseline,
                'event direct parent differs')
        require(state['interval_energy']['records'][:136] == core.restore(origin)[1]['interval_energy']['records'],
                'event parent lost the 20ms cumulative prefix')
        if parent != origin:
            require(state['implementation_source_sha256'] == frozen and state['branch_origin_sha256'] == origin_hash,
                    'approach accepted parent implementation/origin differs')
        batch.make_trial(pa, state, p)
        per_root = max(per_root, event['wall_seconds'])
        if stage == 'halves':
            macro = post.read(BASE/'macro/manifest.json')
            require(macro['status'] == 'MACRO_NUMERICALLY_QUALIFIED_NOT_ACCEPTED' and
                    macro['source_sha256'] == frozen and macro['parent_sha256'] == core.sha(parent) and
                    macro['event_record_sha256'] == core.sha(event_file), 'same-parent qualified macro absent')
            require(abs(macro['candidate_time_s']-(state['time_s']+DT)) < 1e-15,
                    'macro candidate endpoint differs')
            references[checked_file(macro['candidate_checkpoint'], macro['candidate_sha256'])] = macro['candidate_sha256']
            references[checked_file(macro['candidate_evidence'], macro['candidate_evidence_sha256'])] = macro['candidate_evidence_sha256']
            per_root = max(per_root, macro['candidate_wall_seconds'])
    count = MAX_APPROACH_STEPS if stage == 'approach' else 1 if stage == 'macro' else 2
    forecast = 30+1.5*count*per_root
    require(forecast <= 450, 'registered event forecast exceeds 450s')
    return dict(p=p, provenance=provenance, pa=pa, state=state, parent=parent,
                parent_hash=core.sha(parent), origin=origin, origin_hash=origin_hash,
                frozen=frozen, per_root=per_root, forecast=forecast, event=event, references=references)


def check_references(context):
    post.check_sources(context['frozen'])
    require(core.sha(context['parent']) == context['parent_hash'] and
            core.sha(context['origin']) == context['origin_hash'], 'root parent drift')
    for path, digest in context['references'].items():
        require(core.sha(path) == digest, 'bound event evidence drift')


def solve_candidate(trial, dt, directory, deadline, per_root):
    guess = batch.guess_for(trial, dt)
    attempt = directory/'broyden'; attempt.mkdir(); fallback = None
    try:
        result, metrics = batch.solver.solve(trial, dt, guess, attempt, deadline, 24)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        require(time.perf_counter()+per_root+15 <= deadline, 'insufficient budget for same-parent Newton retry')
        fallback = repr(error); attempt = directory/'newton_retry'; attempt.mkdir()
        result, metrics = core.solve(trial, dt, guess, attempt, deadline, 24)
        metrics['controller'] = 'ORIGINAL_MODIFIED_NEWTON_RETRY_SAME_PARENT'
    metrics['fallback_reason'] = fallback
    # Only detached numerical candidate states are made here, never published.
    pa, state, fresh = core.accepted_copy(trial, result)
    return pa, state, fresh, metrics


def audit_candidate(trial, result, directory):
    geometry = post.geometry_audit(trial, result)
    core.write(directory/'geometry_audit.json', geometry)
    shell = audit_contact_shell(trial.contact, trial.x0, result['position'],
        core.pose(trial.p, trial.state['time_s']),
        core.pose(trial.p, trial.state['time_s']+result['dt_s']))
    for entry in shell['hidden_contacts']:
        entry['particle_id'] = int(trial.pa.particle_id[entry['particle_index']])
    core.write(directory/'continuous_shell_audit.json', shell)
    require(shell['status'] == CERTIFIED, 'continuous endpoint-hidden shell audit did not certify')
    require(not geometry['transient_contacts_missed_by_endpoints'], 'sampled hidden shell contact')
    events = []
    for fraction in (0., .25, .5, .75, 1.):
        endpoint = trial.contact.geometry(trial.x0+fraction*result['displacement'],
            core.pose(trial.p, trial.state['time_s']+fraction*result['dt_s']))
        for i, j in np.argwhere(endpoint['active']):
            if j in (0, 1, 6, 7):
                events.append(dict(fraction=fraction, particle_id=int(trial.pa.particle_id[i]),
                                   feature=trial.contact.features[j], gap_m=float(endpoint['gaps_m'][i,j])))
    core.write(directory/'sampled_nonbore_events.json', events)
    return geometry, shell, events


def save_candidate(directory, trial, pa, state, result, frozen):
    core.journal(directory, trial, result, 'FULL_UNCOMMITTED_candidate')
    staged = dict(state, restart_allowed=False, event_candidate_accepted=False,
                  implementation_source_sha256=frozen)
    target = directory/'NOT_RESTARTABLE_candidate.npz'
    core.checkpoint(target, pa, staged); core.verify_snapshot(target, pa, staged)
    saved = core.load(str(target))['solver_data']
    for key in ('interval_energy', 'bore_observation', 'implementation_source_sha256'):
        require(core.clean(saved[key]) == core.clean(staged[key]), 'uncommitted metadata readback differs')
    evidence = directory/'candidate_evidence.npz'
    np.savez_compressed(evidence, particle_id=pa.particle_id, i=result['i'], j=result['j'])
    core.write(directory/'energy.json', dict(mechanical=state['energy'], interval=state['interval_energy'],
                                           thermal=result['thermal']['ledger']))
    return target, evidence


def scalar_equivalence(trial, result, directory):
    old_trial = core.DiscreteContactTrial(trial.pa, trial.state, trial.p)
    old = old_trial.evaluate(result['displacement'], result['dt_s'])
    old['energy_comparison'] = core.energy_diagnostics(old_trial, old)
    candidate = dict(result, energy_comparison=core.energy_diagnostics(trial, result))
    difference = batch.candidate_equivalence(old, candidate)
    # The helper compares all contact arrays. These remaining model/sweep fields
    # also bind the finite-map classification and continuous-solid certificate.
    for key, value in old['contact'].items():
        if not isinstance(value, np.ndarray):
            require(core.clean(value) == core.clean(result['contact'][key]), 'contact metadata differs: '+key)
    core.write(directory/'scalar_batch_equivalence.json', dict(status='PASS_SAME_CANDIDATE_COMPLETE_MAP', differences=difference))
    return difference


def one_hop_ids(seeds, graphs):
    selected = set(int(v) for v in seeds)
    for ids, i, j in graphs:
        left, right = ids[np.asarray(i, int)], ids[np.asarray(j, int)]
        adjacent = np.isin(left, list(seeds)) | np.isin(right, list(seeds))
        selected.update(int(v) for v in np.r_[left[adjacent], right[adjacent]])
    return sorted(selected)


def local_compare(fine, coarse, selected):
    require(np.array_equal(fine.particle_id, coarse.particle_id), 'local comparison ID order differs')
    index = np.flatnonzero(np.isin(fine.particle_id, selected))
    require(len(index) == len(selected) and len(index) > 0, 'event neighborhood IDs incomplete')
    dx = np.column_stack((fine.x-coarse.x, fine.y-coarse.y, fine.z-coarse.z))[index]
    sf, sc = fine.S.reshape(-1,9)[index], coarse.S.reshape(-1,9)[index]
    sr = np.linalg.norm(sf-sc,axis=1)/np.maximum(np.linalg.norm(sf,axis=1),1e6)
    ep = np.abs(fine.ep[index]-coarse.ep[index])/np.maximum(np.abs(fine.ep[index]),1e-5)
    position = np.linalg.norm(dx,axis=1); temperature = np.abs(fine.temperature[index]-coarse.temperature[index])
    rows = [dict(particle_id=int(fine.particle_id[k]), position_error_m=float(position[n]),
        temperature_error_K=float(temperature[n]), stress_relative=float(sr[n]), plastic_relative=float(ep[n]))
        for n,k in enumerate(index)]
    return dict(particle_ids=selected, rows=rows, position_limit_m=5e-6, temperature_limit_K=.5,
                stress_relative_limit=.1, stress_denominator_floor_Pa=1e6,
                plastic_relative_limit=.2, plastic_denominator_floor=1e-5,
                passed=bool(np.all(position<=5e-6) and np.all(temperature<=.5) and np.all(sr<=.1) and np.all(ep<=.2)))


def worker(stage):
    directory = BASE/stage; directory.mkdir(parents=True, exist_ok=False)
    begin = time.perf_counter(); history = []; pending = []
    report = dict(stage=stage, status='RUNNING', accepted_steps=0, accepted_time_s=None,
                  classification=core.CLASSIFICATION, physical_accuracy='ENGINEERING_EVENT_SCREEN_NOT_EXPERIMENTAL_VALIDATION')
    def stop():
        try:
            core.write(directory/'HARD_TIMEOUT.json', dict(accepted_history=history,
                published_checkpoints=[str(x.relative_to(directory)) for x in directory.glob('step_*/accepted.npz')],
                uncommitted_candidates=[str(x.relative_to(directory)) for x in directory.glob('step_*/NOT_RESTARTABLE_candidate.npz')]))
        finally:
            os._exit(124)
    timer = threading.Timer(540., stop); timer.daemon=True; timer.start()
    try:
        c = preflight(stage); pa,state,p = c['pa'],c['state'],c['p']; per_root=c['per_root']
        report.update(source_sha256=c['frozen'], provenance=c['provenance'],
            parent_checkpoint=c['parent'].relative_to(core.ROOT).as_posix(), parent_sha256=c['parent_hash'],
            branch_origin_sha256=c['origin_hash'], starting_time_s=state['time_s'], accepted_time_s=state['time_s'],
            forecast_wall_s=c['forecast'], thread_environment={k:os.environ.get(k) for k in
                ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')})
        if stage != 'approach':
            report['event_record_sha256'] = core.sha(BASE/'approach/event.json')
        core.write(directory/'manifest.json',report); core.write(directory/'accepted_history.json',history)
        initial,initial_state=core.restore(batch.OLD/'quarters/independent_initial.npz')
        core.checkpoint(directory/'independent_initial.npz',initial,initial_state)
        core.checkpoint(directory/'continuation_initial.npz',pa,state)
        direct_file,direct_hash=c['parent'],c['parent_hash']
        count=MAX_APPROACH_STEPS if stage=='approach' else 1 if stage=='macro' else 2
        dt=DT/2 if stage=='halves' else DT
        for number in range(1,count+1):
            require(time.perf_counter()-begin+1.5*per_root+15<=525, 'remaining event budget exhausted')
            stepdir=directory/f'step_{number:03d}';stepdir.mkdir();started=time.perf_counter()
            trial=batch.make_trial(pa,state,p)
            newpa,newstate,fresh,metrics=solve_candidate(trial,dt,stepdir,begin+525,per_root)
            candidate,evidence=save_candidate(stepdir,trial,newpa,newstate,fresh,c['frozen'])
            geometry,shell,events=audit_candidate(trial,fresh,stepdir)
            new_features=sorted(set(e['feature'] for e in events))
            require(new_features==sorted(geometry['new_non_bore_features']), 'sampled event feature accounting differs')
            if stage=='approach' and new_features:
                check_references(c);require(core.sha(direct_file)==direct_hash,'event parent drift')
                event=dict(source_sha256=c['frozen'],parent_checkpoint=direct_file.relative_to(core.ROOT).as_posix(),
                    parent_sha256=direct_hash,parent_time_s=state['time_s'],parent_global_step=state['step'],
                    dt_s=dt,candidate_time_s=newstate['time_s'],new_features=new_features,
                    event_particle_ids=sorted(set(e['particle_id'] for e in events)),events=events,
                    candidate_checkpoint=candidate.relative_to(core.ROOT).as_posix(),candidate_sha256=core.sha(candidate),
                    candidate_evidence=evidence.relative_to(core.ROOT).as_posix(),candidate_evidence_sha256=core.sha(evidence),
                    wall_seconds=time.perf_counter()-started,accepted=False,restart_allowed=False)
                core.write(directory/'event.json',event)
                report.update(status='EVENT_CANDIDATE_NOT_ACCEPTED',event=event,event_record_sha256=core.sha(directory/'event.json'))
                break
            if stage!='approach':
                require(set(new_features)<= {'SHOULDER'}, 'non-shoulder feature outside registered event scope')
                scalar_equivalence(trial,fresh,stepdir)
                pending.append(dict(directory=stepdir,trial=trial,pa=newpa,state=newstate,result=fresh,
                    metrics=metrics,events=events,candidate=candidate,evidence=evidence,
                    wall_seconds=time.perf_counter()-started))
                pa,state=newpa,newstate
                per_root=max(per_root,pending[-1]['wall_seconds'])
                continue
            check_references(c);require(core.sha(direct_file)==direct_hash,'direct parent drift')
            batch.publish(stepdir,newpa,newstate,trial,fresh,c['frozen'],c['provenance'],direct_hash,c['origin_hash'])
            pa,state=newpa,newstate
            row=dict(step=number,global_step=state['step'],time_s=state['time_s'],dt_s=dt,parent_sha256=direct_hash,
                checkpoint=f'step_{number:03d}/accepted.npz',sha256=core.sha(stepdir/'accepted.npz'),
                wall_seconds=time.perf_counter()-started,**metrics)
            history.append(row);core.write(directory/'accepted_history.json',history)
            direct_file=stepdir/'accepted.npz';direct_hash=row['sha256'];per_root=max(per_root,row['wall_seconds'])
            report.update(accepted_steps=len(history),accepted_time_s=state['time_s'],interval_energy=state['interval_energy'],
                          crossing_mass_kg=state['crossing_mass_kg'])
            core.write(directory/'manifest.json',report);print(json.dumps(row),flush=True)
        else:
            if stage=='approach':
                report['status']='BOUNDED_APPROACH_COMPLETE_NO_EVENT'
        if stage=='macro':
            require(any(e['feature']=='SHOULDER' for e in pending[0]['events']), 'registered shoulder event not reproduced')
            item=pending[0]
            report.update(status='MACRO_NUMERICALLY_QUALIFIED_NOT_ACCEPTED',candidate_time_s=state['time_s'],
                candidate_checkpoint=item['candidate'].relative_to(core.ROOT).as_posix(),candidate_sha256=core.sha(item['candidate']),
                candidate_evidence=item['evidence'].relative_to(core.ROOT).as_posix(),candidate_evidence_sha256=core.sha(item['evidence']),
                candidate_wall_seconds=item['wall_seconds'],event_particle_ids=sorted(set(e['particle_id'] for e in item['events'])))
        if stage=='halves':
            macro=post.read(BASE/'macro/manifest.json')
            raw=core.load(str(checked_file(macro['candidate_checkpoint'],macro['candidate_sha256'])))
            coarse,coarse_state=raw['arrays']['solid'],raw['solver_data']
            require(coarse_state['restart_allowed'] is False and abs(coarse_state['time_s']-state['time_s'])<1e-15,
                    'macro is not the same-time uncommitted reference')
            all_events=[e for item in pending for e in item['events']]
            require(any(e['feature']=='SHOULDER' for e in all_events),'halves did not reproduce shoulder contact')
            seeds=set(c['event']['event_particle_ids'])|set(macro['event_particle_ids'])|{e['particle_id'] for e in all_events}
            with np.load(checked_file(macro['candidate_evidence'],macro['candidate_evidence_sha256'])) as data:
                graphs=[(data['particle_id'].copy(),data['i'].copy(),data['j'].copy())]
            graphs += [(item['pa'].particle_id,item['result']['i'],item['result']['j']) for item in pending]
            selected=one_hop_ids(seeds,graphs)
            global_screen=post.compare(pa,coarse,p);local_screen=local_compare(pa,coarse,selected)
            core.write(directory/'temporal_comparison.json',dict(global_screen=global_screen,local_screen=local_screen,
                event_particle_ids=sorted(seeds),all_half_events=all_events))
            report.update(global_time_screen=global_screen,local_time_screen=local_screen)
            require(global_screen['pass'] and local_screen['passed'],'first shoulder global/local temporal screen failed')
            # No restartable half-step file exists before BOTH candidates and all
            # event screens pass. A later publication failure retains a qualified
            # prefix; each checkpoint publication itself is atomic.
            check_references(c)
            for number,item in enumerate(pending,1):
                require(core.sha(direct_file)==direct_hash,'direct parent drift before qualified event commit')
                batch.publish(item['directory'],item['pa'],item['state'],item['trial'],item['result'],
                    c['frozen'],c['provenance'],direct_hash,c['origin_hash'])
                path=item['directory']/'accepted.npz'
                row=dict(step=number,global_step=item['state']['step'],time_s=item['state']['time_s'],dt_s=dt,
                    parent_sha256=direct_hash,checkpoint=f'step_{number:03d}/accepted.npz',sha256=core.sha(path),
                    wall_seconds=item['wall_seconds'],**item['metrics'])
                history.append(row);core.write(directory/'accepted_history.json',history)
                direct_file=path;direct_hash=row['sha256']
                report.update(accepted_steps=len(history),accepted_time_s=item['state']['time_s'])
                core.write(directory/'manifest.json',report)
            report.update(status='PASS_FIRST_SHOULDER_LOCAL_TIME_SCREEN_NOT_TRANSPORT',
                          interval_energy=state['interval_energy'],crossing_mass_kg=state['crossing_mass_kg'])
        check_references(c);code=0
    except Exception as error:
        traceback.print_exc();report.update(status='FAILED_WITH_ACCEPTED_PREFIX' if history else 'FAILED_NO_NEW_ACCEPTED_STATE',error=repr(error));code=2
    finally:
        timer.cancel()
    report.update(wall_seconds=time.perf_counter()-begin,peak_process_working_set_bytes=core.memory_bytes())
    core.write(directory/'manifest.json',report);return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['approach','macro','halves'],required=True)
    parser.add_argument('--worker',action='store_true');parser.add_argument('--preflight-only',action='store_true');args=parser.parse_args()
    if args.preflight_only:
        context=preflight(args.stage)
        print(json.dumps(dict(stage=args.stage,forecast_wall_s=context['forecast'],parent_sha256=context['parent_hash'],
            parent_time_s=context['state']['time_s'],parent_global_step=context['state']['step'],physical_evaluations=0)))
        raise SystemExit(0)
    if args.worker:raise SystemExit(worker(args.stage))
    logs=core.ROOT/'logs/coupled_cost_20261004/batched/first_shoulder_event';logs.mkdir(parents=True,exist_ok=True);log=logs/(args.stage+'.log')
    require(not log.exists() and not (BASE/args.stage).exists(),'refuse evidence overwrite')
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--stage',args.stage,'--worker']
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1');started=time.perf_counter()
    with log.open('w',encoding='utf-8') as stream:
        try:code=subprocess.run(command,cwd=core.ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600).returncode
        except subprocess.TimeoutExpired:code=124
    core.write(logs/(args.stage+'_supervision.json'),dict(command=command,exit_code=code,wall_seconds=time.perf_counter()-started,limit_s=600))
    print(json.dumps(dict(stage=args.stage,exit_code=code,log=str(log))));raise SystemExit(code)
