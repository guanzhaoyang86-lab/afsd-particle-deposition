"""Unpublished ramp candidate with one common step motion for every audit."""
import time
import numpy as np
import run_visible_support_process as old
import ramped_hot_trial as model
from visible_support_path_tight import audit_visible_support_path

core, post, event, batch, solver = old.core, old.post, old.event, old.batch, old.solver
require = old.require


def geometry_audit(trial, r):
    """Same sampled gate as the frozen audit, using the declared step motion."""
    t, dt, contact = trial.state['time_s'], r['dt_s'], trial.contact
    first = contact.geometry(trial.x0, trial.motion.pose(t))
    last = contact.geometry(r['position'], trial.motion.pose(t+dt))
    transient, counts, events = [], [], []
    active = first['active'] | last['active']
    for fraction in (0., .25, .5, .75, 1.):
        middle = contact.geometry(trial.x0+fraction*r['displacement'], trial.motion.pose(t+fraction*dt))
        for i, j in np.argwhere(middle['active']):
            events.append(dict(fraction=fraction, particle_id=int(trial.pa.particle_id[i]), feature=contact.features[j]))
        if fraction in (.25, .5, .75):
            missed = middle['active'] & ~first['active'] & ~last['active']
            transient.extend(dict(fraction=fraction, particle_id=int(trial.pa.particle_id[i]), feature=contact.features[j]) for i,j in np.argwhere(missed))
            counts.append(middle['active'].sum(axis=0).tolist())
            active |= middle['active']
    neighbors = np.bincount(np.r_[r['i'], r['j']], minlength=len(trial.x0))
    report = dict(endpoint_counts=last['active'].sum(axis=0).tolist(), interior_counts=counts,
        features=list(contact.features), transient_contacts_missed_by_endpoints=transient,
        maximum_penetration_m=float(np.max(last['overlap_m'])), min_neighbor_count=int(neighbors.min()),
        maximum_correction_condition=float(np.linalg.cond(r['correction']).max()),
        new_non_bore_features=[contact.features[j] for j in (0,1,6,7) if active[:,j].any()],
        limitation='Samples supplement continuous certificates for the declared step motion only.')
    return report, events


def candidate(pa, state, dt, directory, deadline, per_root, frozen):
    started = time.perf_counter()
    trial = model.make_trial(pa, state, dt)
    guess = batch.guess_for(trial, dt)
    attempt = directory/'broyden'; attempt.mkdir(); fallback = None
    try:
        r, metrics = solver.solve(trial, dt, guess, attempt, deadline, 24)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
        require(time.perf_counter()+per_root+15 <= deadline, 'insufficient budget for same-parent Newton retry')
        fallback = repr(error); attempt = directory/'newton_retry'; attempt.mkdir()
        r, metrics = core.solve(trial, dt, guess, attempt, deadline, 24)
        metrics['controller'] = 'ORIGINAL_MODIFIED_NEWTON_SAME_PARENT_RAMP_FULL_MAP'
    metrics['fallback_reason'] = fallback
    core.journal(directory, trial, r, 'SOLVED_UNCOMMITTED_candidate')
    old.save_graph(directory, trial, r)
    nextpa, nextstate, fresh = model.advance_candidate(trial, r)
    require(np.array_equal(r['i'],fresh['i']) and np.array_equal(r['j'],fresh['j']) and
        r['graph_fingerprint']==fresh['graph_fingerprint'], 'fresh graph changed at identical displacement')
    require(nextstate['interval_energy']['records'][:-1]==state['interval_energy']['records'], 'full prefix changed')
    target, evidence = event.save_candidate(directory, trial, nextpa, nextstate, fresh, frozen)
    old.save_graph(directory, trial, fresh)
    p = trial.p
    scalar = model.make_trial(pa, state, dt)
    scalar.contact = old.ChainRuleDiscreteContact(tool_radius=p['tool_radius'], bore_width=p['bore_width'],
        tool_height=p['tool_top']-p['layer'], particle_radius=p['dx']/2, normal_stiffness=p['normal_stiffness'],
        tangential_stiffness=p['normal_stiffness'], friction=p['friction'], base_z=-p['substrate_depth'])
    sr = scalar.evaluate(fresh['displacement'], dt)
    differences = batch.candidate_equivalence(dict(sr,energy_comparison=core.energy_diagnostics(scalar,sr)),
        dict(fresh,energy_comparison=core.energy_diagnostics(trial,fresh)))
    require(sr['graph_fingerprint']==fresh['graph_fingerprint'] and np.array_equal(sr['g'],fresh['g']) and
        np.array_equal(sr['thermal']['pair_conductance_W_K'],fresh['thermal']['pair_conductance_W_K']), 'scalar graph/heat differs')
    core.write(directory/'scalar_batch_equivalence.json', dict(status='PASS_SAME_RAMP_FILTERED_COMPLETE_MAP', differences=differences))
    geometry, events = geometry_audit(trial, fresh)
    core.write(directory/'geometry_audit.json', geometry)
    require(not geometry['transient_contacts_missed_by_endpoints'] and set(geometry['new_non_bore_features']) <= {'SHOULDER','TOP'},
        'unqualified or endpoint-hidden contact feature')
    t = state['time_s']; pose0,pose1 = trial.motion.pose(t),trial.motion.pose(t+dt)
    shell = old.audit_contact_shell_finite_enclosure(trial.contact,trial.x0,fresh['position'],pose0,pose1)
    core.write(directory/'continuous_shell_audit.json',shell)
    require(shell['status']==event.CERTIFIED, 'continuous shell unresolved')
    start = time.perf_counter()
    path = audit_visible_support_path(trial.x0,fresh['position'],fresh['i'],fresh['j'],pose0,pose1,p,pa.particle_id)
    path.update(wall_seconds=time.perf_counter()-start,candidate_sha256=core.sha(target),graph_sha256=core.sha(directory/'visible_graph.npz'))
    core.write(directory/'visible_path_audit.json',path)
    require(path['status']==old.CERTIFIED, 'visible path unresolved')
    core.write(directory/'motion.json',trial.motion.to_dict())
    core.write(directory/'tool_pose.json',pose1)
    core.write(directory/'sampled_events.json',events)
    core.write(directory/'transport.json',nextstate['bore_observation'])
    core.vtk(directory/'fields.vtk',nextpa,p,trial.bulk,nextstate['time_s'],model.CLASSIFICATION)
    loaded = np.linalg.norm(np.asarray(fresh['contact']['force_by_feature_N'])[:,:7],axis=2).max(axis=1) > fresh['force_limit_N']
    nonbase = np.asarray(fresh['contact']['active'])[:,:7].any(axis=1)
    response = dict(nonbase_loaded_ids=pa.particle_id[loaded & nonbase].astype(int).tolist(),
        new_plastic_ids=pa.particle_id[fresh['ep']>pa.ep].astype(int).tolist(),
        max_plastic_increment=float(np.max(fresh['ep']-pa.ep)), max_plastic_strain=float(fresh['ep'].max()),
        temperature_min_K=float(fresh['temperature'].min()),temperature_max_K=float(fresh['temperature'].max()),
        crossing_ids=nextstate['bore_observation']['crossing_ids'])
    core.write(directory/'response.json',response)
    return dict(pa=nextpa,state=nextstate,trial=trial,result=fresh,metrics=metrics,directory=directory,
        checkpoint=target,evidence=evidence,events=events,response=response,wall_seconds=time.perf_counter()-started)
