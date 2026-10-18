"""Append audited actual frames to the fixed full deposition-observer history.

Pure saved-array and checkpoint reading: no PySPH, NNPS, contact, material or
thermal solver imports. Prior failures are not changed or promoted.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
import traceback
import numpy as np
from deposition_stage_observer import observe_deposition_stages

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT/'outputs/coupled_cost_20261004/batched/readback_pre_occlusion_39_375ms'
BASE = ROOT/'outputs/visible_support_20261004/readback_v1'
OBSERVER = ROOT/'outputs/deposition_observer_20261004/readback_40_375ms_v1/summary.json'
OLD_SHA = '4632c3cfc7430ad00942aeea58edb227607333c5b078610ca66af1838595b600'
BASE_SHA = '54543bdc917afc3b2be431549818bf8b523815a2e36c8d48d12d6b04a95a50d8'
OBSERVER_SHA = '37943c25fe49cb413ffe5f63965222bae9c4af063d0a5697dce9228981d4aa8e'
ACCEPTED = frozenset(('PASS_SAVED_EVIDENCE_READBACK_NOT_FRESH_PHYSICS',
    'PASS_VISIBLE_CONTINUATION_READBACK_NOT_PHYSICAL_VALIDATION',
    'PASS_VISIBLE_TOP_EVENT_READBACK_NOT_PHYSICAL_VALIDATION',
    'PASS_VISIBLE_LOCAL_REFINEMENT_READBACK_NOT_PHYSICAL_VALIDATION',
    'PASS_TIGHT_PATH_PROCESS_READBACK_NOT_PHYSICAL_VALIDATION'))
FRAME_KEYS = ('time_s','position_m','plastic_strain','tool_origin_m','tool_rotation')


def require(value, message):
    if not value: raise ValueError(message)


def sha(path):
    with path.open('rb') as stream:
        result = hashlib.sha256()
        for block in iter(lambda: stream.read(1024*1024), b''): result.update(block)
    return result.hexdigest()


class Reader:
    def __init__(self): self.bound = {}

    def bind(self, path, expected=None):
        path = Path(path).resolve(); require(path.is_relative_to(ROOT), 'input outside repository')
        name = path.relative_to(ROOT).as_posix(); actual = sha(path)
        require(expected is None or expected == actual, 'changed evidence '+name)
        require(name not in self.bound or self.bound[name] == actual, 'input changed during read '+name)
        self.bound[name] = actual
        return path

    def read(self, path, expected=None):
        return json.loads(self.bind(path, expected).read_text(encoding='utf-8'))

    def arrays(self, path, expected=None):
        with np.load(self.bind(path, expected), allow_pickle=False) as data:
            return {key:data[key] for key in data.files}

    def checkpoint(self, path, expected, legacy_audited=False):
        # Only hash-bound locally produced accepted checkpoints are unpickled.
        with np.load(self.bind(path, expected), allow_pickle=True) as data:
            state = data['solver_data'].item(); arrays = data['particles'].item()['solid']['arrays']
        require(state.get('restart_allowed') is True or
                (legacy_audited and 'restart_allowed' not in state), 'unaccepted candidate in observation chain')
        return arrays, state

    def unchanged(self):
        for name, expected in self.bound.items(): require(sha(ROOT/name) == expected, 'input drift '+name)


def check_frame(arrays, state, data, index, ids, mass, components):
    require(np.array_equal(arrays['particle_id'], ids), 'checkpoint IDs differ')
    require(np.array_equal(arrays['m'], mass) and np.array_equal(arrays['component_id'], components), 'mass/component drift')
    require(state['time_s'] == data['time_s'][index], 'checkpoint time differs')
    require(np.array_equal(np.column_stack([arrays[k] for k in ('x','y','z')]), data['position_m'][index]), 'trajectory positions differ')
    require(np.array_equal(arrays['ep'], data['plastic_strain'][index]), 'trajectory plastic strain differs')


def append_readback(reader, summary_path, data, static, frames, previous_summary, previous_sha, state, last_path, last_sha, base=False):
    summary_path = Path(summary_path)
    if summary_path.is_dir(): summary_path = summary_path/'summary.json'
    summary_path = summary_path.resolve()
    summary = reader.read(summary_path, BASE_SHA if base else None)
    require(summary['status'] in ACCEPTED and summary['exit_code'] == 0, 'readback did not pass')
    require(summary['physical_solver_evaluations'] == 0 and summary['nnps_evaluations'] == 0, 'readback is not saved-array evidence')
    require(summary['inherited_interval_never_reset'], 'readback reset inherited energy interval')
    if not base:
        require((ROOT/summary['parent_readback_summary']).resolve() == previous_summary and
                summary['parent_readback_sha256'] == previous_sha, 'wrong readback ancestry')
        require((ROOT/summary['parent_checkpoint']).resolve() == last_path and summary['parent_sha256'] == last_sha,
                'window parent is not the previous actual accepted endpoint')
    trajectory = summary['trajectory']; new = reader.arrays(ROOT/trajectory['path'], trajectory['sha256'])
    history = summary['accepted_history']; require(len(new['time_s']) == len(history)+1, 'nonactual or missing frames')
    require(summary['initial_accepted_global_steps'] == state['step'], 'parent global step differs')
    for key in FRAME_KEYS: require(np.array_equal(data[key][-1], new[key][0]), 'shared true parent differs '+key)
    for key in ('particle_id','mass_kg','component_id'): require(np.array_equal(static[key], new[key]), 'static trajectory identity differs '+key)
    active = []; prior_prefix = list(state['interval_energy']['records']); prior_events = state['bore_observation']['events']
    frozen_eligible = state['bore_observation']['eligible_ids']; parameters = state['parameters']
    evidence = summary['evidence_sha256']
    for number, row in enumerate(history, 1):
        path = (ROOT/row['checkpoint']).resolve(); arrays, saved = reader.checkpoint(path, row['sha256'])
        check_frame(arrays, saved, new, number, static['particle_id'], static['mass_kg'], static['component_id'])
        require(saved['step'] == state['step']+1 and row['global_step'] == saved['step'], 'accepted step gap')
        require(saved['parameters'] == parameters and saved['bore_observation']['eligible_ids'] == frozen_eligible, 'physics or eligibility changed')
        require(saved['interval_energy']['records'][:len(prior_prefix)] == prior_prefix, 'complete energy prefix changed')
        require(saved['bore_observation']['events'][:len(prior_events)] == prior_events, 'crossing history changed')
        require(saved['continuation_parent_sha256'] == last_sha, 'checkpoint parent hash differs')
        journal_path = path.parent/'FULL_UNCOMMITTED_candidate.json'; relative = journal_path.relative_to(ROOT).as_posix()
        require(relative in evidence, 'contact journal missing from passed readback binding')
        journal = reader.read(journal_path, evidence[relative]); contact = journal['contact']
        require(contact['features'][0] == 'SHOULDER', 'contact ordering differs')
        require(np.array_equal(np.asarray(contact['history']), saved['contact_history']), 'contact history differs from accepted checkpoint')
        flags = np.asarray(contact['active'], dtype=bool)
        require(flags.shape == (len(static['particle_id']),8), 'contact activity shape differs')
        active.append(flags[:,0]); frames.append(dict(global_step=saved['step'], time_s=saved['time_s'],
            checkpoint=row['checkpoint'], checkpoint_sha256=row['sha256'], contact_source=relative,
            contact_source_sha256=evidence[relative], readback_summary=summary_path.relative_to(ROOT).as_posix()))
        state, last_path, last_sha = saved, path, row['sha256']
    for key in FRAME_KEYS: data[key] = np.concatenate((data[key], new[key][1:]), axis=0)
    if active:
        data['shoulder_contact'] = np.concatenate((data['shoulder_contact'], np.asarray(active)), axis=0)
        data['shoulder_contact_available'] = np.concatenate((data['shoulder_contact_available'], np.ones(len(active),bool)))
    require(state['step'] == summary['accepted_global_steps'] and state['time_s'] == summary['accepted_time_s'], 'readback endpoint differs')
    return state, last_path, last_sha, summary_path, reader.bound[summary_path.relative_to(ROOT).as_posix()]


def inventory(args, state, report):
    p = state['parameters']; local = (args['position_m'][-1]-args['tool_origin_m'][-1])@args['tool_rotation'][-1]
    ids = args['particle_id']; events = state['bore_observation']
    remaining = np.isin(ids,events['eligible_ids']) & ~np.isin(ids,events['crossing_ids'])
    above = np.flatnonzero(remaining & (local[:,2]>0)); order = above[np.argsort(local[above,2])]
    initial = args['position_m'][0]; substrate = args['component_id']==0; bar = args['component_id']==1
    length = report['eligible_initial_mass_kg']/p['rho0']/p['width']**2
    return dict(scope='GEOMETRY_AND_PRESCRIBED_KINEMATICS_ONLY_NOT_MATERIAL_TRANSPORT_OR_VALIDITY_PREDICTION',
        total_mass_kg=float(args['mass_kg'].sum()), bar_mass_kg=float(args['mass_kg'][bar].sum()),
        initial_substrate_centre_bounds_m=[initial[substrate].min(axis=0).tolist(),initial[substrate].max(axis=0).tolist()],
        nominal_substrate_width_m=p['substrate_width'],tool_shoulder_diameter_m=2*p['tool_radius'],
        full_shoulder_footprint_initially_contained=False,
        prescribed_feed_travel_since_initial_m=p['feed']*state['time_s'],
        prescribed_traverse_travel_since_initial_m=p['traverse']*state['time_s'],
        prescribed_feed_time_per_spacing_s=p['dx']/p['feed'],
        eligible_initial_equivalent_length_over_nominal_bar_area_m=length,
        prescribed_feed_time_for_eligible_equivalent_length_s=length/p['feed'],
        tool_rotation_period_s=2*np.pi/p['omega'],
        nearest_uncrossed_eligible_centres_above_outlet=[dict(particle_id=int(ids[k]),height_m=float(local[k,2]),
            square_projection_margin_m=float(p['bore_width']/2-np.max(np.abs(local[k,:2])))) for k in order[:8]],
        limitations=['The 4 mm coupon does not contain the 5.4 mm shoulder footprint even initially; no full-width stable deposition interval is established.',
            'The eligible inventory is finite and no sustained inlet supplies new material.',
            'Feed-length conversions are actuator kinematics, not free-material arrival predictions, execution estimates or validated runtime limits.',
            'An observed finite-coupon outlet/spreading event cannot establish full-track or steady deposition.'])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--readback', action='append', required=True, help='Chronological extra PASS summary.json or its directory; repeat for later windows.')
    parser.add_argument('--output', required=True)
    options = parser.parse_args(argv); output = Path(options.output).resolve()
    require(output.is_relative_to(ROOT/'outputs'), 'output must be under repository outputs')
    require(not output.exists(), 'refuse to overwrite evidence'); output.mkdir(parents=True)
    start = time.perf_counter(); reader = Reader(); report = {}; code = 2
    try:
        for name in ('scripts/readback_deposition_progress_tight.py','scripts/deposition_stage_observer.py','scripts/check_deposition_stage_observer.py'):
            reader.bind(ROOT/name)
        previous_observation = reader.read(OBSERVER, OBSERVER_SHA)
        prior = reader.read(OLD/'summary.json', OLD_SHA)
        require(prior['status'] == 'PASS_ACTUAL_CHAIN_READBACK_NOT_PHYSICAL_OR_EXPERIMENTAL_VALIDATION', 'old actual chain audit absent')
        data = reader.arrays(OLD/'accepted_trajectory.npz', prior['trajectory_sha256'])
        static = dict(particle_id=data['particle_id'], mass_kg=data['particle_property__m'][0,:,0],
            component_id=data['particle_property__component_id'][0,:,0])
        data = {**{k:data[k] for k in FRAME_KEYS}, 'shoulder_contact':data['contact_active'][:,:,0],
            'shoulder_contact_available':data['contact_activity_available']}
        frames = []; source = prior['source_and_input_sha256']; initial_path = ROOT/'outputs/discrete_contact_process_20261004/v3/quarters/independent_initial.npz'
        last_path = initial_path; last_sha = source[initial_path.relative_to(ROOT).as_posix()]
        arrays, state = reader.checkpoint(last_path, last_sha, legacy_audited=True)
        check_frame(arrays,state,data,0,static['particle_id'],static['mass_kg'],static['component_id'])
        frames.append(dict(global_step=0,time_s=state['time_s'],checkpoint=last_path.relative_to(ROOT).as_posix(),checkpoint_sha256=last_sha,contact_available=False))
        for number,row in enumerate(prior['rows'],1):
            require(row['previous_checkpoint_sha256'] == last_sha, 'old accepted hash chain differs')
            last_path = ROOT/row['checkpoint']; last_sha = row['sha256']; arrays,state = reader.checkpoint(last_path,last_sha,legacy_audited=True)
            check_frame(arrays,state,data,number,static['particle_id'],static['mass_kg'],static['component_id'])
            require(state['step'] == number, 'old actual step gap')
            require(sorted(static['particle_id'][data['shoulder_contact'][number]].tolist()) == sorted(row['shoulder_active_ids']), 'audited old contact activity differs')
            frames.append(dict(global_step=number,time_s=state['time_s'],checkpoint=row['checkpoint'],checkpoint_sha256=last_sha,
                contact_source='outputs/coupled_cost_20261004/batched/readback_pre_occlusion_39_375ms/accepted_trajectory.npz',contact_source_sha256=prior['trajectory_sha256'],contact_available=bool(data['shoulder_contact_available'][number])))
        current = (state,last_path.resolve(),last_sha,(OLD/'summary.json').resolve(),OLD_SHA)
        current = append_readback(reader,BASE/'summary.json',data,static,frames,*current[3:],*current[:3],base=True)
        # Replay the immutable observer baseline before extending it.
        state,last_path,last_sha,previous_summary,previous_sha = current
        args = dict(data,**static,eligible_ids=state['bore_observation']['eligible_ids'],crossing_events=state['bore_observation']['events'],
            dx_m=state['parameters']['dx'],bore_width_m=state['parameters']['bore_width'],tool_radius_m=state['parameters']['tool_radius'])
        baseline = observe_deposition_stages(**args)
        require(baseline['frame_count']==448 and baseline['final']==previous_observation['final'] and
            baseline['crossed_particle_rows']==previous_observation['crossed_particle_rows'],'old observation baseline did not replay')
        for summary in options.readback:
            current = append_readback(reader,summary,data,static,frames,previous_summary,previous_sha,state,last_path,last_sha)
            state,last_path,last_sha,previous_summary,previous_sha = current
        args.update(data,crossing_events=state['bore_observation']['events'])
        require(len(data['time_s']) == state['step']+1 == len(frames), 'complete actual history length differs')
        for value in args.values():
            if isinstance(value,np.ndarray): value.flags.writeable = False
        report = observe_deposition_stages(**args); report['finite_inventory_and_geometric_references'] = inventory(args,state,report)
        np.savez_compressed(output/'observation_inputs.npz',**{k:v for k,v in args.items() if isinstance(v,np.ndarray)})
        (output/'actual_frame_bindings.json').write_text(json.dumps(frames,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        reader.unchanged()
        report.update(status='PASS_PASSIVE_ACCEPTED_ARRAY_OBSERVATION',command=[sys.executable,'-B',str(Path(__file__).resolve()),*(sys.argv[1:] if argv is None else argv)],
            input_and_source_sha256=reader.bound,accepted_checkpoint=last_path.relative_to(ROOT).as_posix(),accepted_checkpoint_sha256=last_sha,
            accepted_global_steps=state['step'],actual_frames_only=True,interpolated_frames=0,added_accepted_time_s=0.,
            geometry_evaluations=0,nnps_evaluations=0,old_448_frame_observer_replayed=True,
            observation_inputs_sha256=sha(output/'observation_inputs.npz'),actual_frame_bindings_sha256=sha(output/'actual_frame_bindings.json'))
        code = 0
    except Exception:
        report.update(status='FAIL_PASSIVE_OBSERVATION',error=traceback.format_exc(),input_and_source_sha256=reader.bound)
    report.update(exit_code=code,wall_seconds=time.perf_counter()-start)
    (output/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:report.get(k) for k in ('status','frame_count','time_end_s','accepted_global_steps','final','crossed_particle_rows','error','wall_seconds')},indent=2))
    return code


if __name__ == '__main__': raise SystemExit(main())
