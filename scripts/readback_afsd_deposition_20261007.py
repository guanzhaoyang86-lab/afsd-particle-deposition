"""Read screened deposition stage prefixes without performing new physics."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import afsd_deposition_setup_20261007 as model
import afsd_process_framework_20261007 as identity_model
import run_afsd_process_framework_20261007 as io
from afsd_deposition_evidence_20261007 import summarize_frames

ROOT = model.ROOT


def checked_screened_rows(meta):
    rows = meta.get('screened_history', [])
    screens = meta.get('pair_screens', [])
    if len(rows) != 2*len(screens):
        raise ValueError('each screened macro must contain exactly two fine states')
    for index, declared_screen in enumerate(screens):
        first, last = rows[2*index:2*index+2]
        screen_path = ROOT/first['screen_path']
        if (first['screen_path'] != last['screen_path'] or
                first['screen_sha256'] != last['screen_sha256'] or
                io.core.sha(screen_path) != first['screen_sha256']):
            raise ValueError('fine states do not share the same intact pair screen')
        screen = io.read(screen_path)
        if (screen != declared_screen or not screen['passed'] or
                not screen['global_screen']['pass'] or not screen['local_screen']['passed'] or
                not screen['crossing_ids_match'] or
                [first['half_index'], last['half_index']] != [1, 2] or
                first['macro_index'] != last['macro_index'] or
                first['macro_index'] != screen['macro_index'] or
                [first['sha256'], last['sha256']] != screen['fine_checkpoint_sha256'] or
                not first['time_s'] < last['time_s'] or last['time_s'] != screen['time_s']):
            raise ValueError('screen must bind these exact passing fine states in order')
        for row in (first, last):
            observer_path = ROOT/row['observer_path']
            if io.core.sha(observer_path) != row['observer_sha256']:
                raise ValueError('saved observer differs')
    return rows


def readback(directory):
    chain = []
    seen = set()
    references = {}
    current = directory / 'manifest.json'
    while current is not None:
        current = current.resolve()
        if current in seen or not current.is_relative_to(ROOT/'outputs'):
            raise ValueError('invalid or cyclic stage evidence chain')
        seen.add(current)
        meta = io.read(current)
        supervision_path = current.parent/'supervision.json'
        supervision = io.read(supervision_path)
        if (not supervision['source_inputs_still_match'] or
                supervision['manifest_sha256'] != io.core.sha(current) or
                supervision['log_sha256'] != io.core.sha(current.parent/'run.log')):
            raise ValueError('saved stage supervision differs')
        if ((meta['status'].startswith('PASS') or chain) and
                (not meta['status'].startswith('PASS') or supervision['exit_code'] != 0)):
            raise ValueError('claimed successful or ancestor stage did not exit successfully')
        io.verify(meta['source_sha256'])
        io.verify(meta['artifacts_sha256'])
        for path in (current, supervision_path, current.parent/'run.log'):
            references[io.rel(path)] = io.core.sha(path)
        references.update(meta['artifacts_sha256'])
        chain.append(meta)
        parent = meta.get('parent_stage_manifest')
        if parent:
            current = ROOT/parent
            if io.core.sha(current) != meta['parent_stage_manifest_sha256']:
                raise ValueError('parent stage manifest differs')
        else:
            current = None
    chain.reverse()
    first = chain[0]
    initial = ROOT/first['initial_checkpoint']
    if io.core.sha(initial) != first['initial_sha256']:
        raise ValueError('independent initial checkpoint differs')
    checkpoint_rows = [dict(checkpoint=io.rel(initial), sha256=first['initial_sha256'])]
    references[io.rel(initial)] = first['initial_sha256']
    for meta in chain:
        if meta['initial_checkpoint'] != first['initial_checkpoint'] or meta['initial_sha256'] != first['initial_sha256']:
            raise ValueError('stages do not share the same independent initial state')
        checkpoint_rows.extend(checked_screened_rows(meta))
    frames = []
    previous_state = None
    p = None
    for row in checkpoint_rows:
        path = ROOT/row['checkpoint']
        if not path.resolve().is_relative_to(ROOT/'outputs') or io.core.sha(path) != row['sha256']:
            raise ValueError('screened candidate checksum differs')
        pa, state = io.load(path)
        program = model.validate_state(state)
        if identity_model.identity_digest(pa) != state['immutable_material_identity_sha256']:
            raise ValueError('immutable particle identity differs')
        if previous_state is not None:
            n = previous_state['step']
            if (state['step'] != n+1 or state['time_s'] <= previous_state['time_s'] or
                    state['machine_history'][:n] != previous_state['machine_history'] or
                    state['interval_energy']['records'][:n] != previous_state['interval_energy']['records']):
                raise ValueError('discontinuous physical history in screened prefix')
        elif state['step'] != 0 or state['time_s'] != 0.:
            raise ValueError('readback must begin at the independent initial state')
        p = state['parameters']
        pose = program.global_pose(state['time_s'])
        bore = state['bore_observation']
        frames.append(dict(time_s=state['time_s'], ids=pa.particle_id.copy(),
            position=np.column_stack((pa.x, pa.y, pa.z)), mass=pa.m.copy(), ep=pa.ep.copy(),
            eligible_ids=bore['eligible_ids'], crossing_ids=bore['crossing_ids'],
            tool_origin=pose['origin'], tool_rotation=pose['rotation']))
        references[io.rel(path)] = row['sha256']
        previous_state = state
    result = summarize_frames(frames, p)
    result.update(status='PASS_SAVED_SCREENED_PREFIX_READBACK_NOT_NEW_PHYSICS',
        stage_count=len(chain), last_stage_status=chain[-1]['status'],
        last_stage_successful=chain[-1]['status'].startswith('PASS'),
        source_sha256={io.rel(Path(__file__)): io.core.sha(Path(__file__)),
            'scripts/afsd_deposition_evidence_20261007.py': io.core.sha(ROOT/'scripts/afsd_deposition_evidence_20261007.py')},
        evidence_sha256=references, final_checkpoint=checkpoint_rows[-1]['checkpoint'],
        final_checkpoint_sha256=checkpoint_rows[-1]['sha256'],
        solver_qualification='LOCAL_STEP_DOUBLING_SCREENED_PREFIX_ENGINEERING_APPROXIMATE',
        formal_accepted_time_advanced_s=0., old_formal_process_connected=False)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    stage, output = (ROOT/args.stage).resolve(), (ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT/'outputs') or output.exists():
        raise ValueError('new readback output inside project outputs required')
    start = time.perf_counter()
    result = readback(stage)
    result['wall_seconds'] = time.perf_counter()-start
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: result[k] for k in ('status', 'frame_count', 'first_supplied_shoulder_time_s',
        'first_supplied_trailing_time_s', 'wall_seconds')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
