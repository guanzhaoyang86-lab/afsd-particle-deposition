"""Bounded component tests and saved-path enclosure; no physical evaluation."""
from pathlib import Path
import json
import os
import subprocess
import sys
import time
import numpy as np
import run_first_shoulder_event as event
from batched_contact import BatchedChainRuleContact
from contact_shell_enclosure import audit_contact_shell_enclosure, CERTIFIED

core, post = event.core, event.post


def main():
    output = core.ROOT/'outputs/coupled_cost_20261004/shoulder_shell_enclosure.json'
    logs = core.ROOT/'logs/coupled_cost_20261004/shell_enclosure'
    if output.exists() or logs.exists():
        raise ValueError('refuse evidence overwrite')
    logs.mkdir(parents=True)
    sources = {p.relative_to(core.ROOT).as_posix(): core.sha(p) for p in
        (Path(__file__).resolve(), core.ROOT/'scripts/contact_shell_enclosure.py',
         core.ROOT/'tests/review/test_contact_shell_enclosure.py',
         core.ROOT/'scripts/contact_shell_audit.py', core.ROOT/'scripts/batched_contact.py')}
    command = [sys.executable, '-B', '-m', 'unittest', 'discover', '-s', 'tests/review',
               '-p', 'test_contact_shell_enclosure.py', '-v']
    begin = time.perf_counter()
    with (logs/'tests.log').open('w', encoding='utf-8') as stream:
        try:
            code = subprocess.run(command, cwd=core.ROOT, stdout=stream, stderr=subprocess.STDOUT,
                env=dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1'), timeout=60).returncode
        except subprocess.TimeoutExpired:
            code = 124
    core.write(logs/'receipt.json', dict(command=command, exit_code=code, source_sha256=sources,
        physical_evaluations=0, physical_time_advanced_s=0., wall_seconds=time.perf_counter()-begin))
    if code:
        return code
    directory = event.BASE/'approach'
    meta, history = post.read(directory/'manifest.json'), post.read(directory/'accepted_history.json')
    if meta['status'] != 'FAILED_WITH_ACCEPTED_PREFIX' or len(history) != 15:
        raise ValueError('registered failed event evidence absent')
    post.check_sources(meta['source_sha256'])
    parent = directory/history[-1]['checkpoint']
    if core.sha(parent) != history[-1]['sha256']:
        raise ValueError('direct parent hash differs')
    pa, state = core.restore(parent)
    candidate = directory/'step_016/NOT_RESTARTABLE_candidate.npz'
    data = core.load(str(candidate)); newpa, newstate = data['arrays']['solid'], data['solver_data']
    if newstate['restart_allowed'] is not False or newstate['step'] != state['step']+1:
        raise ValueError('candidate classification differs')
    if not np.array_equal(pa.particle_id, newpa.particle_id):
        raise ValueError('particle IDs differ')
    p = state['parameters']
    contact = BatchedChainRuleContact(p['tool_radius'], p['bore_width'], p['tool_top']-p['layer'],
        p['dx']/2, p['normal_stiffness'], p['normal_stiffness'], p['friction'], base_z=-p['substrate_depth'])
    first, last = np.column_stack((pa.x, pa.y, pa.z)), np.column_stack((newpa.x, newpa.y, newpa.z))
    started = time.perf_counter()
    result = audit_contact_shell_enclosure(contact, first, last,
        core.pose(p, state['time_s']), core.pose(p, newstate['time_s']))
    old = post.read(directory/'step_016/continuous_shell_audit.json')
    if result['original_audit'] != old:
        raise ValueError('original certificate replay differs')
    for row in result['planar_certificates']:
        row['particle_id'] = int(pa.particle_id[row['particle_index']])
    paths = [parent, candidate, directory/'manifest.json', directory/'accepted_history.json',
        directory/'step_016/continuous_shell_audit.json', logs/'receipt.json']
    input_sha = {path.relative_to(core.ROOT).as_posix(): core.sha(path) for path in paths}
    result.update(source_sha256=dict(meta['source_sha256'], **sources), input_sha256=input_sha,
        parent_time_s=state['time_s'], candidate_time_s=newstate['time_s'], physical_evaluations=0,
        physical_time_advanced_s=0., wall_seconds=time.perf_counter()-started,
        old_failed_stage_reclassified=False, candidate_accepted=False)
    post.check_sources(result['source_sha256'])
    core.write(output, result)
    print(json.dumps({key: result[key] for key in ('status', 'geometry_queries',
        'additional_finite_patch_queries', 'plane_path_evaluations', 'physical_evaluations', 'wall_seconds')}))
    return 0 if result['status'] == CERTIFIED else 3


if __name__ == '__main__':
    raise SystemExit(main())
