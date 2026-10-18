"""Freeze/test the general finite-patch enclosure on a saved failed path."""
from pathlib import Path
import json
import os
import subprocess
import sys
import time
import numpy as np
import run_shoulder_window as window
from batched_contact import BatchedChainRuleContact
from contact_shell_finite_enclosure import audit_contact_shell_finite_enclosure

core, post, event = window.core, window.post, window.event


def main():
    output = core.ROOT/'outputs/coupled_cost_20261004/finite_shell_enclosure_event.json'
    logs = core.ROOT/'logs/coupled_cost_20261004/finite_shell_enclosure'
    if output.exists() or logs.exists():
        raise ValueError('refuse evidence overwrite')
    logs.mkdir(parents=True)
    sources = {p.relative_to(core.ROOT).as_posix():core.sha(p) for p in
        (Path(__file__).resolve(), core.ROOT/'scripts/contact_shell_finite_enclosure.py',
         core.ROOT/'tests/review/test_contact_shell_finite_enclosure.py',
         core.ROOT/'scripts/contact_shell_audit.py', core.ROOT/'scripts/batched_contact.py')}
    command = [sys.executable, '-B', '-m', 'unittest', 'discover', '-s', 'tests/review',
               '-p', 'test_contact_shell_finite_enclosure.py', '-v']
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
    directory = window.BASE/'coarse'
    manifest, history = post.read(directory/'manifest.json'), post.read(directory/'candidate_history.json')
    if manifest['status'] != 'FAILED_NO_NEW_ACCEPTED_STATE' or len(history) != 18 or manifest['accepted_steps'] != 0:
        raise ValueError('registered failed coarse evidence absent')
    post.check_sources(manifest['source_sha256'])
    parent = directory/history[-1]['checkpoint']
    if core.sha(parent) != history[-1]['sha256']:
        raise ValueError('reference candidate hash differs')
    candidate = directory/'step_019/NOT_RESTARTABLE_candidate.npz'
    old = core.load(str(parent)); new = core.load(str(candidate))
    pa, state = old['arrays']['solid'], old['solver_data']
    newpa, newstate = new['arrays']['solid'], new['solver_data']
    if state['restart_allowed'] is not False or newstate['restart_allowed'] is not False:
        raise ValueError('coarse references unexpectedly restartable')
    if state['step'] != 171 or newstate['step'] != 172 or abs(state['time_s']-.0265)>1e-15 or abs(newstate['time_s']-.02675)>1e-15:
        raise ValueError('saved reference interval differs')
    if newstate['interval_energy']['records'][:-1] != state['interval_energy']['records']:
        raise ValueError('candidate energy prefix differs')
    if not np.array_equal(pa.particle_id, newpa.particle_id):
        raise ValueError('particle IDs differ')
    p = state['parameters']
    contact = BatchedChainRuleContact(p['tool_radius'], p['bore_width'], p['tool_top']-p['layer'],
        p['dx']/2, p['normal_stiffness'], p['normal_stiffness'], p['friction'], base_z=-p['substrate_depth'])
    first, last = np.column_stack((pa.x,pa.y,pa.z)), np.column_stack((newpa.x,newpa.y,newpa.z))
    started = time.perf_counter()
    result = audit_contact_shell_finite_enclosure(contact, first, last,
        core.pose(p,state['time_s']), core.pose(p,newstate['time_s']))
    old_audit_file = directory/'step_019/continuous_shell_audit.json'
    old_audit = post.read(old_audit_file)
    if result['original_audit'] != old_audit['original_audit']:
        raise ValueError('original Lipschitz certificate replay differs')
    paths = [parent,candidate,directory/'manifest.json',directory/'candidate_history.json',
             old_audit_file,logs/'receipt.json']
    result.update(source_sha256=dict(manifest['source_sha256'], **sources),
        input_sha256={p.relative_to(core.ROOT).as_posix():core.sha(p) for p in paths},
        parent_time_s=state['time_s'],candidate_time_s=newstate['time_s'],physical_evaluations=0,
        physical_time_advanced_s=0.,wall_seconds=time.perf_counter()-started,
        old_failed_stage_reclassified=False,candidate_accepted=False)
    post.check_sources(result['source_sha256'])
    core.write(output,result)
    print(json.dumps({key: result[key] for key in ('status','geometry_queries',
        'additional_finite_patch_queries','physical_evaluations','wall_seconds')}))
    return 0 if result['status'] == event.CERTIFIED else 3


if __name__ == '__main__':
    raise SystemExit(main())
