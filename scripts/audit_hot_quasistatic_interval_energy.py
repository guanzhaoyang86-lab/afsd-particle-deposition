"""Additional whole-interval audit; retains the original per-step decisions."""
import json
import math
from pathlib import Path
from run_hot_quasistatic_probe import BASE, ROOT, write, sha


def main():
    directory = BASE/'refinement'
    history = json.loads((directory/'accepted_history.json').read_text())
    output = BASE/'readback/interval_energy_audit.json'
    if output.exists(): raise ValueError('refuse audit overwrite')
    raw_terms, scale_terms, rows = [], [], []
    for entry in history:
        energy = json.loads((directory/entry['checkpoint']).with_name('energy.json').read_text())['mechanical']
        raw = energy['global_raw_residual_J']; scale = energy['absolute_work_scale_J']
        raw_terms.append(raw); scale_terms.append(scale)
        accumulated = math.fsum(raw_terms); accumulated_scale = math.fsum(scale_terms)
        limit = 1e-7+.05*accumulated_scale
        rows.append(dict(time_s=entry['time_s'], per_step_raw_J=raw,
            per_step_absolute_work_scale_J=scale, per_step_raw_to_scale=raw/scale,
            per_step_registered_acceptance=abs(raw)<=energy['registered_energy_limit_J'],
            accumulated_raw_J=accumulated, accumulated_absolute_work_scale_J=accumulated_scale,
            single_interval_absolute_allowance_J=1e-7, single_interval_limit_J=limit,
            single_interval_gate_pass=abs(accumulated)<=limit))
    result = dict(status='FAIL_WHOLE_INTERVAL_ENERGY_SCREEN', physical_evaluations=0,
        original_per_step_acceptance_records_unchanged=True,
        interval_screen_is_additional_not_claimed_original_driver_gate=True,
        per_step_pass_does_not_certify_interval=True,
        accepted_time_s=history[-1]['time_s'],
        sum_raw_J=math.fsum(raw_terms), sum_absolute_work_scale_J=math.fsum(scale_terms),
        accumulated_raw_fraction_of_scale=math.fsum(raw_terms)/math.fsum(scale_terms),
        interval_limit_with_one_absolute_allowance_J=1e-7+.05*math.fsum(scale_terms),
        first_prefix_fail_time_s=next(r['time_s'] for r in rows if not r['single_interval_gate_pass']),
        process_continuation_qualified=False,
        reason='Absolute tolerance was applied once per step; it is not a whole-process energy certificate.',
        rows=rows, source_sha256={str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))})
    assert all(row['per_step_registered_acceptance'] for row in rows)
    assert not rows[-1]['single_interval_gate_pass']
    write(output, result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','source_sha256')},indent=2))


if __name__ == '__main__': main()
