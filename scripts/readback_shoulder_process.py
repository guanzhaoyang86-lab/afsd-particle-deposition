"""Read-only 22 ms shoulder-chain adapter over the frozen batched reader.

Only the explicit stage list is adapted. All checkpoint, field, parent, source,
energy, shell, scalar-equivalence, and passive-crossing audits remain in the
original package function. Future shoulder windows are deliberately unsupported.
"""
from pathlib import Path
import argparse
import json
import readback_batched_process as base

# Capture the original function before replacement: the adapter must never call
# itself through base._default_stages after installing the controlled list.
_ORIGINAL_DEFAULT_STAGES = base._default_stages


def shoulder_stages():
    stages = list(_ORIGINAL_DEFAULT_STAGES())
    base.require(len(stages) == 5, 'unexpected frozen base stage layout')
    return stages+[
        ('first_shoulder_approach', base.BASE/'first_shoulder_event/approach', 15, None, 'event'),
        ('first_shoulder_enclosed_halves', base.BASE/'first_shoulder_enclosed/halves', 2, None, 'event')]


def package(output):
    output = Path(output).resolve()
    base.require(not output.exists(), 'refuse to overwrite shoulder readback')
    source_paths = (Path(__file__).resolve(), Path(base.__file__).resolve(),
                    Path(base.reader.__file__).resolve())
    source_hashes = {path.relative_to(base.ROOT).as_posix(): base.sha(path) for path in source_paths}
    old = base.BASE/'readback_21_875ms'
    old_hashes = {path.relative_to(base.ROOT).as_posix(): base.sha(path)
                  for path in (old/'summary.json', old/'accepted_trajectory.npz')}
    stages = shoulder_stages()
    base.require(base._default_stages is _ORIGINAL_DEFAULT_STAGES,
                 'underlying reader already has another stage adapter installed')
    base._default_stages = shoulder_stages
    try:
        # No append argument: event behavior is selected by the explicit modes,
        # so the original directory whitelist is neither changed nor bypassed
        # through an unvalidated arbitrary CLI path.
        summary = base.package(output)
    finally:
        base._default_stages = _ORIGINAL_DEFAULT_STAGES
    base.require(summary['accepted_steps'] == 153 and summary['samples'] == 154,
                 'unexpected 22 ms accepted chain length')
    base.require(abs(summary['accepted_time_s']-.022) <= 1e-15,
                 'unexpected shoulder readback endpoint')
    for name, digest in dict(source_hashes, **old_hashes).items():
        base.require(base.sha(base.ROOT/name) == digest, 'source or earlier readback changed: '+name)
    summary.update(
        shoulder_reader_adapter_sha256=source_hashes[Path(__file__).resolve().relative_to(base.ROOT).as_posix()],
        unchanged_batched_reader_sha256=source_hashes[Path(base.__file__).resolve().relative_to(base.ROOT).as_posix()],
        unchanged_discrete_reader_helper_sha256=source_hashes[Path(base.reader.__file__).resolve().relative_to(base.ROOT).as_posix()],
        controlled_stage_list=[dict(name=name, directory=directory.relative_to(base.ROOT).as_posix(),
            accepted_count=count, required_status=status, audit_mode=mode)
            for name, directory, count, status, mode in stages],
        original_default_stages_restored=True, append_argument_used=False,
        future_shoulder_window_schema_supported=False,
        earlier_readback_21_875ms_sha256=old_hashes,
        event_mode_preserves_continuous_shell_and_scalar_equivalence_checks=True,
        shoulder_global_and_local_time_screens_required=True)
    (output/'summary.json').write_text(json.dumps(base.clean(summary), indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(base.clean(dict(adapter_status='PASS_22MS_ACTUAL_CHAIN_ADAPTER',
        physical_evaluations=0, accepted_steps=summary['accepted_steps'], samples=summary['samples'],
        time_s=summary['accepted_time_s'], shoulder_active_count=summary['final']['shoulder_active_count'],
        crossing_mass_kg=summary['final']['crossing_mass_kg']))))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=base.BASE/'readback_22ms')
    arguments = parser.parse_args()
    package(arguments.output)
