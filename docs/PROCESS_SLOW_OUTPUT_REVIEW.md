# Independent saved-output audit of the slow-loading S0 prefix

2026-10-03. Status: **PASS_OUTPUT_INTEGRITY_RAW_PREFIX_COMPARISON_NOT_PROCESS_ACCURACY**.
This is a read-only audit, not a trajectory run, restart, convergence result or
20-microsecond feasibility claim. The original scout remains
`STOPPED_LOCAL_SCOUT_NOT_PROCESS_FEASIBLE`. No source, original output, physical
state or prior failure was modified.

## Verified publication and field identity

Audited `outputs/process_slow_20261003/pilot_initial`, `S0_initial`, and
`logs/process_slow_20261003/reference_repair1` against their recorded receipts.
The preflight's 77 source and 19 artifact hashes, each run's 49 source hashes,
the pilot's 8 artifact hashes, S0's 38 artifact hashes, and the independent
reference's 4 source and 7 artifact hashes all matched. Both runs bind the
same preflight, reference history and physical configuration.

S0 contains exactly eight published stages and nine frames including the
initial state, at 250 ns intervals through 2 microseconds. Every checkpoint's
publication hash, independently recomputed semantic state hash, parameter
hash, full 40-property particle hash and saved particle-property snapshot
matched. Published parent hashes form the expected chain. IDs 1..27, literal
trapezoid masses totaling 0.0027 kg and passive labels remain intact.

All eight VTK files read back exactly to the corresponding checkpoint:
positions, velocity, canonical HHT acceleration, IDs, mass, tags, density,
temperature, pressure, equivalent stress, plastic strain and all six published
full-stress components. The JSON frames independently match checkpoint full
stress, force, acceleration, temperature, density, velocity and position;
per-face reactions and cumulative ledgers agree with their saved records.
The pilot's first accepted state equals S0's first accepted state bitwise.

The last published time is `2e-6 s`, with semantic state hash
`8c94bbdc6d474545a3eb58988dc10b9bfe4339fb2d79a33977a97778140f9283`.
The rejected 2–2.5 microsecond candidate has no ninth stage, frame or accepted
history entry. Its record explicitly forbids acceptance and restart. It records
one nonlinear evaluation, sampled Coulomb defect `0.018714892092337464` versus
`0.01`, `whole_path_passive=true` and `requires_event_split=false`. These are
saved-candidate facts; this output audit does not infer a physical stop or the
cause of the candidate's trajectory error.

All 4067 entries in the preceding HHT freeze manifest matched before and after
the audit: zero unauthorized changes and zero index-document changes at the
time of this check.

## Raw errors at identical saved times

The independent reference certifies common snapshots through
`5.9374999999999995e-6 s`; its event status remains `EVENT_APPROACH_BLOCKED` and
its complete-20-microsecond flag is false. The audit compares all nine S0 frames
against their existing independent snapshots at the same times. No interpolation,
ODE integration or tested HHT force map generates the comparison. Full stress
uses saved deviatoric stress and density; face force is independently reconstructed
from saved oracle positions/velocities and the literal slow wall programme.

These are raw errors for one coarse, incomplete schedule. The registered
finest-grid complete-trajectory bounds and refinement gates are not applied
to reclassify this prefix as accurate or converged.

| Quantity | Maximum over 0–2 us | At 2 us |
|---|---:|---:|
| Volume-weighted RMS position, m | 1.46648747e-9 | 1.24368048e-9 |
| Volume-weighted RMS velocity, m/s | 2.97608739e-3 | 2.97608739e-3 |
| Volume-weighted RMS full stress, Pa | 26104.56495 | 15915.23503 |
| Volume-weighted RMS deviatoric stress, Pa | 21823.92221 | 14010.42469 |
| Maximum density error, kg/m3 | 1.47468561e-3 | 1.18375782e-3 |
| Maximum temperature error, K | 2.84842372e-8 | 2.84842372e-8 |
| Maximum face-force component error, N | 0.141803686 | 0.100499057 |
| Maximum face-impulse component error, N s | 1.20610149e-7 | 1.20610149e-7 |
| Maximum face wall-work error, J | 4.25216326e-9 | 4.25216326e-9 |
| Maximum face raw-friction-heat error, J | 8.71401694e-9 | 8.71401694e-9 |
| Maximum signed face-reservoir-heat error, J | 6.00761758e-13 | 6.00761758e-13 |

Per-time rows, exact hashes and checks are retained in the JSON evidence.
No 20-us completion, temporal convergence, post-event branch or transport
result follows from these numbers.

## Reproduction evidence

Command: `.envs/pysph/python.exe -B logs/process_slow_20261003/output_readback.py`.
External PowerShell supervisor imposed a 60 s limit. Exit 0; external elapsed
3.0497596 s, in-script elapsed 2.6927744 s. No repair was needed. The script
refuses to overwrite its JSON evidence.

- Script: `logs/process_slow_20261003/output_readback.py`;
  SHA256 `d894a6a5b935f869b6105dc417416bc613e555e74f190600efc8f0777f0fef07`.
- Results: `logs/process_slow_20261003/output_readback.json`;
  SHA256 `5080c30d2709dcbf3aa2f5019dd81c29ffe3d03caadbf23dfad3e9b137828289`.
- Raw logs: `output_readback_initial.log`, `output_readback_initial_stderr.log`;
  supervision: `output_readback_initial_supervision.log`, in the same directory.
