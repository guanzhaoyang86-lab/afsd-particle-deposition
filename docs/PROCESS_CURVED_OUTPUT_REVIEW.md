# Curved-contact output readback and raw reference comparison

2026-10-03. Independent read-only audit of `outputs/process_curved_20261003/`
`saved_red_initial`, `pilot_initial`, and `S0_initial`. No trajectory, ODE,
NNPS solve, physical event, AFSD step, or production-source edit was performed.

Result: **PASS_OUTPUT_INTEGRITY_RAW_PREFIX_COMPARISON_NOT_PROCESS_ACCURACY**.
This certifies the saved outputs and the stated readback identities only. It
does not certify temporal accuracy, convergence, classical HHT stability,
post-stop contact, the complete 20 us scout, a 1 ms process, or AFSD transport.

## Evidence and bounds

The new reproducible reader is
`logs/process_curved_20261003/output_readback.py`. Its initial execution and
machine-readable results are `output_readback_initial.log`,
`output_readback_initial_stderr.log`, and `output_readback.json` in that
directory. `output_readback_supervision.json` records 3.0099023 s internal,
3.3647168 s external, exit 0, and a 60 s external ceiling. No repair was needed.

Before execution the reader declared exact equality for stored arrays/hashes
and 1e-12 absolute SI plus 1e-8 relative agreement for independently rebuilt
integrals, work and identities. These are integrity tolerances, not newly
chosen physical-accuracy allowances. Reference differences are raw numbers
with no acceptance threshold. An additional direct read-only NPZ inspection
is transcribed in `rejected_candidate_readback.json`.

Both before and after the audit, all **4,221** files in
`docs/PROCESS_SLOW_FREEZE.json` were present and unchanged. The three index
documents were allowed by the preservation policy but had not changed either.

## Published states and receipts

The 189-test preflight is PASS. Its 87 source hashes and 22 artifact hashes
match. The independent reference's four source hashes and seven artifact
hashes match; it remains `PASS_PREFIX_REFERENCE_EVENT_BLOCKED`, certified
through 5.9375 us, with `EVENT_APPROACH_BLOCKED` thereafter. Each run has 56
matching source hashes. Artifact counts checked were two for saved RED,
eight for pilot, and 59 for S0. Preflight/reference/history/regression hashes
bind the runs to the same gated inputs; pilot and S0 resolved configurations
and first published state hashes agree.

S0 contains exactly **13** published stage directories, checkpoints and VTK
files, plus **14** frames including the initial condition. The last accepted
time is 5 us. There is no stage 014. Checks include:

- The new `HHT_CURVED_CONTACT_IMPULSE_ELASTIC_BE_STATE_V1` schema, complete
  recursive semantic hashes, parent chains, parameters and particle hashes.
- All saved particle properties, persistent IDs 1–27, literal unscaled
  trapezoidal physical masses, passive labels and zero elastic plastic strain.
- Exact NPZ/frame/VTK positions, velocities, density, temperature, acceleration,
  all six full-stress components and the derived pressure/equivalent stress.
- The complete `last_interval`, including actual contact impulse/work/raw Q,
  DeltaJ, DeltaW, quadrature, signed filter/root work and corrected identities.
- The actual accumulated face impulse, wall work, raw friction heat,
  reservoir exchange and numerical filter/root work at each published stage.

These are solved saved fields, not synthetic visualization-only values.

## Independent interval reconstruction

The reader calls neither the tested boundary/contact integrator nor its
identity helper. It reconstructs the compatible HHT cubic position and
quadratic velocity from saved endpoint histories, then applies the literal
5 Hz moving-plane geometry, penalty stiffness 7e6 N/m and vector Coulomb
coefficient 0.2. Each saved Gauss fraction/weight, actual normal/tangential
force, impulse, material work, wall work and raw friction heat is checked.
The reader separately integrates the reconstructed curve using **64-point
Gauss quadrature**, rather than the frozen old two-node chord approximation.
All 13 published intervals used the new 16-point rule.

Across the published intervals, pilot, private saved-RED candidate and rejected
6 us candidate, the largest 64-point versus saved-rule component difference
is **1.95496e-18 N s/J**. The largest independently reconstructed per-particle
momentum and kinetic identity residuals are **4.26495e-21 N s** and
**8.34730e-22 J**. DeltaJ/DeltaW are retained explicitly; no balance residual
is assigned to heat. Normal work is separated from tangential raw Q.

The squared norm of the complete quadratic tangential-slip vector is
independently minimized over endpoints and all real stationary points.
No accepted path has a full-vector zero; the minimum over accepted paths is
0.0437809716 m/s. Scalar sign changes alone are not rejected as stop events.

Thermal raw Q equals the mechanical raw Q without clipping; the declared
one-half partition is preserved. Reservoir heat independently equals
`h * .025 W/K * (Twall - T1)` at each contacted particle. Internal conductive
heat cancels globally. Saved `m * cp * DeltaT` agrees with total heat input to
**8.48951e-15 J** per particle. This checks source attribution and the saved
heat balance, not independent reassembly of the thermal SPH graph.

## Raw same-time reference errors

The 14 new frames are matched to independently saved reference states at
identical times through 5 us, with no interpolation. Vector/tensor RMS errors
use the fixed reference-volume weights; scalar temperature/density and face
force/ledger errors use maximum absolute components.

| Quantity | Maximum over 0–5 us | At 5 us |
|---|---:|---:|
| Position RMS, m | 2.11582939e-8 | 2.11582939e-8 |
| Velocity RMS, m/s | 1.33871128e-2 | 1.32493766e-2 |
| Full stress RMS, Pa | 4.08995291e5 | 4.08995291e5 |
| Deviatoric stress RMS, Pa | 3.64583921e5 | 3.64583921e5 |
| Density maximum, kg/m3 | 2.48724090e-2 | 2.48724090e-2 |
| Temperature maximum, K | 4.96245093e-7 | 4.96245093e-7 |
| Face-force maximum, N | 4.15425495e-1 | 3.27516689e-1 |
| Face-impulse maximum, N s | 3.05680022e-7 | 3.05680022e-7 |
| Wall-work maximum, J | 5.14178426e-8 | 5.14178426e-8 |
| Raw-friction maximum, J | 8.26192942e-8 | 8.26192942e-8 |
| Reservoir-heat maximum, J | 2.80374342e-12 | 2.80374342e-12 |

The frozen old S0 has only nine frames through 2 us. Comparing the old and new
methods against the same reference at those same nine times is a separate
baseline comparison, **not temporal refinement or convergence**:

| Maximum over the common 0–2 us prefix | Old endpoint chord | New curved contact |
|---|---:|---:|
| Position RMS, m | 1.46648747e-9 | 1.37365564e-9 |
| Velocity RMS, m/s | 2.97608739e-3 | 2.92403599e-3 |
| Full stress RMS, Pa | 2.61045649e4 | 2.49625804e4 |
| Temperature maximum, K | 2.84842372e-8 | 2.99140197e-8 |
| Face-force maximum, N | 1.41803686e-1 | 1.31628083e-1 |
| Face-impulse maximum, N s | 1.20610149e-7 | 1.54263697e-7 |
| Wall-work maximum, J | 4.25216326e-9 | 6.06646536e-9 |
| Raw-friction maximum, J | 8.71401694e-9 | 8.93455400e-9 |

Several errors decrease and several increase. Neither improved local contact
integration nor a passing output audit establishes global process accuracy.
The JSON retains every per-frame error, including the omitted table fields.

## Private saved RED and actual stop reason

The retained 2–2.5 us, 500 ns old-method regression reproduces
`contact turning-traction quadrature law defect` with defect
0.018714892092337464. The explicitly migrated new numerical variant solves
that private interval in 12 evaluations. The physical parent and numerical
acceleration history are unchanged by migration. Raw interval friction heat
is 2.4823136327698588e-6 J. Its checkpoint wrapper remains private,
unpublished and non-restartable; it is not an extra accepted S0 frame.

S0's next [5,6] us candidate numerically converged, but publication was denied
because 6 us lies beyond the certified reference prefix ending at 5.9375 us.
`REJECTED.json` says `UNREFERENCED_INTERVAL`, with no contact-event rejection
diagnostic. `REJECTED_UNREFERENCED.npz` is a separate rejected wrapper bound to
the last published 5 us parent. Its stored corrected-root force residual
implies a velocity residual of **3.43364e-10 m/s**, below the unchanged 1e-9
root bound; its raw friction heat is 4.726110714482718e-6 J.

At the rejected endpoint, particle ID 4 has material velocity
`(0.0796679797, 9.74010e-18, 0.1773477918) m/s`. The wall's tangential velocity
is `(0,0,0.2) m/s`; hence tangential slip is
`(0, 9.74010e-18, -0.0226522082) m/s`, with norm **0.0226522082 m/s**. This is
also the minimum full-vector slip over that candidate path. It is not zero.
The supplemental JSON records its x/v/S/rho/T, canonical acceleration, true
endpoint force and correction sizes directly from the saved NPZ.

The nested numerical `candidate_state` retains its internal restart flag,
but the enclosing saved wrapper explicitly has `accepted=false` and
`restart_allowed=false`; it is absent from accepted history. Continuation
authority must come from the publication history and wrapper, not from
extracting an unpublished nested state.

The independent reference approaches a physical full-vector stop near
5.97466 us. This output audit neither reruns nor resolves that event. The
new candidate's nonzero 6 us slip does not validate crossing the reference
stop. Physical stop/static-mode handling remains unimplemented, and no
6 us reference-accuracy, complete-20-us, convergence or AFSD PASS is claimed.
