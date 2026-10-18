# Stop-scout output audit and raw same-time errors

2026-10-03. Independent read-only review of
`outputs/process_stop_20261003/scout_initial`. Result:
**PASS_OUTPUT_INTEGRITY_RAW_SAME_TIME_ERRORS_NOT_EVENT_OR_PROCESS_ACCURACY**.
No solver/ODE advancement, 3D event solve, AFSD step, production-source change
or existing-evidence edit was performed by this audit.

## Registration and evidence

`logs/process_stop_20261003/output_readback_registration.json` preregistered
exact stored-field/hash equality and an integral/identity reconstruction
bound of **1e-12 absolute SI + 1e-8 relative** before execution. Independent
reference comparisons have no accuracy PASS threshold and are not refinement
or convergence tests.

The new script `logs/process_stop_20261003/output_readback.py` completed on its
initial run: exit 0, 2.7864601 s internal / 3.1869995 s external, under a
60 s external ceiling. No repair was needed. That directory contains the
full JSON report, initial stdout/stderr, supervision JSON and
`OUTPUT_READBACK_COMMANDS.md`.

All **4,365** files in `PROCESS_CURVED_FREEZE.json` were unchanged before and
after the audit. The three authorized index-document exceptions had not
changed at audit time either; later authorized index updates are outside
this historical check.

## Publication and receipt integrity

The scout is `STOPPED_LOCAL_EVENT_SCOUT_NOT_PROCESS`. Its 90 source hashes
and 19 artifact hashes match. The same-source **233-test PASS** preflight's
90 sources and 27 artifacts match. The independent `reference_endpoint2`
receipt's six sources and 13 artifacts match, including the same physical
configuration and saved reference history. Parent report/checkpoint/state
hashes bind migration to the previously published 5 us state.

Exactly two NPZ/VTK stages and two frames are published:

| Record | Time | Meaning |
|---|---:|---|
| `stage_000_transition` | 5 us | Same physical fields and numerical histories, new event-capable schema |
| `stage_002` | 6 us | One actual accepted 5–6 us interval |

Both use `HHT_LOCATED_STOP_ELASTIC_BE_STATE_V1`. Their semantic hashes,
particle hashes, parameter hashes, parent chain and publication receipts
agree. All **40 particle properties**, literal unscaled masses, persistent
IDs, passive labels, saved x/v/S/rho/T, acceleration and force histories were
checked. NPZ/frame/VTK fields match exactly, including all six full-stress
components, pressure, equivalent stress and vector velocities/accelerations.

Migration preserves every particle property, canonical acceleration, actual
force, boundary snapshot, previous interval record and cumulative integral.
It has zero physical-energy jump; it is not an additional time step. The
accepted 6 us state's last interval is bound to the strict contact diagnostic
and the 5 us parent. There is **no published HHT stop/event-left/right restart**,
and no `event.json` receipt. Last published time remains 6 us.

## Independent impulse/work/heat reconstruction

The audit reuses only frozen pure readers and independent formulas from the
previous output reader, not the tested contact or kinetic-identity helpers.
From saved histories it reconstructs the compatible cubic position/quadratic
velocity, literal 5 Hz moving-plane geometry, 7e6 N/m penalty force and
coefficient-0.2 vector Coulomb law. Saved quadrature nodes/weights and each
force, impulse, material work, wall work and raw Q contribution are checked.
A separate **64-point Gauss** calculation cross-checks the saved 16-point rule.

Largest reconstructed differences/residuals for the accepted interval:

| Check | Maximum absolute value |
|---|---:|
| Independent Gauss64 versus saved integral component | 1.91599e-18 N s/J |
| Per-particle momentum identity | 4.43369e-21 N s |
| Per-particle kinetic/work identity | 7.10929e-22 J |
| Saved thermal-energy balance | 6.98534e-15 J |
| Cumulative face impulse | 1.08420e-19 N s |

The interval's total raw friction heat is **4.726110714482718e-6 J**. Mechanical
raw Q exactly supplies the thermal raw-Q ledger, with the declared one-half
partition. Reservoir exchange is independently rebuilt from the endpoint
temperature and 0.025 W/K per contacted particle. Internal conductive heat
cancels globally. DeltaJ/DeltaW, actual material work and numerical filter/root
work remain separate; no identity residual becomes heat. This checks the
saved thermal source/balance, not an independent SPH graph reassembly.

Independent full-vector polynomial minimization gives minimum slip
**0.0226522082 m/s** over the accepted 5–6 us HHT path. That certifies the
nonzero branch of its own numerical path, not agreement with the physical
reference's event time.

## Raw comparison at exactly 6 us

The independent `reference_endpoint2` history includes exactly 6 us; no time
interpolation is used. It has located two sequential full-vector-zero events,
PID 4 at 5.974656587951323 us and PID 22 at 5.974656587951857 us, with static
infeasibility and reverse sliding. Its certified common prefix ends at 6.5 us;
the complete 8 us reference is still not available. The JSON also records the
matching 5 us migration-frame errors.

| Error at 6 us | Raw value |
|---|---:|
| Position volume-weighted RMS, m | 1.95202326e-8 |
| Velocity volume-weighted RMS, m/s | 4.18508307e-2 |
| Full stress volume-weighted RMS, Pa | 3.30632460e5 |
| Deviatoric stress volume-weighted RMS, Pa | 2.96691453e5 |
| Density maximum, kg/m3 | 1.90998852e-2 |
| Temperature maximum, K | 3.29998443e-7 |
| Face-force maximum, N | 6.93016430 |
| Face-impulse maximum, N s | 8.56263086e-7 |
| Wall-work maximum, J | 1.71252617e-7 |
| Raw-friction maximum, J | 1.21072094e-7 |
| Reservoir-heat maximum, J | 3.88116460e-12 |

The branch-timing discrepancy is explicit: at 6 us HHT PID 4 has tangential
z-slip **-0.0226522082 m/s**, while the reference has **+0.000559327090 m/s**.
PID 22 gives the same opposing signs. HHT remains incoming while the reference
has reversed. Publication integrity and consistent heat/work do not eliminate
this temporal error. No accuracy, convergence or process-feasibility PASS is
inferred from these raw comparisons.

## PRIVATE/rejected attempts cannot continue the published trajectory

The private NPZ semantic hashes and particle snapshots match their diagnostic
receipts. Explicit read-only restore attempts reject trials 001, 003 and 004
with `private/rejected candidate is not a restartable state`; no advance was
called. Trial 001's physical fields equal the later certified 6 us state, but
the original private file itself remains unaccepted/non-restartable.

| Attempt | Endpoint / target | Saved outcome |
|---|---:|---|
| 002 | 6.5 us | Initial full-vector-zero rejection, no NPZ |
| 003 | 6.3016715285 us | PRIVATE incoming continuation, restore rejected |
| 004 | 6.2755694698 us | PRIVATE incoming continuation, restore rejected |
| 005 | 6.2725476233 us | Bounded line-search failure, no NPZ |

Trials 003/004 retain six/two negative local raw-Q node contributions, with
minima -3.74638e-11/-1.32802e-12 J respectively. They explicitly lack physical
contact-path certification, are absent from accepted history, and were not
compared as accepted physical states. Positive integrated heat cannot promote
them. Trial 005 stopped after 18 evaluations; its last velocity residual is
3.83375e-6 m/s versus the unchanged 1e-9 m/s limit. Its rejected JSON is bound
to the still-published 6 us parent and has no continuation checkpoint.

Therefore this run advances one accepted interval to 6 us, but does **not**
complete a coupled HHT stop/reversal. The independent reference's successful
branch handling cannot be transferred to the failed HHT event search. The
8 us/20 us/1 ms HHT horizons, temporal convergence and real AFSD remain NOT_RUN
or unresolved as recorded by their respective contracts.
