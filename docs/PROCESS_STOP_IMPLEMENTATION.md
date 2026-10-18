# Located stop implementation scope

2026-10-03. This is a local27-node elastic3D benchmark, NOT AFSD transport.
The old4365-file inventory, all old failures, material and paired corrected
spatial operator are retained. No GPU, new material law or production job.

## Code path and provenance

`hht_stop_stage.advance_private` is an explicit variant copied from the frozen
`hht_curved_coupled.advance`, retaining its24-evaluation nonlinear root, current
NNPS, midpoint Jaumann/FKS elastic branch, implicit density and BE temperature.
Only the contact integration dispatch/state lifecycle differs. It returns a
private state with accepted=false and restart_allowed=false, even if its
mechanical/density/thermal residual gates pass. A private incoming continuation
is a root-search equation, NOT physical reversal or accepted heating.

`hht_stop_contact.integrate` retains refined8/16/32 Gauss integration of actual
forces along the same HHT dense path. A selected naturally-collinear incoming
branch can continue mathematically beyond zero for private root search; signed
negative node heat remains visible. Negative integrated contact heat rejects.
Ordinary contacts retain actual full-vector Coulomb. Registered endpoint modes
must have complete vector zero at roundoff, the correct actual polynomial
one-sided direction, no resolved earlier zero and the whole-path gap certificate.
There is no velocity snap, traction freeze for a face or residual-to-heat transfer.

`hht_stop_flow.locate_duration` uses its own coupled candidates to solve duration.
Each trial completely resolves the mechanical/density/temperature state from
the same accepted parent. Its initial predictor comes from an actual rejected
dense path, never the independent-reference event time. Maximum12 duration
trials, each maximum24 residual evaluations. The coupled certificate is raw
full-vector kinematic zero plus the unchanged coupled residual gates. It does
NOT claim the independent oracle's5e-19s bracket accuracy or continuous event
accuracy: its final bracket and later raw time error must be reported separately.
The initial coupled driver still supports only a jointly certified event group;
separated events must reject. Sequential oracle support does not silently grant
that capability to the coupled solver.

`certify_interval` re-evaluates under the strict endpoint rules and requires
the solved actual impulse/work/heat arrays to remain exactly equal before a
candidate can become a local accepted state. Publication additionally requires
the physical domain, output readback, source/parent/reference receipts and
certified reference horizon in `probe_hht_stop.py`.

`restart_sliding` evaluates the actual event-state static traction requirement.
The free disjoint contact can restart sliding only with an unambiguous infeasible
static disk. General3D sticking, shared-contact mobility and actuator reactions
remain unsupported. It preserves every particle array and records old canonical
HHT acceleration separately from instantaneous force/m. Right acceleration is
explicitly initialized from actual right-limit total force. Physical impulse,
work and heat of this zero-duration operation are zero. The signed acceleration
history-storage change is separately retained and is NOT material heat or proof
of the full nonlinear algorithmic energy.

## Publication and budget

Only the frozen curved S0 accepted_history's13th checkpoint at5us is a valid
parent. The previously rejected6us candidate is never read as a starting state.
The explicit same-fields numerical variant transition is itself checkpointed;
all later parent hashes remain traceable. Zero-duration left/right event records
share physical time and are marked separately from integrated intervals.

At most17 complete coupled solves in this single scout, including failed
predictors and duration trials. The forecast uses the measured previous
15-evaluation root cost, scales to24, doubles it and adds4s. This is a conservative
estimate, not measured complete I/O cost. Forecast must be<=50s; internal55s,
external60s via `run_process_stop_command.py`. Both are below user's600s ceiling.
No1ms or process run is authorized. See the result report for what actually ran.

Before any coupled scout, the independent extension encountered a later real
zero near6.565us in its coarse path. A certified common prefix>=6.5us may thus
be used instead of an8us horizon. Targets retain6/7/8us where available and add
the certified final prefix boundary. This is a declared horizon truncation,
not a modified physical wall/load or a convergence claim. Reference event time
is still never used to select the coupled stop time. If that common prefix is
not independently certified, the actual coupled scout remains NOT_RUN.

## Independent references and retained risks

The independent reference imports only literal inputs and the separate
all-pair spatial formulas, never these contact/material/root implementations.
Its initial simultaneous-event certificate failed; original sources and complete
rejected candidates are preserved. The separately preregistered sequential
oracle change is in PROCESS_STOP_SEQUENTIAL_ADDENDUM; it changes event ordering,
not physical parameters or numerical bounds.

An elastic fixture does not solve positive-plastic coupling, finite open-tool
entry, shared faces/actuators, long-time temporal accuracy or actual material
transport. General thermodynamic material consistency remains an open risk.
Current AFSD accepted time is still850.393632910904ns and transport crossings0.
