# Saved slow-reference event: independent causal readback

2026-10-03. READ_ONLY_SAVED_REFERENCE_EVENT_CAUSAL_AUDIT_NOT_AFSD.
No HHT step, independent ODE step, post-event branch or physical trajectory was
executed by this audit. The reference inspected here is explicitly
`logs/process_slow_20261003/reference_initial`; a later reference-driver
metadata/bracket repair is not silently substituted for that recorded source.
All saved artifact hashes match its report. The imported independent all-pair
formula and literal-input sources match their recorded hashes. The changed
slow-reference driver is NOT imported; its old/current hash difference is
recorded. The initial report's event brackets are inherited observations, not
an independently re-localized event certificate from this readback.

## Method and scope

Reader: `logs/process_slow_20261003/event_audit_initial.py`. It reads the last
accepted coarse/fine NPZ arrays and uses `reference_hht_block.pair_fields`:
independent current all-pair corrected-SPH algebra, not tested HHT/backend
mechanics, material, contact, heat or integration calls. Penalty/Coulomb wall
forces, rigid-wall acceleration, FKS elastic margin, static friction disk and
work identities are evaluated directly in the audit. Neither a rejected
candidate nor a dense root state becomes an accepted physical state.

Command:

```text
.envs/pysph/python.exe -B logs/process_slow_20261003/event_audit_initial.py --reference logs/process_slow_20261003/reference_initial --output logs/process_slow_20261003/event_audit_initial_run/audit.json
```

Exit0; internal0.111281s, external0.316612s, within60s. JSON/logs and exact
command record are in `logs/process_slow_20261003/event_audit_initial_run/`.
Readback algebra bounds were fixed before execution: work/heat identity1e-11J,
momentum1e-12Ns, net internal force1e-10N and spatial power1e-10W. These are
readback/roundoff checks, not relaxed temporal-accuracy or process gates.

## Decisive force evidence

Fine last accepted time is5.974655151367073us. The reference separately
reports full-vector zero-slip approach near5.97466us, face0/PID4. At the saved
state, PID4 has mass8.4375e-5kg and relative tangential velocity
`(0,2.6939e-18,-1.26233461e-7)m/s`: genuinely approaching full-vector zero,
not merely an arbitrary signed-z threshold. Its transverse component is
roundoff-small. PID22 on the opposite face approaches the same event to the
recorded mechanical symmetry accuracy; its speed is1.26233508e-7m/s.

For the fixed normals and disjoint single-face cohorts here, the instantaneous
static traction required in the zero-slip limit is

`f_required = P_t [m*a_wall - f_internal - f_other_contacts]`.

There are no actuators, no other contact on either particle and zero
tangential wall acceleration. This is an instantaneous physical force check,
NOT a finite-step coupled effective-mass or sticking-solver certification.

|Fine saved PID4 quantity|Value|
|---|---:|
|Internal force z|+4.6748433418N|
|Normal penalty force magnitude|13.6961965337N|
|Coulomb capacity mu*Fn|2.7392393067N|
|Current negative-slip friction z|+2.7392393067N|
|Static traction required z|-4.6748433418N|
|Capacity minus required norm|-1.9356040351N|
|Current physical acceleration z|+87870.6092m/s2|
|Algebraic acceleration z if friction reverses at same state|+22940.4923m/s2|

Static sticking is **infeasible at the saved limiting state**, with a1.936N
margin: this is not rounding. Coarse PID4 gives margin-1.9356040352N;
fine/coarse disagreement is below1e-10N. PID22 gives the same negative margin.
Internal tangential force alone exceeds friction capacity and continues
accelerating the material beyond wall speed even after friction changes sign.
Consequently **restart as reverse sliding is strongly supported**, not a
reason to clamp material to the wall. The exact event-state evaluation and
post-event trajectory remain NOT_RUN, so this is a near-event causal inference,
not a claim that a reverse-sliding integrator has passed.

The correct next event implementation must examine both nearly simultaneous
PID4/PID22 roots. Shared timing must be certified from brackets/symmetry, not
chosen by which contact ID happens to be iterated first. The fixture still
does not test shared-particle multi-face or actuator constraints.

## Other causes and balances

The following fine last-accepted checks PASS:

- Minimum overlap1.6188183347um: no observed normal opening.
- Minimum independently evaluated FKS elastic margin313.363378MPa: no yield.
- Correction condition3.38586, density relative change5.8419e-5 and
  temperature400.000025–402.000095K: no observed invalid material/geometry domain.
- Mass0.0027kg and all persistent IDs match the literal input.
- Net internal force maximum2.14e-14N; net torque maximum1.39e-16Nm.
- Momentum change minus accumulated external impulse maximum1.63e-19Ns.
- Paired stress/kinetic power residual0W at the inspected state.
- Mechanical work identity residual1.56e-17J; thermal balance residual-9.21e-13J.

Mechanical identity is
`deltaK + delta_normal_potential + integral(sigma:D*V) - Wwall + Qfriction`.
It uses actual integrated stress power, **not** an assumed hyperelastic
stored-energy formula. It does not certify the old Jaumann closed-cycle
free-energy problem. Thermal identity is
`sum(m*cp*deltaT) - Qreservoir - 0.5*Qfriction`; no residual/algorithmic heat is
introduced. Plastic work remains zero on this elastic path, not a positive
plastic-work validation.

Adding observational metadata and reassigning IDs in the readback
data leaves current independent force, gradient, stress, conductance and
correction matrices bitwise identical. This checks readback equations only;
the separate actual-PA runner label/reorder tests establish the corresponding
initial boundary wiring. Neither test is a post-event checkpoint-restart run.

## Evidence classification and next dependency

- CONFIRMED: saved-state static disk infeasibility, nonzero internal force
  driving the approach, open mechanical free DOFs and unchanged physical arrays.
- SUPPORTED: a genuine negative-to-positive tangential slip reversal at the
  located approach, involving the mechanically symmetric opposite-face point.
- RULED_OUT at the inspected states: opening, yielding, zeroed free velocity,
  label force, or merely0.05m/s guard failure.
- BLOCKED/NOT_RUN: exact accepted event state, event restart, post-event
  coupled trajectory, shared-face/actuator effective mass and1ms accuracy/cost.

The next substantive feature is located stop/reverse-sliding event handling
with right-limit acceleration/work bookkeeping while retaining physical
position, velocity, stress, density and temperature. Additional equivalent
microsteps before the event do not supply that missing branch. Nothing in
this local readback advances the real AFSD checkpoint or establishes transport.

## Actual S0 driver failure is earlier and different

Readback of `outputs/process_slow_20261003/S0_initial` establishes eight
published local stages ending at2us. The next registered interval is
2–2.5us, the first500ns interval. It fails at residual evaluation1 on
face0/PID3 with `contact turning-traction quadrature law defect`. The
checkpoint and all semantic arrays still represent the published2us state;
the2.5us predictor is not accepted. This is distinct from the independent
physical PID4/PID22 stop/reversal approach near5.97466us.

Independent saved-polynomial reader:
`logs/process_slow_20261003/s0_failure_readback.py`; result:
`logs/process_slow_20261003/s0_failure_readback.json`. Exit0, internal0.031216s,
external0.403888s under60s. No integration is called. Its only PySPH operation
is checkpoint loading; force directions and five-point Gauss work are direct
algebra on recorded boundary/slip polynomials. All run/reference artifact
hashes and the last accepted checkpoint fields match. It uses the completed
`reference_repair1` report, not an unverified future reference.

|Rejected predictor check|Independent reconstruction|
|---|---:|
|Maximum sampled Coulomb defect|0.018714892092337464|
|Unchanged acceptance bound|0.01|
|Minimum full tangential slip|0.04600241847m/s|
|Maximum friction power|-0.08971677436W|
|Raw friction heat|5.58822741254e-8J|
|Slip direction turn across interval|15.83514352degrees|

The saved whole-path passivity flag is true, interior-zero list empty and
zero-speed roundoff bound4.98e-15m/s. Independent reference PID3 slip norms
are0.06408515190m/s at2us and0.04890221780m/s at2.5us: ordinary curved sliding,
not an actual stop. No negative friction heat, clipping, wall jump or
physical zero-slip event explains this first driver rejection.

The actual call order matters: `advance` initially sets `a1=a0`, calls
`evaluate`, then `_contact_work` rejects the interpolated traction before
the initial residual returns. Consequently no Newton correction was attempted.
This **confirms an infeasible initial predictor under the acceptance guard**;
it does not prove that a converged500ns coupled root cannot meet that guard.
Nor does it prove that allowing root iteration would pass: that remains
NOT_RUN and must not be inferred from the reference's smooth solution.

A minimal next numerical test may distinguish provisional residual-evaluation
traction defect from final accepted-path law acceptance. The0.01 final bound,
whole-path event/passivity checks, raw heat law and finite-domain checks must
remain unchanged. Independently demonstrate that a defective converged
candidate still rejects; merely suppressing the exception is not a valid fix.
No such patch or rerun was made in this read-only review. The later physical
stop/reverse-sliding branch remains a separate required feature even if an
initial-predictor handling repair succeeds.

### Subsequent private diagnostic rules out a sufficient predictor-only fix

The main agent subsequently preregistered and executed exactly one private
500ns root diagnostic from the same stage008. I reviewed its saved contract,
source and report, not a rerun:
`docs/PROCESS_SLOW_PREDICTOR_DIAGNOSTIC.md`,
`logs/process_slow_20261003/diagnose_predictor_veto.py`, and
`logs/process_slow_20261003/predictor_diagnostic_initial/report.json`.

Only the exploratory sampled-direction veto was deferred inside that private
process. The original unwrapped contact-work function was reapplied with the
unchanged0.01 bound after the root solve. Source and parent artifacts remained
unchanged; the private candidate is non-restartable and no accepted journal
entry was added. The diagnostic is not a physical trajectory continuation.

The12-evaluation root converged: velocity residual6.176e-10m/s, density
relative residual3.369e-16, temperature residual0K. But its maximum final
direction defect is0.01252632750448, still above0.01. The original final guard
rejects contact0:3. Full-path passivity remains true; minimum slip0.04920452m/s
and rawQ5.78006247e-8J remain positive. Outcome is therefore
**FINAL_ORIGINAL_GUARD_FAIL**, not a green repair or an event.

This supersedes the earlier uncertainty about whether merely reaching the
root could remove the defect: premature predictor rejection is real, but it
is **not a sufficient explanation/fix** for this registered500ns interval.
The converged endpoint-linear traction representation itself fails the
original sampled Coulomb law accuracy gate on a rotating slip direction.
Do not implement only a delayed exception and call the issue solved. A next
variant needs contact-path traction accurate to its actual curved slip, with
consistent impulse/work/heat and independent quadrature/direction benchmarks,
while preserving the old failure and final acceptance requirements. The later
5.97466us physical stop/reversal is still a different missing event branch.
