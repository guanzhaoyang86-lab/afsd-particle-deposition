# Curved-slip contact coupled consistently: scoped result

2026-10-03. ENGINEERING_APPROXIMATE development, not project completion.
PySPH1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5 remains the only backend.
The real AFSD accepted state remains850.393632910904ns, crossings0. New AFSD
steps0. All outputs here are a27particle,81free-DOF manufactured3D elastic
startup block, not a feedstock/substrate deposition result.

## What changed and why it is not the previous failed fix

The original saved[2,2.5]us interval failed an endpoint-linear friction
interpolation check. Merely postponing that check had already converged to
another inadmissible root. The new named numerical variant evaluates actual
Coulomb traction along the same candidate HHT dense position/velocity path,
integrating normal/tangent impulse, material/wall work and signed friction heat.

New files `hht_curved_contact.py`, `hht_curved_identity.py` and
`hht_curved_coupled.py` reuse old PA/NNPS, material, corrected operator, bounded
Schur iteration and BE thermal facilities. The advance body retains the old
material/root workflow with the contact impulse correction; it is not a new
backend. The original HHT code, thresholds, files and FAILs are preserved.

With dJ and dW denoting actual minus endpoint-linear contact impulse/work,
the new mechanical residual includes dJ/(gamma*h). The signed algorithmic work
includes vbar dot dJ - dW. Raw Q alone enters the existing heat partition/BE
solver; neither correction nor root residual becomes heat. F1 is still the
true physical endpoint force, A1 the numerical canonical acceleration. This
is not a claim of inherited classical-HHT stability; see independent derivation.

Whole-vector zero-slip paths still reject. Whole-path contact gap and wall
tangential kinematics are source-bound certificates, not quadrature-only
guesses. Direct exact Coulomb force is used; a callback inconsistent even by
0.5% rejects rather than being silently overwritten. Gauss8/16, then32 if
needed, controls an explicit quadrature-error ESTIMATE (not a universal bound).
No negative heat, temperature, yield, mass/time/modulus or free-velocity clipping.

## Actually run

|Check/run|Result|Scope/cost|
|---|---|---|
|Prior freeze check|PASS4221|Integrity only|
|Independent impulse/work identity suite|7 PASS|Nonzero history, a=1, exact old reduction, Galilean and analytic work|
|Independent contact suite|8 PASS|Analytic rotating-vector and full HHT turning-slip integrals, zero/sign-only cases, bad-force rejection|
|Coupled transaction suite|6 PASS|Actual3D short step, unchanged labels, heat, migration, new checkpoint hashes and retained failure details|
|Complete selected preflight|189 PASS|168 old +21 new;11.5068s|
|Saved RED regression|PASS_PRIVATE,not published|Original500ns FAIL reproduced; new root12evaluations,external1.1159s|
|Cost pilot|PASS|One250ns step;measured trajectory.41663s|
|Fixed original S0|STOPPED_UNREFERENCED|13published steps to5us;external9.29646s|
|Published4–5us checkpoint replay|PASS_BITWISE|Exactly one recomputation;40properties/full state/interval/ledgers equal;external1.4764s|
|S1/S2,20us convergence,1ms,AFSD production|NOT_RUN|No completed common process/reference horizon|

Analytic contact impulse/heat agree near1e-16 in the independent SI benchmarks.
The saved-RED root passes unchanged raw velocity residual2.78861e-10m/s,
actual Coulomb defect0 and kinetic identity residual7.15e-22J. No physical
array or original checkpoint was changed by its explicit numerical migration.
That private regression is not an accepted trajectory or restartable state.

Independent saved-output audit verifies all13NPZ/VTK and14same-time frames.
A separate literal path plus64point quadrature reconstruction agrees with
saved contact impulse/work/heat within1.96e-18 in their SI units; maximum saved
temperature energy balance8.49e-15J. The one published-interval replay preserves
every physical property, canonical force/acceleration and complete interval
history bitwise, with no new published state. Neither check proves convergence
or restart through an event. See OUTPUT_REVIEW and REPLAY_REPORT.

Source review caught and fixed diagnostics, ledger and gate wiring before the
final preflight: preserve contact exceptions/candidates, bind interval dJ/dW/
quadrature into the new checkpoint hash, require saved-RED PASS before pilot,
and remove unregistered S1/S2 launch modes. No empirical physical repair or
tolerance adjustment was needed after the first coupled trajectory attempt.

## Actual S0 outcome and first stopping condition

The physical configuration hash equals the previous slow fixture. Fixed steps
remain250ns to2us,500ns to4us,then1us. The run reaches5us in13steps; laststep
15root evaluations,.769345s root cost. This is not a sequence obtained by
shrinking dt until it passes. All accepted contact quadratures use actual
forces, and the5us step's normal potential/work residual is2.42e-17J.

The candidate[5,6]us also converges, but6us exceeds the independent oracle's
certified common prefix5.9375us. The driver therefore saves it as
REJECTED_UNREFERENCED, never publishes stage014 and never authorizes restart.
This is a reference-domain stop, NOT negative heat, a root failure or a solved
physical zero-slip event. The separate oracle approaches genuine stop near
5.97466us; the rejected coarse numerical candidate has PID4/22 slip still
negative at about-.0226522m/s. Its path has not reproduced that event yet.

Same-time5us coarse errors against the oracle are RMS position2.1158e-8m,
RMS velocity.0132494m/s, full-stress RMS408995Pa and maximumT4.96245e-7K.
These are real temporal discrepancies, not algebraic residual noise. They
are reported raw, not assessed as a completed finest-grid sweep or convergence.
More details and independently recomputed ledgers: PROCESS_CURVED_OUTPUT_REVIEW
and PROCESS_CURVED_CAUSAL_REVIEW. No paper trend/area-ratio comparison is claimed.

## Cost and scope limits

The launch forecast for28steps was33.3313s, below50s; actual stopped S0 took
8.71057s internally,9.29646s externally,peak working set81,080,320bytes.
The full run includes13accepted intervals plus one rejected/unreferenced
candidate and output. Root costs are not mislabeled total step costs.

Naively maintaining these local startup costs would make1ms roughly769–1742s,
already above600s. This is a warning from27nodes, NOT a production forecast;
neither process-scale step growth nor full AFSD economics is established.
No long run, GPU, parallelism, mass scaling or paid computation was launched.

## Next mandatory work

The earlier contact-chord numerical obstruction is repaired for the saved
case and fixed startup prefix. That does not settle event timing or the later
physical stop/reversal. Extend an independent event/reference solution first,
then implement located stop, static feasibility and reversed sliding with
consistent right-limit acceleration/history. Compare only equal physical
times; do not promote the6us rejected state. Keep inertia and the existing
physical boundary unchanged. See PROCESS_CURVED_NEXT_WORK_PACKAGE.

Positive plastic coupling, shared-face/actuator response, actual finite open
square-bore/shoulder/inlet integration, material transport and affordable
process-scale progression remain incomplete. A local prefix is not their proxy.
Commands, failed/private/accepted records and new freeze remain inspectable;
the freeze proves preservation, not production/scientific acceptance.
