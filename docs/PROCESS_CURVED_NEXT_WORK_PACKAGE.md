# Next: physical stop/reversal reference and coupled event handling

2026-10-03. Current sources and evidence: PROCESS_CURVED_REPORT, DERIVATION,
OUTPUT_REVIEW, CAUSAL_REVIEW and REPLAY_REPORT. All old FAILs remain.

Do not retry the old chord interpolation or tune its1% bound. Its saved RED
now passes using actual curved Coulomb force and impulse-consistent mechanics.
The new accepted prefix reaches5us, not the oracle's true stop5.97466us.
The rejected6us numerical candidate is still slipping the old way; do not
claim an event merely because the reference stopped at a similar time.

## One substantive dependency

1. Extend the independent slow reference with a located full-vector-zero event
   and actual static-force feasibility. PID4/PID22 are nearly simultaneous;
   determine joint event timing from bracket/symmetry, not enumeration order.
   Saved near-event internal tangential force4.675N exceeds capacity2.739N,
   supporting reversal rather than sticking. Confirm at the actual event state.
   Preserve last accepted/pre-event/rejected states and exact right-limit data.
2. Before integrating the3D fixture after that event, independently verify the
   discontinuous branch with the original analytic reversal case:
   m1kg,u0+1m/s,Fother-2N,muFn1N,total1s; tstop1/3s,ufinal-2/3m/s,
   Q7/18J,Wother1/9J,deltaK-5/18J. No abs/max clipping or automatic failure on
   genuine reversal. Add feasible-stick and nonzero-transverse component-turn
   cases so scalar roots cannot masquerade as full-vector zero.
3. Add located event splitting/restart to the named curved HHT variant, keeping
   x/v/S/rho/T continuous and recording the right-limit force/acceleration and
   signed numerical storage/work jump. The general shared-face/actuator case
   needs coupled mobility/reaction; a disjoint-plane force test is not enough.
   Do not blindly carry left acceleration across the force discontinuity,
   silently project velocities, or integrate a smooth polynomial across a
   static interval it cannot represent.
4. Only after independent branch tests, use a registered same-time prefix or
   event-completed horizon to test temporal accuracy/cost. The current coarse
   errors and delayed stopping are not just root tolerance noise. Keep all
   unsuccessful paths; no different-failure-time convergence comparison.

Interfaces likely touched: new event helper alongside
`hht_curved_contact.integrate` (currently rejects full zero),
`hht_curved_coupled.advance/initialize` (new event history definition),
`reference_hht_slow` in a NEW preserved variant, and a source-bound driver.
Re-use the existing exactCoulomb/zero-path/impulse/work/heat helpers and PySPH
storage; no backend comparison, material rewrite or complex frontend.

Register physical event tolerances and evaluation/wall budgets BEFORE runs.
Single actual run<=600s; prefer short tests. Do not launch a1ms reference or
process run while its measured forecast exceeds budget. Need fewer accurately
resolved process steps or another explicitly justified numerical strategy,
not undisclosed scaling/removal of inertia. If changed physical friction,
compliance, quasi-static assumptions or long compute is required, ask direction.

Positive plasticity and actual AFSD finite-tool/inlet/transport remain later
mandatory integration; do not redefine completion around the elastic fixture.
