# Same-branch material response: minimal numerical integration design

2026-10-04. DESIGN, not implemented/validated integration. This review changes
no physical source, builds no matrix and calls no material update or trial.
The actual old driver is `scripts/run_process_global.py`; there is no current
`scripts/probe_process_global.py` file. Old source and FAIL evidence stay frozen.

## Recommendation and causal limit

Add one named **candidate-refreshed, fixed-midpoint, rate-aware material plus
pressure/continuity Schur proposal**. Keep the exact existing
ProcessThermalTrial, FKS update, contact, actuator, BE sources and raw residual
controller. This is a nonsymmetric approximate search matrix, not a new
constitutive law or a full process Jacobian.

First validate its action on the already saved 650 ns positive-plastic parent /
50 ns bridge. Do **not** presume this will fix the older 500 ns failure:
that run's final candidate was elastic at every particle and its remaining
velocity defect was concentrated at stop-slip contacts. These are different
candidate states and different causal evidence. The new material block is
supported for the bridge, not established as the cause of the 500 ns tail.
The 500 ns rerun remains NOT_RUN; no contact extension is part of this design.
The agreed scope of this package ends with saved-journal matrix construction
and action: **zero new physical evaluations or roots**. The optional two-call
integration below is a future separately registered task, not current authority.

## Fixed-midpoint algebra; retain the factor one half

Use the already returned candidate midpoint `i,j,g,B,mass,rho,L,T` journal,
not a new NNPS search. Construct its identical linear gradient map D such that
`vec(L)=D*vbar`, volume V=m/rhomid and E=trace(D). Check this reproduction
against the saved L. All particle/axis ordering follows persistent IDs and
the actual arrays; labels never enter the assembly.

The new component supplies `C_alg = d(vec(S1))/d(vec(L))` [Pa s] at fixed
S0, ep0, Tmid, dt, on its certified fixed material/rate branch. It must include
the rate derivative required by the independent material tests. Its stated
spin/prestress scope must be preserved: if their derivatives are not supplied,
call this a fixed-spin/prestress approximation, not an exact update tangent.
The physical update still computes its original rotation and 48 return steps.

Power-adjoint assembly at this frozen midpoint gives

```text
F = -D^T V vec(sigma)
sigma = (S0+S1)/2 - pressure(rhomid) I
J_v = dF/dvbar|rho,T,geometry = -1/2 D^T V C_alg D       [kg/s]
J_r = dF/drho1|pressure only = 1/2 E^T diag(V K/rho0)
J_T = dF/dT1|material only = -1/4 D^T V (dS1/dTmid)
kappa = -dt*rho_parent/(1+dt*trace(L)/2)^2
A = M - dt/2 J_v - dt/2 J_r diag(kappa) E               [kg]
```

`C_alg` already contains the material update's dt. Do not multiply it by an
extra dt or add the old full elastic stiffness on top of it. Temperature has
**two** factors one half: stress midpoint and candidate temperature midpoint.
If temperature response is not included, record `J_T_OMITTED` rather than
silently supplying a zero certified derivative. Its branch certificate is
required when included. In the unstressed, zero-spin elastic limit
`C_alg=dt*C_dev`, with div(L)=0 and rho_parent=rho0, A reduces exactly to
`M+dt^2*(K_dev+K_bulk)/4`, where `K_bulk=E^T diag(V K) E`.

Use the same pressure closure and exact **candidate** midpoint-continuity
derivative in a newly named numerical Schur block:

```text
rhs = M*dbar + dt/2 J_r*drho + dt/2 J_T*dT
delta_bar_c = dbar_c
A_ff delta_bar_f = rhs_f - A_fc delta_bar_c
delta_rho = drho + diag(kappa)*E*delta_bar
delta_x = dx + dt*(delta_bar-dbar)
delta_v = dv + 2*(delta_bar-dbar)
delta_T = dT
```

Here d=actual_map(C)-C, rho0=p[rho0]. Pressure derivative with respect to
candidate rho1 contains its separate midpoint factor 1/2. No density or
pressure formula in the physical map is changed. This deliberately replaces
only the old proposal's reference-density coefficient `-dt*rho0` with kappa;
it is not described as the unchanged old preconditioner. Its general block
derivation is `(M-dt*J_v/2) delta_bar -dt*J_r*delta_rho/2
-dt*J_T*delta_T/2=M*dbar`, `delta_rho-kappa*E*delta_bar=drho`,
`delta_T=dT`. The eliminated matrix and RHS above follow directly.
J_r omits the other density effects through V, D, B and material geometry;
J_T omits BE feedback. These omissions are explicit, not physical changes.
On constrained axes use `delta_bar_c=dbar_c` and retain their original columns
in A_fc and E: density back-substitution includes the full solved delta_bar.
Do not project away an actuator column before constructing the Schur block.
Likewise the factor two endpoint lift is not asserted for real stop events;
the independent event-average field and original full map remain authoritative.

Omitted from this limited matrix: derivative of changing geometry/volumes/B,
contact/event and reservoir response and BE temperature feedback. No claim of full
Newton, symmetry, positive definiteness or process accuracy follows.

## Smallest interface, without modifying frozen root logic

`ConstrainedResponse` rejects nonsymmetric stiffness and uses a positive-
definite elastic A. `StagePreconditioner._install` also requires that exact
class type. The new rate/prestress response need not satisfy either premise.
Do not symmetrize it, forge a ConstrainedResponse instance, clip eigenvalues,
add diagonal stiffness or use a pseudoinverse to force it into that interface.

Add a separate small `GenericConstrainedResponse(A, mass_per_dof, dt, free)`
with detached finite arrays and actual positive mass. Solve the free block
using `numpy.linalg.solve`; report the unchanged full algebraic reaction
`A*delta-rhs`, its free residual, nonsymmetry and conditioning. Require the
existing scaled linear backward error <=1e-12. Singular/nonfinite failure
stops the candidate. No arbitrary condition cutoff or hidden matrix scaling
is proposed; this matrix condition is not the SPH moment condition <=1e6.
These reactions remain search diagnostics, never contact impulse, work or heat.

Add a proposal class with the same external `.apply(candidate, mapping)`
result fields, but its own explicit installation path. It must not subclass
the old class merely to bypass its validated symmetric-PSD assumptions.

The old root order is `trial.evaluate -> raw residual -> pre.apply -> journal`.
Its `on_evaluation` callback is **after** `apply`, so it cannot refresh the
current candidate matrix. Recommended binding, avoiding a second physical call:

1. A thin `CapturingTrial` delegates `.candidate()` unchanged. Its
   `.evaluate(C,dt)` clears stale context, calls the existing physical trial
   exactly once, then binds that completed evaluation to the proposal provider.
2. The provider's `.apply(C,map)` checks the bound parent hash, full five-field
   candidate hash, map hash, dt, parameter/source signature and ID inventory.
   It builds the new matrix lazily from that evaluation's saved midpoint and
   material state, only when the old root actually requests a proposal.
3. Each subsequent physical evaluation invalidates the previous context.
   A raw-growth-rejected or already-converged candidate builds no matrix.
   Failed/stale/cross-candidate context is a hard error, not a reused P.
4. Reuse `stage_structural_root.solve_stage` unchanged for its raw residuals,
   Anderson beta .75/history6, growth test, evaluation cap and failure journal.
   A named new driver/delegator records the new proposal classification while
   retaining the old controller method as a separate provenance field; do not
   misreport its matrix as `FROZEN_ELASTIC`.

No cache is inserted into the physical particles or checkpoint. The wrapper
binds data, not an alternative force or candidate map. No `evaluate` call is
allowed inside `.apply`. Return full branch margins and component scope, source
hashes and timing; no hardcoded PID269, expected force direction or tracer use.

## Required independent tests before any integration call — NOT_RUN here

- Existing independently owned material-point tests: rate/strain/temperature
  directions, fixed elastic/plastic and below/above-floor branches; explicit
  yield/floor switch rejection. The continuous implicit-return derivative is
  distinguished from the finite 48-bisection map uncertainty.
- Nonsymmetric direct solve: literal A=[[2,1],[0,3]], rhs=(3,3) -> (1,1).
  Prescribing second increment .5 gives first 1.25 and second reaction -1.5.
  Embed in the real three-axis inventory; retain nonzero cross-block coupling.
  Singular, NaN, prescribed-on-free and stale-context cases fail explicitly.
- Fully coupled manufactured two-particle block (six physical velocity DOFs,
  with each particle's actual mass repeated on its axes): choose selected x
  masses (2,3), dt=2, J_v=[[-1,1],[0,-2]], one active divergence row
  E=(1,-1), kappa=-1, J_r=(1,-1)^T, J_T=(1,-1)^T. Other velocity and density
  directions are direct/uncoupled. J_r is proportional to E^T, so the pressure
  sign remains power-adjoint; this is a manufactured linear block, not SPH.
  With dbar=(1,2), drho=1/2, dT=1/4, literal Schur A=[[4,-2],[-1,6]],
  rhs=(11/4,21/4), delta_bar=(27/22,95/88), delta_rho=31/88.
  Constraining the second selected increment to its residual2 gives
  delta_bar=(27/16,2), delta_rho=13/16 and algebraic second reaction81/16.
  Test these fractions independently, including the original uneliminated
  block equations; a missing Tmid factor, omitted constrained column or wrong
  continuity sign then fails. No tested assembler generates the reference.
- Elastic-limit and Schur signs: recover the old literal pressure benchmark
  delta_bar=.32 m/s, delta_rho=-1.2 kg/m3. A separate two-particle analytic
  material force action catches missing/doubled midpoint stress factors.
- Nonparallel nonzero initial stress versus imposed deviatoric rate: include
  an independently derived nonsymmetric material response, not only coaxial
  loading where accidental symmetry can hide the missing term.
- Zero residual yields zero proposal; a zero/tiny proposal cannot accept a
  nonzero raw residual. Event-average velocity remains independent. The matrix
  reaction is not entered into physical work/heat.
- Fake delegated trial count and exact map equality; context hashes reject
  stale C/map/parent/dt. Sorting transforms vectors/matrices consistently;
  label-only changes do not affect values. No accepted history mutation.

Then one offline action comparison using saved base/plus/minus, their fixed D
and epsilon, unchanged physical journals. Compare `delta_C-P_new*delta_R`
against the old two P actions, reporting all fields and support, not choosing
only an improved norm. Their contact branch differs, so this is a finite
signed action measurement, never a certified derivative. No real evaluation
is needed for this first integration check.

## Budget and actual-run boundary

Old evidence has 32 evaluations, not eight. The stored evaluation+proposal
timings exclude output callback IO: mean5.9658480 s, sum190.9071365 s; the
complete old run was265.0000492 s internal /265.5593475 s external. Old probe
build .7774535 s and complete call7.4581122 s produced forecast389.1555658 s.
Those costs do not automatically bind a new per-candidate matrix refresh.
No 500 ns retry is recommended by the new material evidence alone.

If component and offline action checks pass, register **at most two** actual
50 ns evaluations from the same saved650.393633ns positive-plastic parent:
first the fixed bridge base C; second the single original-beta numerical
proposal `C1=C+.75*P_new*(F(C)-C)`. No line search, epsilon sweep, restart from
REJECTED500ns, dt reduction or third evaluation. Report raw residual change,
actual branches and independent ledgers even on deterioration. This probes
integration of a search proposal; neither result is automatically accepted.

Before those calls, measure the new P build/apply offline from saved journals
(no physical solve). With cost b and the bridge's largest complete call
9.2744745 s, proposed launch forecast is `1.5*2*(9.2744745+b)+5 <=50 s`,
internal55/external60 s. If not affordable, mark NOT_RUN and propose a new
explicit budget <=600 s before execution; do not silently extend the timer.
After the first full evaluation+audit+archive, update the estimate and skip
the second if the registered budget cannot hold. Full32-step root and500ns
accuracy remain separately NOT_RUN; material action improvement is insufficient
to authorize a blind repeat of the265s contact-dominated failure.

## Registered bounded read-only check of the old500ns evidence

Before executing the new reader: one invocation, internal55/external60 s.
Bind the old manifest/source/artifact hashes, original stage004 parent NPZ,
latest trial and final candidate032. Reconstruct plastic increment from saved
ep1-ep0, infer qtrial=qout+3G*dep, and evaluate the FKS formula literally on
saved ep0/rate/Tmid, without calling `material.update` or `yield_stress`.
Report its yield and rate-floor margins as readback values, not a new physical
test or derivative oracle. The identity-based inferred qtrial may carry
ordinary floating-point error; no new acceptance threshold is introduced.
Also summarize the32 existing costs and final raw residual. Execution status
before this reader: NOT_RUN; results will be appended below, no old evidence
overwritten.

### Executed old500ns readback

Command `.envs/pysph/python.exe -B logs/material_response_20261004/read_saved500_supervisor.py`:
exit0; external1.1334166 s, internal0.6584298 s. Report:
`logs/material_response_20261004/saved500_readback_initial.json`; log and
supervision receipt are adjacent. All95 bound files (55 sources,39 artifacts
and the parent checkpoint) matched before/after. No matrix/update/trial ran.

The final candidate hash is
`4e71a8f7639adb4a4050a8581045dcb0b773b5611b55843f07a34d1b825ee7b4`,
identical in latest_trial, candidate032 and the REJECTED archive; parent NPZ
hash `18d97b64c5bcb960e5a12a7491a849287811551ae89c67ae0fdc0583507966b8`.
All342 dep are exactly zero; saved plastic work sum is exactly0 J. The
independently evaluated smallest yield margin is2915884.21525942 Pa at PID341;
minimum rate-floor distance is1.2947483201684484/s. Saved rate range is
2.2947483201684484–111.20723749171252/s. These values verify that the final
500ns candidate is not the positive-plastic650ns bridge state.

Its mean-velocity residual-squared fractions are0.9953570732433161 at particles
with recorded stop events and0.7266643249854793 at multi-contact particles.
This is localization, not a proof of a unique contact derivative cause.
Its final scaled raw residual is72.55309545164955, with last five
528.4221,236.4473,153.9794,107.0168,72.5531: still declining, not a proven
stationary plateau. The32-evaluation failure remains valid and nonrestartable.
**Do not advertise the new FKS response as its demonstrated repair.**

## Executed manufactured search-block tests

The implementation owner chose a small independent `MaterialStageResponse`
class, with the nonsymmetric direct solve implemented inside it rather than
extracting an additional generic class. This satisfies the required separation
from the old symmetric ConstrainedResponse without an unnecessary abstraction.
The physical map and old sources remain unchanged.

Independent contract: `docs/MATERIAL_STAGE_RESPONSE_TEST_CONTRACT.md`, unchanged
after its first run. Test: `tests/review/test_material_stage_response.py`.
Command `.envs/pysph/python.exe -B logs/material_response_20261004/run_stage_tests_initial.py`.
First execution: **12/12 PASS**, unittest0.058 s, external0.2780245 s, exit0.
No repair, actual342 matrix construction or physical evaluation occurred.
Source/test/contract snapshots and before/after hashes:
`logs/material_response_20261004/stage_initial_sources` and
`stage_tests_initial_supervision.json`; complete output `stage_tests_initial.log`.
Tested module SHA256:
`18b32815934a688885c59438aa41438588306cdf77642dc206bffa1d68836161`.

Actually covered: literal scalar and nonsymmetric two-particle uneliminated
block references, exact candidate-continuity kappa, both temperature midpoint
factors, nonzero prescribed-column feedback and algebraic reaction, all-prescribed
empty free block, input/output detachment, zero-defect independent event mean,
unclipped numerical proposals, invalid domains, singular rejection without
pinv, particle/ID covariance, and manufactured dimensional unit conversion.

Scope not run by this suite: `from_journal`, real342 assembly/action, material
derivative accuracy itself (separate suite), refreshed root integration and
physical convergence. Those do not inherit PASS from these12 algebra tests.
