# Next single work package: inertial process-step accuracy and algorithmic energy

2026-10-03. Entry point after the current derivative experiment, irrespective
of whether its full root runs or converges. Status: DESIGN / NOT_IMPLEMENTED /
NOT_RUN. This document creates no new accepted physical state and authorizes
no production run. Preserve all previous source, checkpoints and failures.

## Decision and bounded outcome

Develop and independently test ONE named inertia-retaining time-integration
variant, `PROCESS_HHT_EVENT_THERMAL_APPROXIMATE`, on one small three-dimensional
loaded/contact/thermal benchmark. The deliverable is a runnable same-history
step, an independently checked algorithmic-energy account, common-end-time
accuracy evidence, and a measured feasible-step/cost decision. More accepted
500 ns windows are not the deliverable. Do not simultaneously develop another
backend, a quasistatic solver, new constitutive law or production AFSD case.

Select HHT with explicit future-force weight a=0.9, gamma=0.6, beta=0.3025;
a=1, gamma=0.5, beta=0.25 is the nondissipative comparison within this family.
These are numerical choices, not fitted aluminum or physical damping values.
For the smooth force-only limit, define the convention directly:

    x1 = x0 + h v0 + h^2 [(1/2-beta) acc0 + beta acc1]
    v1 = v0 + h [(1-gamma) acc0 + gamma acc1]
    M acc1 = a F(x1,v1,t1) + (1-a) F(x0,v0,t0).

The contact extension must be derived, not inferred from these smooth formulas.
The official [OpenSees HHT documentation](https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/integrator/HHT.html)
was checked on 2026-10-03: it gives these gamma/beta formulas and limits linear
stability claims to linear problems. We use the displayed weight convention
without importing a package or inferring nonlinear/contact stability.
[MOOSE Dynamics](https://mooseframework.inl.gov/modules/solid_mechanics/Dynamics.html),
also checked, distinguishes inertia, physical damping and HHT time weighting.
Neither source validates this project's contact, SPH, FKS or process result.

## Why this work precedes another physical continuation campaign

- `docs/PROCESS_GLOBAL_STEP_AUDIT.md`: at 500 ns the frozen elastic
  omega_max*h is approximately10; root convergence does not resolve that mode.
- `src/afsd_backend/pysph/process_thermal_trial.py`: material uses one L/T
  per candidate interval; the local kick freezes force, normals and wall speed.
  Its actual event-average velocity is independent of endpoint mean velocity.
- `src/afsd_backend/pysph/material.py`: noncoaxial/thermal evolution within a
  larger interval is not resolved by the fixed48-iteration return calculation.
  Existing FKS assumptions and Jaumann free-energy failures remain evidence.
- `logs/process_event_20261003/event_response_initial.log`: continuous-ODE
  response accuracy failed at4.06763264e-5 versus2e-5 s/kg. The separate discrete
  derivative PASS explains that failure; it does not restore continuous accuracy.
- `logs/process_event_20261003/preflight_stable/report.json`:138 passing tests
  and1 retained continuous-reference failure; this is a qualified numerical gate.
- `outputs/process_event_20261003/response_probe/manifest.json`: the initial
  cloud probe stopped at PID1's one-sided derivative check, with0 new accepted
  states. Its failure is retained even if a later registered variant succeeds.

## Implementation scope and compulsory identities

Reuse the fixed PySPH particle/NNPS/ID facilities, physical masses/moduli,
corrected paired spatial force, current material update and conservative BE
thermal component. Add the named variant alongside the current implementation.
Keep candidate material/thermal histories private and restart every nonlinear
evaluation from the same accepted state; only a fully accepted step commits.

Specify acceleration history explicitly. An interval-average acceleration in
an old journal is not automatically the endpoint acc0 needed by HHT. Derive
its consistent initialization from the selected problem and contact constraints;
store it in the new accepted schema and test its restart equivalence. Do not
zero old velocity/stress or borrow history from any rejected candidate.

Derive a compatible event path, its displacement/actual mean velocity, and
contact/actuator impulse and work identities. The old kick is a reference and
reusable local facility, not proof that its frozen-force path can be inserted
unchanged into HHT. Resolve contact/control switches at their time locations;
update geometry and material time layers consistently at registered subintervals.
Any altered physical map or quadrature must be named and tested explicitly.

Derive the algorithmic-energy identity before running benchmarks. First prove
it for the elastic smooth limit with known stored energy. Report its raw signed
terms; do not force positivity by clipping. In the contact/thermal case, keep
wall/actuator work, raw friction heat, plastic heat, return diagnostic and
algorithmic energy separate. A leftover balance residual is not automatically
algorithmic dissipation, and algorithmic loss must never become material heat.
The old hypoelastic energy remains diagnostic where no free energy is certified.

## One benchmark package and preregistration gates

Use two inexpensive independent references supporting ONE coupled3D fixture:

1. A two-mode elastic oscillator with omega=1 and100 rad/s, nonzero initial
   displacement/velocity and no physical damping. Use common final time10 s
   with h=.1,.05,.025 s. Compare low-mode displacement/velocity to the analytic
   solution; report high-mode attenuation, phase and the independent discrete
   amplification/energy identity for a=.9 and the a=1 control.
2. The already documented two-capacity analytic thermal problem in
   `PROCESS_COMPLETION_TIME_ASSESSMENT.md`, including signed reservoir work.
3. One small3D loaded block, reusing an existing independently tested fixture
   and its manufactured material provenance, with nonzero initial stress and
   velocity, two contact faces, loading/unloading and a slip reversal. Couple
   mechanics and BE heat in the same accepted transaction. Keep all free
   material degrees of freedom force driven.

Before implementing the fixture driver, save its exact geometry, load/control
program, material values, common final time, h/h2/h4 and independent reference
construction. Fix numerical/physical error bounds before its first execution;
do not choose them from observed coarse-step errors. The independent reference
must not call the new step or generate its own oracle through the tested map.
Require low-frequency response, reaction/work, cumulative Q, plastic work and
temperature accuracy, not only a small nonlinear residual. Keep original raw
residual/physical-domain tolerances; any new acceleration or energy gate needs
its stated units and independent rationale in that preregistration.

Each component/reference suite has an external60 s bound. Allow at most two
evidence-based implementation repair rounds, keeping every initial failure.
Measure complete build/solve/event/output cost. A later actual AFSD accuracy
probe requires its own fixed common-end-time registration and a forecast under
the existing600 s ceiling; it is not launched automatically by component PASS.
If accuracy-acceptable and budget-acceptable steps do not overlap, report that
result and stop this work package. Do not respond by extending iterations,
loosening physical gates or accumulating more equivalent microsteps.

## Transport dependency and authority boundary

The current accepted configuration has maximum displacement50 um while the
nearest original particle center is250 um from the shoulder. See
`PROCESS_COMPLETION_TIME_ASSESSMENT.md` and the continuation resolved config.
The nominal250 um/feed scale is0.1179245283 s, not a predicted crossing time.
Use it only to show the required cost gap. A feasible integrator still needs
finite-edge work-consistent contact (`PROCESS_FINITE_TOOL_COMPONENT.md`) and
accepted-history-bound inlet/ID/source transactions (`PROCESS_INLET_COMPONENT.md`).
Those remain separate subsequent integration dependencies, not completed here.

Existing engineering-approximation authority permits this explicitly named,
independently tested numerical variant while retaining inertia and old failures.
Deleting inertia, omitting the physical preload/startup interval, replacing it
with a separately equilibrated hot initial state, or starting long/paid/GPU/
production runs requires new explicit user direction. Material/friction changes
are outside this work package; any later alternative needs its own named scope
and independent benchmarks under ENGINEERING_APPROX_CONTRACT, not a silent swap.
No mass/time/modulus scaling, free-particle scripted trajectory, velocity reset,
thermal clipping, residual heat or unmeasured crossing may bridge the gap.
