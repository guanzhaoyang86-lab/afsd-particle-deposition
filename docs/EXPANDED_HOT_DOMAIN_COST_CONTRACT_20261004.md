# Independent 8 mm visible hot-domain cost and optional root diagnostic

Registered before any new physical evaluations. This controller is prepared
for explicit scheduling alongside the separately continuing 342-particle
branch; preparation and syntax checking do not authorize concurrent execution.
It never reads a 342-particle accepted checkpoint, never extends its particle
array, never calls an old accepted-copy function, and never publishes acceptance.

Use the frozen `expanded_hot_domain.py` factory/trial adapter and its independent
8 mm visible initialization receipt/checkpoint. Require 918 particles (768
substrate, 150 bar), 1961 free DOF, time/step zero, 723.15 K, zero initial stress,
plasticity and velocity, empty interval-energy history, `accepted=false`, and
`restart_allowed=false`. Every particle/property is loaded from that explicitly
independent snapshot and its source/parameter/geometry hashes are checked.
Record initial mass 0.000309825 kg and initial sensible heat 119.902275 J relative
to 293.15 K separately from zero initial elastic/plastic/kinetic storage. No
claim of initial mechanical coercivity, thermal equilibrium, or process accuracy
comes from the initialization receipt. Initial elastic rank remains NOT_RUN
in this minimum-cost probe and prevents promotion to a process-qualified root.

The physical model and 1.0 mass/time/modulus scales are unchanged. Width is the
explicit new-domain assumption, not a continuation from the small substrate.
The 8 mm initial half-width equals the 2.7 mm shoulder radius plus the 1.3 mm
compact-kernel support radius (2h). This is only initial nominal coverage;
translation immediately reduces that buffer. It does not establish boundary
independence or an unlimited valid process duration. Finite bar inventory and
absence of inlet remain limitations.

## Scheduled stages

Each CLI invocation selects exactly one independent interval, dt=62.5 us or
31.25 us, and one stage (`cost` or `root`). Both dt choices begin from the same
independent time-zero snapshot. They are not successive steps and do not by
themselves constitute an equal-end-time temporal comparison. The guess uses
`batch.guess_for`: inherited particle velocity times dt, with only the registered
fixed axes replaced by `trial.prescribed(dt)`, using the new domain parameters.

`cost`: call `expanded_hot_domain.make_trial`, then exactly one real complete
`evaluate(guess,dt)`. This includes real 3D NNPS, visible graph, density/correction,
material/thermal feedback, paired force, contact, and raw work/heat ledgers.
Record the full arrays/graph and the elapsed time before doing anything else.
Then call the unchanged `core.newton_direction`, which builds one actual frozen
tangent and solves the physical free-axis matrix with the inherited backward
error check. Save the direction and linear diagnostics as unaccepted; do not
take a step. Record stage durations, peak process working set, all map counts,
observed moment/correction conditioning, dense storage estimates, and complete
source/input/artifact hashes. Cost may fail and must retain available evidence.

`root` is a separate optional invocation, allowed only after that dt's complete,
matching cost receipt. Use the unchanged modified-Newton `core.solve`, no
regularization, no reset, and no second root/retry. Let e be measured full-map
time and t be measured tangent-plus-linear time. Choose the largest integer
M<=24 satisfying `30+1.5*((1+8*M+2)*e+M*t)<=450 s`; require M>=2.
The two extra maps reserve one fresh evaluation and one independent scalar
contact full map. The 30 s reserve covers recorded geometry/journal overhead;
it is a forecast, not a guarantee. Enforce a soft deadline at 525 s and actual
remaining-budget checks before fresh maps; an external timer still dominates.
Do not run an unforecasted root when the useful M>=2 budget does not fit.

After a numerical root, recompute the complete map at the same displacement.
Check fresh force, inherited raw and all-prefix energy, sweep certification,
0.1dx penetration, unchanged graph at identical displacement, prescribed axes,
mass/identity preservation, material/density/thermal finite/domain conditions
already enforced by the full map, scalar-versus-batched complete-map identity,
and independent continuous contact-shell and retained-graph visibility paths.
Only SHOULDER/TOP non-bore tool features are in the existing event-screen scope;
unexpected or transient finite contacts fail this diagnostic. Visibility uses
the separately qualified tight certificate, still 256 queries and depth 32.
Store crossing observation only as an unaccepted diagnostic; never update the
initial accepted/history time. Initial elastic rank, temporal/domain convergence,
inlet, sustained transport and deposition remain NOT_RUN even if these root
residual checks pass. A passed root is not a process-acceptance authority.

## Artifacts and budgets

One single-thread worker (`OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=OMP_NUM_THREADS=1`)
per invocation; 540 s internal kill timer and 600 s external timeout, no physical
run during script preparation. New evidence tag/directories only; refuse all
overwrites. Protocol and hashes are written before spawning. Persist the
initial snapshot reference, parameters/provenance/initial energy, explicit guess,
full-map particle/material/contact/thermal/graph arrays and ledgers, tangent
direction/metrics, iterations, failures, wall times, peak memory and supervision
receipt. Every candidate artifact has `accepted=false`, `restart_allowed=false`.
No `accepted.npz` is ever written. A partial/full failure cannot generate a PASS
cost or root receipt. Source hashes are checked before and after each stage.

Exact commands, to be scheduled by the parent controller:

```
.envs/pysph/python.exe -B scripts/probe_expanded_hot_domain_cost.py --stage cost --dt-us 62.5 --tag v1
.envs/pysph/python.exe -B scripts/probe_expanded_hot_domain_cost.py --stage cost --dt-us 31.25 --tag v1
```

Optional roots use the same command with `--stage root`, only when each matching
cost receipt forecasts a root within 450 s. Outputs are under
`outputs/expanded_hot_domain_cost_20261004/v1/dt_62_5us/` and `dt_31_25us/`.
