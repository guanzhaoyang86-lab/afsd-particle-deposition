# HHT local 3D fixture design and reuse audit

2026-10-03. Classification: `LOCAL_BENCHMARK_PROCESS_HHT_3D_CONTACT_BE`.
Status: SUPERSEDED_BY_PROCESS_HHT_FIXTURE_CONTRACT / COMMON_END_TIME_NOT_RUN.
The exact current-geometry fixture requested subsequently is registered in
`PROCESS_HHT_FIXTURE_CONTRACT.md`. The earlier reference-configuration proposal
below is retained as design evidence only; its parameters are not the selected
coupled problem. Its read-only cost measurements remain real but are not a
whole-driver forecast for the selected current-geometry fixture.
This file is a preregistered candidate problem, not an accepted physical state,
an AFSD run, or a certification of the new HHT contact extension. No reference
CSV or fabricated numerical trajectory has been created.

## Reuse findings

- `scripts/probe_engineering_block.py::Block` and
  `configs/benchmarks/engineering_block_e2.json` already provide the required
  genuine 3D PySPH particle array, persistent IDs, NNPS neighbors, corrected
  gradient and paired internal force. `n=3` is supported and independently
  checked by `tests/review/test_engineering_block_e2.py`. Reuse the geometry,
  volumes and operator setup, not its initial zero state or its integrator.
- E2 is explicitly a SMALL_STRAIN_REFERENCE_CONFIGURATION benchmark. Its
  reference geometry and reference density are fixed. Reusing that convention
  here does not exercise moving NNPS, finite deformation or the production EOS.
- `PerfectPlasticBenchmark` in `probe_engineering_material.py` explicitly
  rejects changing temperature. It cannot be used unchanged for coupled BE
  heat, and removing that guard would be a material change. Instead use the
  existing unmodified `material.FKS` and actual `material.update` in a path
  that must remain elastic. The E2 moduli below remain manufactured benchmark
  values; they are not a newly calibrated AA6061 material.
- `tests/review/test_process_thermal_trial.py` already supplies a 27-point
  actual-FKS thermal fixture and an independently enumerated scalar radial
  heat graph. Its equilibrium and single-candidate tests are not a coupled
  HHT trajectory or an independent dynamic reference.
- Q3/Q4/Q6/Q7 scripts supply useful trial/commit, face-ledger and history
  patterns, but their quasistatic equilibrium paths cannot supply an inertial
  reference. The 64-point rotating block supplies actual FKS and contact, but
  uses prescribed normal loads/clamps and has no solved thermal field.

## Exact proposed problem

All coordinates, times and constants below are SI. The following values are
fixed before any coupled fixture execution; selecting a different problem
requires a separately named revision rather than changing a failed run.

Geometry: reference cube `[0,L]^3`, `L=0.01 m`, `n=3`, `dx=0.005 m`.
Enumerate `(ix,iy,iz)` lexicographically in `{0,1,2}^3`, with IDs `1..27`.
`X=dx*(ix,iy,iz)` and `V=dx^3*w(ix)*w(iy)*w(iz)`, where
`w(0)=w(2)=1/2`, `w(1)=1`. Thus `sum(V)=1e-6 m^3` and mass is
`rho0*V`, with `rho0=2700 kg/m^3`, total `0.0027 kg`.
Smoothing length is `1.3*dx`. Use the original E2 reference NNPS/corrected
operator; do not invoke the projected SBP prototype.

Elastic constants are inherited E2 manufactured values:
`E=70e9 Pa`, `nu=0`, `G=35e9 Pa`, `K=E/3`. Use current FKS values
`a1=317.6e6 Pa`, `a2=1e6 Pa`, `a3=70`, `b1=12.44`, `b2=5.65`,
`c1=13`, `c2=0.01`, `room=293.15 K`, `melt=855.15 K`,
`reference_rate=rate_floor=1 /s`. Those FKS values retain their existing
APPROXIMATE/unresolved provenance; their presence is not material validation.
Actual `material.update` must process every candidate history. Never replace
the preload and evolving stress with a frozen elastic force in the coupled
driver. The frozen linear operator is useful only for the separately named
smooth HHT component and cost/spectral diagnostics.

Thermal/friction values are inherited from `configs/debug/starved.json`:
`conductivity=167 W/(m K)`, `cp=900 J/(kg K)`, `mu=0.3`,
`contact_conductance=1000 W/(m^2 K)`, material friction partition `0.5`,
and plastic heat fraction `1`. These remain assumed constants, not fitted
benchmark results. Initial temperature is `400+2*Xz/L K`, initial `ep=0`.

Initial total stress is `diag(-p0,0,0)`, `p0=700000 Pa`. Set
`u0=(-p0/E*(Xx-L/2),0,0)` and
`S0=diag(-2*p0/3,p0/3,p0/3)` so that the E2 linear volumetric term
`K*trace(grad(u0))*I` reconstructs that stress. Initial velocity is
`(0.001*cos(pi*Xx/L),0,0.02) m/s`; it is not subsequently prescribed.
Initial acceleration must be obtained from the actual corrected internal
force plus the initial contact force, divided by each physical mass. It is
not zero and is not borrowed from an old interval-average journal.

Two opposing infinite rigid planes contact the nine `ix=0` and nine `ix=2`
particles. Reference face area per particle is `Ai=dx^2*w(iy)*w(iz)`;
each face has total area `L^2`. Let `delta0=1e-7 m`,
`kn_area=E/L=7e12 N/m^3`, and `ki=Ai*kn_area N/m`. This is an explicit
LOCAL_BENCHMARK linear penalty stiffness, not a change to AFSD tool geometry
or its stiffness calibration. The initial pressure is `kn_area*delta0=p0`.

Set common final time `tf=2e-6 s`, normal stroke
`d(t)=2.5e-8*(1-cos(2*pi*t/tf)) m`, and plane positions

```
left(t)  = min(x0_x) + delta0 + d(t)
right(t) = max(x0_x) - delta0 - d(t).
```

The inward normal force is `ki*max(overlap,0)` on each face, with exact
normal potential `0.5*ki*max(overlap,0)^2`. This maximum defines unilateral
contact, not clipping of an inadmissible computed reaction. Rigid-plane normal
velocities are the derivatives of these positions. Both planes have prescribed
tangential velocity `(0,U)` in their yz tangent coordinates, with
`U=+0.2 m/s` for `0<=t<tf/2` and `U=-0.2 m/s` for `tf/2<=t<=tf`.
Each contact uses the full circular Coulomb disk with its solved normal force.
The left reservoir is `395 K`, the right reservoir `405 K`.

All 81 material degrees of freedom are force driven. There is no roller,
clamp, body-force cancellation of discrete residuals, forced free-particle
motion, or prescribed material velocity after initialization. Normal stroke
loads and unloads; the material law must remain elastic. The tangential switch
creates an exact-time, control-induced reversal of relative slip without
resetting material velocity, stress or temperature. This minimum fixture does
not claim a continuously traversed zero-slip/sticking event or simultaneous
two-face contact on one particle. If either is required, this fixture alone
does not meet that stronger requirement.

The designed wall switch is a discontinuity of force/acceleration, not of
material position or velocity. Split exactly at `t=1e-6 s`. The HHT extension
must derive and record its left/right acceleration history at this switch;
silently carrying a pre-switch endpoint force into an unrelated post-switch
HHT interval is not certified. Detect any earlier contact creation/loss. A
reference or driver that assumes a persistent branch must verify that branch
and stop on its first violation. It may not omit a detected event.

Nominal step grids are `h=250 ns`, `h/2=125 ns`, `h/4=62.5 ns`, giving
`8,16,32` nominal steps to the same `tf`. The designed switch lies on all
three grids. Record actual event subintervals separately from nominal steps.
Run `a=0.9` and the `a=1` comparison only after the event/contact extension
can advance this problem faithfully; the smooth HHT component alone cannot.

## Independent semidiscrete reference

Use one small independent reference function, not another solver/backend.
It must never call the tested HHT step, its trial/force-history map, the local
kick, `material.update`, `internal_force`, `corrected_operators`,
`HyperbolicSpline.gradient`, or either tested heat-advance routine.

1. Enumerate all 351 unordered pairs from the literal reference positions.
   Include a pair if `q=r/h<2`; evaluate the scalar radial derivative as
   `15/(62*pi*h^4)*(3*q^2-6)` for `q<1`, or
   `15/(62*pi*h^4)*[-3*(2-q)^2]` for `1<=q<2`.
   Form radial vector `g=derivative*(Xi-Xj)/r`.
2. Construct each moment `Mi=sum_j Vj*(Xj-Xi) tensor g_ij`, invert it,
   and independently assemble `D`. Form forces either by independent paired
   accumulation `Vi*Vj*(sigma_i*Bi^T+sigma_j*Bj^T)*g`, or by the adjoint
   `-D^T V sigma`. Check the two independent constructions agree before
   integration. Reading a stiffness matrix or gradients out of the tested
   HHT step and integrating those is not sufficient independence.
3. Integrate the small continuous semidiscrete equations for `u,v,S,T`:
   `u'=v`, `m*v'=Fint+Fcontact`,
   `S'=2*G*dev(sym(Dv))+W*S-S*W`, `W=skew(Dv)`, and
   `sigma=S+K*trace(Du)*I`. This retains the continuous corotational elastic
   law; replacing it by `sigma=E*sym(Du)` silently would omit the actual
   material update's spin contribution. Evaluate the literal FKS yield
   expression independently and require the entire reference path to remain
   elastic. Plastic work and `ep` then have the independent reference zero.
4. Independently evaluate the two plane forces and Coulomb slip formula on
   the persistent slipping branches. Integrate both smooth halves separately
   with a high-order adaptive method, e.g. existing SciPy DOP853. At an
   unexpected zero-slip or contact event, stop and report the missing branch
   rather than regularizing friction or extrapolating across it.
5. Independently construct the heat graph with
   `gij=Vi*Vj*2*k*(-radial_derivative/r)`, capacities `mi*cp`, face
   conductances `1000*Ai` while contact is active, and signed reservoir heat
   `Hi*(Twall-Ti)`. Integrate `C*T'=-(L+H)*T+H*Twall+0.5*Qfriction_rate`.
   Include separate integral states for normal/tangent wall work, face
   impulses, signed reservoir heat, raw friction heat, and plastic work.
   The driver uses BE; the reference solves the continuous heat problem, so
   agreement checks temporal accuracy rather than merely copying BE algebra.

Use nondimensional reference state scaling: displacement `1e-7 m`, velocity
`0.02 m/s`, stress `p0`, temperature increment `1 K`, and energy `1e-5 J`.
Run reference tolerances `rtol=1e-10`, scaled `atol=1e-12`, then one independent
tightening to `rtol=1e-11`, scaled `atol=1e-13`. The two results must agree
within one percent of every registered fixture error allowance. This is an
oracle self-check, not a generated HHT truth trajectory. Keep both raw results
and runtime; if the 60 s reference budget is exceeded, report NOT_VERIFIED.

An independent semidiscrete oracle certifies temporal response of this spatial
model only. E2's reported finest-grid continuum stress error was 27.53%, above
its 20% bound; it remains FAIL. Constant stress plus continuum face traction
does not automatically balance each corrected-SPH particle at a free surface.
Consequently, a claimed exact rigid translation obtained by cancelling each
particle's tested force with an applied body force would be tautological and
would violate the force-driven fixture requirement. Do not use that shortcut.

## Preregistered accuracy and accounting gates

These proposed bounds precede any coupled run. They are numerical error
budgets for this local problem, not relaxed AFSD acceptance thresholds.
Compare every grid at all common output times `t=k*250 ns`, including both
sides of the prescribed control switch when force is discontinuous.

- Finest-grid volume-weighted RMS displacement error <= `1e-9 m`, velocity
  <= `2e-4 m/s`, full stress <= `7000 Pa`, and each face reaction component
  <= `0.7 N`. Also report center-of-mass velocity/displacement and the first
  nonrigid eigenmode projection, so local energy closure cannot substitute
  for low-frequency response accuracy.
- Each face impulse error <= `1.4e-7 N s`; each face wall-work and cumulative
  raw-friction-Q error <= `2e-7 J`. Compare signed reservoir heat separately
  per face with absolute error <= `2e-8 J`; do not cancel their opposite signs.
- Maximum temperature error <= `1e-7 K`, measured on `T(t)-T(0)` to expose
  the small thermal response. Also report the reference increment magnitude.
  If the actual increment is below 10 times numerical resolution, call the
  temperature-accuracy part inconclusive; do not claim heat coupling passed
  solely because absolute 400 K temperatures agree.
- `ep` and plastic work must remain zero within `1e-15` and `1e-14 J`,
  respectively; record the yield margin on both paths. This verifies a
  reversible elastic path and zero plastic heating only. Positive plastic
  work accuracy is NOT_TESTED; a plastic benchmark would be a separate scope.
- Require decreasing errors under refinement; whenever a coarse error is
  larger than ten times its oracle uncertainty, the fine/coarse error ratio
  must be <= `0.7`. Do not demand second-order temperature convergence from
  BE or infer spatial convergence from this fixed geometry.
- Retain mass relative `1e-12`, ID/tag independence, correction condition
  `1e6`, finite state, supported solid-temperature domain, per-particle BE
  heat balance, signed contact/impulse/work identities and existing physical
  guards. Boundary force and work must be recorded separately on each face.
  HHT algorithmic energy, return diagnostic, friction heat and plastic heat
  remain distinct. A mechanical balance remainder is not algorithmic loss
  and must never be deposited into material temperature.

The inherited E2 displacement domain is retained:
`max(norm(u))/L < 1e-4`; do not increase it after observing a failure.
The nominal initial displacement already includes the nonzero preload.

## Measured cost and bounded launch decision

A read-only component-cost audit was actually executed using
`.envs/pysph/python.exe -B -c <in-memory audit>` on 2026-10-03. It constructed
E2 `n=3`, inspected its elastic spectrum, and timed 100 calls each to the
existing actual FKS update and BE thermal component. It did not advance a
coupled trajectory or save a checkpoint. Exit code: 0.

| Quantity | Measured value |
|---|---:|
| Particles / reference NNPS pairs | 27 / 305 |
| Total mass | 0.0027000000000000006 kg |
| Maximum correction condition | 3.385745485548915 |
| Rigid null eigenvalues in the 81-DOF elastic matrix | 6 |
| First positive elastic angular frequency | 575220.73849 rad/s |
| Largest elastic angular frequency | 1560692.005997 rad/s |
| E2 setup wall time | 0.0033949 s |
| Mean material-plus-BE call wall time | 0.001787699 s |
| Entire in-process audit wall time | 0.184324 s |
| External process wall time | 0.625952 s |

The timing probe used `k=180 W/(m K)` solely as a timing input; the proposed
fixture uses the inherited `167 W/(m K)` specified above. It used actual
unmodified FKS and physical masses/moduli. The spectrum is a reference elastic
diagnostic, not a substitute for the evolving material/contact path.
The coarse `omega_max*h` is about `0.390173`. This does not certify contact
stability or HHT accuracy.

For 56 nominal steps, eight material/thermal evaluations per step would cost
about 0.80 s in those two kernels alone. This excludes HHT residual assembly,
linear solves, event localization, independent reference, imports and output;
it is NOT a verified whole-driver forecast. Before the first common-end-time
run, measure one complete private candidate including event/ledger/output
work and register a whole-run estimate. Require each component/reference
suite externally <=60 s, with a conservative launch forecast <=45 s.
Stop after at most the authorized two evidence-based repair rounds.

No actual AFSD probe is proposed here. Its separate existing 600 s ceiling
and own common-end-time/cost registration remain in force. Successful smooth
HHT algebra, BE heat, or this cost audit cannot authorize a physical run.
If the HHT endpoint path currently accepts only event-free passive intervals,
the honest joint status is BLOCKED_ON_CONTACT_EVENT_GLUE: preserve the design,
report smooth components separately, and do not label a frozen-stress or
force-only substitute as completion of this coupled fixture.

## Selected current-geometry reference: actual evidence

The subsequent selected contract is `PROCESS_HHT_FIXTURE_CONTRACT.md`; its
exact 2 micrometre preload/current geometry parameters supersede the earlier
candidate above. The independent implementation is
`scripts/reference_hht_block.py`. Its shared `fixture_data()` returns only
physical arrays/parameters; the coupled HHT driver must implement its own
force, contact, material and heat calls.

First launch and tests both exited 1 before executing an ODE because SciPy
is absent from the pinned PySPH environment. Logs are preserved as
`logs/process_hht_20261003/reference_initial.log` and
`reference_tests_initial.log`. No installation was attempted. Before the first
actual integration, the contract was amended to independent classical RK4
with 128 and 256 steps per smooth half (implementation repair round 1).

Actual independent reference run: exit 0,
`logs/process_hht_20261003/reference_round1/report.json` and
`reference_history.npz`. Total in-process wall time 1.333157 s; external
process wall time 1.660182 s, below the external 60 s bound. The two grids used
1024 and 2048 RHS evaluations. Four literal fixture/operator/reference tests
passed in 0.055 s, external 0.429901 s, recorded in
`reference_tests_round1.log`. No physical AFSD run was launched.

| Reference refinement difference | Measured maximum |
|---|---:|
| Volume-weighted displacement RMS | 1.98788e-17 m |
| Volume-weighted velocity RMS | 3.61540e-11 m/s |
| Volume-weighted full stress RMS | 2.98708e-4 Pa |
| Temperature | 4.71800e-12 K |
| Face force | 1.54219e-9 N |
| Face impulse | 1.54721e-15 N s |
| Face wall work | 3.01019e-15 J |
| Face raw friction heat | 8.79999e-17 J |
| Face reservoir heat | 1.35525e-20 J |

Every refinement difference is below the preregistered one-percent oracle
uncertainty allowance. Recorded minimum signed slip is 0.116063 m/s, minimum
overlap 1.507035 micrometres, and minimum independently evaluated FKS yield
margin 313.302 MPa. The RHS also checked each RK stage against the persistent
contact, signed-slip >.05 m/s and elastic-domain guards.

An independent raw-ledger readback is saved as
`logs/process_hht_20261003/reference_readback.json`: maximum
`Delta K + Delta U_normal + stress_work - wall_work + Qraw` defect
1.13620e-15 J; total momentum/face-impulse defect 1.05032e-19 N s;
particle-versus-face cumulative raw heat defect 1.35525e-20 J.
The maximum thermal-energy balance defect is 3.96895e-12 J. Final raw
friction heat is about 9.765585e-6 J per face, reservoir heat into material is
4.499996e-6/8.099995e-6 J, and maximum temperature increment is 3.656319e-5 K.
These are reference-ODE quantities, not results of the tested HHT step.

Status: independent semidiscrete elastic/contact/thermal reference PASS;
coupled HHT common-end-time accuracy and cost remain NOT_RUN in this audit.
The preload, elastic stress-rate history, current correction operators and
signed ledgers were retained. Positive plastic work, continuous zero-slip
sticking, shared-particle multi-face contact, continuum boundary accuracy and
physical AFSD process validation remain outside this minimum reference.
