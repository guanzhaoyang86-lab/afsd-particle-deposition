# Bounded current-geometry coupled HHT component

2026-10-03. Preregistered before its unit execution. One new component in the
existing PySPH backend; no old solver, accepted AFSD state, material parameters,
or physical mass is replaced. Classification:
`PROCESS_HHT_EVENT_THERMAL_APPROXIMATE_COUPLED_ELASTIC_BENCHMARK`.
The forthcoming independent loaded/contact/thermal 27-node common-final-time
fixture is a separate required accuracy/cost gate. These component tests do not
complete that gate or authorize AFSD continuation.

## API and force convention

`initialize(pa,p,boundary,time=0)` returns accepted numerical history.
`advance(pa,state,p,boundary,dt)` returns `(detached_pa,new_state,diagnostics)`.
`restart_at_control_event(pa,state,p,new_boundary)` returns
`(new_state,reset_diagnostics)` while preserving every physical particle field.
`restore_checkpoint(template_pa,state,p)` restores all saved particle properties
into an equivalent ParticleArray layout and checks the complete semantic hash.

The boundary callback receives a detached ParticleArray and time. It must return
`force_definition='NONCONTACT_APPLIED_FORCE_ONLY'`, an Nx3 `force_N` containing
only noncontact body/applied forces, and a contact list. The component adds each
normal and tangential contact force exactly once. Every row supplies
`particle_index`, stable `contact_id`, unit `normal`, `normal_force_N`,
`tangential_force_N`, `wall_velocity_m_s`, `wall_temperature_K`,
`conductance_W_K`, `penetration_m`, `normal_potential_J`, and
`normal_stiffness_N_m`. Penalty force/potential must match the supplied physical
stiffness and penetration. No free material velocity is prescribed by this API.

Parameters include `rho0`, `shear`, `bulk`, `cp`, `conductivity`, `room`, `melt`,
`friction`, `plastic_heat_fraction`, `heat_partition`. Optional FKS fields are
passed unchanged to the existing FKS class; omitted fields retain its existing
defaults. `hht_weight` defaults to .9; gamma=.6 and beta=.3025 follow the
registered family. The named current-geometry fixture's manufactured E, nu,
density and initial conditions are fixture provenance, not calibrated AA6061.

## Coupled map and bounded nonlinear solve

Accepted history stores true endpoint total physical force F0 and acceleration
A0. Every candidate starts from the same old x/v/S/ep/rho/T. Newmark provides
endpoint x/v from candidate A1; `hht_path.sample(...,.5)` supplies the actual
quadratic-velocity/cubic-position midpoint. At that position use the existing
corrected paired SPH gradient, candidate midpoint density and temperature, and
the existing 48-iteration Jaumann/J2/FKS update. This component explicitly
rejects any plastic increment; the actual FKS elastic branch remains active.
Density uses implicit midpoint continuity's Cayley factor. Endpoint force uses
the actual endpoint geometry, density and updated stress, not frozen K.

Endpoint-linear contact force and affine wall velocity use the existing new
HHT path audit: exact Gauss raw friction work, whole-polynomial passivity,
zero-slip event detection, stable identities and fixed normals. Interior
events or unsupported changes reject the private candidate. In addition,
the capacity-relative defect from instantaneous Coulomb sliding traction is
checked at endpoints, both Gauss nodes, .25/.5/.75 and all friction-power
stationary points. `contact_coulomb_relative_tolerance` defaults to .001 and
may be preregistered up to .01. This sampled error check is explicitly a
quadrature approximation, not a continuous Coulomb-law proof. It is separate
from the exact polynomial passivity guard and the independent fixture reference.

BE heat uses raw radial pairs at the dense midpoint, actual Gauss friction Q,
and end-state reservoir conductances/temperatures. The same nonlinear solve
checks candidate T1 against BE's mapped T1. Numerical force-filter work,
normal temporal residual, and mechanical return diagnostics never become heat.

The initial numerical correction used `diag(m)+a*beta*h^2*K_frozen` but updated
density separately. The retained coarse step-6 failure exposed missing pressure/
density cross terms. Repair 1 uses the same existing frozen structural operator
in the coupled Schur preconditioner derived below. It only proposes numerical
corrections; every trial recomputes real current-geometry material/boundary/
thermal forces. The physical residual equations are unchanged.
Use at most 24 full residual evaluations, including at most seven halvings per
line search. No finite-difference global Jacobian, increased budget, or physical
scaling is allowed. Acceptance thresholds, fixed before execution, are
`h*max|A1-Amap| <= 1e-9 m/s`, associated Newmark position defect <=1e-12 m,
density defect/rho0 <=1e-10, and temperature defect <=5e-12 K. These retain
the prior mechanical scales and a tightly solved BE candidate; they are
numerical-root criteria, not physical accuracy thresholds.

All state, force/acceleration history, particle IDs, complete particle-property
snapshot, physical/parameter hashes, parent hash, timestamps and control-reset
record participate in the semantic state hash. Rejected evaluations cannot be
restarted or committed. Save/load must preserve every array value and dtype.

## Work and reset accounting

The exact signed `hht_path.kinetic_identity` records physical force-interpolation
impulse/work, algorithmic force-filter impulse/work, root residual impulse/work,
and their kinetic identity separately. Spatial conjugacy and BE balances are
checked with existing tolerances. Normal potential/work mismatch remains a
signed temporal diagnostic. The actual hypoelastic model has no certified
free energy here; no positive total HHT dissipation claim is made for it.

An explicit control restart recomputes the real right-limit force and `A=F/m`,
preserving x,v,S,ep,rho,T and all labels. It records both histories and the
change in the clearly defined last-step acceleration storage
`.5*h_last^2*(beta-gamma/2)*sum(m*A^2)`. That partial numerical-storage jump is
computable; full nonlinear HHT algorithmic-energy jump is explicitly
NOT_CERTIFIED. Neither is thermal energy. This all-free API does not implement
static contact constraints or generic stop/stick reinitialization.

## Unit execution budget and independent references

One suite, external wall-time cap 60 s, at most the work package's two retained
repair rounds. Planned tests: 27-node zero equilibrium/full semantic checkpoint;
literal constant-body-acceleration motion; independently assembled all-pairs
thermal graph; loaded passive contact with nonzero initial stress/velocity,
raw heat, same-history repeatability and passive-tag invariance; exact control
force restart and continuation equivalence; ambiguous-force and plastic-scope
rejection without input mutation. Analytic motion uses <=2e-18 m and
<=2e-14 m/s roundoff checks; independent thermal graph uses <=2e-12 K.
The contact test checks transaction/work properties, not coupled accuracy.

Command: `.envs/pysph/python.exe -B tests/review/test_hht_coupled.py`.
Results and every initial failure must be saved under a new dedicated log
directory. Production and the 27-node common-final-time physical comparison
remain NOT_RUN until their separate driver executes and reports evidence.

Initial component execution: **6 PASS**, 0.174 s test time / 0.680323 s external
command time, exit 0 under an explicit 60 s subprocess bound. No repair needed.
Evidence: `logs/process_hht_coupled_20261003/COMMANDS.md` and `initial.log`.
This result covers only the six stated component tests.

## Evidence-based numerical block repair1 after the first complete-sweep attempt

The above initial preconditioner is retained as historical design, not the
current correction. Coarse step6 of the registered driver failed bounded line
search. Saved evidence and original source are under
`logs/process_hht_coupled_20261003/repair1/`. A standalone replay establishes
that the previous block omitted the pressure-density residual term and the
density response to the acceleration correction. Actual physical residuals
and all24-evaluation/error/grid/physical limits are unchanged by this repair.

Let H be the frozen divergence map, V=m/rho, C=H^T diag(V K/rho0),
kappa=gamma/2+1.5*(beta-gamma/2), and
L=-diag(rho_old*h^2*kappa/(1+h*div(vmid)/2)^2)H.
With Kdev=Kfrozen-H^T diag(V K)H, A=M+a*h^2*kappa*Kdev,
the numerical block is A*da-a*C*drho=-RA, -L*da+drho=-Rrho.
Use (A-a*C*L)*da=-RA-a*C*Rrho and drho=L*da-Rrho.
It remains only a correction proposal; candidate actual force/current geometry,
material, continuity and heat are all reevaluated before acceptance.

Independent uneliminated block and C=0 tests pass. The exact saved failed step
now accepts in6evaluations (velocity residual7.08412e-11m/s). Component suite
9PASS, external.873s. Added failure candidates keep raw A/rho/T/x/v/force and
residuals, separately marking last line-search base and last evaluated trial
REJECTED, accepted=false/restart_allowed=false. No bad state becomes a restart.
This one-step regression does not itself certify the completed time sweep.

## Repair 1: pressure/continuity block, independently diagnosed

The initial 8-step grid failed its sixth step at saved time 1.25 us after 19
evaluations. Retain `sweep_initial/grid_008/REJECTED.json`, all five accepted
stages, and `logs/process_hht_coupled_20261003/repair1/source_before.py.txt`.
At the line-search stall the acceleration force residual norm was .001816 N,
but the omitted density-induced pressure correction was .003485 N. An actual
1% old direction increased merit 1050.0136 -> 1050.4333. The analytically coupled
direction decreased it to 1039.1241; its full step decreased it to 60.2470.
These saved-state diagnostics took .3898 s and changed no source or history.

Let H be the accepted frozen corrected divergence map, V=m/rho, and
`kappa=gamma/2+1.5*(beta-gamma/2)=.30375`, the exact dense midpoint velocity
response to h*A1. Let `R_A=M*A1-a*F1-(1-a)*F0` and
`R_rho=rho1-rho_map`. A consistent frozen block uses

    C = H.T diag(V*Kbulk/rho0)
    L = -diag(rho_old*h^2*kappa/(1+h*div(vmid)/2)^2) H
    Kdev = Kfrozen-H.T diag(V*Kbulk) H
    Adev = M+a*kappa*h^2*Kdev
    Adev*dA-a*C*drho = -R_A
    -L*dA+drho = -R_rho.

Eliminate density without dropping its residual:

    (Adev-a*C*L)*dA = -R_A-a*C*R_rho
    drho = L*dA-R_rho.

The old numerical update omitted both `-a*C*R_rho` and `L*dA`, even though
its acceleration matrix contained bulk stiffness. The corrected block's
measured algebraic residual was 2.52e-18. Frozen geometry, prestress and
density-dependence omissions remain numerical-preconditioner approximations;
the unchanged full residual decides every candidate's acceptance. No damping
coefficient, tolerance, physical value, evaluation cap, or physical grid changes.

New regressions compare against an independent uneliminated joint block,
check the C=0 case, rerun only the saved failed step under the same root gates,
and ensure rejected diagnostics retain acceleration/rho/T, all raw residuals,
and the physical candidate hash with `restart_allowed=False`. No diagnostic
run advances or publishes a whole trajectory. Results are retained under
`logs/process_hht_coupled_20261003/repair1/`.
