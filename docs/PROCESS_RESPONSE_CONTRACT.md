# Numerical contact-response work package, registered before execution

2026-10-03. The previous turn was progress:233 components passed, an independent
reversal reference reached6.5us, and a newly solved local interval reached6us.
The AFSD state did NOT advance;850.393633ns and zero crossings remain.
Preserve the4516-entry PROCESS_STOP_FREEZE and every old failure. Only HANDOFF,
IMPLEMENTATION_PLAN and OUTPUT_SCHEMA may be updated as continuation indices.

## Target and explicit change in hypothesis

Test whether the implemented finite-quadrature nonlinear residual benefits
from its contact response in the numerical preconditioner. The old hypothesis
of obtaining a classical continuous cross-zero transverse Jacobian is NOT
assumed true: independent analysis shows a PRIVATE impulse |epsilon| cusp.
The continuous integral and its finite-node approximation are different maps.
The cusp and endpoint singularity must be independently demonstrated, not
hidden by central finite differences or a smooth-looking Gauss result.

New component name: FIXED_QUADRATURE_DISCRETE_CONTACT_PRECONDITIONER_EXPERIMENT.
It is not a new friction law, an implicit contact-model certification or a
zero-event solver. It changes only the trial search direction. The actual
force, adaptive contact quadrature, material, density, thermal solve and ALL
old acceptance/collinearity/zero/passivity/domain gates remain untouched.
The experiment may fail and must then stop with full rejected data.

## Independent checks before any private coupled call

1. Literal analytic point Kx/Kv, physical and PRIVATE smooth-side signs,
   normal gap perturbations, rotations and finite-difference force tests.
2. Complete endpoint-plus-DeltaJ FORCE residual against independent fixed-rule
   Gauss construction/differences, including material-point position/velocity
   acceleration responses. Also compare with independent higher quadrature on
   nonzero smooth paths. Bounds are separately registered in RESPONSE_TEST_CONTRACT.
3. Continuous PRIVATE crossing cusp/one-sided limits and failure of a classical
   transverse Jacobian. This is a PASS counterexample, NOT a successful event.
4. Retain scalar reversal and true-vector guards; test coupled wrappers for
   no physics/parameter/parent mutation and nonrestartable output.

Every matrix has declared units: Kx N/m, Kv N*s/m, D_A F kg,
D_A J kg*s, D_A R(force units) kg. No zero-speed floor, clipping, pinv,
mass/time/modulus scaling or derivative residual assigned to heat.

## Single saved-state diagnostic, conditional on the checks

Use ONLY published6us checkpoint from process_stop/scout_initial, not private
trial3/4 or rejected trial5. Restore its full physical/canonical history.
Duration is the unchanged saved trial005 value2.7254762326048973e-7s.
Selected contacts/directions come from its saved predictor, not desired tags
or independent oracle time. This is one same-parent, same-equations numerical
experiment; no trajectory, publication, event restart, or output to old dirs.

At each evaluable current candidate, assemble a block response from the exact
saved Gauss nodes of that candidate. The derivative holds node selection fixed;
it does NOT differentiate adaptive rule switching or claim a continuous
cross-zero tangent. PRIVATE negative-side derivative has its correct opposite
sign. Zero-override samples, nonfinite blocks, changed normals/kinematics or
out-of-scope constraints reject. The next candidate still goes through the
unmodified true nonlinear residual and all old physical gates.

Add only D_A Rcontact to the old elastic/density-pressure Schur preconditioner.
Never add the mass twice, alter the material residual, silently symmetrize the
contact block, or assume a PRIVATE block is positive-definite. Record its norms,
condition/spectrum diagnostics where meaningful, and all line-search outcomes.
At most24 complete residual evaluations, seven backtracks, as before.

Forecast from the old18-evaluation failed diagnostic and measured point cost;
require estimated full call<=50s. Internal55/external60s, both below600s.
At most one PRIVATE saved diagnostic this package. No repeated real scouts.
Whether or not the root converges, save accepted=false/restart_allowed=false.
Compare raw residuals to old1e-9m/s,1e-12m,1e-10density,5e-12K gates; do not
reclassify the old failure. Convergence is numerical-only, not event completion.

## Keep the end-to-end target visible

The goal remains finite open-tool transport in a real3D deformable AFSD model,
not extending an elastic benchmark forever. An independent source/evidence
review must list actual integration gaps and feasible next milestones before
this experiment can justify another trajectory. Fixed-normal elastic walls
do not certify a rotating bore, plastic deposition or a usable process cost.
New physical approximations or long/paid jobs require explicit direction.
