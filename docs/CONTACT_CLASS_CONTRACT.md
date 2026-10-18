# Contact class certification contract — 2026-09-30

Status: PRE_TEST. This stage replaces particle-by-particle progress with a fixed
local-method matrix. No AFSD real trial, GPU, production, backend switch, material
replacement or spatial-operator edit in this stage. Keep every previous freeze,
failure record and physical parameter. Prior CONTACT454_REPAIR_FREEZE: 360 files
verified PASS before these edits. PySPH pin unchanged.

## Mathematical problem and independent acceptance

Frozen normal penalty/internal forces/geometry and actuator mobility P/m remain
the model. Tangent columns T_a are orthonormal. For each local midpoint/static
solve, minimize phi(f)=b:f + (c/2)(sum T_a f_a)^T M(sum T_a f_a),
subject to ||f_a||<=limit_a. M=P/m is PSD; H_ab=c*T_a^T M T_b.
The force decomposition need not be unique: rank(H)<=3 while f has 2N entries.
No new physical regularization/penalty is introduced to pick a preferred answer.

Independent convex optimality certificate: g=b+Hf. Feasibility plus
gap=sum_a(g_a dot f_a + limit_a*||g_a||) bounds phi(f)-phi(opt).
This follows by minimizing each linear g_a dot y_a over its disk; it does not
call the tested solver or require its iteration path. Also independently check
the projected stationarity map and known optimal objective. Keep raw signed gap
within floating point error; do not use it to alter heat or force.
For nonunique solutions compare M times resultant force, objective and KKT,
not arbitrary per-face components. Unique analytic solutions retain original
1e-11+1e-10 relative force checks (454/439).

Numerical acceptance (double precision, registered BEFORE new tests):
- force disk feasibility normalized by largest capacity: <=2e-12;
- independent dimensionless KKT gap and stationarity: <=2e-10;
- manufactured known objective error <=2e-10 after scaling by common force
  and velocity/acceleration scales; unique endpoint/resultant <=2e-9 scaled;
- solver requested tolerance remains2e-12; certificate internally targets
  min(tolerance,2e-14) to support strict old unique-force tests, not a loosening;
- contact time/heat scalar exact benchmarks keep their old tolerances;
  family endpoint/heat/work checks use1e-10 absolute+relative where analytic;
- large common-velocity invariance has additional explicit input-representation
  allowance 64eps*boost. It must not hide material-relative O(1) errors.
- old physical heat sign checks, event accuracy and iteration budget stay in
  force. No abs(Q), max(Q,0), artificial heat, velocity clipping or mass scaling.

## Fixed matrix and budget

Cone manufactured KKT: 1/2/3 faces; free mobility ranks0/1/2/3;
all-interior/mixed/all-sliding; orthogonal and nearly parallel normals;
capacities spanning1e-3..1, plus zero capacity; rotation and face permutation;
common force/velocity numerical rescalings1e-6,1,1e6. Scales in tests change
representation/problem magnitudes only, NOT production parameters or mass.
Construct f* and nonnegative boundary multipliers FIRST, then independently
set b=-H f*-lambda f*. This is a known convex KKT reference, not a solver output.
Persist every individual case, including failures. Historical439/454 are fixed
holdout cases; no parameter tuning against their force values.

Kick family:1/2/3 parallel faces sharing total friction; free/constrained axes;
stop/stick/reversal; stationary/moving and different wall velocities;
rigid coordinate rotations, common velocity boosts including a boost orthogonal
to the slip plane. Existing corotation/actuator/work/curved-slip tests remain.
Separate frozen cone validity from event integration validity.

First run the new matrix against current code, preserving baseline failures.
At most TWO bounded local-solver implementations/repair rounds in this stage;
do not grow scope to new constitutive/geometry algorithms. Each subprocess <=60s.
If certification cannot be obtained, record the unresolved class and stop local
patching. No real2us trial even if this matrix passes: material/boundary decisions
and a separate small coupled benchmark must be resolved first.

## Allowed implementation and non-contact decisions

Only contact_kick.py/local cone helper may change physics numerics, with no law
change. Retain old projected iteration callable as a historical reference. A
scaled semismooth projected-KKT solve may use an EXPLICIT rank-revealing least
squares Newton correction for a nonunique force decomposition, always verified
by the full original convex certificate; this is not pinv(M_i), stiffness
regularization or a spatial-neighborhood bypass. Singular spatial M_i still fails.
Tests may correct the invalid unique-force oracle, preserving its snapshot and
adding a regression accepting two distinct certified equilibria.

Material decision: retain declared Jaumann approximate path unchanged; do not
certify finite-strain conservative elastic energy. Free-energy alternative stays
DESIGN_ONLY until independently approved/verified, never a silent paper substitute.
Boundary decision: retain corrected/adjoint pair and old nodal FAIL. Inspect
continuum traction vs weak work under fixed physical cube refinement n=3,5,7;
do not define external traction as minus tested internal force to force a PASS.
Use independent analytic affine/cubic surface-work integrals; report convergence
or lack thereof without moving physical boundaries.

Pre-implementation fixture audit: initial baseline exposed a zero-objective
normalization overflow in the independent oracle; baseline_valid_oracle exposed
normalization of a theoretically zero rotated tangent block by its own rounding
noise. Both logs remain. Final reference scale is max(||b||,c||M||*force_scale),
with unit scale only for an identically zero objective. This is the full physical
mobility scale, not division by numerical null-space noise. The matrix also
checks c*M*sum(force-reference)/velocity_scale <=2e-9. No UUT change preceded
these fixture corrections; final pre-implementation baseline is baseline_final.

Round1 outcome: all432 cone cases and439 pass; existing curved-stop fraction.5
rejects Q=-3.91654e-19J at a residual relative velocity about1e-16m/s.
The instantaneous-only rounding bound shrank to6.11e-31J after subtracting the
wall frame and stopping. It omitted accumulated velocity error from an initial
relative speed0.728011m/s. Round2 adds this SAME invariant characteristic scale
to the existing256eps impulse-times-speed rounding bound. No raw heat is changed;
no physical negative heat is made positive. Store the bound per accepted interval
and test it against an independent recomputation. Existing exact reversal and
large-boost -0.5J counterexample must still be distinguished from rounding.
