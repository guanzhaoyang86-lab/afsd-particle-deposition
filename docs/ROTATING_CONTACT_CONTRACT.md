# Two bounded rotating-contact units, preregistered 2026-10-02

Prior archive PROCESS_PREDICTOR_SENSITIVITY_FREEZE:1776 entries PASS. Preserve
all prior physics/config/FAIL/output; only HANDOFF/IMPLEMENTATION_PLAN/OUTPUT_SCHEMA
may be updated. No backend switch, GPU, performance work or production run.
This is ENGINEERING_APPROXIMATE contact-component development, NOT AFSD deposition.

## Unit 1: kinematics and work interface, no old physics changes

Reuse geometry.tool_pose/local_xy/surfaces and stable body-fixed faces:
0 BASE,1 SHOULDER,2 BORE_X_MINUS,3 BORE_X_PLUS,4 BORE_Y_MINUS,5 BORE_Y_PLUS.
Contact identity=(persistent particle_id,feature_id); retain two corner contacts.
Store force_application_point x and finite-face projection c=x-gap*n separately.
The latter is NOT a complete closest point on a closed solid mesh. Top/outer-side/
true edge contacts remain unsupported; keep original buried-centre guard.
Compute wall velocity at BOTH points: V+omega cross (point-origin).
Existing force DOFs are particle translations at x. Run CENTRE_COLLOCATED_APPROXIMATE:
material v_i and wall v(x) are paired. Do not swap only wall velocity to v(c).
An optional L*(c-x) material reconstruction is diagnostic ONLY. Production use
would require its transpose force distribution to all affected DOFs, not fictitious
particle spin. Record delta torque=(c-x) cross F and delta power=F dot(vwall(c)-vwall(x));
neither is heat. Raw friction power=-F_t dot u_t; normal relative work separate.

Independent literal geometry (m,s): width1,dx.2,layer1,top2,radius2,traverse.3,omega2.
At t0,x=(.45,.2,1.5): X_PLUS gap.05,n=(-1,0,0),c=(.5,.2,1.5),
wall(x)=(-.1,.9,0),wall(c)=(-.1,1,0). For F=(-10,3,4)N,
delta torque=(0,-.2,.15)Nm,delta power=.3W. With matching rigid v_i and
L=[[0,-2,0],[2,0,0],[0,0,0]], reconstructed material and wall slip both zero.
Hole centre remains open above/below shoulder; (.45,.45,1.5) has two faces.
Test independent 90-degree pose, finite-difference wall point velocity, general
3D omega, passive tags/ID permutation, nonfinite/duplicate input guards, signed work.
Tolerances: positions/velocity/power1e-12 SI for unit-scale literals; centred
finite-difference velocity1e-8 m/s (subtractive/truncation allowance, no fitted value).

Actual old kick at rotating-shoulder frozen sample: x=(.5,0,.95),wall=(.3,1,0),
n=(0,0,-1),m=Fn=mu=1,other=(0,0,1),v0=(1.3,1,0),dt=.2:
v1=(1.1,1,0),Jt=(-.2,0,0)Ns,Q=.18J,Wwall=-.06J,deltaK=-.24J.
Reference constants analytic; tolerances1e-10 SI. Do not call it a rotating trajectory.

## Unit 2: evolving shoulder contact, not a deformable-material solver

Reuse PySPH ParticleArray, pinned HostLeapfrog drift/kick/drift, existing coupled
contact_kick.kick, and checkpoint/VTK facilities. Four point masses, FREE tangential
coordinates, constrained z velocity0. Prescribed normal load1N per contact has
recorded actuator reaction; no normal penalty calibration or new material law.
No continuum internal force/material/NNPS needed in this isolated contact benchmark;
do not fabricate temperature/stress/plasticity/strain/density fields as solved values.
This does not rigidify or replace substrate/feedstock in the AFSD solver.

SI parameters and sources: configs/benchmarks/rotating_shoulder_component.json.
Initial xy=(.01,0),(0,.012),(-.01,0),(0,-.012); z=.0013m; zero velocity.
Tool layer.0015m,top.02m,radius.02m,open square width.008m,dx.0005m,
base depth.01m,traverse0,omega20rad/s; particle masses.001kg,mu.2,Fn1N.
ALL manufactured component values, not fitted AA6061 or paper configuration.
Centre/projection difference is parallel to omega for shoulder, so their wall
velocities coincide; bore discrepancy is NOT resolved by this special case.
Physical end.0008s, dt2e-5/1e-5s (40/80 steps), same geometry/loads/end;
no CFL/continuum convergence claim. Wall force determines all xy motion, no prescribed
circular trajectory or hidden centripetal force. Contact feature must stay SHOULDER.

Independent sliding ODE: r'=v; v'=(mu Fn/m)*(omega cross r-v)/|omega cross r-v|
in xy, z fixed. Simultaneously integrate Q'=mu Fn*|slip| and
Wwall'=F_t dot(omega cross r), using plain RK4, not kick/cone/kinematics functions.
Reference steps1e-7 and5e-8s. Reference must stay |slip|>=.03m/s; otherwise this
sliding-only reference is invalid/BLOCKED, never force continuation through zero.
Reference refinement differences <=1e-11m,1e-8m/s,1e-11J for Q/W.
Component vs finer reference: position max vector error<=1e-8m, velocity<=1e-5m/s;
Q/W difference <=1e-10J+2e-4*|reference|. Engineering component criteria, not paper accuracy.
Record signed energy residual K-K0-Wwall-Wconstraint+Q (normal work zero), <=5e-12J;
mass/ID exact, z drift<=1e-14m, finite arrays, signed subinterval Q unchanged and
existing local-kick rejection bounds unmodified. Negative raw values never abs/clipped.
Additional independent rotation ledger (added before any full trajectory/reference
run): delta Lz=sum[(r_mid cross (Jn+Jt+Jconstraint))_z], tolerance2.5e-13 N m s;
Wwall=omega*delta Lz, tolerance5e-12J. Midpoint r is actual application point
after the first drift, not geometric face projection or final r. This follows
discrete impulse/angular momentum algebra for drift/kick/drift and the fixed
origin, translating DOFs of this benchmark; not proof of general contact torque.
No broadening old near-zero/large-geometry solver certification from this case.

Gate sequence: unit1 literals/frozen kick PASS -> unit2 short component and reference
tests PASS -> ONE run each40/80steps (no tuning retries), independent saved-array readback.
Initial cost estimate<30s per path and<30s per reference pair; hard external60s each
(all below user600s), singlethread. Estimate from tiny4-body/40–80-kick workload,
not measured claim. Report actual times. Failure stops next full-path attempt;
preserve accepted/rejected state. Additional runs require a new declared purpose.
No source/screen changes in response to failed physical results; retain all logs.

Deliver history with time,x,v,ID/tags,tool pose,contact points/normal/wall velocity,
impulses,normal reaction,rawQ/W/K,hashes,checkpoint plus NOT_SOLVED field list.
Sort/restart tests compare same ID. Benchmark restart metadata must not authorize
AFSD production restart. Rotating bore/multiface geometry integration, deformable
block under rotation, feeding, normal opening, thermal/finite strain, production NOT_RUN.
Existing full historical suite NOT_RUN; selected old contact regression may run.
All config parameter entries inherit the top-level source and their explicit
unit/status; none are unspecified zeros. Required PySPH scratch rho/e/h properties
in checkpoints are NOT_SOLVED, not exported as physical VTK fields. Dynamics of
this manufactured component are planar xy under z constraint despite 3D interfaces.
