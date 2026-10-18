# Coupling failure investigation v1 (2026-09-30)
User-authorized diagnostic/minimal-contact revision; not production approval.
Prior SCIENTIFIC_FREEZE.json remains immutable, verified PASS before editing.
New sources intentionally invalidate its current-source entries; retain old snapshots.

## Predeclared sequence and thresholds
1. Diagnostics only: legacy law/material/corrected spatial operator unchanged.
   Max40 full kicks, <=600s/run; predicted <30s with logging. Save preceding accepted
   kick and rejected candidate, distinguish virtual impulse decomposition from
   actually applied stages. Compare observer-off/on particle arrays and ledger
   bit-for-bit at every accepted state and same first failure.
2. A independent scalar reversal: m=1,u0=1,Fother=-2,muFn=1,dt=1.
   Exact: tstop=1/3,ufinal=-2/3,Q=7/18,Wother=1/9,deltaK=-5/18.
   Event-aware reference comparison abs1e-10 SI; work identity abs1e-12+rel1e-10.
   Retain legacy algorithm RED; true reversal is allowed.
3. B actual material.update, elastic-only test law (infinite yield, NOT production).
   F=exp(e*A), A=diag(1,-.5,-.5), e=0.001; Jc=.99; V0=1/rho0,m=1.
   Four exact exponential path segments, 16/32/64/128 substeps per segment.
   Log sigma:D*V work with midpoint volume/stress, final F/S, Es and e scratch status.
   Closed-path geometry/stress tolerances 1e-11 relative; work noise1e-9J.
   Determine limiting work rather than assume zero; compare independent formula.
4. Contact changes only after identified evidence, separately named variant.
   Local contact work rel1e-10+abs1e-12J; endpoint analytic rel/abs1e-10.
   No heat clipping/abs; 600s real run cap. Preserve legacy code/tests.
5. Spatial probe tolerances: affine/objectivity/power/force/torque 1e-10 relative,
   constant-stress patch with independently defined surface traction (do not
   mistake nonzero free-boundary traction for failure), distorted neighbors,
   and explicit zero-mode measurement. cond>1e6 rejects, no pinv.
6. Only after all user preconditions run dx=.5mm,CFL=.1 to2us, then .05,.025;
   only equal-time success allows comparison, then dx=.25mm same physics.
   Original N08/N09/N10 thresholds unchanged. Any absent run is NOT_RUN.
No GPU/parallel/mass/time/modulus scaling; no material replacement in this task.

## Contact-variant pre-run addendum (after phase 2, before integration)
Frozen geometry/normal/internal force within each upstream kick, collocated wall
and force velocity retained. Per-particle free-DOF projector P; face coupling
C_ab=R_a P R_b/m, R_a=I-n_a*n_a^T. Disk projection solves Coulomb admissibility;
this is not clipping negative work. Straight-slip zero event is split exactly;
curved slip uses midpoint step-doubling (rtol2e-7, atol1e-12 m/s, max4096
intervals) and a 0.2 direction-step guard. Cone iterations max2000, relative
force update tolerance2e-12; nonconvergence rejects, never falls back silently.
Independent curved-slip RK4: 8192 steps, endpoint/work abs2e-6 SI.
First frozen failing particle: independent polar-coordinate analytic integral
of tangent ODE; default heat/stop-time relative1e-3, tight solve relative1e-5. This is a
local frozen-load benchmark, NOT an independent complete AFSD reference.
Straight analytic cases retain abs1e-10 SI, kick identity abs1e-12+rel1e-10.
Closed loop demonstrated nonzero limiting work; no constitutive substitution.
Nodal boundary patch compares independent trapezoidal surface tractions on
[-1,1]^3; weak affine balance is separately reported, not substituted for nodal
equilibrium. Single-cloud zero-mode SVD threshold1e-10 largest singular value.

The initially planned DOP853 check was NOT_RUN: isolated environment has no
SciPy (attempt2 log retains import failure). No dependency was installed.
Use stronger closed-form independent reference instead, with unchanged tolerances.
For constant tangent force magnitude f<c=mu*Fn, r0=|u0|,
b=(1+cos(theta0))/2, d=(1-cos(theta0))/2:
tstop=m*r0*(b/(c-f)+d/(c+f));
Q=m*c*r0^2*(b^2/[2(c-f)]+b*d/c+d^2/[2(c+f)]).
Derivation: r(theta)=C*tan(theta/2)^(c/f)/sin(theta);
integrate dt=-m*r/(f*sin(theta)) dtheta and Q=integral c*r dt to theta=0.
