# PROCESS_HHT_EVENT_THERMAL_APPROXIMATE: smooth linear component

Preregistered 2026-10-03, before executing this component's tests. This new
component does not modify or validate old physical states, contact, material,
SPH, thermal or AFSD calculations. Its independent references are algebra and
analytic linear ODE solutions. All quantities below use SI units.

## Equations, acceleration history and API

For constant symmetric positive-semidefinite K [N/m], positive diagonal M [kg],
external load f [N], displacement x [m], velocity v [m/s], and acceleration z
[m/s^2], set w=0.9, d=1-w, gamma=1/2+d, beta=(1+d)^2/4. The control is w=1,
gamma=1/2, beta=1/4. This is the future-force convention in
PROCESS_EVENT_NEXT_WORK_PACKAGE.md, not a physical damping coefficient.

    x1=x0+h*v0+h^2*((1/2-beta)*z0+beta*z1)
    v1=v0+h*((1-gamma)*z0+gamma*z1)
    M*z1=w*(f1-K*x1)+d*(f0-K*x0).

The solved matrix is M+w*beta*h^2*K. No mass/modulus scaling, velocity reset,
stress reset, clipping or pseudoinverse occurs. All input arrays are copied.
The smooth start is z0=M^-1*(f0-K*x0); subsequent z0 is the saved previous
endpoint acceleration, not an interval-average acceleration.

`kinematic_prediction(x0,v0,acc0,dt,weight=.9)` returns a dictionary with
`x_predictor`, `v_predictor`, `beta`, `gamma`. `initialize_acceleration` performs
the smooth start. `step(mass,stiffness,x0,v0,acc0,f0,f1,dt,weight=.9,
previous_displacement_increment=None,previous_load_increment=None)` returns
a detached result with `x`, `v`, `acc`, `displacement_increment`,
`load_increment`, and a SI-labelled `ledger`. Missing increments mean an
explicit virtual start with zero previous increments; they are not a restart.
A restart stores x, v, acc, the two previous increments, h, and weight.
No result from this component is a contact or physical-state acceptance.

## Independently derived algorithmic-energy equality

Let dx=x1-x0, dz=z1-z0, df=f1-f0; let b and g denote the previous displacement
and load increments. Define p=f-K*x, c=beta-gamma/2=d^2/4, and q=d*(1/2-d).
The force-history defect is

    e0=M*z0-p0+d*(g-K*b) [N].

Valid saved HHT history has e0=0. The virtual smooth start b=g=0 with the
initialization above also has e0=0, including nonzero x0 and v0. Let bars mean
arithmetic endpoint means. Direct subtraction of the kinematics gives

    dx=h*vbar+h^2*c*dz,
    dv=h*zbar+h*d*dz.

Set r_i=M*z_i-p_i. The equilibrium and history definitions imply

    r1=-d*(p1-p0),  r0=-d*(g-K*b)+e0,
    d*M*dz+rbar=q*((p1-p0)-(g-K*b))+(1/2-d)*e0.

Multiply the two kinematic identities by their conjugate M/K terms. With
E=(v^T*M*v+x^T*K*x)/2, this yields first

    Delta(E+h^2*c*z^T*M*z/2)+h^2*d*c*dz^T*M*dz
       =dx^T*fbar+dx^T*(d*M*dz+rbar).

Using dx^T*K*(dx-b)=(dx^T*K*dx-b^T*K*b+(dx-b)^T*K*(dx-b))/2
then proves the final identity, without defining any term from its residual:

    H0=E0+h^2*c*z0^T*M*z0/2+q*b^T*K*b/2
    H1=E1+h^2*c*z1^T*M*z1/2+q*dx^T*K*dx/2
    Da=h^2*d*c*dz^T*M*dz
    Dx=q*(dx-b)^T*K*(dx-b)/2
    W=dx^T*(f0+f1)/2                         [endpoint-load quadrature]
    Wf=q*dx^T*(df-g)
    Wh=(1/2-d)*dx^T*e0
    H1-H0+Da+Dx=W+Wf+Wh.                         [all terms J]

For the registered weights and exact PSD K, H>=0, Da>=0 and Dx>=0. For
consistent unforced history and CONSTANT h, H is nonincreasing. Physical E
alone has no asserted per-step monotonicity. For w=1 and consistent history,
H=E and the equality reduces to trapezoidal-load work conservation. Wf and Wh
are raw signed terms; neither is heat. Wh must be reported if history is
inconsistent. The equality is still algebraically exact with that defect.
W is not generally exact continuous-path external work. For the explicitly
registered cubic Hermite displacement interpolation whose velocity is
v(s)=v0+s*(v1-v0)+6*s*(1-s)*(dx/h-(v0+v1)/2), 0<=s<=1, and the linear force
interpolation f(s)=f0+s*df, direct polynomial integration gives
Wexact=W+h*df^T*(v1-v0)/12. The ledger reports this separately as a conditional
path-integral value, together with Wexact-W; it does not replace W in the
algorithmic identity or certify this interpolation for contact.
For changing h, each local identity holds, but H0/H1 use that step's h, so an
extra change in the acceleration-energy coefficient must be included before
telescoping. No adaptive-step monotonicity or contact energy theorem is claimed.

## Fixed tests, references, bounds and budget

Independent amplification: solve the three simultaneous scalar equations for
[x1,v1,z1], then advance by the resulting 3x3 matrix. This reference must not
call the component or its prediction, coefficients, ledger or initialization.
The energy reference must evaluate the displayed quadratic forms independently.

Two-mode fixture: M=diag(2,3) kg, frequencies (1,100) rad/s,
K=diag(2,30000) N/m, x0=(.3,.002) m, v0=(-.2,.12) m/s, f=0;
z0=-M^-1*K*x0. For each w in (.9,1), run to common T=10 s with
h in (.1,.05,.025) s. Report both modal amplitudes and wrapped phase errors,
spectral radii, energies, runtime, and low-mode analytic displacement/velocity
errors. The analytic solution is x=x0*cos(omega*t)+v0*sin(omega*t)/omega,
v=-omega*x0*sin(omega*t)+v0*cos(omega*t).

Acceptance fixed before tests:

- Independent trajectory maximum errors <=2e-12 m, 2e-11 m/s, 2e-9 m/s^2.
- Each raw energy identity <=5e-11 J; cumulative identity <=2e-9 J;
  consistent-history defect <=1e-10 N. These absolute roundoff allowances are
  generous for O(1 J) energy, O(100 N) force and at most 400 steps, and do not
  replace a physical error test. Independently recomputed ledger terms must
  agree to 2e-12 J. Nonnegative terms have 1e-14 J roundoff allowance only.
- Low-mode absolute x and v error <=0.6*h^2 m and m/s, respectively; the joint
  x/v error decreases by a factor between 3 and 5 at each halving. These use
  a conservative second-order coefficient relative to unit-frequency motion
  over 10 s, fixed here before results. They do not certify the high mode's
  phase accuracy; the high mode is intentionally unresolved at coarse h.
- Amplification spectral radius <=1+1e-12; high-mode final attenuation for
  w=.9 is below the w=1 control. Unforced H may increase by at most 5e-11 J
  per step; w=1 final physical energy error <=2e-9 J.
- Free flight (K=f=0) with nonzero x/v reproduces its analytic path within
  2e-12 m and 2e-12 m/s. Ramp force f(t)=.7+.2*t N on m=2 kg, K=0,
  x0=.3 m, v0=-.2 m/s, T=10 s uses the same h, weight cases and absolute
  0.6*h^2 bounds for x/v against its cubic/quadratic analytic solution.
- Moving spring target y(t)=.2*sin(.4*t) m, m=2 kg, K=2 N/m,
  x0=.3 m, v0=-.2 m/s, f(t)=K*y(t), T=10 s uses the same h, weights.
  Its analytic forced oscillator gives x/v bounds 0.6*h^2 m,m/s and final
  physical energy/work error <=0.8*h^2 J. Independent direct discrete
  forcing remains a separate roundoff-level implementation test.
- Restart results must be bitwise equal when the full saved history is copied;
  inputs remain unchanged and output mutation must not change inputs.
  A deliberately wrong initial acceleration must produce a nonzero explicit
  history defect and history-work term, while satisfying the signed equality.

All component runs have an external 60 s limit. Retain every log, including
failures; allow at most two evidence-based implementation repairs. A passing
component is not permission or evidence for a physical AFSD run.

## Initial execution evidence (added after the preregistered run)

`PYTHONPATH=src .envs/pysph/python.exe -B -m unittest discover -s tests/review
-p test_hht_linear.py -v`: 5 PASS, 0 FAIL, exit 0, 0.885 s unittest time,
1.080343 s external elapsed time, external ceiling 60 s. No implementation
repair was needed. Standard output and standard error are preserved in
`logs/process_hht_20261003/hht_linear_initial.log` and
`hht_linear_initial_stderr.log`; supervision metadata is in
`hht_linear_initial_supervision.log`.

For w=.9 and h=.1/.05/.025 s, low-mode displacement errors at T=10 s were
.0034493938826/.00086782001183/.00021744795449 m. Joint displacement/velocity
refinement ratios were 3.9908287575 and 3.9971356740; the control ratios were
3.9954975361 and 3.9988748446. Independent amplification, signed identities,
consistent history, analytic free/ramp/forced motion, and conditional exact
path work all met their fixed bounds. The copied full-history restart was
bitwise equal. Deliberately inconsistent acceleration produced force defects
(2,3) N, signed history work -0.00918916076252 J and identity residual
-6.07153216592e-17 J; the defect was exposed, not converted to heat.

The HHT two-mode physical energy increased on 31/73/146 steps while the
proved augmented H satisfied its nonincrease bound. This is direct evidence
against imposing a physical-energy per-step monotonicity gate. High-mode final
amplitude ratios were 2.43472654e-8/5.76595531e-12/7.84237815e-11; these show
attenuation, not resolved high-frequency phase accuracy. No AFSD run occurred.
