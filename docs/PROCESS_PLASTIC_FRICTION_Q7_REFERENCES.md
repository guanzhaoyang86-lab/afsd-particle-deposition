# Q7 independent plasticity/friction coupling references

Status at definition: **NOT_RUN**. This document fixes the independent reference
before executing Q7 tests. Constants are explicitly manufactured SI data, not
AA6061, FKS, calibrated values, or paper parameters. No evaluated solver output
is used to construct an expected answer below. Existing material, spatial and
real-AFSD failures are not superseded by this interface benchmark.

## Why a single perfectly-plastic spring in series is insufficient

A single rate-independent perfectly-plastic tangent spring with constant yield
force Y in series with a constant friction slider C generally activates only the
lower force threshold. If Y<C the material yields while the wall sticks; if
C<Y the wall slips while the spring remains elastic. The equality Y=C is
nonunique and is not a useful simultaneous-mechanism verification. Adding
hardening just to avoid this issue would introduce new unrequested physics.

The minimum reference here therefore uses an actual J2 normal material point
to generate the current normal reaction, which determines the circular tangent
friction capacity of an independent elastic two-component support. This is a
one-way coupled manufactured assembly: normal plasticity changes the capacity
used by the tangent solve. Both irreversible mechanisms operate during the
same accepted intervals. It does **not** certify combined normal/shear yielding
in one material point, the full 3D block discretization, or arbitrary finite
strain. Those require the separate actual-body checks and limitations.

## Manufactured assembly and boundary program

The normal bar has reference length L=1 m, area A=1 m^2, volume V=1 m^3,
shear modulus G=1500 Pa, bulk modulus K=1000 Pa, hence E=3000 Pa and nu=0.
J2 yield stress is a fixed 3 Pa, temperature is prescribed 293.15 K, and there
is no hardening, thermal solution, inertia, or physical heat conversion.
These low moduli/yield/large size belong only to a toy reference; they are not
rescaled versions of the actual 3D block and cannot replace its parameters.

Define positive normal compression e and the additive small strain
epsilon=diag(-e,0,0). Its spin vanishes. The actual existing material updater
must receive this gradient history; the benchmark must not use a scalar
replacement to generate its tested normal force. Pressure closure is
sigma=S+K*trace(epsilon)*I. Force on material is -N in the normal direction.
The independently obtained scalar J2 solution is

    q = 3000*e for 0<=e<=.001, then q=3 Pa,
    S = diag(-2*q/3,q/3,q/3),
    ep = 0 for e<=.001, then (2/3)*(e-.001),
    N = -A*sigma_xx = 1000*e + 2*q/3.

The free tangent displacement is a two-vector u. The elastic tangent support
has stiffness k*I with k=1000 N/m. Its force balance is tau=k*u. Circular
Coulomb friction has mu=.2 and current capacity C=mu*N. Wall tangent motion is
along a=(3/5,4/5); the friction direction must be obtained from actual relative
displacement, not simply copied from a.

The two-stage program uses seconds only as a quasi-static load parameter:

1. t=0 to 1 s: e=.001*t and tangent wall w=(0,0).
2. t=1 to 2 s: x=.002*(t-1), e=.001+x, and w=2*x*a m.

At t=1 the material is exactly at first yield, N=3 N, C=.6 N, u=tau=0.
Throughout stage 2, N=3+1000*x and C=.6+200*x. Sticking gives u=2*x*a and
norm(tau)=2000*x. Equating this norm with C gives the first slip event

    x*=1/3000, t*=7/6 s,
    w*=u*=(1/1500)*a m, tau*=(2/3)*a N.

After this event the exact branch is

    u=(.0006+.2*x)*a m,
    tau=(.6+200*x)*a N,
    d(u-w)=-1.8*a*dx,

so slip and friction oppose each other without any absolute-value repair.
Normal plastic strain continues increasing while the tangent contact slips.

## Literal terminal states and separate work ledgers

All terminal values are at the common t=2 s, not different failure times.

| quantity | independent reference | unit |
| --- | --- | --- |
| normal strain e | .003 | 1 |
| deviatoric stress S | diag(-2,1,1) | Pa |
| total stress sigma | diag(-5,-2,-2) | Pa |
| equivalent plastic strain ep | 1/750 | 1 |
| normal reaction N | 5 | N |
| tangent wall w | (.0024,.0032) | m |
| material tangent u | (.0006,.0008) | m |
| tangent force tau | (.6,.8) | N |
| accumulated relative tangent displacement u-w | (-.0018,-.0024) | m |
| normal external/stress work Wn | .0095 | J |
| material plastic work P | .004 | J |
| normal elastic diagnostic En | .0055 | J |
| tangent support energy Ut | .0005 | J |
| tangent wall work Wt | .003 | J |
| raw friction dissipation Q | .0025 | J |
| total external work Wn+Wt | .0125 | J |

These follow from independent integrals:

    Wn = integral_0^.001 3000*e de
         + integral_0^.002 (3+1000*x) dx = .0015+.008 = .0095 J,
    P = 3*ep*V = .004 J,
    En = V*(S:S/(4*G)+K*e^2/2) = .001+.0045 = .0055 J,
    Ut = k*(u dot u)/2 = .0005 J,
    Q = integral_(1/3000)^.002 (.6+200*x)*1.8 dx = .0025 J,
    Wt = Ut+Q = .003 J,
    Wn+Wt = En+Ut+P+Q = .0125 J.

All these force/strain branches are linear or constant between their explicitly
located events. Trapezoidal stress/contact work is exact for this aligned path.
The reference is not a claim that general nonlinear Q7 block work must equal
its current elastic diagnostic plus plastic work at finite increments.

The return diagnostic sum is separate from physical dissipation. For exactly
three accepted intervals ending at t=1, 7/6, 2 s, the two plastic increments are
dep_1=1/4500 and dep_2=1/900. Thus the explicit algorithmic reference is

    sum (1.5*G*dep^2*V) = 13/4500 J.

This diagnostic changes under subdivision: it must not be added to P, Q, or a
thermal source merely to improve energy closure. If more intervals are used,
13/4500 J is no longer that run's reference; the force/state/work references
above remain valid for event-aligned branches.

### Additional unload with a stationary tangent wall

For t=2 to 3 s hold w=.004*a m and let e=.003-.001*(t-2). The normal material
unloads elastically: S goes from diag(-2,1,1) to zero, ep stays 1/750 and N
decreases linearly from 5 to 2 N. Capacity decreases from 1 to .4 N, so the
tangent spring releases while the wall remains stationary: u changes from
.001*a to .0004*a m. The signed friction work is still positive,

    delta Q = -.7*(-.0006) = .00042 J,
    delta Wn = .5*(5+2)*(-.001) = -.0035 J,
    delta Wt = 0.

The final S=0, sigma=-2*I Pa, u=(.00024,.00032) m, tau=(.24,.32) N,
En=.002 J, Ut=.00008 J, cumulative P=.004 J and Q=.00292 J. Total external
work is .009 J and equals .002+.00008+.004+.00292. The negative normal work
is elastic energy release, not a negative dissipation to be clipped. The
return diagnostic adds zero on this elastic unloading segment.

## Tests and preregistered tolerances

The independent test assembly may use the existing actual material updater and
the actual cone-response/row/event functions as the system under test. Literal
expected values in assertions must stay independent of those functions. An
accepted state is copied into pure candidates; probes/Newton/bisection must
not accumulate ep or plastic work. Only a successful commit changes either
material or contact state. Repeated rejected/probe evaluations must leave the
accepted snapshot unchanged.

- Force/stress/displacement/work references: absolute 1e-10 in their stated SI
  unit, rtol=0. Ep absolute 1e-12. Event time absolute 1e-10 s. These allowances
  are much smaller than the manufactured .0025/.004 J mechanism signals.
- Fixture root target 1e-12 N, slip target 1e-12 m, <=48 event bisections to a
  <=1e-14 s bracket, finer than the 1e-10 s reference tolerance so that the
  3.6 N/s cone-margin derivative leaves <1e-12 N boundary error. Numerical
  budgets are not physical time or mass scaling.
- Zero tangent wall motion control: P remains .004 J and N remains 5 N while
  Q=Wt=Ut=0. This catches conflating the two dissipation ledgers.
- Basis representations 0 and 37 degrees must leave global terminal u/tau,
  event time and all work unchanged. This is not physical tool rotation.
- Ep/force/work references are checked after repeated candidate evaluations
  so double accumulation cannot pass merely by finding a force root.
- Raw Q and material P retain their signs; neither abs nor clipping is allowed.
  No fictitious temperature, thermal heat, or mass/time/modulus scaling exists.

The actual 27-point Q7 block must separately pass predeclared equilibrium,
contact-law, work-identity, purity/rollback and passive-label tests. It need not
match this toy bar/spring terminal table: geometry, stiffness and boundary
conditions differ. Space refinement, thermal conversion, rotated normal/tool,
production restart, full AFSD and exact paper reproduction remain NOT_RUN.

Sources: independent force balance and integrals in this document; the tested
material update is src/afsd_backend/pysph/material.py, and the explicit constant
law is scripts/probe_engineering_material.py. The existing E1/Q2 references and
Q6 disk references remain unmodified. No new constitutive model is introduced.
