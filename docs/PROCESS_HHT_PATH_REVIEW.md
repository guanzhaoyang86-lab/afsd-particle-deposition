# Independent review of the initial HHT path component

2026-10-03. Read-only review of `hht_path.py`, `test_hht_path.py`,
PROCESS_HHT_PATH_CONTRACT and PROCESS_HHT_CONTACT_DESIGN. The component's
initial PASS evidence remains intact. This review changes no implementation,
test, material, contact law, or accepted state. Findings below refer to the
initial reviewed source; later repair evidence must be recorded separately.

## Confirmed algebra

The Newmark-compatible quadratic velocity integrates to its stated endpoint
displacement. Two-point Gauss correctly integrates a linear force times that
quadratic velocity. The exact conditional physical path work is

    Wpath = Fbar dot dx + h*DeltaF dot DeltaV/12.

For delta=1-weight, q=delta*(1/2-delta), r0=M*a0-F0, and
e1=M*a1-weight*F1-(1-weight)*F0, the returned impulse identity correctly uses

    Jphysical=h*Fbar,
    Ialgorithmic=h*((1-gamma)*r0+q*DeltaF),
    Iroot=h*gamma*e1.

Its kinetic-work expression is the corresponding exact signed identity; it
does not pretend that algorithmic work is positive or thermal heat. This
agrees with the independently derived smooth energy identity when its distinct
endpoint-load quadrature and exact linear-load path work are kept separate.

## P1: an interior zero-slip event is missed at a double zero

Initial `contact_audit` locates zero speed only by finding stationary roots of
the quartic squared-speed polynomial. At a double zero of a slip component,
the derivative of its squared speed has a triple root, which is ill-conditioned
under an ordinary companion-matrix polynomial root solve.

Reproducer: normal `(0,0,1)`, h=1 s, weight=.9, zero wall velocity,
v0=v1=`(.25,0,0)` m/s, a0=`(40,0,0)` m/s^2,
a1=`(-80/3,0,0)` m/s^2, f0=f1=`(-1,0,0)` N and both capacities 1 N.
The exact specified path is `u_x(s)=(s-.5)^2` up to arithmetic roundoff; the
midpoint is zero slip. The initial implementation returned:

    whole_path_passive=True
    requires_event_split=False
    interior_zero_fractions=[]
    minimum_checked_slip_m_s=1.3505752072262567e-11
    zero_speed_bound_m_s=1.2789769243682066e-13

This is a false event-free certificate. Add this regression before repair.
Component quadratic roots, followed by the original full-vector zero-speed
check, provide a small independent candidate source and are already allowed
by the design. Retain stationary-point checks for vector-speed minima and
whole-path passivity. Do not fix the failure by enlarging the zero-speed gate.

## P2: invalid interval/data can return an apparent certificate

The initial path-work/contact functions do not check finite positive dt.
With the same normal/force/capacities, v0=v1=`(.25,0,0)` and a0=a1=0,
`contact_audit(dt=-1)` returns `whole_path_passive=True` and raw heat -0.25 J.
dt=0 also returns a certificate. These invalid intervals should fail before
building polynomial or work records.

`linear_work` checks f1 and path arrays but omits f0 from its finiteness check;
NaN f0 returns NaN work. `kinetic_identity` checks mass positivity but omits
mass finiteness; NaN mass returns NaN identity values. Add rejection regressions
for these concrete holes. A shared finite-array/positive-dt validation helper
is sufficient; no new solver or tolerance is needed.

## Required scope regression: passive interpolation is not full Coulomb

This is already correctly disclosed by
`interior_coulomb_direction='NOT_CERTIFIED_BY_PASSIVITY'`; it is not an algebra
failure. Preserve that distinction with an explicit regression:

    h=1, v0=(1,0,0), v1=(0,2,0), a0=a1=(-1,2,0)
    normal=(0,0,1), wall0=wall1=0
    f0=(-1,0,0), f1=(0,-1,0), capacity0=capacity1=1.

Both endpoint forces are exact saturated Coulomb forces and interpolation is
passive. At s=.5, however, interpolated force is `(-.5,-.5,0)` N while the
continuous Coulomb force is `(-1/sqrt(5),-2/sqrt(5),0)` N. The interpolated
magnitude is 1/sqrt(2) N despite 1 N capacity. A coupled caller therefore
needs the registered independent direction/time-error gate, not only this
path's passivity boolean. No extra full-Coulomb claim should be added here.

The initial source also evaluates reported power through expanded polynomial
coefficients although the design specifies direct original-vector evaluation
at candidate locations. Its existing scale bound can cover ordinary evaluation
roundoff, but using original f(s) dot u(s) is a small alignment with the stated
contract; no demonstrated false-passivity example was found in this audit.

## Probe evidence

The read-only in-memory probe exited 0 in 0.2831634 s external tool elapsed
time (0.0057168 s inside Python), below the 60 s component limit. It reproduced
the event and input-validation findings without modifying source/tests. Exact
observations are retained in
`logs/process_hht_20261003/hht_path_review_probe.log`. No physical run occurred.
