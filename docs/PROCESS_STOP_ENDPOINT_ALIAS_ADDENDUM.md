# Endpoint-roundoff alias certificate — numerical repair 2, final

Registered 2026-10-03 BEFORE implementation, tests or the last reference run.
Repair1 source, report, all candidate and accepted files are preserved.
This correction is NOT a change to the full-zero bound, physical state,
Coulomb law, raw heat, event time or any shared-reference uncertainty bound.

Saved fine-grid evidence: after separately locating both initial stops,
PID4 has positive outgoing z-slip1.2268e-14m/s at the second event state. A
y-component root in the first dense interval has positive z-slip2.093e-14m/s,
positive norm-squared derivative9.60e-10m2/s3 and positive outgoing projected
acceleration22940.3948m/s2. The full cubic has no z-root or norm-stationary
root in (0,1). This is a false full-zero report from testing a component
root against the full-vector roundoff bound near its already located zero.

For a candidate at fraction s with full vector norm <=the ORIGINAL512eps
bound, construct the scalar outgoing projection of the entire cubic and
its quadratic time derivative on the dense-start-to-candidate interval.
Evaluate endpoints and every real derivative/extremum root in that interval.
Require nonnegative minimum projected slip, strictly positive projected
slip at the candidate, and strictly positive minimum projected time
derivative a_min. Thus this interval cannot contain a new full-vector zero.
Also require that elapsed time since this contact's independently located
event is nonnegative and <=original_velocity_bound/a_min. This is the
kinematic extent of the original numerical-zero uncertainty, NOT a new
physical slip/time tolerance. Save the complete certificate in alias logs.

For the observed case, elapsed9.1225e-19s is below the derivative-derived
4.95575e-18s window. Keep raw slip and signed power unchanged. A negative or
nonmonotone projected branch, a zero outside that local window, an unlocated
contact or an uncertifiable polynomial remains rejected. In particular the
coarse later candidate near6.56576us has static requirement.49424N below
capacity2.70526N: this is a separate feasible-stick branch and stays BLOCKED.

Before the final invocation add literal cubic tests: tiny-positive-z/y-root
alias PASS certificate; actual z-root, decreasing return, negative excursion,
outside-window and large-transverse cases receive no alias certificate.
Original scalar/algebra1e-12 SI and all nine reference bounds stay fixed.
One final combined reference invocation <=55s internal /60s external,
predicted<20s; stop at next unsupported branch. A common6.5 or6.5625us prefix
may PASS while8us remains BLOCKED. No further source repair this turn.
