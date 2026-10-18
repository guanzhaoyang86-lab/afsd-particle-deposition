# Same-endpoint time check, registered after the 500 ns root

2026-10-04. The 500 ns root accepted in 15 evaluations /154.4919222 external
seconds. Its numerical convergence is insufficient evidence of time accuracy.
Run exactly one two-stage 250+250 ns branch, from the SAME original accepted
850.393632910904 ns parent to the SAME 1350.393632910904 ns endpoint. This does
not extend the latest physical time and is not a second deposition experiment.

Unchanged mass/modulus/time/FKS/geometry/thermal map/acceptance. At each stage,
reuse the 500 ns accepted journal only as a frozen numerical force/geometry
anchor, recomputing force secants with the current branch's actual v0, mass,
250 ns interval and actuator command. Its forces/normals are not used by the
physical trial. Every physical trial rebuilds its usual candidate NNPS/material/
contact/heat. Report this frozen-anchor approximation; no Jacobian claim.

At most14 actual root evaluations per stage,28 total; no retries or dt changes.
One process inner570/external600 s, initial launch screen540 s. Measured upper
planning cost from the 500 ns package is 9.0645269 s/full trial+proposal+journal
and42.5737043 s/all-contact response. Forecast
1.5*(28*9.0645269+2*42.5737043)+30 =538.6312427 s (<540 s).
Before each root, use its actually measured secant build cost and remaining
budget; stop if a conservative remaining estimate exceeds the internal budget.
No continuation after either stage fails. Accepted partial states are separate
branch outputs and do not overwrite the already accepted full-step evidence.

Save both stage checkpoints, complete journals, work/heat, trajectory fields,
tool poses and iteration histories. Readback compares persistent-ID-aligned
endpoint position, velocity, density, temperature, deviatoric stress, plastic
strain/work, friction heat, tool/actuator work and normal potential. Use absolute
differences plus characteristic scales (dx,feed,rho0,1 K,shear,unit plastic
strain); do not invent a post-hoc PASS tolerance. Two levels provide sensitivity,
not order-of-convergence, process accuracy or experimental validation. Negative
raw friction entries within original roundoff remain in the record.

The 50 um guard remains. Crossing/deposition cannot be inferred. No full project
regression, production, GPU, or uninterrupted long-run qualification is claimed.
