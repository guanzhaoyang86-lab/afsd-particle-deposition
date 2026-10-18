# Cancellation-aware visibility certificate, independent version 1

Registered before implementation or qualification on 2026-10-04. This is a new
geometry certificate, not a change to material, contact, thermal, support-graph,
or trajectory mechanics. Frozen `visible_support_path.py`, its failure at
50 ms, and every existing accepted state remain unchanged. A saved-candidate
replay is diagnostic only; a new process controller must re-solve from a real
accepted parent and bind this contract, helper, qualification, and dependencies.

## Declared path and exact bound

Retain the inherited linear world particle/origin chords and shortest SO(3)
tool path. This contract was refined BEFORE implementation and test execution
to differentiate the frozen Rodrigues polynomial directly, avoiding any
assumption that its floating-point axis or R0 is exactly orthogonal.
Write row vectors `r(s)=r0+s*v`, `v=(x1-x0)-(o1-o0)`. Reconstruct the frozen
axis extraction without renormalizing it, and the identical stored `K2=K@K`.
Treat their stored entries as coefficients of the real reference polynomial:
`R(s)=R0*(I+sin(theta*s)*K+(1-cos(theta*s))*K2)`.
Then `R'=R0*theta*(cos(theta*s)*K+sin(theta*s)*K2)` and
`R''=R0*theta^2*(-sin(theta*s)*K+cos(theta*s)*K2)`.
For `q=rR`, `q'=vR+rR'`, `q''=2vR'+rR''`.
Set `C >= ||R0||F*(||K||F+||K2||F)` and
`rho=max(||r0||,||r1||)`. Submultiplicativity and `|sin|,|cos|<=1` imply
`B=C*(2*theta*||v||+theta^2*rho) >= sup ||q''||`. The speed bound is
`L=||q'(1/2)||+B/2 >= sup_{0<=s<=1} ||q'(s)||`.
This keeps cancellation inside `q'`; its error allowance must not be relative
to the possibly tiny result of that cancellation. An edge uses the maximum of
the endpoint bounds. An interval of width h stays within `h*L_edge/2` of its
mid-time straight segment for every edge parameter.

Use this new formula only for `1e-6 <= theta <= pi/2`. Outside this range use
the inherited bound, including its original arithmetic margin. The exact-pi
ambiguity rejection is inherited. No new angle or geometry acceptance is added.

## Floating-point engineering domain and conservative enlargement

The inherited input validation (finite arrays, proper rotations with defect
at most 1e-12, bounded geometric arithmetic margin) and all partition numerical
domain rejections remain mandatory. This remains an explicitly limited IEEE
binary64 engineering enclosure, **not a formal interval-arithmetic proof**.
Assume round-to-nearest binary64 arithmetic and NumPy sin/cos/norm evaluation
error enclosed by a generous 8192-epsilon allowance on absolute input scale;
no assumption of relative accuracy is made for subtracting large terms.

Inflate each computed Frobenius/vector norm by `1+eta`, where
`eta=8192*eps`, and then `nextafter(+inf)`. Inflate C after multiplication.
This directly covers permitted defects in R0 and in K/K2; no exact SO(3)
identity for the stored matrices is assumed. Nonfinite computations reject.
For rounding estimates also use the outward bound
`A=||R0||F*(sqrt(3)+||K||F+2*||K2||F)` for `sup ||R(s)||`.
With S equal to the inherited absolute world/geometry coordinate scale,
bound the evaluated path's uniform positional error by
`E=eta*A*(S+max(rho))`. Add **twice** E to the inherited displacement margin,
covering the evaluated moving point and its midpoint reference. This can only
enlarge the old arithmetic clearance. No rotated endpoint gets snapped.

Evaluate q' at the midpoint using the separate terms `v@Rmid` and
`rmid@Rprime_mid`. Its added absolute derivative allowance is
`eta*(A*(S+||v||)+theta*C*(S+rho))`. This bounds cancellation errors using
both large multiplication scales, not the residual derivative norm. Inflate
the derivative-plus-curvature-plus-error result again by `1+eta` and
`nextafter`. The scale includes absolute world positions and origins, so
subtraction cancellation is not hidden. Curvature B uses outward norms and
the same outward multiplication allowance. The selected speed
is `min(inherited_speed,new_speed)` only when the above domain is met, with
the enlarged displacement margin always retained in that case. The former
speed is already independently valid; choosing the smaller of two bounds is
valid with the common enlarged margin. Every nonfinite bound rejects.

Do not silently relax the fixed limits: 256 combined strict/enclosure batch
queries, subdivision depth 32, zero strict positive-intersection tolerance,
all retained edges, no extra graph deletions, identical degenerate-edge and
partition numerical-domain fail-closed behavior. A certificate cannot be
obtained by sampling alone. FOUND, UNRESOLVED, and CERTIFIED remain distinct.

## Qualification registered before results

New helper and independent test file must statically adapt the frozen files;
no monkeypatch of imported global functions. Reuse all old path fixtures and
the numerical-domain fixtures with their calls explicitly targeting the new
helper. Add analytic co-rotating-world-chord, translation/stationary, generic
rotation/world-coordinate covariance, small orthogonality defect, midpoint
cancellation/error-scale, fallback-angle, and purity/finite-bound checks.
Sampled analytic derivative comparisons are regression checks accompanying
the derivation above, not a proof obtained from dense sampling. Existing hidden
intersections, insufficient query/depth budgets, degeneracy, and numerical
domain failures must still reject. Do not claim material/physics validation.

After components pass, replay the original 49.9375--50 ms failed candidate's
entire supplied graph (7789 edges) against the new certificate, at the same
256/32 limits. Verify the original input-array, graph, candidate, and source
hashes and preserve `restart_allowed=false`. Output an independent diagnostic
JSON whether CERTIFIED, FOUND, or UNRESOLVED. Do not split or drop edges, modify
the original result, or promote the failed candidate to an accepted checkpoint.

Each component or saved-geometry worker has an external limit of 120 seconds;
protocol, exact command, exit code, log, elapsed time, before/after hashes and
summary must be recorded. Qualification binds this contract, new helper and
test/runner source, frozen fixture sources, partitioner, numerical guard,
rotation source, and frozen path source. No force, NNPS, material, thermal, or
real process evaluations are permitted in this stage. The summary is suitable
for a subsequent independent process preflight to bind, not automatic permission
to increase dt or claim global temporal/physical convergence.
