# Angle-aware endpoint enclosure for the unchanged Rodrigues path

This is an explicit supplemental geometric certificate. All frozen sources and
failed candidates remain unchanged. It does not change material, mechanics,
thermal coupling, particle positions, support masks, or event timing. A saved
chord replay is diagnostic only and cannot become a continuation parent.

The inherited `visible_support_prism_20261010.py` certificate uses the exact
reference polynomial with stored coefficients:

`R(s) = R0 [I + sin(theta*s) K + (1-cos(theta*s)) K2]`, `K2 = K@K`.

For the inherited domain `1e-6 <= theta <= pi/2` and `0 <= s <= 1`,
`|sin(theta*s)| <= theta`, `|1-cos(theta*s)| <= theta^2/2`. Consequently

`A_theta = ||R0||F (sqrt(3) + theta ||K||F + theta^2 ||K2||F / 2)`

bounds the same rotation polynomial. The inherited angle-independent bound
uses 1 and 2 in place of theta and theta squared over 2. The derivative
coefficient `C = ||R0||F (||K||F + ||K2||F)` and the inherited second derivative
bound are retained unchanged; no exact orthogonality identity is assumed.

Retain the inherited binary64 engineering assumption `eta = 8192*eps`; this
is not formal interval arithmetic. With inherited absolute world scale S and
`rho = max(||x0-o0||, ||x1-o1||)` over particles, the supplemental uniform
evaluation error is

`E = eta*A_theta*(S+rho) + eta*C*(1+theta)*rho`.

The first term covers world-coordinate arithmetic, dot/matrix evaluation,
and subtraction in `1-cos`. The second term explicitly covers absolute
sin/cos errors of eta and argument error at most eta*theta, since sine and
cosine are 1-Lipschitz. Thus no smaller transcendental-error assumption is
introduced. Norms, sums, products and the final result are outward enlarged by
`1+eta` and `nextafter(+inf)`, also covering mixed rounding terms.

The supplemental padding is the original displacement margin plus `2*E`.
The original margin retains pose-orthogonality and endpoint-polynomial
mismatch allowances. On each interval, inherited curvature `B*h^2/8` plus
this padding encloses every moving endpoint around its endpoint chord. The
four-corner convex hull therefore encloses every edge parameter. The original
four signed-horizontal/z wall rectangles and outward SAT projections prove
separation from a superset of the tool solid. Exact touching does not pass.

The old prism is attempted first. Only if it cannot prove separation may the
supplemental enclosure be tried, within the same prism batch. Strict checks,
partition numerical guards, old inflated-solid checks, degenerate-edge
rejection, 256 total batches, depth 32, and every original outcome remain
unchanged. Outside the inherited angle/curvature domain this proof is
unavailable. No sampled geometry alone can authorize a path.

Qualification requires all original prism tests retargeted to the new helper,
explicit arithmetic/angle-domain checks, and exact input-hash replay of the
saved 2026-10-15 v4 after-event failure. A real process run must re-solve the
complete coarse/fine pair from a previously verified complete checkpoint;
the original physical and pair gates remain required.
