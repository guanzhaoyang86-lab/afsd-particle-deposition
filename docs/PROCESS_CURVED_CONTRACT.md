# Curved-slip impulse-corrected HHT variant — registered before tests/runs

2026-10-03. Scope: named numerical variant within the existing PySPH backend,
not new physics, production or an AFSD result. Preserve every4221 previous
freeze entry except the three handoff/index documents. Prior failures stay.

Use identical PROCESS_SLOW_CONTRACT27particle physical problem and original
S0 schedule, material/temperature/geometry/5Hz walls. No mass/time/modulus
scaling, no clipping, label forces or free velocity prescription. Existing
material, corrected operator, PA/NNPS/checkpoint and implicit heat are reused.

## Numerical definition and non-negotiable accounting

Use the existing compatible HHT dense position/velocity path. At each actual
path point evaluate physical normal penalty and full-vector Coulomb traction,
not a linear chord between endpoint friction vectors. Define dJ and dW as
actual path-contact impulse and material work minus the endpoint-linear
contact impulse/work. Correct the acceleration equation by dJ/(gamma*h):

    M*A1 = a*F1 + (1-a)*F0 + dJ/(gamma*h) + R.

F0/F1 remain TRUE endpoint total physical forces; canonical A remains numerical
HHT acceleration, not relabeled F/M. Physical total impulse is Jlinear+dJ;
physical force work is Wlinear+dW. Old force-filter impulse stays unchanged;
its work acquires vbar dot dJ - dW. Root residual impulse is gamma*h*R.
Derivation and independent algebra tests must pass before any coupled run.
No correction/root/filter term becomes material heat. Actual friction heat is
the sum of signed -Jt_node dot u_node only; normal/contact wall work is separate.

This is NOT classical HHT stability certification. For smooth force dJ=O(h^3)
gives an O(h^2) acceleration correction; temporal convergence and cost still
need empirical proof. It reduces to the old method when contact force is
linear along the same path. Persist a NEW schema/integration-definition, never
silently reinterpret an old state. Any diagnostic migration explicitly preserves
all physical arrays and records old/new force/acceleration/history hashes.

## Contact path scope and verification

Fixed normals and stable face/particle identities only. The caller supplies a
source-bound conservative whole-path gap certificate and affine tangential
wall-velocity contract; endpoints/nodes must agree. For the slow planes, bound
the cubic material x extrema and analytic normal stroke range0..2e-7m.
No path with zero full-vector slip, opening or a changing contact normal is
accepted. Use full-vector polynomial zero check, not signed-z or.05m/s.
These events still need a later located event/restart method.

Quadrature8/16points, then32only if necessary; convergence requires each
impulse/work/heat component difference<=1e-14(in its SI unit)+1e-8relative.
Finite quadrature difference is an estimator, not a rigorous continuum bound.
Unconverged32point result rejects. The actual callback traction must satisfy
the unchanged1% Coulomb-law bound at nodes, positive normal force and consistent
penalty potential; no energy clipping. Actual constitutive traction has no
endpoint-chord interpolation; do not count its removal as a relaxed threshold.
For the claimed exact-Coulomb branch, additionally require callback traction
to equal the direct formula within512eps force-scale roundoff; otherwise
reject, never silently overwrite a different supplied law. The direct formula
is the force actually integrated. This supplies analytic pointwise passivity
on the certified positive-gap/nonzero-slip path, not a claim that arbitrary
callbacks are proven passive merely by finite sampling.

Independent tests: analytic straight and turning slip, moving/stationary walls,
Galilean work, full-vector zero vs component reversal, rejected bad traction,
full work/momentum identities with nonzero force history, and old-method
reduction. Root raw residual gates and maximum24 evaluations remain unchanged.
Temperature/domain/ID/mass/paired-power/overlap guards remain unchanged.

## Bounded execution order

1. New algebra/contact tests and all168old selected regressions. Record source
   hashes and no mid-test changes. Independent oracle remains untouched.
2. One private migration of savedS0stage008 and one500ns candidate, not a new
   accepted trajectory. Reproduce the original FAIL first; apply new full
   physical-law/work/thermal gates. Estimate<10s, external60/internal55s.
3. If that passes, fresh-source250ns cost pilot. Use the same original S0
   schedule only if2*28*measured_whole_step+10<=50s. Internal55/external60s.
   No smaller-dt retries, and never publish beyond independent5.9375us prefix.
4. Stop on first new branch/failure. S1/S2 and20us convergence require a fully
   certified common horizon;1ms/production/actualAFSD remain NOT_RUN here.

At most two evidence-based numerical implementation repairs; preserve first
sources/failures. A successful local prefix is not transport or project end.
