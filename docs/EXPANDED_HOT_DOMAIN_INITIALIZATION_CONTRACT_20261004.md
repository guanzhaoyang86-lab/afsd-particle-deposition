# Independent 6/8 mm hot-domain initialization

Registered before initialization checks or any new physical evaluation. This
adapter only builds independent, initially unloaded 450 C particle states and
inspects parameters, arrays and finite-tool distances. It does not evaluate
forces, material response, thermal increments, real NNPS, tangents or a root.
No process time is accepted. Prior 342-particle trajectories and failures are
neither changed nor used as parents; adding particles to their checkpoints is
explicitly outside this contract.

Use the pinned PySPH commit
`1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5`, existing ParticleArray construction,
and existing independent discrete-contact initial-state factory. Parameters are
read from the immutable contact-material resolved JSON, SHA256
`37a153e2df917070ed707a93abcda8709cc8468f82393fc6cc90e07071bbdfd7`.
Only the already declared unloaded square bore width of 2.504 mm and the new
substrate width (6 or 8 mm) override that file. dx=0.5 mm, bar width/length
2.5/3 mm, substrate depth 1.5 mm, physical velocities, material, thermal and
contact constants remain unchanged. Every parameter retains its value, unit
and source/assumption status. The inherited FKS dimension/low-rate and
Jaumann closed-cycle limitations are not resolved by initialization.

The complete three-dimensional bar and substrate are deformable. Rigid tool
geometry retains the finite 5.4 mm diameter shoulder, 1 mm tool height, finite
open rotating square bore and rigid base. All bottom-substrate axes are fixed;
only the z axis of the 25 initially topmost bar particles is prescribed by the
existing feed rule. x/y remain mechanical unknowns, including at the bar top.
The independent hot quasistatic state starts at rest, not by resetting any old
dynamic state. Passive material/tracer/whole-bar IDs use the original factory;
IDs are unique inside each independent domain, not shared trajectory identities
between different substrate widths. No inlet/replenishment is implemented.

`independent_initial_state` is a new named factory. The visible variant has
classification `EXPANDED_HOT_QS_VISIBLE_SUPPORT_APPROXIMATE_V1` and an explicit
`independent_domain_initialization` record. It has no parent checkpoint hash,
branch origin or fabricated `spatial_model_transition`. Its new trial adapter
initializes the existing discrete-contact base, then reuses the frozen visible
residual/tangent methods without calling the old transition function. The
original operator classification remains separately identified. A discrete
variant is also available as an independently initialized operator comparison;
it retains the existing discrete model classification and the new domain scope.

Both variants start at step/time zero, with empty interval-energy records,
zero stress/plastic/contact history and an independent bore observer. Only bar
centres initially above the shoulder are eligible; initially lower material is
not counted as throughput. Initial mass and sensible heat relative to 293.15 K
are recorded separately from initial zero kinetic/elastic/plastic storage.
Finite-patch distances permit a separately reported analytic normal-potential
calculation; no contact-force function is called. Floating-point shell overlap
at exact substrate/base tangency is reported and bounded by the explicit
64-epsilon geometry roundoff bound, never silently removed from that potential.

Initialization snapshots are `NOT_RESTARTABLE`, `accepted=False`, and lack
physical-process qualification. Uncomputed rate/scratch fields are identified
as unevaluated. The later cost/root controller must separately qualify actual
NNPS, filtered graph, moment rank/condition, scalar/batched full map, contact
events, complete force/thermal/material/energy response and time comparisons.
The old visible accepted-copy transition gate is inapplicable: a new controller
must explicitly bind this independent initialization and its own qualifications.
Constructing the diagnostic trial object alone confers none of those passes.

Expected populations: 6 mm has 432 substrate +150 bar particles and 1289 free
displacement DOF; 8 mm has 768+150 particles and 1961 free DOF. Report actual
arrays and float64 theoretical dense allocation sizes for D(9N,3N), full
J(3N,3N), free J/inverse and thermal(N,N), without allocating those matrices.
These sizes are not measured peak memory or runtime estimates.

Initial shoulder-edge margins are 0.3/1.3 mm; compact-kernel support radius is
1.3 mm. The 6 mm domain does not establish low boundary sensitivity. Advancing
tool coverage obeys W/2 >= R + traverse*t + buffer for the undeformed nominal
plate; it is only a geometric observation, not a validated duration. A paired
6/8 mm physical short-run domain comparison remains NOT_RUN, as do inlet,
spatial convergence, stable deposition, three feed cases and paper validation.
