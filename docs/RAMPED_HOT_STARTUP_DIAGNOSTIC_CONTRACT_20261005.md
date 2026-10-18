# Independent smooth-start feasibility experiment

This is an ENGINEERING_APPROXIMATE boundary-history sensitivity experiment,
not a restart of any accepted process, not deposition acceptance and not a
paper startup reproduction. The sole backend remains pinned PySPH. Independently
create the same unloaded 918-particle, 8 mm substrate, 723.15 K state at time zero.
No old evolved particle state is consumed. Material, low-rate floor, thermal,
contact, free axes, mass, modulus and nominal final speeds stay unchanged.

Nominal speeds are multiplied by 3*s^2-2*s^3 for s=t/0.020 in [0,1], then by one.
Endpoint translation, rotation and prescribed top feed use the exact integral
R(t)=T*(s^3-s^4/2), then t-T/2. This changes physical loading history; it is not
time scaling. Within each numerical step the model explicitly uses linear
translation, constant-angle geodesic rotation and linear actuator displacement
between those integrated endpoints. This is a piecewise approximation to the
smooth path. Contact, tangent, thermal contact, visibility, continuous path
certificates, sampled geometry and crossing observer use that identical step
motion. No old module globals are rebound. The inherited material map uses
total equivalent strain rate, retains its 1/s floor and unresolved provenance.

One predeclared schedule only: four chunks of 1 ms; each coarse chunk has two
500 us steps, each fine chunk four 250 us steps. Both paths start independently
at the same time-zero state and continue their OWN complete histories. The
comparison is COMMON_INITIAL_CONTINUOUS_TWO_RESOLUTION, not a claim that later
chunks share the same immediate parent. At every shared 500 us endpoint use
the unchanged global limits and the unchanged local pointwise limits, applied
conservatively to all 918 particle IDs. Compare crossing ID sets exactly.
Stop on any physical or temporal failure. No automatic halving or fourth path.

At most 24 roots plus the initial cost/qualification stage. Every actual stage
has forecast <=450 s, internal watchdog 540 s and external timeout 600 s.
Single BLAS/OpenMP threads. Before each root stage use
30+1.5*N*(max historical/experiment complete-candidate or stage-per-candidate
time + 1 s), retaining the actual expanded-window maximum as a cost floor.
The cost stage uses existing map/tangent/rank measurements plus an explicit
120 s allowance, not an unmeasured root forecast. Before every stage require
elapsed supervised experiment time plus its forecast <=1200 s total.
Independent saved-only readback must pass each completed chunk before the next.
The initial geometry and elastic rank are qualified at t=0, where the old and
new poses coincide exactly; this does not qualify later deformed states.

Every candidate stays accepted=false, restart_allowed=false and diagnostic_only.
Complete prefix energy gates, fresh force evaluation, scalar/batched equivalence,
penetration, density, material/temperature domain and continuous geometry
certificates remain unchanged. No stress/temperature/dissipation clipping or
residual heating. All failures and partial files remain evidence.

A loaded feasibility PASS additionally requires both paths to have real non-BASE
contact force above their unchanged force tolerance and new plastic strain,
then at least one further complete shared 500 us interval with physical and
temporal gates passing. Merely reaching a low-load startup endpoint is not a
loaded PASS. A completed 4 ms study still cannot establish a full-speed usable
step, deposition, mesh convergence, experimental agreement or production cost.
This experiment neither advances nor adds time to the accepted 342/918 states.
