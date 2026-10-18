# Independent 8 mm initial elastic rank and equal-map qualification

This is a new bounded diagnostic for the independently initialized 918-particle,
8 mm, 450 C engineering approximation. It accepts no process time and creates
no restartable checkpoint. Neither a 342-particle parent nor an invented spatial
transition is admissible. The two immutable inputs are the visible and discrete
`initial_8mm_*_v2` snapshots recorded by the separate initialization contract.
Their arrays must agree exactly; only their declared operator classifications
may differ. All inherited material/thermal/contact assumptions remain explicit.

Reference interval: 15.625 us. The initial rank check uses the actual time-zero
PySPH NNPS graph and its real kernel gradients, followed by the existing finite
tool visibility filter at time zero. This registered initial comparison requires
zero removed edges. Thus one corrected-gradient map D, in 1/m, is shared by the
two operators; the test fails if that identity is absent. Local moment matrices
are dimensionless. Every moment must have positive eigenvalues and condition
number at most 1e6, the existing correction-domain requirement. Save all local
moments, eigenvalues, corrections, the graph and physical free-axis mask.

At zero stress/strain rate the unchanged material tangent must equal the
isotropic elastic derivative. With volume V in m^3, shear G and bulk B in Pa,
the strain energy is sum V*(G*dev(sym(Du)):dev(sym(Du)) + B*tr(Du)^2/2).
Its Hessian K in N/m is the negative original frozen internal-force derivative.
Restrict K to the physical free axes: all bottom substrate axes fixed, and only
bar-top z prescribed. No inertia, mass shift, pseudoinverse, extra axis fixing,
contact spring regularizer, or matrix symmetrization is permitted. Verify
symmetry to 1e-12 of its maximum entry before the symmetric eigenproblem; save
the unmodified free matrix, full eigenvalue spectrum and least eigenvector.
The existing rank criterion is lambda_min > 1e-12*lambda_max > 0. Verify the
least eigenpair residual relative to the matrix infinity norm at 1e-10 and
its quadratic form against the explicit strain energy at 1e-9 relative scale.
These are numerical consistency tests, not material stability over future flow.

Check the unloaded initial contact at fixed time-zero tool poses using the
original scalar and batched contact routines: max force <=1e-10 N and normal
potential <=1e-20 J, as in the old initial equilibrium audit. This does not
claim thermal equilibrium. There is no artificial zero-temperature or rate
substitution. The rank derivative uses the real zero-rate elastic branch.

For the full-map identity, use the actual unchanged prescribed-axis displacement
guess at 15.625 us, with the real rotating/translating tool endpoint. Evaluate
independent discrete-batched and visible-batched trials. Require exact graph,
gradients, corrections, volume, material rate, internal force, material feedback
temperature, pair heat conductance/transfer and heat ledger. Apply the existing
`run_batched_process.candidate_equivalence` field and raw-energy tolerances,
and the existing nested-contact/material/thermal comparison. Compare their
actual unmodified search tangents exactly. A third scalar-contact visible map
must match the batched visible map under the same criteria. Save complete raw
fields, contact/energy/thermal journals and graph arrays, including failures.
This is a same-input map check, not a solved step or temporal convergence test.

Both initial states and all initial particle arrays must remain unmodified.
Tracer/whole-bar labels remain the frozen passive identity fields; no tags enter
the new rank, force or thermal expressions. The qualification binds all backend
sources, adapter/helper/contract, both initial snapshots and initialization
receipts, and completed 15.625 us cost evidence. The source map is checked before
and after the worker. Root, pairwise temporal, contact-path, domain, inlet and
paper qualifications are separate; this receipt never replaces their gates.

One single-thread worker has a 540 s internal / 600 s external hard bound. The
pre-run estimate is 30 + 1.5*(3*measured_full_map_s + 2*measured_tangent_linear_s)
+120 seconds for the 1961-DOF elastic assembly/eigenproblem and evidence I/O;
this conservative planning allowance is not a measured eigenproblem runtime.
It must be <=450 s. The theoretical D allocation is 182028384 bytes, full J
60676128 bytes, each free K/eigenvector matrix 30764168 bytes, and each 918^2
thermal matrix 6741792 bytes. Peak memory and actual timing are separately
reported; these sizes do not promise a particular implementation peak.

Success is `PASS_EXPANDED_INITIAL_ELASTIC_RANK_AND_SAME_MAP_NOT_PHYSICAL_VALIDATION`
with `qualification_passed=true`, `accepted=false`, `restart_allowed=false`, and
zero added accepted time. A failed test preserves all available diagnostics and
does not issue the success status. The later `expanded_hot_acceptance` adapter
requires a hash-bound successful receipt and calls the unchanged
`discrete_process.accepted_copy`, including its fresh physical evaluation,
raw/all-prefix energy and passive observer. Its returned candidate stays
unpublished/nonrestartable; only a separate process controller can publish after
all its root, geometry, source and same-parent coarse/fine checks pass.
