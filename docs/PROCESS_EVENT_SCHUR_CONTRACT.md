# Event-response candidate Schur contract (2026-10-03)

Scope: a new unaccepted-candidate preconditioner only. The existing physical
contact kick, five-field map, raw tolerances, material/thermal operators and
accepted-state commit remain unchanged. This component supplies no new
accepted state, restart permission, real AFSD advancement or full contact
validation. ENGINEERING_APPROX_CONTRACT and old failures remain applicable.

## Frozen linearization and elimination

All particle arrays retain accepted-state persistent-ID order. Let D be the
accepted-geometry corrected velocity-gradient operator, E=trace(D), V=m/rho
the accepted volume, rho_ref the specified reference density, and kappa the
unchanged physical bulk modulus. Define

    H = (1/2) E^T diag(V kappa/rho_ref),
    K = Kdev + E^T diag(V kappa) E.

Kdev is the frozen physical deviatoric elastic stiffness. The linearized
midpoint internal-force change before density elimination is

    deltaF = -(dt/2) Kdev delta_bar + H delta_rho.
    delta_rho = d_rho - dt diag(rho_ref) E delta_bar.

Thus deltaF=H d_rho-(dt/2)K delta_bar. In these expressions d=F(y)-y is
the direct defect of the actual five-field map, never a preconditioned defect.
Supplied per-particle responses are independently B=d(vbar_actual)/dF_internal
and C=d(v1_actual)/dF_internal, both in s/kg, at a documented contact branch.
They are not interchangeable, and B=C/2 is not assumed. The eliminated block is

    A delta_bar = rhs,
    A = I + (dt/2) B K,
    rhs = d_bar + B H d_rho.

Original actuator-prescribed OUTPUT rows of B and C must be exactly zero;
inconsistent input is rejected. Original prescribed rows of A are identity
and delta_bar_c=d_bar_c. The solver sets these rows exactly, shifts their
elastic contribution into the free RHS, and uses numpy.linalg.solve on the
actual generally nonsymmetric free block. B is never inverted; singular B
is valid. Singular A fails explicitly, without pseudoinverse, clipping,
regularization or an assumption of positive definiteness. Contact sticking
does not remove axes; unconstrained rigid modes remain represented.

Back-substitution uses the independently supplied B and C:

    delta_x   = d_x   + dt B deltaF,
    delta_v1  = d_v1  + C deltaF,
    delta_bar = d_bar + B deltaF,
    delta_rho = d_rho - dt diag(rho_ref) E delta_bar,
    delta_T   = d_T.

The solved delta_bar is returned, with its independent back-substitution
defect recorded. Position feedback is evaluated explicitly with B, including
its zero prescribed rows. A new candidate is y+delta; it requires evaluation
of the original physical map and all original acceptance conditions.

For free inertial response B=dt M^-1/2, C=dt M^-1, premultiplying the
Schur equation by M yields the old continuity preconditioner's
(M+dt^2 K/4) delta_bar=M d_bar+(dt/4)E^T diag(V kappa/rho_ref)d_rho.
The prescribed-row shift and all five corrections reduce to that old block.

Geometry, prestress, objective rotation, plasticity, thermal coupling,
contact geometry and changes of event response are explicitly omitted
derivatives. This is a local frozen linear preconditioner, not the full map
Jacobian. Branch validity belongs to the response producer and is preserved
in copied provenance.

## API and auditable linear evidence

`StageEventPreconditioner(accepted_pa,p,dt,time,B,C,provenance)` builds from
accepted geometry. `from_response_blocks(elastic,B,C,dt,free_axes,rho_reference,
provenance=None)` supports a detached FrozenElastic fixture or reuse of a
previously built accepted geometry. B and C have shape (N,3,3), with both
axes in Cartesian order and particles in FrozenElastic.particle_ids order.
Optional provenance.particle_id is checked against that order. Inputs and
outputs are copied, the physical ParticleArray is neither retained nor changed.

The full nonsymmetric A is exposed as readonly `.matrix`. Its SHA256 over
C-order float64 bytes, shape and infinity norm are in each diagnostic; a
driver must save that matrix once per construction when archiving a run.
Full d, rhs, solved delta_bar, deltaF, full/free residuals, prescribed-row
defect and back-substitution defect are returned for every application.
No algebraic constrained defect is claimed to be a physical reaction force.

Preregistered algebra gate: both full and free dimensionless infinity-norm
backward errors must be <=1e-12. The full error is

    ||A delta_bar-rhs||inf / (||A||inf ||delta_bar||inf + ||rhs||inf).

The zero-denominator case is zero only for the zero system RHS/solution.
Invalid/nonfinite inputs and singular systems fail explicitly. No candidate
density or temperature is clipped to meet a physical acceptance test.

## Independent baselines registered before execution

New tests: literal reversal B=19/54, C=7/9; a two-particle elastic spring with
nonsymmetric A and cross-particle response; nonzero original prescribed-axis
defect; independent coupled density/elastic affine-map root; old free-response
limit; singular B; singular A rejection; rigid translation/rotation retention;
input/provenance purity; passive-label invariance; persistent-ID permutation;
and exact preservation of zero five-field defects. These are algebra and
frozen linear-map baselines, not full nonlinear contact tests. They are not
generated by rerunning the new component as their oracle.

Execute only the new test file, with an external 60-second timeout. Maximum
two repair rounds; preserve first-run output and each later run separately
under logs/process_event_20261003/event_schur*. No real AFSD run is authorized
by this component contract. Test outcome is NOT_RUN until separately recorded.

## Initial execution evidence

2026-10-03: all 15 new tests PASS on the initial run, exit 0, unittest
0.046 s and externally measured 0.4042099 s under the registered 60 s timeout.
No repair round was used. Exact command and stdout/stderr are preserved in
`logs/process_event_20261003/event_schur_initial.log`. These outcomes cover
only the stated algebra/frozen linear-map scope; real nonlinear contact or
AFSD process advancement is NOT_RUN by this component's test execution.
