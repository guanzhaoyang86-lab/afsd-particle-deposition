# Actual joint HHT local benchmark execution gates — before driver runs

The selected exact geometry/physics and finest accuracy bounds are frozen in
PROCESS_HHT_FIXTURE_CONTRACT. Shared reference fixture_data() supplies ONLY
physical input arrays; driver wall/contact computation is independently written.
The wall normal stroke is the registered cosine (not the earlier informal
triangular suggestion). This is a 27-particle, all81-free-DOF LOCAL_BENCHMARK,
not the existing AFSD checkpoint, deposition or a plastic-flow demonstration.

Driver uses current PySPH 3D NNPS/corrected paired mechanics and unmodified FKS
material update, current density and BE heat. Only this benchmark's elastic
branch is admitted; positive-plastic-work accuracy is NOT_TESTED. Endpoint-linear
contact traction is named quadrature approximation. The full quadratic path
must be passive and free of interior zero slip; additionally sampled Coulomb
traction defect <=.01 of capacity. Never discard negative Q or add residual heat.
The rigid-wall velocity jump at1us splits every grid and explicitly initializes
right-limit acceleration/force history while preserving physical fields.

Prerequisites: independent reference PASS, smooth linear HHT tests PASS, dense
path tests PASS (including independently found double-touch zero regression),
and actual coupled component tests PASS. All source hashes must remain fixed
during each run. Preserve initial missing-SciPy and every other failure.

Execution sequence:

1. One coarsest250ns pilot step. Prior 27-point material+heat evaluation measured
   .001788s, but fulljoint cost unmeasured. Allow60s internal/70s external.
   This is a cost/integration probe; not a common-end-time accuracy result.
2. Only if pilot accepted and 2*(8+16+32)*max(pilotstep, .01)+10<=50s,
   run all3 fixed grids (250/125/62.5ns) to identical2us, internal55/external60s.
   Otherwise mark full sweep NOT_RUN, no timing extension or retuned grid.
3. Save checkpoint and actual3D VTK at every accepted stage. State contains
   acceleration and true force/boundary history. Source/config and complete
   NPZ roundtrip identity must pass before publication. Nonlinear candidates
   and REJECTED diagnostics never become physical output or restart points.
4. Compare all common times to independent RK4 reference, never different end
   times. Apply finest bounds and decrease criteria from FIXTURE_CONTRACT.
   No tuning physical values/thresholds after errors are observed.

Local extra physical-domain checks registered here: displacement from initial
position <=5um, absolute relative density change<=.005, normal overlap<=5um,
temperature in original solid domain, finite every output field; mass/ID exact.
These are benchmark validity limits, not relaxation of any AFSD continuation
guard. Node kinetic identity <=1e-12+1e-9 sum(|physical work|+|filter work|) J;
node momentum identity <=1e-12+1e-9 sum(|physical impulse|) Ns; thermal ledger
<=1e-12+1e-9 sum(|sources|) J. Nonzero algorithmic force-filter terms are signed
numerical terms and not wall reactions/material heat. Current hypoelastic
energy and normal temporal work residual remain explicitly diagnostic.

One same-source resumed fine-grid trace from its stored halfway checkpoint
may be recomputed only after a successful sweep and measured forecast<=50s,
internal55/external60s. Compare full solved physical and acceleration/force
history, not a storage-only roundtrip; raw bitwise or max differences recorded.
This does not authorize external AFSD restart or process-scale steps.

Before the first pilot, independent read-only driver audit added enforcement
of already registered bounds: compare BOTH1us wall-force limits; analytically
minimize signed quadratic slip; bound cubic particle normal motion against
the complete cosine stroke range; recompute BE balance from the actual saved
candidate temperature. Checkpoint readback and unchanged source hashes precede
atomic accepted_history publication; unpublished candidate files never count
as accepted restart points. Report the low-rate FKS yield margin and nonzero
temperature increment. These are diagnostic/acceptance checks, no physical
parameter change. Coupled a=1 control remains NOT_RUN; smooth a=1 is tested.
