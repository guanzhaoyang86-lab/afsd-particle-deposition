# Onset negative-normal diagnostic result (2026-10-02)

Preregistration: PROCESS_PLASTIC_FRICTION_ONSET_FAILURE_DIAGNOSTIC.md, unchanged
after execution; its hash is stored in diagnostic.json. New diagnostic source:
scripts/diagnose_process_plastic_friction_onset.py. No solver/source repair or
second physical path run occurred in this diagnostic.

## Actual execution: diagnostic PASS, physical path still BLOCKED

Command, wrapped in single-thread subprocess.run(timeout=60):

```
.envs/pysph/python.exe -B scripts/diagnose_process_plastic_friction_onset.py --run outputs/process_plastic_friction_onset_20261002/n3_i4 --output outputs/process_plastic_friction_onset_20261002/negative_normal_diagnostic
```

Exit0; outer .4162345 s, inner .0173560 s. Exactly one branch and one actual
material trial were evaluated. Accepted record reconstructed exactly; every
accepted particle/reference/array/metadata signature is identical before/after.
All saved run source and input artifact hashes match. Production restart remains
forbidden; both commit and advance were replaced by diagnostic-only exceptions.

The first negative-normal raw candidate occurs at TARGET .125 s, branch source
line182, iteration0, before the first Newton residual/iteration record or any
line search. Accepted state remains step16 at .1125 s. Contact particle ID19,
feature1, has normal force -6.9721105498611475 N; its accepted normal was
+212.14481166443747 N. This is not a floating-point tail. All nine prescribed
normal gaps equal exactly zero. The arithmetic assigned-contact free mechanical
balance is .1853146131 relative to the force scale; this is not a cone root or
an accepted equilibrated state. The cone response rejects its domain before a
contact-law residual can be evaluated. No friction work or accepted heat follows.

## What is established, and what is not

- CONFIRMED: the normal-domain guard is hit by the initial target displacement
  guess (all forced-stick tangent displacements imposed at target, remaining
  free displacements still at the accepted state), not by a converged root.
- CONFIRMED: no Newton correction was attempted; increasing its iteration budget
  alone cannot affect this particular first failure.
- SUPPORTED: an admissible initial-guess/continuation issue must be examined
  before declaring physical opening. This does not prove a positive-normal root
  exists at the target or that a smaller request interval solves the path.
- NOT_ESTABLISHED: true normal separation, release/recontact feasibility, or
  an exact physical detachment time. Fixed normal gap here is a constraint.
- NOT_RUN: alternative predictor, interval subdivision, normal active-set repair,
  further accepted steps, 8-increment path, thermal or actual AFSD process.

Evidence: outputs/process_plastic_friction_onset_20261002/negative_normal_diagnostic/
diagnostic.json and REJECTED_first_negative_normal.npz; log:
logs/process_plastic_friction_onset_20261002/negative_normal_diagnostic.log.
The candidate carries REJECTED_DIAGNOSTIC_NOT_RESTART, original particle IDs,
raw material/internal-force arrays, signed contact normal/tangent forces, mask,
prescribed gaps and explicitly invalid-candidate force-balance bookkeeping.
It must not be plotted/exported as an accepted physical endpoint.
