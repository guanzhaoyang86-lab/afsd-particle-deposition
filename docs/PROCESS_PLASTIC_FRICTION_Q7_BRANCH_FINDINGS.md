# Q7 first rejected branch: bounded independent diagnosis

Evidence: `outputs/process_plastic_friction_q7_20261001/branch_diagnostic_12_14_layout/diagnostic.json`
and `logs/process_plastic_friction_q7_20261001/branch_diagnostic_12_14_layout.log`.
Preregistration: `PROCESS_PLASTIC_FRICTION_Q7_BRANCH_DIAGNOSTIC.md`.

The input remains the original BLOCKED `n3_i4` run, accepted t=.05 s/step 6,
target candidate t=.0625 s. No state was accepted or physical path advanced.
The independent saved-data reader passed; every source and artifact hash was
verified before rebuilding the diagnostic instance. Every raw particle and
reference property, block array, work ledger and metadata was unchanged after
each pure branch solve. No production restart guard was weakened or called.

## Actual observations

| fixed Newton budget | result | final relative root | wall time |
| --- | --- | ---: | ---: |
| 12, original | BLOCKED; exact replay of all 19 rejected arrays and full trace | 1.1528835562823164e-6 | .21993 s |
| 14, preregistered | root returned at iteration 13; NOT_ACCEPTED | 1.637588335220718e-10 | see diagnostic.json |

The original relative root threshold is unchanged at 1e-8. The conditioned
reduced tangent was about 15.53 to 27.79, far below the existing 1e12 rejection
limit. The final additional correction used a full line-search step. The full
diagnostic process took .92902 s; the budget was 60 s.

The converged branch is **not an admissible contact interval**. Particle ID 21,
feature 1, remains forced sticking on this branch, but its current circular
friction margin is -5.691232925791425 N. Its normal force is 349.3139740022227 N,
tangent force (38.288546133237766,-13.57148958688364) N, and slip is only
(4.24e-22,-8.47e-22) m because of the imposed sticking branch. Its all-contact
natural-map residual is 5.691232925791427 N, well above 9e-5 N. Existing event
splitting must locate the contact transition before anything is committed.

For released rows, the maximum natural-map residual is 1.4737061416745039e-6 N,
below 9e-5 N; vector flow opposition passes. No moving row begins off its old
cone boundary or has old/new tangent forces in opposite hemispheres. The
minimum raw candidate friction work is -6.402578620333036e-20 J, retained as
roundoff-scale data. That tiny negative value is not the 5.69 N cone violation
and is not the first stopping condition. It has not been clipped or turned into
heat. Normal detachment was not observed in this candidate.

## Root-cause classification

- **CONFIRMED:** the first thrown exception is exhaustion of the original
  12-correction branch budget. One extra correction finds that branch root at
  unchanged tolerance, physical parameters, state and initial guess.
- **CONFIRMED:** at the target endpoint the observed forced-stick branch is
  outside ID 21's friction disk. Root success must return to event handling;
  it must not be promoted directly to an accepted interval or path completion.
- **RULED_OUT for this first exception:** a singular reduced tangent, NaN,
  negative friction-work stopping condition, material-state accumulation during
  pure probes, or loss of the accepted checkpoint/particle correspondence.
- **SUPPORTED:** an explicitly bounded Q7-only branch budget adjustment may
  allow the existing event search to continue. This is not a demonstrated full
  remedy: the event search and later branches have not been executed here.
- **POSSIBLE/NOT_RUN:** later event-location, changing slip direction, plastic
  tangent or path-dependent material-work issues. This diagnosis cannot rule
  them out. Old finite-strain-energy and spatial failures remain evidence.
- **BLOCKED:** full .2 s Q7 completion and same-terminal increment comparison
  remain unsupported. No three-condition AFSD production run is authorized.

## Minimal next change to consider, not implemented by this diagnosis

Keep Q6 and production code untouched. In a separately named/evidenced Q7
companion variant, make the branch correction budget an explicit 14 rather
than changing force/slip/energy tolerances or material parameters. Preserve the
original 12-correction failure and source snapshot before any source change.
Add a regression requiring the same accepted-state replay, a converged branch
root, and rejection of direct commit while ID 21 violates the cone. Only then
exercise the existing event search with a bounded actual short attempt; do not
increase iteration limits again in a loop until a desired result appears.

The first diagnostic invocation stopped before any Newton solve because raw
history restoration changed Fortran-layout contact/support arrays to C layout.
That changed reductions by at most 4.55e-13 N and broke strict record equality.
The diagnostic loader now restores the actual `zeros_like(X)` force layout;
the complete accepted record and original rejected state then replay exactly.
Both the initial error log and read-only inspection log remain preserved. This
was an evidence-loader correction, not a physics or tolerance change.
