# Feasible predictor independent audit (2026-10-02)

Status: PASS_INDEPENDENT_READBACK_NOT_PHYSICAL_ACCEPTANCE.
The returned target branch remains REJECTED: eight forced-stick contacts violate
the cone. Accepted state remains .1125 s / step 16, not .125 s.

## Scope and evidence

Read-only audit of the one completed old/new comparison in
`outputs/process_feasible_predictor_20261002/target_0125/` and the frozen accepted
input in `outputs/process_plastic_friction_onset_20261002/n3_i4/`.
The audit imported only Python stdlib and NumPy. It did not import or call the
material, contact, branch, assess_candidate, restore, or solver implementation.
No additional branch was run, and no physical state was changed.

Detailed independent output: 
`logs/process_feasible_predictor_20261002/independent_readback.log`.
Actual comparison log: `logs/process_feasible_predictor_20261002/actual_comparison.log`.
Primary audit inputs: `diagnostic.json`, `REJECTED_old_raw_0000.npz`, both
`REJECTED_*_final_candidate.npz`, and all six predictor raw NPZ candidates.
The reported input/source/diagnostic-source/output artifact hashes were freshly
recomputed and all matched. This audit does not independently reconstruct the
solver's full constitutive force from scratch.

## Independent calculations and results

Contact rows were located by persistent ID matching, not by assuming ID equals
an array index. IDs 1 through 27 are unique and unchanged; contact IDs 19 through
27 and feature IDs equal the frozen history. The only released row is ID 27.

For every raw candidate, independently read normal force as
`N = internal_force[contact_rows, 0]` and tangent force as
`tau = -internal_force[contact_rows, 1:]`. All prescribed support components are
exactly zero and all nine normal gaps `-4e-5 - u_x` are exactly zero. The predictor
step preserves all contact displacement components bit-for-bit, including the
released row held temporarily during initialization.

| Candidate | Minimum raw normal force (N) | Negative IDs |
| --- | ---: | --- |
| Old initial/final and predictor initial | -6.9721105498611475 | 19 |
| After one interior predictor correction | 50.70480928224773 | none |
| After first cone correction | 13.559468417211061 | none |
| After second cone correction | 14.579225350706638 | none |
| After third cone correction | 14.59915777178539 | none |
| Final after fourth cone correction | 14.599258321039724 | none |

The trace independently contains exactly one predictor and four cone correction
entries with an accepted line-search factor, all equal to 1. Total corrections
are 5, below the shared 14 limit; six raw trials are retained for this branch.
Old initial/final and new initial match in eight mechanics/material arrays;
new final matches its last raw trial in the same arrays: displacement, internal
force, stress, S, ep, plastic increment, return diagnostic increment and strain
rate. There is no array evidence of substituting a different material history.

The final force margins `0.1*N - norm(tau)` for IDs 19 through 27, in N, are:

```
-72.353127215975, -102.6905416911407, -20.997320686048447,
-122.35992557255031, -138.83385336851438, -125.08959422934903,
-16.979194021816305, -101.46422685038048, -0.0000014068161391378453
```

Thus IDs 19 through 26 all violate the forced-stick cone by much more than the
unchanged 9e-5 N tolerance. ID 27's margin is within that tolerance. This is not
an acceptable physical endpoint even though the selected branch equations have
a small root residual.

### Signed contact work

Wall displacement was evaluated by the declared linear reversal segment
`a(t)=2e-5-4e-5*(t-.1)/.05`, `wall=(-4e-5,.6*a,.8*a)`.
With old displacement/force taken from frozen history, independently computed
`slip = u_new-u_old-(wall_new-wall_old)` in the tangent plane and
`Q_i = -0.5*(tau_old+tau_new) dot slip_i`. No abs, max, clipping or residual heat
was applied. Results for IDs 19 through 27, in J:

```
-2.516069288894096e-20, -4.8472077452881056e-20,
-1.2095558691504962e-21, -1.3872744137579997e-19,
-0.0, -0.0, -0.0, -1.151519985255503e-19,
 1.6616717278334899e-4
```

Sum = 1.6616717278334866e-4 J; minimum = -1.3872744137579997e-19 J.
The small signed negative terms are retained. Agreement with the saved diagnostic
is within 1e-18 J. None of these rejected-candidate values is accepted heat.
Nonnegative total work is not sufficient to accept over-cone sticking rows.

### Independent residual construction

Constructed contact force from raw N/tau, support reaction from the fixed
component mask, and residual `F_internal + F_contact + F_support`. For the
released row only, replaced tangent equations by the independently implemented
disk natural-map relation: `q=tau-7e8*slip`, radial projection of q onto radius
`0.1*N`, minus tau. This projection evaluates a root equation; it does not
modify any actual force. The unknown mask excludes support, all contact normals,
and the eight forced-stick tangent rows. Force scale is `90e6*.01^2 = 9000 N`.

- Mechanical free residual / scale: 3.797257399141952e-11.
- Interior residual / scale: 3.7972573991419525e-11.
- Selected coupled-branch residual / scale: 2.924011096486816e-10,
  below the existing 1e-8 threshold and matching the saved assessment.
- Released-row law residual norm: 2.60932479703236e-6 N, below 9e-5 N.
- Maximum law residual over ALL contact rows: 138.83385336851438 N,
  grossly above 9e-5 N because the forced-stick assumptions are no longer valid.

This independently confirms ROOT_RETURNED_EVENT_SPLIT_REQUIRED_NOT_ACCEPTED,
not an accepted contact solution or completed loading path.

## Accepted state and limitations

For both comparisons, the complete captured accepted-state fingerprint
dictionaries before/after are exactly equal, as are the accepted records.
Displacement, stress, normal force and tangent force fingerprint hashes, dtypes
and shapes were independently recomputed from accepted history and matched.
All frozen input files, including checkpoint and accepted journal, retain their
recorded hashes. The postprocessor cannot retrospectively inspect process RAM;
its in-process immutability evidence is the saved before/after fingerprints,
supported by unchanged input bytes and selected independently hashed arrays.

CONFIRMED: negative normal at the original target guess was avoidable through
an unchanged-material interior predictor; a positive-normal selected branch
exists at .125 s. NOT_ESTABLISHED: a valid full contact endpoint at .125 s,
event locations, event-resolved path continuation, normal separation/recontact,
or full .2 s completion. These checks are NOT_RUN in this audit. Existing J2
energy warnings and prior failed tests are not resolved or reclassified here.
