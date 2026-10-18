# Independent saved material-response action readback

2026-10-04. Scope: read already saved arrays and receipts only. No matrix build,
preconditioner application, material/contact/thermal evaluation, solver or
trajectory is permitted. This is an output/accounting review, not a new action.

## Registered checks before execution

- Source, artifact, persistent-ID, test-log and receipt hashes must match exactly;
  the same watched files must retain their hashes after readback.
- Verify all five saved identities E=dC+dcorr exactly in stored float64 arithmetic.
  Independently reconstruct dC and dR from the original plus/minus gzip journals.
- Recompute every field/group input, new-material, old-structural and
  old-continuity infinity, Euclidean and scaled-infinity norm. Euclidean norms
  use scalar `math.fsum` of squared elements, independently of the writer's
  `np.linalg.norm`. Numerical comparisons allow only512 machine-epsilon times
  the compared magnitude; exact zero remains zero. This is a floating-point
  serialization/algebra gate, not a physical-error allowance.
- Independently recompute both recorded L2 ratios; a zero denominator must yield
  null, never an invented reduction. Preserve all worse groups/fields, including
  cross-response in fields with zero input direction. Do not impose an
  improvement threshold or certify nonlinear convergence.
- The base proposal must remain UNEVALUATED, accepted/restart_allowed/can_publish
  false. Its stored proposal must equal base candidate plus stored correction.
- Validate both component test receipts (11 material-point and12 algebraic-Schur
  tests), without rerunning them. The `from_journal` factory has not yet received
  independent synthetic/domain/fake-context tests; the actual saved input's
  successful factory call is not full root-wiring certification.

The single readback is externally capped at60s, uses only Python standard
library and NumPy, and writes fresh `output_readback*` evidence under
`logs/process_material_response_20261004`. Initial status: NOT_RUN.

## Executed result

PASS_READBACK_ONLY, exit0, internal0.722366s/external0.883457s, no timeout.
Command: `.envs/pysph/python.exe -B logs/process_material_response_20261004/output_readback_run.py`.
The reader imported only standard-library modules and NumPy. No production
module was imported, and no build/apply/evaluate was called. All175 watched
source/artifact/receipt files retained identical hashes. The11+12 component
test receipts and complete OK logs bind the current tested source hashes.

All21 NPZ arrays have the declared inventory,342 unique persistent IDs and
finite values. Five E=dC+dcorr identities, original-gzip dC/dR, old-NPZ dC/dR,
and five base-proposal addition identities are exact. All300 reported norm
values and50 ratios (including11 null zero-denominator cases) agree with the
independent recomputation under the preregistered floating-point-only bound.

Whole-population raw action-error comparison follows. These are differences
under the existing saved C+/C- direction, not errors against a physical oracle.

| Field (SI) | Old continuity infinity norm | New material infinity norm | New/old L2 |
| --- | ---: | ---: | ---: |
| Position (m) | 3.5773849601e-15 | 1.5387806815e-16 | .0842077792 |
| Velocity (m/s) | 1.4308974618e-7 | 6.1551548946e-9 | .0913714318 |
| Density (kg/m3) | 6.5303902716e-9 | 4.4781904311e-10 | .0812734642 |
| Temperature (K) | 3.2605385059e-10 | 3.2605385059e-10 | 1 |
| Mean velocity (m/s) | 7.1544873074e-8 | 3.0775935247e-9 | .0841838178 |

The complete JSON also retains old-structural norms and all plastic,
no-contact, multiple-contact and zero-event subsets. The temperature action
array is bitwise unchanged from the old-continuity array. Density and
temperature have exactly zero dC, so their nonzero cross-response is recorded
without a spurious new/input reduction ratio. Temperature's multiple-contact
subset has zero old and new L2, so its new/old ratio also remains null.

No deterioration is hidden: compared with old continuity, the multiple-contact
group infinity norms are larger by factors1.0108395001 (position),1.0108410912
(velocity), and1.0108426800 (mean velocity). Their scaled-infinity counterparts
are the same three deteriorations in root units: six explicitly saved metric
rows. The corresponding L2 ratios are about.9901711, so a slightly improved L2
does not imply uniformly improved components. No other listed group/field
norm worsens against either saved old baseline.

Maximum new velocity/mean-velocity/position component is PID221,y, with two
contact faces and an observed zero event; it is not plastic. Maximum density
component is PID198 (no contact), and temperature is PID269 (the sole plastic
point, no contact). The new velocity no-contact L2/old-continuity ratio is
.0019029788; the plastic-point ratio is.0001107495. This supports the local
material-response repair while leaving contact/event-dominated behavior as a
limitation, not proving a unique causal contribution from the FKS rate term.

## Permissions and remaining limits

`UNEVALUATED_BASE_PROPOSAL.json.gz` retains its explicit UNEVALUATED label and
accepted/restart_allowed/can_publish=false. Its diagnostics also retain
physical_map/full_tangent=false, and its algebraic reaction is explicitly not
physical reaction/work. It is not a checkpoint or evaluated physical endpoint.
The source archive remains NONSMOOTH_NOT_CERTIFIED; no root convergence,
physical acceptance, new step or transport event was established here.

The12 independent Schur tests cover the numeric constructor and apply map;
none calls `from_journal`. That new factory's independent domain, deliberately
inconsistent journal and fake-context tests remain NOT_RUN. Its one successful
real saved-journal use and this hash/algebra audit do not certify complete root
wiring. Matrix backward-error values are bound saved diagnostics, not a second
independent matrix solve in this readback. Geometry/volume, contact response and
BE heat-feedback derivative omissions remain, and midpoint geometry/continuity
also differ from old P, preventing unique attribution to one material term.

Evidence: `logs/process_material_response_20261004/output_readback.py`,
`output_readback_run.py`, `output_readback.json`, `output_readback.log`, and
`output_readback_supervision.json`. The JSON retains all raw comparison values,
worsening/zero-denominator lists, checks and175 before/after hashes.
