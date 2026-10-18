# Real342 bridge output readback contract

2026-10-03. Registered before executing the new readback script. This is a
read-only review of the existing four saved candidate evaluations, not another
map evaluation, root solve, acceptance test or trajectory. The reviewer does
not call `afsd_bridge_audit.audit`, any material/contact/thermal function, or
the bridge driver as an oracle. NumPy algebra, standard JSON/hash operations
and the upstream PySPH checkpoint reader are sufficient.

Scope: verify all report-bound source/input/artifact hashes and the 4600-file
prior freeze (only its three already-authorized continuation documents may
change); read all parent particle properties and semantic state from the
original accepted NPZ and compare with the archived parent; verify latest
published history/checkpoint has not changed. Reconstruct the fixed geometric
direction, actuator free-axis mask, exact C/C+/C-/C input arithmetic, signed
residual quotient, and complete first/fourth numerical evaluation equality.
Report branch differences without declaring a smooth Jacobian.

Independently reconstruct, for all four saved evaluations, every recorded
face/subinterval signed friction identity and its existing roundoff bound,
interval normal/constraint impulse and work, positive plastic work and return
dissipation, and mapped-temperature BE source-once/pair/wall balances. This is
not an independent re-derivation of the stress/contact solution itself. Keep
negative raw Q and distinguish a balanced map from a converged fixed point.

Bounds copied unchanged from the registered bridge/audit contract:

- Exact hashes, property/semantic state identity, branch serialization and
  repeat output. Independently rebuilt D: absolute 2e-15 dimensionless, from
  double-precision arithmetic; candidate arithmetic is exact using saved D.
- Friction identity: 1e-25 J + 1e-12 relative. Roundoff certificate identity:
  1e-30 J + 1e-12 relative; raw Q >= minus its existing bound. No clipping.
- Momentum: m*1e-9 m/s + 1e-17 Ns. Free-axis reaction: 1e-17 Ns.
  Work: 1e-15 J + 1e-9 times sum of absolute work terms.
- Plastic/return work: 1e-18 J + 1e-10 relative. Plastic strain monotone.
- BE source/pair identity: 1e-25 J + 1e-12 relative. Reconstruct radial
  conductance in its recorded scalar expression; its comparison uses
  1e-15 W/K + 1e-12 relative. The original physical BE bounds use
  scale=sum(abs(total_heat_input_i)), particle 1e-12 J+1e-10*scale,
  global 1e-12 J+1e-9*scale, pair net 1e-15 J+1e-10*sum(abs(pair heat)).
  Tiny reassociated diagnostic residuals are reported, not independently
  treated as heat sources or used to widen the physical balance threshold.

The script records maximum errors and contextual examples, not another huge
per-identity journal. The actual full evaluation gzip remains authoritative.
Candidate residual magnitudes, temperature/rho/penetration/conditioning
bounds, authority flags and budget receipts are explicitly inspected.
`NONSMOOTH_NOT_CERTIFIED` must not hide an accounting FAIL or imply physical
acceptance. Even a successful readback cannot establish time/space accuracy,
macroscopic transport, a smooth contact Jacobian, or production permission.

Readback budget: internal 55 s, independent external 60 s, one execution,
no new evaluations and no writes to any original output/history. The actual
bridge's pre-run budget amendment was 110/120 s; report it as such, not the
earlier provisional 50/60 s estimate. Its four actual measurements took about
40 s. Status of this new readback before execution: NOT_RUN.

## Initial readback tooling failure; bounded repair1 proposal

The initial readback exited FAIL in 2.647 s external after 5093 checks and zero
failed assertions: `close()` tried to report `max(limit)` for a zero-contact
empty array. This is an empty-array reporting bug, not evidence of a physical
failure. Its source and initial FAIL output remain untouched. Repair1 is a
new wrapper that treats matching empty shapes as the vacuous identity, with
the original nonempty comparator and all tolerances unchanged. It writes a
fresh output directory and still performs no evaluation. The replacement
must be explicitly authorized by the main task and uses the same 55/60 s
budget; no further readback retries are planned. The main task authorized this
single repair. Empty comparisons are counted separately as
`EMPTY_NOT_APPLICABLE_shape_finite`, never as nonempty numerical PASS; only
their equal shape and finite contents are checked.
