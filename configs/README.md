# AFSD case configuration contract

`base.json` contains parameters shared by every paper case. Files in `cases/`
contain only the case identity and the parameters that differ.  Consumers must
resolve `extends` before use; `afsd_scaffold.config.load_resolved_case()` is the
reference loader.

Rules:

- Calculation values are SI.  `display_units` are presentation preferences and
  must never be fed back into the solver without an explicit conversion.
- Every parameter is an object with exactly `value`, `unit`, `source`, and
  `status`.
- An unresolved value is JSON `null`, never numeric zero.  `TBD` belongs in the
  provenance text, not in a numeric field.
- A resolved parameter with a physical unit must be a finite JSON number;
  booleans, strings, `NaN`, and infinities are rejected. Unit `1` additionally
  permits non-empty enum strings and booleans; any numeric value is still
  required to be finite.
- These are incomplete paper-target configurations.  Loading and validation do
  not make them runnable physics cases.  The solver launcher must reject a
  physics run while a required closure or geometry value is unresolved.
- Case labels and `tracer_mode` are passive metadata.  They must not select
  material properties, forces, stress updates, heat transfer, or kinematics.

Validate all three cases with:

```text
python tests/run_tests.py
```

Every solver launcher must call the fail-closed gate after resolution:

```python
from afsd_scaffold.config import assert_physics_ready, load_resolved_case

case = load_resolved_case("configs/cases/starved.json")
assert_physics_ready(case)  # PAPER_COMPARISON: blocks on any unresolved value
```

For a `DEBUG` or `APPROXIMATE` run, a backend may pass its reviewed, non-empty
required path set and the actual run classification. Only those paths are then
allowed to gate capability; unknown/group paths and unresolved required values
still fail. `PAPER_COMPARISON` always checks the full tree even when a path set
is supplied.

