# Action inspection: serialization-only repair 1

2026-10-04. Registered before this repair execution. The initial invocation in
logs/process_bridge_action_20261003 failed before either matrix was built:
TypeError: Object of type ndarray is not JSON serializable. Its scripts, report,
log and supervision are retained unchanged; zero physical calls occurred.

This single repair copies the checker and supervisor to a fresh directory,
changes only canonical JSON state/parameter comparisons to use the existing
recursive plain() ndarray/NumPy-scalar conversion, and binds this repair note.
All state fields remain compared; none is removed. Original algebra, physical
inputs, proposed plus/minus actions, group definitions and tolerances are unchanged.
No evaluate, contact kick, material update, new epsilon or physical step is added.

Initial mechanical matrix actions: NOT_RUN. This repair permits one bounded
execution, internal55s/external60s, forecast below15s. Preserve FAIL or partial
outputs if it fails; no further repair in this package. The original one-attempt
rule is amended only for this pre-algebra serialization defect, not for numerical
results or physics. This file must not be interpreted as evidence that an action
check or physical root has passed.

