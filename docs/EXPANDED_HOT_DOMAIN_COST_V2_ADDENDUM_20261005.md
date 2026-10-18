# Cost probe v2: source-map preflight repair only

The v1 CLI failed before spawning a worker or evaluating any physical map:
`dict()` received duplicate keyword names from overlapping dependency maps.
Preserve the original driver and the preflight failure receipt. No physical
result or budget estimate exists from that attempt.

The new static driver `probe_expanded_hot_domain_cost_v2.py` merges each source
map explicitly, requiring every repeated path to have the same SHA256. This
repairs source bookkeeping only. The independent 8 mm v2 initial snapshot,
62.5/31.25 us intervals, physical factory/map, tangent, forecast, full root
diagnostics, and 540/600 s limits are unchanged. Bind both the original contract
and this addendum. Use new evidence tag `v2`, and keep all v1 files unchanged.

```
.envs/pysph/python.exe -B scripts/probe_expanded_hot_domain_cost_v2.py --stage cost --dt-us 62.5 --tag v2
.envs/pysph/python.exe -B scripts/probe_expanded_hot_domain_cost_v2.py --stage cost --dt-us 31.25 --tag v2
```

No root is scheduled by this repair. Each root still needs its matching measured
cost receipt and a separate parent-controller schedule. The frozen visible
full-map implementation sets `wall_seconds`; the scalar diagnostic may read
that field without assuming a missing timing result.
