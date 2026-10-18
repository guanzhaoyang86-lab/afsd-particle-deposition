# Independent expanded startup refinement v3

Register before v3 cost/root work. Preserve v1 bookkeeping failure, both v2
cost passes, and both v2 root failures. The 918-particle independent branch has
**zero accepted steps/time**. No previously failed candidate is a restart.

The v2 62.5 us root converged in 6 iterations/10 full-map evaluations, followed
by one fresh map. Force residual 1.3954685902650699e-7 N passed 1.01e-5 N, but
raw energy 3.8382123964101473e-7 J exceeded 2.1402874845967925e-7 J. The material
work/storage gap was 3.838207501325283e-7 J, including BE elastic deviatoric
3.148187755165156e-7 J and BE volumetric 6.900324589603698e-8 J.

The v2 31.25 us root converged in 5 iterations/6 full-map evaluations plus one
fresh map. Force residual 3.810886791544635e-6 N passed 1.01e-5 N; raw energy
2.889192566095835e-7 J exceeded 1.4349068327744203e-7 J. The material gap was
2.8891910133194426e-7 J, including BE elastic deviatoric 2.386392706067829e-7 J
and BE volumetric 5.028060231539308e-8 J. Normal-contact work error was zero
for both roots. The raw gap must not be subtracted, compensated, or converted
to heat. Scalar/full-path checks were NOT_RUN after the fresh energy failures.
Different final times are not an equal-time convergence comparison.

The new static `probe_expanded_hot_domain_cost_v3.py` allows **only** 15.625 us
and 7.8125 us, each from the same independent 8 mm v2 time-zero snapshot. All
physical routines, guesses, absolute/relative energy thresholds, whole-path
limits, fresh/scalar checks, 540/600 s limits and no-acceptance policy are exactly
those in the original cost contract and v2 bookkeeping addendum. Default/new
evidence tag is v3. Every protocol additionally binds both original v2 cost
receipts and failed-root manifests, fresh gates and complete mechanical/thermal
ledgers. No failure is overwritten or retrospectively promoted.

Schedule `cost` then `root` at 15.625 us, followed by `cost` then `root` at
7.8125 us. Root remains conditional on its measured <=450 s forecast. This is
four bounded diagnostic stages, not an accepted startup sequence. Independent
elastic-rank qualification and a separately registered same-parent equal-time
comparison (15.625 us against two 7.8125 us increments) are still required before
a real 0.125 ms startup window can be accepted. This addendum does not run or
approve that process controller and does not authorize larger time steps.

```
.envs/pysph/python.exe -B scripts/probe_expanded_hot_domain_cost_v3.py --stage cost --dt-us 15.625 --tag v3
.envs/pysph/python.exe -B scripts/probe_expanded_hot_domain_cost_v3.py --stage root --dt-us 15.625 --tag v3
.envs/pysph/python.exe -B scripts/probe_expanded_hot_domain_cost_v3.py --stage cost --dt-us 7.8125 --tag v3
.envs/pysph/python.exe -B scripts/probe_expanded_hot_domain_cost_v3.py --stage root --dt-us 7.8125 --tag v3
```
