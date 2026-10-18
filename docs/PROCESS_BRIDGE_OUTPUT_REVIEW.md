# Independent saved-output review: real342 process bridge

2026-10-04. Result: **PASS_INDEPENDENT_READBACK_NOT_PROCESS_ACCEPTANCE**.
No physical evaluation, advance, thermal solve or acceptance function was
called by this review. The review used NumPy algebra on saved journals, JSON
and file hashes, and upstream PySPH's existing checkpoint reader. It did not
reuse `afsd_bridge_audit.audit` as its reference.

This establishes reproducible accounting for the four existing real342 map
responses. It does **not** make their candidates converged, their directional
quotient a Jacobian, or their outputs accepted AFSD states.

## Execution and preserved failure

Scope and tolerances were registered in
`docs/PROCESS_BRIDGE_OUTPUT_REVIEW_CONTRACT.md` before execution. Initial
command ` .envs/pysph/python.exe -B logs/process_bridge_20261003/output_readback_supervisor.py`
failed after 2.647360 s external: a statistics helper attempted `max(empty)`
for an actual zero-contact array. It had completed 5093 checks with no failed
assertion. This was a readback-tool failure; the initial source and FAIL
receipt remain untouched.

The main task authorized exactly one empty-array-only repair. New wrapper
`output_readback_repair1.py` preserves the original nonempty calculations and
tolerances. Empty shapes are recorded as `EMPTY_NOT_APPLICABLE_shape_finite`,
not counted as nonempty numerical PASS. Executed command:

```text
.envs/pysph/python.exe -B logs/process_bridge_20261003/output_readback_repair1_supervisor.py
```

Exit 0; external 10.8231223 s, internal 10.4777769 s, limits 60/55 s.
471275 check records: 467339 nonempty/identity assertions and 3936 separately
identified empty shape/finite checks; zero failures. These are repeated
accounting assertions, **not** 471275 independent scientific benchmarks.
No original report, output, source or accepted history was written.

Primary evidence:

- `logs/process_bridge_20261003/output_readback_repair1/report.json`:
  aggregate errors, independent physical totals, branch differences and bounds.
- `logs/process_bridge_20261003/output_readback_repair1_supervision.json`:
  actual command, elapsed time, exit and three readback-source hashes.
- `logs/process_bridge_20261003/output_readback_initial/report.json` and
  `output_readback_initial_supervision.json`: retained tooling FAIL.
- `outputs/process_bridge_20261003/real342_initial/*`: unchanged original data.

Readback source SHA256:
`71bfb7337d613eea44c44b4632e63efb146f8d3de36cf3c9abf6b8a339893f6f`;
repair wrapper:
`a57a30aeae3d0acccdd906d69b3eb70fa4a0309ae8cb6546644a9ba7337089ea`.
Reviewed original report:
`436085fa1694a2bcb45a1c615123e293d146f9e84474ad021574baad041939b0`.

## Identity, provenance and publication: PASS

All 128 bound sources, 159 inputs and 17 artifacts matched before and after
readback. The 4600-file prior freeze matched except its already-authorized
continuation changes to IMPLEMENTATION_PLAN and OUTPUT_SCHEMA. The actual
273-test preflight source/artifact receipt and actual-run supervisor were
also bound by their hashes.

All **41** particle properties and the entire semantic state matched the
original accepted NPZ; not only x/v/stress. Parent hash:
`ff94e2a7ae4c4a3861e2d6600f322cdaa44c42318fccba8b0a74cd5f997b11ea`.
The real342 historical parent remains 650.3936329109041 ns. The diagnostic
endpoint is 700.3936329109042 ns; it is not a new accepted state. Public latest
history remains four stages ending at 850.393632910904 ns, its existing NPZ
and manifest unchanged. The historical parent was not substituted for latest.

The fixed affine direction was independently reconstructed from geometry,
mass and persistent IDs; the 25 actuator-particle constrained z components
remain zero. All five fields of C/C+/C-/C obey the registered epsilon=1e-7 m/s
and dt=5e-8 s arithmetic exactly using the saved direction. Density and
temperature perturbations are exactly zero. Complete base/repeat candidate,
evaluation and branch outputs match, not just a selected norm. The saved
signed five-field residual quotient equals `(Rplus-Rminus)/(2*epsilon)`.

Exactly four actual map evaluations were recorded. All candidate and report
accepted/restart/publish/production authority flags remain false; new physical
steps and transport events remain zero. The actual measurement used the
pre-run **110 s internal / 120 s external** budget amendment and took
39.4136653 / 39.9551181 s. This review does not mislabel that watchdog as the
earlier provisional 50/60 s budget; both actual and review remain below the
user's 600 s single-run limit.

## Independent accounting: PASS within unchanged bounds

Reconstructed all four outputs' face/subinterval Q, original roundoff
certificate, normal impulses, constraint impulses and work, mapped momentum,
positive plastic work, and mapped-temperature BE sources. No tested physical
function generated the reference. Base examples:

| Quantity | Independent reconstruction |
| --- | ---: |
| Contact faces / multi-contact particles | 104 / 8 |
| Contact subintervals / recorded zero events | 7581 / 28 |
| Signed friction total | 5.241676700179341e-9 J |
| Actuator work | -5.595200466721312e-10 J |
| Largest mapped momentum defect | 2.497429856843555e-22 Ns |
| Positive plastic work, PID 269 | 6.484598586830537e-9 J |
| Separate return-dissipation diagnostic, PID 269 | 4.529236322291872e-9 J |
| BE plastic heat | 6.484598586830537e-9 J |
| BE material share of friction heat | 2.620838350089671e-9 J |
| Signed reservoir heat into material | -1.100000972505583e-8 J |
| Net reconstructed thermal input | -1.894572788135624e-9 J |
| Measured thermal increment | -1.894573417189349e-9 J |
| Global / maximum-particle BE balance error | -6.290537250843286e-16 / 1.724841862641914e-17 J |

PID269 has ep increment 3.03067126444021e-5 and equivalent stress
1711764.2603377479 Pa. The independent identity is `P=(m/rhomid)*q*delta_ep`.
Its nonzero return-dissipation diagnostic was **not** added to heat. Pair heat
was rebuilt on 8131 undirected edges using raw radial gradients; base pair net
is -1.542386533011108e-28 J. The original BE physical balance limits were
retained (base particle/global approximately 1.000003e-12 / 1.000030e-12 J).
The much smaller observed balances do not certify constitutive energy accuracy.

Negative raw Q was not hidden: base has 179 negative face/subinterval records,
minimum -1.304139436554811e-25 J at PID221, face221:2, interval274, versus its
original roundoff bound 7.064299832617929e-23 J. Their signed sum is
-2.156233703060425e-25 J. Plus/minus contain 147 each; all remain in the
source journals and below their unchanged rounding certificates. Independent
lab-vector reconstruction differs from recorded tangent-basis arithmetic by
at most 5.140164e-26 J across the four maps, inside the preregistered identity
bound. There is no absolute-value heat replacement or zero clipping here.

All saved candidate and map arrays meet their existing T/rho/speed/displacement
domain limits. Base mapped T is 723.14998899–723.15006315 K; maximum relative
density change is 7.380022e-5; maximum speed 0.257973 m/s. Recorded pre/mid/post
normal penetration is at most 1.975901e-6 m versus the 5e-5 m limit. Recorded
operator condition is 2.215058 versus 1e6. These last geometric/condition
checks inspect their saved diagnostics; the readback did not rerun geometry
or NNPS and does not grant an independent spatial-accuracy certificate.

## Why NONSMOOTH remains, and what it does not mean

Relative to base, plus changed interval-kind sequences for PIDs
202,209,257,281 and event signatures for PID281. Minus changed interval-kind
sequences for 202,209,257,286,329,334 and event signatures for 209,286,329.
Face inventories, free-axis masks, initial/final active normal-contact sets,
plastic-active flags and rate-floor flags did not change. Both perturbed maps
still contain 28 zero-event records and one positive-plastic particle.

Thus `NONSMOOTH_NOT_CERTIFIED` correctly refuses a smooth-response certificate
when the implemented numerical contact/event partitions change. It is not,
by itself, proof that the exact continuous contact response lacks a derivative,
and it does not conceal a discovered bookkeeping failure: all independently
reconstructed identities passed. No epsilon sweep was run.

Crucially, base and repeated base have maximum scaled coupled residual
**178.32917365 > 1**; plus/minus are 197.61448126 / 210.93455871. The base
temperature residual is about 8.916459e-10 K under the named BE map, even
though its input was accepted by a different historical explicit-thermal map.
Accounting correctness therefore does not make this candidate a converged
BE-mechanical stage. No physical acceptance call was made or implied.

Status: provenance/reentry/accounting PASS; smooth Jacobian certification
BLOCKED; full-stage convergence NOT_ESTABLISHED; new time advance, transport,
resolution/time-step convergence and production NOT_RUN. A safe next task may
use these real-map responses to assess a registered numerical search method,
but must preserve these branches, physical laws, signed ledger and guards.
It cannot simply publish the present four candidates or treat this review as
completion of AFSD reconstruction.
