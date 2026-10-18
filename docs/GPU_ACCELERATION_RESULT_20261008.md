# GPU acceleration of the deposition process

User explicitly requested GPU acceleration. New execution provider: CPU PySPH
particle/neighbor/physics operations with CuPy CUDA **float64** dense mechanical
solve/inverse and implicit thermal solve. Original CPU files/checkpoints stay
unchanged. The new static adapters retain all force, energy, contact, geometry,
path and coarse/fine gates. This does not change the model or prescribe particle
trajectories; it is a numerical execution change.

## Environment and installation

- NVIDIA GeForce RTX 5090 Laptop GPU, device memory 25,650,855,936 bytes.
- Driver 595.97, compute capability 12.0; installed CUDA toolkit12.8.
- Existing PySPH Python3.11.16 and NumPy2.4.6 remain unchanged.
- Isolated `.envs/gpu_overlay_20261007`: cupy-cuda12x14.2.0,
  cuda-pathfinder1.8.3, fastrlock0.8.3. Install package hashes are retained in
  `outputs/gpu_acceleration_20261007/install_report.json` and repair receipt.
- Initial cuda-pathfinder1.3.1 was too old for this CuPy version, corrected to1.8.3.
- The desktop inherited an unrelated Anaconda CUDA12.1. The GPU backend now
  preloads and verifies absolute toolkit12.8 DLL paths in its own process.
  Existing conflicting loaded DLLs are refused. System settings were not edited.
- Actual NVRTC12.8 verified. Driver-reported CUDA13.2 is not the installed toolkit.
- No silent CPU fallback, lower precision, system-wide install, or paid compute.

Official installation reference: https://docs.cupy.dev/en/stable/install.html

## Component and actual-state verification

`outputs/gpu_acceleration_20261007/components_v1/manifest.json` is PASS with
all four executed commands exit0. It binds source, dependency and log hashes.
The suites include17 thermal,8 search-adapter and15 candidate/receipt tests,
plus real FP64 solve/inverse comparisons and singular/float32 rejection.

Two exploratory thermal test receipts remain FAIL (`thermal_gpu_components_v1`
and`v2`). Those tests incorrectly demanded nearly cancelled CPU/GPU residuals
and global sums agree to a fixed1e-12 J. The revised **test** checks each residual
against the original heat-balance bound, each ledger against its own particle
fsum, and ledger differences against propagated strictly checked particle
differences. Production equations, heat gates and particle/temperature comparison
tolerances were not changed. Large-dt/permutation tests and deliberate corruption
negative tests passed in v3 and the final frozen aggregate.

Actual benchmark uses the immutable stage012 parent of CPU stage013,1381
particles and3062 free mechanical DOFs. Three repetitions include CPU/GPU transfer
and synchronization; the first GPU full map includes remaining cold overhead.

| Measured operation | CPU median | GPU median | Ratio |
|---|---:|---:|---:|
| Mechanical solve |0.366973 s|0.053941 s|6.80x|
| Mechanical inverse |1.269084 s|0.175904 s|7.21x|
| Complete residual map |0.485271 s|0.414390 s|1.17x|

Mechanical CPU/GPU solution relative error4.3643e-13; inverse5.5152e-13.
GPU solve posterior residual8.3775e-15, inverse infinity residual1.9957e-9.
All three actual full maps pass the unchanged force/stress/plasticity/temperature/
density/contact/raw-energy equivalence checks with identical support graphs.
Benchmark stage wall25.3929s; this is not an accepted physical time advance.

## Complete-stage evidence

`outputs/gpu_acceleration_20261007/replay_v1`: PASS, exit0. Same stage012 parent,
same474.609375us macro interval, one coarse and two half steps. Original physical
candidate gates and global/all-particle temporal comparison pass.

| Complete-stage wall time | CPU reference | GPU replay |
|---|---:|---:|
| Three candidates and checks |215.426698 s|158.386309 s|
| Coarse candidate |98.732286 s|72.645661 s|
| First half |65.4787 s|47.366579 s|
| Second half |39.6206 s|25.704829 s|

Overall1.3601x, wall reduction26.48%. This is one measured interval, not a
promise for larger particles, later plastic flow, different hardware or two days.
The CPU comparison is the preserved previously executed stage, not a concurrent
new CPU run. The actual same-parent full-map benchmark above is interleaved CPU/GPU.

Across the three endpoints, max CPU/GPU position difference1.1372e-18 m,
temperature2.2738e-13 K, relative stress6.7047e-13, relative plasticity3.8388e-15.
Command histories, masses and IDs agree. GPU calls:619 thermal solves,
88 mechanical solves,9 inverses; no CUDA failure or CPU fallback. Aggregate
GPU call wall16.1601s includes transfers:10.20% of stage wall. The remaining
142.2262s is outside these GPU calls (CPU physics/assembly/diagnostics, validation,
I/O and startup; no complete per-function profile was taken). This integration
is partial. CPU/GPU iteration counts, evaluations, rebuilds and line-search alpha
sequences were independently verified identical: the speedup used the same work.

Replay is a duplicate interval ending11.5625ms; it does not advance beyond
the original CPU result.

New continuation `continue_20261008_001414_714` is PASS, exit0:145.144699s,
three candidates51.144333/57.756127/23.552285s, original temporal checks passed.
It advances the screened fine history to **12.037109375ms /98fine steps**.
New checkpoint SHA256:
`1814c35b94e2482db4195aa5d433eaf1435dd6d03a09db3e82bef7643b68deab`.
New supplied crossing count0, shoulder count0: **deposition is not established**.
The one-interval verification job has ended; no new GPU background job remains.
All states remain diagnostic-only with formal accepted time advance0.

Readback receipt: `outputs/gpu_acceleration_20261007/readback_continue_v1.json`.
It performs no solver evaluations and checks the GPU receipts/parent/screened
history and unchanged old formal checkpoint. Next one-pair forecast289.9026s
is below the existing450s gate; this is a forecast, not an executed future run.

## Usage

Main entry: `scripts/run_afsd_gpu_20261007.py` (benchmark/replay/continue).
Convenience entry after a successful frozen replay:

```powershell
& '.\scripts\continue_gpu_deposition_20261007.ps1' -Stages 1
```

The launcher chooses the latest passing GPU interval; the Python runner independently
verifies every source, reference, parent and screened half step. Each interval has
one coarse and two half solutions. Any physical or cost failure stops and keeps
the evidence. Default is one bounded interval; no long-running job is implied.
Each interval forecasts <=450s, with540s internal/600s external limits. Current
implementation is explicitly bound to this Windows CUDA12.8 environment and
stationary filling; it is not yet a portable Linux/server or traverse launcher.
