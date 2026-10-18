# PySPH fallback baseline evaluation

Evaluation date: 2026-09-24 (America/Chicago)  
Scope: CPU fallback only, after WeldForm exhausted its two allowed build rounds
and WeldFormGPU stopped at its evidenced runtime device-selection blocker.

## Outcome

**WORKING_FALLBACK_BASELINE / NOT_AFSD_READY.**  The pinned PySPH source built
in the project-local Python 3.11 Conda prefix on the first install/build round.
An import probe, one official general example, one official 2D solid-mechanics
example, and one step through the official 3D impact path all ran and produced
readable NPZ outputs.  This establishes an executable upstream baseline; it
does not reproduce the paper and does not establish correct AFSD physics.

Provenance and license are recorded in `third_party/PySPH.UPSTREAM.md`.
Commands, working directories, exit codes, stdout, and stderr are under
`logs/upstream/pysph/`.

## Build and environment evidence

- Pin: `1cd3dea9cfaac7d9f1bda2a2491ce4b42a31dbe5`, detached from
  `https://github.com/pypr/pysph.git`.
- License retained: BSD 3-Clause `third_party/PySPH/LICENSE.txt`.
- Environment: `.envs/pysph`, Python 3.11.16.  It does not reuse global
  Python, Conda `base`, or `torch5090`.
- Procedure followed the cloned documentation for a tracked source checkout:
  install `requirements.txt`, then `python setup.py develop`.
- Build result: exit 0 on round 1; native cp311 Windows `.pyd` modules were
  produced.  Import resolves PySPH to the pinned checkout.
- Saved package inventories:
  `12_conda_list_explicit.stdout.log` and `13_pip_freeze.stdout.log`.
- Non-fatal evidence retained: setuptools warns that `setup.py develop` is
  deprecated; Cython reports integer constant redeclarations; runtime MSVC
  reports that upstream `-O3` is ignored.  These did not change exit status.
- Runtime JIT compilation writes to the upstream hard-coded default
  `<LOCAL_USER_PROFILE>\.pysph\source\py3.11-win-amd64`.  Dependencies remain in the
  project prefix; the JIT location must be considered in future clean-room
  reproducibility work.

Exactly one install/build round was used.  The rejected `--output-dir` smoke
command was a CLI invocation error, corrected to the source-defined
`--directory`; it was not a build repair and both logs are retained.

## Executed cases

All three runs below are **DEBUG_SHORTENED**, not benchmark or paper-comparison
runs.  `logs/upstream/pysph/14_validate_outputs.stdout.log` shows that PySPH
itself reloaded each final NPZ and confirmed the listed state.

| Case | Actual command delta | Result | Evidence boundary |
|---|---|---|---|
| `elliptical_drop_simple.py` | Added `--max-steps 2 --pfreq 1`; retained source `dt=5e-6 s`, `tf=0.0076 s` | PASS; 5,025 particles; 2 steps; `t=1e-5 s`; run CPU 0.0181 s | General 2D WCSPH/JIT/output smoke only |
| `solid_mech/rings.py` | Added `--max-steps 2 --pfreq 1`; retained source `dt=1e-8 s`, `tf=5e-5 s` | PASS; 17,592 particles; 2 steps; `t=2e-8 s`; run CPU 0.187 s | 2D elastic solid path only; cannot prove 3D or plasticity |
| `solid_mech/impact3d.py` | Added `--max-steps 1 --pfreq 1`; retained source `dt=1e-9 s`, `tf=8e-6 s` | PASS; 36,962 particles; 1 step; `t=1e-9 s`; run CPU 2.193 s | Executable 3D call path only; too short for physical validation |

The saved `.info` files say `completed: true`, which means the bounded command
ended normally.  It must not be read as completion to each source case's
original final time.

An optional upstream 3D NNPS unit-test target was considered, but `pytest` is
not installed in the core prefix (`11_nnps_pytest_availability`, exit 1).
No additional dependency was installed merely to expand this evaluation, so
that unit test is **NOT_RUN**.  The one-step 3D impact execution is evidence of
the integrated neighbor/equation path, not brute-force NNPS correctness.

## Capability audit: code versus execution versus AFSD need

| Capability | Existing upstream code | Actually verified here | AFSD status |
|---|---|---|---|
| 3D neighbor path | `LinkedListNNPS` accepts `dim`; three-axis cell counts and a `dim == 3` path exist in `pysph/base/linked_list_nnps.pyx` and `nnps_base.pyx` | The 36,962-particle `impact3d.py` case executed one step with `Gaussian(dim=3)` and `Solver(dim=3)` | Integrated path PASS; neighbor correctness benchmark NOT_RUN |
| 3D elastic stress rate | `HookesDeviatoricStressRate` updates six symmetric stress components using all nine velocity-gradient components | `impact3d.py` called `VelocityGradient3D`, stress momentum/energy equations, and the stress rate for one step | Call path PASS; accuracy/stability NOT_VALIDATED |
| Plasticity | Only `VonMisesPlasticity2D` is defined under `pysph/sph/solid_mech/hvi.py` | The 3D example invoked that class without crashing | **NEEDS_DEVELOPMENT**: its yield invariant uses only `s00`, `s01`, and `s11` even though it rescales six components; this is not verified genuine 3D plasticity |
| Solid integrator | `EPECIntegrator` plus `SolidMechStep` updates x/y/z, u/v/w, density, energy, and all six stress components in predictor/corrector stages | 2D rings and 3D impact bounded runs both passed | Executable; convergence, conservation, and AFSD time-step suitability NOT_VALIDATED |
| Heat conduction | Gas-dynamics code has internal/artificial thermal-energy terms; `EnergyEquationWithStress` adds pressure/stress work to specific energy | No official solid heat-conduction case was found or run | Solid Fourier/temperature conduction case **NOT_AVAILABLE**; must not equate artificial gas conduction or mechanical energy work with AFSD heat transfer |
| Rigid contact | `RigidBodyCollision` and official bouncing/collision examples implement particle-overlap DEM-style rigid-body contact | Inventory only; no contact case was run in this fallback evaluation | Generic rigid-body DEM example **AVAILABLE_NOT_RUN**.  Upstream README calls these proofs-of-concept and says parameters/conservation are unvalidated.  AFSD deformable-solid/rigid-tool contact and thermal contact **NOT_AVAILABLE** |
| Paper AFSD setup | No paper-specific square-bore rotating tool, feedstock/substrate setup, material law, thermal coupling, or validation targets | None | **NEEDS_DEVELOPMENT** |

The 2D `ElasticSolidsScheme` itself imports `VelocityGradient2D`; the 3D impact
example bypasses that scheme and assembles equations manually.  Therefore a
3D plot or the presence of z coordinates must not be used to claim a general
3D thermo-elastoplastic solver.

## Selection recommendation and next-stage blockers

Given the evidenced blockers in the two preferred WeldForm paths, PySPH is the
only evaluated upstream that actually built and ran on this host.  It is a
reasonable CPU fallback foundation for interface/configuration tests and for
incremental solver work, provided the project labels it **NOT_AFSD_READY**.

Before an AFSD paper case can be attempted, the next stage needs at least:

1. a reviewed three-dimensional yield/plastic update with unit tests; do not
   silently reuse `VonMisesPlasticity2D`;
2. solid temperature state, physical heat conduction, constitutive/thermal
   coupling, and energy-balance tests;
3. deformable workpiece/feedstock contact against prescribed rigid tool/base,
   including friction and heat-transfer interfaces, with contact benchmarks;
4. paper geometry, boundary/inflow/rotation handling, persistent particle IDs,
   checkpoint/output adapters, and validation tests wired to the existing
   project schemas;
5. a clean-room policy for the user-profile JIT cache and a pinned environment
   creation file derived from the saved explicit package inventories.

No AFSD result, paper agreement, contact validation, heat-transfer validation,
or full-duration solid benchmark is claimed by this evaluation.
