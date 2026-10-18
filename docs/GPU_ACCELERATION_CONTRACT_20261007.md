# Explicit GPU acceleration, 2026-10-07

User request: “能否接入gpu，让我进程跑快点”. This supersedes the old
no-GPU restriction only. PySPH particle/NNPS backend, material/contact/heat
equations, FP64 precision, force residual, energy, geometry and time screens
remain unchanged. No old source or receipt is edited.

Explicit CuPy dense solve/inverse providers operate in a project-local overlay;
host NumPy arrays, assembly and physics remain on CPU. Each call includes
transfers and synchronization in measured wall time. CPU fallback and float32
are forbidden. Original Broyden/Newton budgets and retry rules are retained.

Before AFSD roots: pass AST provider whitelist, GPU linear algebra and independent
thermal fixtures; then compare actual 1381-particle full maps and 3062-DOF linear
systems from the exact CPU stage013 parent. Full-map comparison uses unchanged
batch.candidate_equivalence; mechanical solve posterior residual <=1e-9,
CPU/GPU relative solution and inverse differences <=1e-8, condition <=1e12.

Then replay stage013 coarse + two halves from its immutable stage012 parent.
Each candidate must pass every original physical gate, and the original global
and all-particle local temporal screens. CPU/GPU endpoint differences must have
position <=0.001 dx, temperature <=1e-4 K, global stress/plastic relative <=1e-3;
IDs, mass, command history and crossed IDs must agree. These endpoint tolerances
allow different floating-point paths through a finite-tolerance nonlinear root;
they are stricter than the original temporal screen and cannot relax it.

Only after replay PASS may its last screened fine state seed new GPU diagnostic
intervals. Explicit parent hash/history and provider hashes accompany every
stage; CPU history remains verified unchanged. Each new interval again compares
one coarse and two half solves. Formal acceptance/restart remain false. Failed
checks stop rather than delete evidence or relabel a CPU run as a GPU run.

Benchmark forecast90s. Replay forecast30+1.5*CPU measured complete-stage time.
Continuation forecast30+4.5*last maximum GPU candidate time. Each must <=450s;
internal540s/external600s. Initial integration is limited to stationary fill;
traverse requires the existing new-material spreading evidence. A GPU speedup
does not validate deposition, constitutive assumptions, or a future time step.

CuPy installation reference: https://docs.cupy.dev/en/stable/install.html
Use existing CUDA12.8; no system installation and no changes to frozen NumPy ABI.
