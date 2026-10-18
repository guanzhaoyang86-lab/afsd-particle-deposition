# Independently initialized full-shoulder domain: bounded startup

Registered before a process step is accepted. This is a new 8 mm square
substrate / 918-particle / 1961-free-DOF, 450 C hot quasi-static branch. Its
initial snapshot is `outputs/expanded_hot_domain_20261004/initial_8mm_visible_v2/INDEPENDENT_INITIAL_NOT_RESTARTABLE.npz`,
SHA256 `31a8ae1910f8189ab3f437f26f663b523fbf0d8195bd41e1a78cc000696fcfa3`.
It starts at time and step zero with its own inventory, unloaded mechanical
fields, initial thermal energy and eligible observer. No old 342-particle state
is enlarged or restarted; its 50.9375 ms cannot be added to this branch.

The fixed PySPH commit, density, modulus, speeds, temperature law and all
coupled operators are unchanged. The independent initialization contract lists
the explicit substrate-width change and every parameter's source and unit.
This model omits inertia and starts hot; it is ENGINEERING_APPROXIMATE, not the
paper's startup or a continuation of the inertial branch. Initial thermal
equilibrium is not claimed. Initial footprint coverage is not a spatial-domain
convergence result or evidence of sustained deposition.

The 62.5 us and 31.25 us diagnostic roots from the new initial state remain
FAILED: their fresh free forces converged but original raw and interval energy
gates rejected the backward-Euler material work/storage defect. This defect is
neither subtracted from the residual nor converted to compensating heat.

Before process execution require: source-bound new-domain initial geometry,
unloaded mechanical equilibrium, common visible/discrete initial graph and
same-displacement complete-map equivalence, full free elastic rank without
mass shifts, pseudoinverse or artificial constraints; qualified tight-path
components; and measured affordable 15.625 / 7.8125 us roots with every fresh
physical, scalar, finite-shell and fixed-retained-graph path gate. Diagnostic
roots alone do not publish process checkpoints.

Register exactly 0.125 ms of startup: eight nonrestartable 15.625 us coarse
increments and sixteen 7.8125 us fine increments from the same independent
initial state. Each pair must retain all physical/history gates and pass the
existing global and raw/retained one-hop local comparisons: 5 um position,
0.5 K temperature, 10% stress (1 MPa local floor), 20% plastic strain (1e-5 local
floor), identical strict crossing IDs. Local seeds include all actuator IDs,
contact events, blocked-edge endpoints and crossed IDs, on both branches.
The initial/new eligible population is never inferred from initially lower
material. No new inlet particles are inserted in this startup experiment.

Only paired fine candidates that pass all gates may be published. New accepted
states explicitly bind the independent origin and qualification. A failed
coarse stage can supply only its already qualified contiguous prefix; the FAIL
is not changed. Any incomplete fine pair remains nonrestartable. All interval
records and actual fields are saved, without filling missing or future frames.
Independent readback must check complete physics journals, saved graph and
certificate arithmetic, ancestry, exact publication, temporal comparisons,
energy ledger and particle identities before continuation is permitted.

Measured same-domain root time plus two map evaluations and one second of
certificate/output allowance gives the first candidate forecast. A complete
stage forecast of 30 s + 1.5 * planned increments * measured candidate cost must
be <=450 s. Each stage has a 540 s internal and 600 s external limit, one CPU
thread and a new immutable output directory. No GPU, mass/time/modulus scaling,
temperature/yield clipping, negative-dissipation clipping or imposed flowlines.
Larger steps, longer windows, replenishment and domain sensitivity require
separately evidenced follow-on qualification; this startup is not deposition.
