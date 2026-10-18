# Saved response output readback: scope and bounds

2026-10-03, before the dedicated readback validation. No solver/contact/material
function or advance is called. Load existing NPZ with the PySPH serializer only;
use independent NumPy algebra and hashes on recorded data. One supervised call,
internal 55 / external 60 s; all writes are new review artifacts.

Exact checks: all recorded source/artifact hashes, both parent and candidate
40-property PA inventories versus their canonical snapshots, semantic hashes,
parent linkage, private/nonrestartable flags, unchanged physical constants and
passive identities. Re-check the old 4516-file freeze, allowing only the three
already-authorized continuation index documents.

Algebraic readback bounds are 2048 machine eps times the largest recorded term
scale plus 1e-24 (corresponding SI units), for node impulse/work/Q reductions,
complete fixed-node response block identities, and kinetic/momentum bookkeeping.
These bounds verify arithmetic, not an AFSD physical acceptance tolerance.
The original raw nonlinear residual gates remain 1e-9 m/s, 1e-12 m,
1e-10 relative density and 5e-12 K. Thermal balance defects are recorded, not
silently converted to heat or judged by a new physical energy tolerance.

For selected PRIVATE contacts, independently find the polynomial longitudinal
zero and check the complete tangential vector there. Integrate the unchanged
declared PRIVATE signed power on [0,s0] and [s0,1] using 96/192 Gauss nodes,
the saved dense kinematics, physical wall program and penalty law. Require
96/192 agreement within 1e-20 J; this is a bounded diagnostic quadrature,
not additional particle advancement. Record negative-tail work even if none
of the original 16 nodes sampled it. Compare, but do not overwrite, the saved
net Q. A PRIVATE root PASS is not whole-path passivity or a located event.
