# Selected HHT 27-particle local reference contract

2026-10-03. `LOCAL_BENCHMARK_HHT_CURRENT_3D_TWO_FACE_THERMAL`.
Registered before the first reference execution. This supersedes the candidate
reference-configuration geometry/material programme in
`PROCESS_HHT_FIXTURE_DESIGN.md`, retaining that document and its cost evidence.
No physical AFSD state, deposition result, production permission or alloy
validation is created. Coupled HHT common-end-time status remains NOT_RUN until
an actual joint driver passes; an independent reference is not that driver.

## Exact inputs and provenance

- E2 manufactured geometry: `[0,0.01]^3 m`, three points per axis, lexicographic
  IDs `1..27`, spacing `0.005 m`, `h=1.3*dx`, tensor trapezoid volumes
  `dx^3*wi*wj*wk`, with endpoint weights one half and midpoint weight one.
  Physical mass is `2700*V kg`, total `0.0027 kg`; no mass/time/modulus scaling.
- E2 manufactured moduli: `E=70 GPa`, `nu=0`, `G=35 GPa`, `K=E/3`.
  Actual FKS is unchanged: a1=317.6 MPa, a2=1 MPa, a3=70, b1=12.44,
  b2=5.65, c1=13, c2=.01, room=293.15 K, melt=855.15 K,
  reference_rate=rate_floor=1/s. The FKS unresolved/APPROXIMATE provenance is
  inherited; this is not a synthetic high-yield replacement or AA6061 fit.
- Initial current position equals reference position, density `2700 kg/m^3`,
  velocity `(0,0,.02) m/s`, `S=diag(20000,-10000,-10000) Pa`, and `ep=0`.
  This is explicitly prescribed nonzero initial deviatoric stress, not an
  unstated prior loading simulation. Its initial energy is recorded.
- Initial temperature `400+2*Xx/.01 K`, cp=900 J/(kg K), conductivity=167
  W/(m K), supported solid domain as above. These thermal constants come from
  the existing debug configuration and remain approximate assumptions.
- Two fixed contact cohorts are the initial x=0 and x=.01 faces, nine points
  each. Normals into material are `(1,0,0)` and `(-1,0,0)`. Particle contact
  radius is `dx/2`; contact planes start at `-dx/2+2e-6 m` and
  `.01+dx/2-2e-6 m`, giving exactly 2 micrometres initial overlap.
- Normal stiffness is `7e6 N/m` per boundary particle, independent of its
  trapezoid weight; initial force is therefore 14 N per contact, 126 N per
  face. This is a specified LOCAL_BENCHMARK penalty input, not an AFSD
  calibration. The Coulomb coefficient `.2` is the already-used manufactured
  rotating-contact coefficient, explicitly selected here rather than the
  debug AFSD value `.3`. Full tangential Coulomb disk, no regularization.
- Common final time `tf=2e-6 s`. Inward normal plane stroke is
  `d(t)=1e-7*(1-cos(2*pi*t/tf)) m`, peak .2 micrometres at 1 microsecond;
  the stroke returns to zero at tf. Use its actual derivative for wall power.
  Both planes have tangential z velocity +.2 m/s before 1 microsecond and
  -.2 m/s from that exact switch onwards. Only rigid walls are controlled.
- Left/right reservoir temperatures are 410/420 K. Each active contact has
  conductance `1000*dx^2=.025 W/K`; no trapezoid area correction is inferred.
  Material gets half the raw friction heat; plastic heat fraction is one.
  All 81 material DOFs are free and force driven. No actuator or clamp.

Current geometry, density, correction moments and raw radial conduction pairs
evolve in the reference and coupled driver. The reference must check moment
conditioning <=1e6, positive density, supported temperature and independently
evaluated FKS elastic margin. Its elastic Jaumann stress-rate equation retains
spin. Replacing this with a frozen linear stiffness omits the tested material
history and is not the selected coupled problem.

Nominal HHT grids are 250,125,62.5 ns, giving 8,16,32 nominal steps to tf;
the prescribed reversal lies on every grid. It is a control-induced slip
reversal, not proof of continuous zero-slip/sticking or a shared two-face
corner constraint. Actual material velocity is continuous at the switch.
The HHT algorithm must separately derive the force/acceleration history on
the two sides and keep raw event work/impulse/heat. The signed z slip margin
must exceed .05 m/s throughout each smooth half. Contact opening, low slip
or yielding encountered by this minimal reference are explicit failures of
its declared smooth branches, not occasions for clipping or changing inputs.

## Independent reference and acceptance

`scripts/reference_hht_block.py` uses only NumPy plus standard-library
I/O. It enumerates all unordered current pairs, evaluates the scalar radial
polynomial, forms current correction moments and paired forces, integrates
the elastic Jaumann stress rate and continuity, and evaluates literal FKS
yield independently. It does not import tested backend force, material, heat,
contact or HHT code. Production uses PySPH NNPS; an independent all-pair graph
is appropriate for 27 nodes and is not a second simulation backend.

The first attempted launch failed during import, before evaluating the ODE:
the existing PySPH environment has no SciPy. Initial import-failure logs are
retained. Before the first actual ODE execution, the method was amended to
independently coded classical RK4 in two exact smooth halves, with 128 steps
per half (7.8125 ns) and 256 steps per half (3.90625 ns). There is no HHT,
trial-map, backend integrator or force reuse. This is implementation repair
round 1; physical inputs and accuracy allowances are unchanged. No library
was installed. Save both references at all 62.5 ns common times,
preserving the left/right force limit at the exact switch.

Finest-grid local HHT accuracy budgets fixed before execution:
volume-weighted RMS x <=1e-9 m, v<=2e-4 m/s, full stress<=7000 Pa;
each face force component<=.7 N; each face impulse<=1.4e-7 N s;
each face wall work and raw friction Q<=2e-7 J; each signed reservoir
heat<=2e-8 J; maximum T error<=1e-7 K. The two independent references must
agree to <=1% of each allowance. Plastic work must remain <=1e-14 J and
ep<=1e-15; positive plastic work accuracy is NOT_TESTED. Require decreasing
temporal error and fine/coarse ratio<=.7 when error exceeds ten times oracle
uncertainty. BE heat is first-order; do not require second-order heat accuracy.

Keep moment, finite state, density, temperature, penetration, ID, mass,
per-particle BE heat, impulse and signed work guards; no existing production
threshold is relaxed. Raw continuum work/thermal integrals are reported
separately. Current-volume hypoelastic total energy remains a diagnostic;
its residual is not HHT algorithmic dissipation and must not become heat.

This reference establishes temporal accuracy of the existing spatial model.
It does not erase E2's 27.53% continuum stress error or turn continuum
traction into exact per-node equilibrium. Force cancellation at every free
DOF, oracle calls into the tested HHT/trial map, or a finer copy of that same
map would be tautological and are excluded.

## Cost, evidence and stop rules

Run the independent reference/test suite under an external 60 s timeout.
The reference also has a 55 s internal total deadline across both tolerance
settings. Record whole execution time, RHS evaluations, source/contract hash,
configuration hash, raw histories and reference-refinement differences.
Any reference or coupled-driver failure is retained. At most two authorized
evidence-based repair rounds; no tolerance adjustment using coarse results.
No actual AFSD run is included; its separate forecast and 600 s ceiling remain.

Initial status at registration: reference NOT_RUN; joint HHT NOT_RUN.
Actual results are to be written under `logs/process_hht_20261003/reference*`.
