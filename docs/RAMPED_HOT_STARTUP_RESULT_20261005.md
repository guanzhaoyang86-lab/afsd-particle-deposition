# Independent smooth-start diagnostic

The full deposition project is not complete. The registered experiment has
terminated with a temporal FAIL at 3.5 ms. Smooth startup passed the first saved
tool-contact endpoint at 3 ms, but did not qualify the subsequent shared
500 us interval. No further refinement or continuation was launched.

This experiment changes the
startup boundary history, independently initializes at 723.15 K and time zero,
and never advances either previously accepted branch. It retains the pinned
PySPH backend and existing material/contact/thermal laws and tolerances.

The immediate motivation is the expanded-domain window002 temporal failure:
PID893 stress difference 10.3881% and plastic-strain difference 21.9901%.
The accepted expanded state remains 68.359375 microseconds / 20 steps with zero
crossings. A saved contact-tangent scale audit did not support a 1 nm contact
branch-switch explanation. Details remain in the current progress report.

New code: `scripts/ramped_hot_motion.py`, `ramped_hot_trial.py`,
`ramped_hot_candidate.py` and `run_ramped_hot_startup.py`. Startup velocity
smoothly reaches nominal values at 20 ms. Exact integrated endpoints and one
explicit piecewise linear/geodesic step path are shared by the full map,
tangent, path certificates and passive crossing observer. This is a new
engineering boundary assumption, not recovery of the paper's startup schedule.

Eight material-free motion/interface checks passed in
`outputs/ramped_hot_motion_20261005/v1/receipt.json`, SHA
`b0110de0d92d0ed364862bbe8a5acb8b2eaafc3b1e2e16131d7956cd5e6e41f8`.
They include a crossing counterexample detecting misuse of the global smooth
path inside a piecewise step. The separate final static wiring review passed.
These checks do not demonstrate coupled loaded feasibility or deposition.

The registered experiment is at most four 1 ms chunks with separate continuous
coarse 500 us and fine 250 us paths from common independent initial arrays.
Every shared endpoint retains original global/pointwise thresholds, applied
to all 918 IDs. No automatic refinement follows failure. Per-stage caps remain
540/600 seconds and forecast <=450 seconds; measured cumulative cost is bounded
by 1200 seconds. The exact contract is
`docs/RAMPED_HOT_STARTUP_DIAGNOSTIC_CONTRACT_20261005.md`.

Actual output root: `outputs/ramped_hot_startup_20261005/v1`.
The initial elastic-rank/same-map/cost stage passed: 918 particles / 1961 free
DOF; full map 0.154499 s and tangent 0.495926 s. Stable loaded timestep remains
unestablished. The first 1 ms chunk completed two coarse and four fine roots,
with stage wall times 7.436401 s and 15.229837 s. Both shared-endpoint temporal
screens passed. Tool contact load and new plastic strain were both absent;
this is only an elastic startup result before tool loading; prescribed axial
feed still performs work. The second chunk also passed, still before tool
contact. At 2.5 ms the coarse path developed plasticity from prescribed axial
feed, without tool contact. At 3 ms both paths had actual tool contact load and
new plastic strain. The 3 ms local stress/plastic differences were 1.66454% /
5.00470%, within the unchanged 10% / 20% limits.

The last coarse stage completed to 4 ms. The fine stage stopped immediately
after the failed shared 3.5 ms endpoint; its final two planned steps were
not run. At that endpoint:

| Field | Measured difference | Original limit |
|---|---:|---:|
| Global stress | 17.3082% | 10% |
| Global plastic strain | 30.8171% | 20% |
| Largest local stress (PID791) | 69.2380% | 10% |
| Largest local plastic strain (PID791) | 67.1944% | 20% |
| Largest position difference | 0.345834 micrometres | 5 micrometres |
| Largest temperature difference | 0.00128471 K | 0.5 K |

This is a time-resolution failure, not a failed force or energy gate. The
coarse 3.0–3.5 ms candidate exhausted its Broyden budget and then passed with
the original same-parent Newton retry. The two fine candidates passed without
that retry. Their complete-candidate costs were 33.5370 s versus 15.6897 s and
10.2687 s. At 3.5 ms cumulative absolute raw energy defects were 291.669 nJ
(coarse) and 272.525 nJ (fine), below 551.939 nJ and 527.341 nJ respectively.
All complete time-zero histories were retained. No residual was removed or
converted into heat. In contrast, the first 1 ms passed chiefly because of the
absolute energy allowance: raw / absolute-work-scale was about 28% coarse and
18% fine. Those early passes must not be called good relative energy accuracy.

The actual run completed 22 candidates (8 coarse, 14 fine), plus the initial
cost/rank stage. Total supervised elapsed time was 192.793469 s. Each stage
stayed below its registered single-run bound and the cumulative experiment
remained below 1200 s. The independent saved-only reader passed all four chunk
audits, including faithful reconstruction of the terminal failure. Its final
`temporal_screens_all_pass` and `loaded_feasibility_screen_passed` are both false.
Readback PASS never changes this physical-study outcome.
Final independent reader: `chunk_004/readback/summary.json`, SHA
`25746d994a5b2d98c68c3e6929dee005b507e3f3527f52ae23ed154b815a50be`.

The failed final fine state contains actual fields: maximum net displacement
0.874291 micrometres, maximum plastic strain 0.0002751054, maximum equivalent
stress 1.517378 MPa, temperature 723.15–723.155383 K, and crossing mass zero.
It is `chunk_004/fine/step_002/NOT_RESTARTABLE_candidate.npz`, SHA
`c8b406f6542785abbbbd2f6d4acbf0d2ad24bd8f51da40f6ec487d58d7efd4f4`.
These are diagnostic observations, not accepted deposition results.

The failure does not establish a unique cause. A separate saved-field
PID791 review examines material and contact branches without rerunning force,
thermal, NNPS or root calculations. It found that PID791 has no direct tool
contact in either path over 2.5–3.5 ms. The fine path adds plastic strain at
3.25 ms and unloads without new plastic strain at 3.5 ms; the coarse 3–3.5 ms
step adds none for that particle. The newly loaded neighboring IDs796/797
belong to its unchanged support graph. Crucially, every particle in both paths
remains below the FKS 1/s rate floor through the failed 3.5 ms endpoint; no
floor crossing occurs. The prior constant-speed PID797/893 rate-transition
explanation cannot be transferred to this failure. The evidence supports
further review of contact-onset time integration and the coupled loading path,
not a unique proven cause or an unqualified material substitution.
The FKS rate floor remains applied throughout that interval; its multiplier
is one. PID791's comparison floors are also active and reduce the reported
relative differences. Its absolute stress/plastic differences are 692379.57 Pa
and 6.71944e-6, exceeding the associated 100000 Pa / 2e-6 allowances.
In total, 228 particle stress comparisons and 92 plastic comparisons fail.
Saved-only forensic receipts are `method_review/summary.json` (SHA
`d27a4a84dfe5c7fd716cfc64d48020ce705afdbf92e9d2320378f80d4aedef54`)
and `method_review/rate_floor_inventory/summary.json` (SHA
`6a947d5ca628e2f9797c2f5555e33c287e0284b633445df63710c79142b175dd`).

The evidence rejects the proposition
that this 20 ms smooth startup alone makes the registered 500/250 us schedule
usable through the loaded onset. It neither excludes every startup schedule
nor proves a narrower timestep would remain affordable at full speed.

No already-qualified material-energy replacement was found in the existing
project source. The contact discrete-gradient implementation covers normal
contact work only; the energy breakdown is diagnostic; the existing midpoint
thermal/dynamic variants explicitly lack material-energy certification. Old
counterexamples to average-stress plastic work remain. None was silently
substituted to turn this temporal failure into a pass.

All generated states remain diagnostic, accepted=false and restart_allowed=false.
No old FAIL was removed and no current full regression or experimental
validation is claimed.

Actual saved-field figure (PID791, both paths through 3.5 ms only):
[PNG](../outputs/ramped_hot_startup_20261005/v1/visualization/pid791_ramp_contact_failure_v2.png)
and [PDF](../outputs/ramped_hot_startup_20261005/v1/visualization/pid791_ramp_contact_failure_v2.pdf).
The figure distinguishes the first saved system tool-load endpoint from an
exact contact-onset time. PID791 itself remains out of direct contact. It
omits the unevaluated initial strain-rate scratch field and joins actual
saved points only; it does not illustrate an accepted material-flow trajectory.
