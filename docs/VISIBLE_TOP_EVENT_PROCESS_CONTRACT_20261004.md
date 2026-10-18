# Visible-support TOP event and bounded continuation contract

Registered before the first physical TOP qualification run. ENGINEERING_APPROXIMATE,
342 particles, independently initialized 450 C hot quasistatic branch. The true
42.875 ms / step 487 parent is
`outputs/visible_deposition_20261004/window_001/fine/step_040/accepted.npz`, SHA256
`97eb0b8ed76892a8dd6c85100033c083a650ca1f0c3c8c854d32c877cdadcf1d`.
Its independent readback SHA256 is
`9e0866f0433c4321470a0a7dcb25a26c55e9ad6e70aec46ec123398892f2f3c0`.
Original coarse step 022 and fine step 042 failures, and uncommitted fine step 041,
remain unchanged and cannot serve as restart states.

The frozen physical residual already evaluates all finite tool patches. New
`visible_top_candidate.py` extends the qualified non-bore event domain from
SHOULDER to SHOULDER+TOP, with no equation, geometry, source, material parameter,
initialization, integration rule or tolerance change. OUTER_CYLINDER and BASE
remain excluded. This is an explicit event qualification, not a model transition.
The existing spatial-model transition, branch origin, all 487 energy records and
all particle and contact histories are preserved exactly.

Prerequisite: the 15 independently referenced TOP component tests recorded in
`outputs/visible_top_components_20261004/v1/summary.json` passed. Its SHA256 is
`7fbb5fdca371335ff528f4132fb840948dc3239de8eb50e420b3533014505cab`.
That receipt does not qualify a real coupled candidate. The first actual window
must observe TOP in both the coarse and accepted fine histories before claiming
the TOP event passed a process screen.

Each numbered window has 16 x 125 us coarse candidates and 32 x 62.5 us fine
increments from the same true accepted parent (2 ms). Coarse states are never
restartable. Only a contiguous fully qualified coarse prefix may be referenced
if coarse fails; original FAIL is retained. A pair of fine increments is published
only after both complete all gates and pass equal-time coarse/fine global and
cumulative event-neighborhood comparisons. TOP particle IDs join the existing
one-hop seeds alongside shoulder, blocked-edge and strict crossing IDs.

Unchanged global/local thresholds: position 5 um, temperature 0.5 K, stress 10%
with 1 MPa local denominator floor, plastic strain 20% with 1e-5 floor. Strict
crossing ID lists must match. Actual filtered graph, full scalar/batch map,
fresh equilibrium/thermal/material/energy acceptance, continuous finite-shell
certificate and whole retained-edge path certificate remain mandatory. No
negative-work, temperature, constitutive or conservation threshold is relaxed.

The first stage cost uses measured *whole-stage* costs from the preceding
42.875 ms continuation, including its failed event work. Forecast is
30 s + 1.5 x planned candidate count x max(measured full cost, slowest root).
It must be <=450 s; each physical stage has 540 s internal and 600 s external
limits. Threads remain 1, fixed PySPH commit, no GPU or scaling. Subsequent
windows reuse their immediately preceding measured full cost.

After each window, an independent saved-evidence readback is required. Only a
full 32-step accepted window with PASS readback can parent the next number;
any event/domain/numerical failure stops automatic continuation. This registers
at most 100 bounded windows (maximum 200 ms), each separately forecast and
verified; it is not a declaration of full deposition or a three-condition run.
The user has repeatedly requested continued short-process development. Once
actual newly crossed material exhibits sustained plastic spreading outside the
bore under the shoulder, stop this continuation series for assessment and
sensitivity planning. The observer never influences force or thermal response.

Independent readback checks hashes, original failure receipts, parent chain,
unchanged fields, complete energy prefix, saved geometry and graph certificates,
contact event qualification, pair comparisons and true accepted trajectory.
Any partial prefix is reported as partial. Component PASS, process numerical
screen, transport observation, spatial sensitivity and experimental validation
remain separate. Existing FKS/Jaumann/free-surface and limited-domain failures
and limitations are not resolved by passing this contract.
