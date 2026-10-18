# Bounded predictor-veto diagnostic (registered before execution)

2026-10-03. One private numerical root from the saved S0 stage008 at2us,
dt500ns, original physical inputs and24-evaluation budget. This is NOT an
accepted trajectory, not a new physical model, and not a retry at a smaller dt.

Question: the first predictor fails the sampled1% Coulomb interpolation check
before its equilibrium residual is formed. Does the unchanged discrete root
also violate that bound? This must be measured, not inferred from the name of
the exception. No old source or published output is edited.

The diagnostic temporarily wraps `_contact_work` inside its own process:
candidate exploration records the original1% defect but defers that one
sampled-direction veto. Exact power/passivity, event-zero, cone, geometry,
thermal and material checks remain active. Raw linear-force work and all
residuals remain identical. No friction formula, heat clipping, masses,
material, time step,24-evaluation limit or final acceptance bound changes.

If the root converges, the original unwrapped `_contact_work` is called again
with the original1% limit and all publication checks are evaluated. Regardless
of outcome, every returned state is saved ONLY as a non-restartable private
diagnostic candidate; no accepted history or accepted checkpoint is written.
The wrapper is removed even if the evaluation raises. The parent state and
source/artifact hashes are checked unchanged.

Possible results: FINAL_ORIGINAL_GUARDS_PASS (predictor veto was premature for
this interval); FINAL_ORIGINAL_GUARD_FAIL (root is also inadmissible); or
ROOT_NOT_CONVERGED (this diagnostic cannot decide). A final failure remains a
failure. No additional dt attempt, relaxed final gate or post-event advance.

Budget: one candidate interval, internal55s/external60s. Prior accepted stages
cost<.15s; conservative estimate<10s, within the registered limit. This narrow
numerical diagnosis does not authorize S1/S2,20us convergence,1ms or AFSD.
