# Independent array audit of the real-process candidate map

2026-10-03, before executing the new checker or the four-candidate bridge.
API: `scripts.afsd_bridge_audit.audit(pa, p, t0, dt, evaluation)` returns a
strict-JSON-compatible dictionary with `passed`, `status`, aggregated `checks`,
`summary`, and always-false accepted/restart_allowed/can_publish. It reads
only parent arrays and the returned journal. No contact, material, thermal,
mechanics, driver or acceptance function is imported as an oracle. No input
is mutated; malformed/missing/nonfinite fields produce an explicit FAIL.

This is a MAP-internal accounting audit, not fixed-point acceptance. In
particular it uses `evaluation.map.temperature_K` for BE fluxes, not the
unconverged input C temperature, and does not borrow the status of
`certify_process_thermal`. Reentry of first/fourth C and archive/source/parent
hash checks remain the surrounding driver's responsibility.

The following SI bounds are copied from PROCESS_RESPONSE_NEXT_WORK_PACKAGE,
not selected after seeing a new failure:

- Momentum per component: m_i*1e-9 m/s + 1e-17 Ns. Free-axis reaction:
  1e-17 Ns; constrained mapped velocity: 1e-12 m/s.
- Contact/constraint/internal work: 1e-15 J + 1e-9 times the appropriate
  sum of absolute work terms. Individual raw friction identity:
  1e-25 J + 1e-12 relative. Negative Q is allowed only within its original,
  independently reconstructed per-event rounding certificate, never clipped.
- Plastic-work and return-dissipation algebra: 1e-18 J + 1e-10 relative.
  ep is monotone; reconstruct q=sqrt(1.5*S:S), P=(m/rhomid)*q*(ep1-ep0).
  Return dissipation 1.5*G*(m/rhomid)*(ep1-ep0)^2 stays separate from heat.
- BE balance: particle 1e-12 J + 1e-10*sum(abs(total_heat_input_i));
  global 1e-12 J + 1e-9 times that same scale. Pair net:
  1e-15 J + 1e-10*sum(abs(pair heat)). Existing exact-once ledger identities
  use 1e-25 J + 1e-12 relative; thermal source sums use stable scalar sums.
- Existing contact journal consistency: interval coverage 128 eps*dt,
  before/after velocity 1e-15 + 1e-13 relative; supplied tangent midpoint
  budget 128 eps*(speed_scale+|v0-frame|+|v1-frame|). Reconstruct the exact
  original 256 eps friction-roundoff certificate; matching tolerance is
  1e-30 J + 1e-12 relative, never a widened dissipation threshold.

For each trace and accepted interval, independently reconstruct Jn, Jt,
constraint impulse, midpoint work and signed Q per actual contact ID/feature.
Then sum to mapped momentum/work/heat; retain negative individual events even
when aggregates are positive. Handle zero-contact traces and prescribed axes
explicitly; labels do not enter any formula.

For BE, reconstruct the undirected radial pair conductance from recorded
geometry/raw gradient/mass/rhomid/conductivity, then q_ij=dt*g_ij*(Tj-Ti),
with opposite signs at its ends. Reconstruct reservoir coefficients from
the actual contact faces, source chi*P and partition*Q once, wall flux at
mapped T1, and each particle's m*cp*(T1-T0). Compare the individual journals,
total ledger, and raw balance. No return dissipation, diagnostic energy defect,
normal penalty error or root residual is inserted into heat.

`checks` aggregates count/max error/worst normalized error and first failure
context by category. Full journals remain in the surrounding driver's gzip;
the checker records negative raw events, interval kinds, contact IDs/features,
plastic IDs and FKS rate-floor branches in the JSON summary. Positive plastic,
multi-contact and actuator activity are reported, not manufactured as fixture
requirements: a zero-physics synthetic checker test may pass. The real bridge
driver must separately require the intended mechanisms to be exercised.

Planned independent synthetic checks: balanced no-contact free acceleration,
fixed-axis reaction/work, positive-plastic BE source and signed wall friction;
corrupt one source/work/impulse/ID at a time and require FAIL. No genuine
342-particle evaluation may be run by this checker subtask. Current new
checker execution status before implementation: NOT_RUN.

## Pre-real-evaluation refinement: cancellation journal only

The initial ten synthetic tests passed in 0.2284862 s; their exact source and
contract versions remain in `logs/process_bridge_20261003/audit_initial_sources`.
Before any real342 evaluation, independent review identified a distinct
arithmetic issue: `mc_delta_T_residual_J` is the difference of large terms,
not an independently supplied heat source. A relative bound based on that
tiny residual itself can falsely reject a legally reassociated sum.

Only this diagnostic journal comparison (per particle and its scalar ledger)
therefore uses the standard finite-arithmetic forward-error bound
`gamma_k=k*eps/(1-k*eps)` with a conservative explicit operation inventory:

```text
k_i = 64 + 8*(number of face/subinterval terms at i)
         + 32*(undirected thermal degree at i) + 8*(active faces at i)
scale_i = |m cp (Tmap-T0)| + |chi P_journal| + |eta Q|
          + dt*(|sum H Twall|+|sum H Tmap|)
          + sum incident |pair heat|
bound_i = 1e-30 J + gamma_k_i*scale_i
global_bound = sum(bound_i) + gamma_(2N+4)*sum|particle residual|
```

The counts cover scalar reductions of per-face friction, face reservoirs,
radial pair coefficients/nodal reductions and final finite source additions;
the formula rejects unsupported k*eps>=1. Stored source and map arrays are
the inputs, not uncertain continuum values. The small floor is only for an
identically zero/subnormal diagnostic. Every individual heat-source identity,
the total-input identity, and the physical particle/global BE balance gates
above are unchanged. This bound is not a way to inject residual heat or excuse
a physically unbalanced map. The conductor expression's multiplication/
division order was also checked against the existing source; different
friction/pair accumulation orders remain the reason not to assume exactness.

A new independent mixed-source fixture initially used T0=300 K, T1=300.425 K,
chi*P=0.3 J and eta*Q=0.125 J. Reassociating `measured-(Pheat+Qheat)` as
`(measured-Pheat)-Qheat` must pass within the derived bound, while adding
1e-6 J to that diagnostic must fail. The physical temperature and every
heat-source field stay fixed in this test. This is a preregistered audit-only
refinement, not a reaction to a real342 result or a physical tolerance change.

As with material P, exact-once thermal journal checks use the original
boundary raw-Q vector only AFTER it has passed the independent per-face,
per-interval and aggregate Q checks. Otherwise a legal difference between
sum(interval then face) and sum(face then interval) can leak into a tiny
diagnostic balance. Independent physical BE balance still uses the separately
reconstructed P and Q, not the copied input vector. This separates original
source-use verification from independent source reconstruction; no face-wise
negative event or forged source is hidden by the aggregate.

### Recorded fixture failure and bounded repair1, before real evaluation

`audit_forward_tests_initial` retained ten PASS and one FAIL: the new test's
assertion that the two reassociations differ was false for 0.3+0.125. It failed
before calling the audit on that case, so it supplied no audit/physics failure.
No checker tolerance or code is changed by this repair. Only the synthetic
mixed case now has eta=0.4, Q=0.25 J -> Qheat=0.1 J, Pheat=0.3 J and
T1=300.4 K. Independent scalar arithmetic before the repair run gives
`measured-(0.3+0.1)=-2.275957200481571e-14 J` and
`(measured-0.3)-0.1=-2.273181642920008e-14 J`, difference
2.7755575615628914e-17 J. The physical reference remains exactly balanced
in real arithmetic. The original source/contract/FAIL log are archived;
real process parameters and all physical bounds remain unchanged.

Repair1's eleven tests passed (external 0.2235798 s). Before any real
evaluation, a final source-level consistency review fixed the physical
balance scale to the original `certify_process_thermal` definition
`sum(abs(total_heat_input_i))`, not the sum of absolute individual heat
sources. The earlier work-package phrase "absolute sources" was ambiguous
under cancellation; the original source definition is authoritative. The
uncancelled magnitudes remain only in the independently derived diagnostic
roundoff bound and a separately named summary field. No physical gate is
widened. This final small consistency change is to be checked by the main
registered preflight; no additional real or standalone trajectory is run.
