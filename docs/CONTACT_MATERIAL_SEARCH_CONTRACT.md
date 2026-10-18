# Contact/material search integration, 2026-10-04

Preregistered ENGINEERING_APPROXIMATE numerical experiment. The current user
request authorizes continued bounded integration; older zero-evaluation work
package limits are historical, not a restriction on this new package.

No physical model is replaced. The unchanged ProcessThermalTrial, kick, FKS,
mass, modulus, control, thermal sources, and acceptance decide every candidate.
Old evidence and source files remain intact. No GPU, production, or three cases.

## Search approximation

CONTACT_FORCE_SECANT_APPROXIMATE_NOT_JACOBIAN uses two actual kicks per free
force axis, at a fixed +/-1e-3 characteristic force. The scale is the maximum
of internal-force norm, total friction capacity, and m*relative_speed/dt.
No-contact blocks are analytic. Prescribed columns and rows are zero.
Both event mean B and endpoint C are measured separately; B=C/2 is not assumed.
The saved baseline is replayed and compared. Each perturbed kick retains its
signed heat certificate and branch diagnostic. A changed branch is reported,
not called a differentiable Jacobian. There is no epsilon search/refinement,
force smoothing, eigenvalue clipping, physical damping, or response-as-heat.

The contact blocks are frozen at the saved rejected 500 ns numerical candidate,
whose physical parent is the published 850.393632910904 ns accepted state.
The rejected candidate is ONLY a search-matrix anchor, never an initial state.
Each actual trial refreshes the existing rate-aware material/pressure blocks
from its own midpoint journal. With Jv,Jr,JT,kappa,E from MaterialStageResponse:

    H = Jv + Jr diag(kappa) E
    A = I - B H
    rhs = dbar + B (Jr drho + JT dT)
    A delta_bar = rhs
    delta_rho = drho + kappa E delta_bar
    delta_force = Jv delta_bar + Jr delta_rho + JT dT
    delta_x = dx + dt B delta_force
    delta_v = dv + C delta_force

All constrained columns of H and E remain. Only physical free output axes
respond to force. Constrained increments equal the raw prescribed defect.
Temperature correction stays dT. Geometry, normal-force, wall-motion and BE
source derivatives remain omitted; neither full Newton nor accuracy is claimed.
Binding rejects stale candidate/map/parent/dt and invalidates on failed trial.
No physical evaluate inside apply. Matrix solves retain backward error <=1e-12.

## Independent gates and budget

Before real evaluation: analytic no-contact response; sliding then sticking
(B nonzero, C zero); reversal; two faces and prescribed-axis reaction; explicit
branch crossing; literal uneliminated coupled algebra; ordering, pure state and
stale context rejection. Analytic smooth secants: rtol 2e-5, atol 2e-8 s/kg
(central truncation at fixed h); algebra: rtol 2e-12, atol 2e-13 in SI test units.
Original material, process thermal, root controller and event tests are rerun.
Tests are components, not AFSD validation.

1. Offline contact response on saved 342-particle journal: external <=240 s,
   no new full trial. Save response blocks, per-particle evidence and timings.
2. One uncommitted full 500 ns cost trial from the latest accepted parent,
   including refreshed material assembly, proposal and complete journal IO:
   external <=90 s. Report memory, particles, dt and measured cost.
3. If 1.5*(32*measured trial/proposal/journal cost)+30 <=450 s, at most one
   500 ns root: <=32 physical evaluations, inner540/external600 s. Original
   beta .75/history6/raw-growth2 controller and physical acceptance unchanged.
   No dt shrinking or retry from rejected states. Failed attempts stay on disk.
   Passed root produces checkpoint, fields, tool pose, transport and energy.

This experiment cannot itself produce passage: the unchanged 50 um displacement
guard and initial 250 um shoulder gap remain. A root PASS is only a numerical
integration milestone. Report the cost gap to 0.1179 s nominal feed length;
do not extrapolate a root into a deposition claim. Finite-tool process integration
and an independently initialized slow hot model remain separate work.

## Preflight correction before any new trial

The first supervised suite reproduced the old test_event_response continuous
ODE mean-response failure (10 PASS, 1 FAIL), while the other five suites passed
64 tests. The first failed preflight is retained in outputs/.../tests and logs.
This is the exact failure already qualified by check_process_event.suite_result;
use that unchanged narrow parser, not an unrestricted nonzero-exit exemption.
The new preflight is tests_review1, with status explicitly containing
RETAINED_ACCURACY_FAIL. Original tolerance, test, contact source, and physical
acceptance are unchanged. Continued work is only numerical discrete-root
integration under the existing engineering approximation contract, not ODE
accuracy or process validation. No new full physical trial preceded this note.
