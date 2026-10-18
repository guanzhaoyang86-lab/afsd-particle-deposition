# Small-angle continuous-path proof and explicit continuation

The user requested continuation after the v4 event-after-path budget stop. This
revision resumes the complete v4 first piece only, at 0.07971112060546842 s and
step 457, checkpoint SHA
e93edb7ca310b551f7a76c2f3636ce5ef8b3172bc64cf7506622001bbf094793.
All 207 existing sources and all failed candidates remain unchanged. The failed
second piece and its event probes cannot be promoted as continuation parents.

The original failure is on the declared continuous retained-edge path after
activation of 1091--1173. It is unresolved within geometric enclosure/query
limits, not evidence of detected penetration. Its old arithmetic position bound
uses a rotation-matrix norm estimate uniform over a much larger angular range.
The supplemental proof uses the actual complete path angle to bound the same
Rodrigues polynomial. For 0 <= t <= theta, abs(sin(t)) <= theta and
abs(1-cos(t)) <= theta^2/2. Hence the inherited Frobenius-norm coefficient can be
bounded by ||R0||F (sqrt(3) + theta ||K||F + theta^2 ||K^2||F / 2).
All products and bounds must retain outward engineering rounding allowances.
This is a tighter derived enclosure, not permission to ignore geometry errors.
The small value of (1-cos(t)) does not justify scaling its absolute arithmetic
error to zero. With eta = 8192 binary64 epsilon and inherited derivative
coefficient C = ||R0||F (||K||F + ||K^2||F), use an explicit position allowance
E = eta A_theta (S + rho) + eta C (1 + theta) rho. The second term separately
retains absolute trigonometric and angle-argument uncertainty. The supplemental
position margin includes the original endpoint/rotation-defect margin plus 2E.

The new proof provider retains the original strict intersection, numerical
domain, degenerate-edge checks, curvature bounds, fixed edge inventory and
continuous path. The added enclosure is sufficient only when strictly separated;
uncertain cases remain unresolved. Original query and depth limits remain.
No coordinate, material state, tool pose, force, thermal graph or constitutive
parameter is adjusted to obtain a certificate. This remains an engineering
binary64 enclosure rather than a formal interval-arithmetic proof.

The new CPU driver changes only the explicitly named path provider and path
protocol. The new receipt reader changes only its explicit driver import.
All nonlinear solvers, event locator and blocked-prefix rules, activation-delay
tolerances, original coarse/fine and physical checks remain. Source AST tests
must establish this boundary. The saved failing chord must be rechecked with
identical inputs and the new proof; this read-only check does not advance physics.
Manufactured intersection/domain/degeneracy and conservative enclosure tests,
plus independent review of the error bound, precede any new physical launch.

The new administrative entry restores the three remaining absolute boundaries
from the saved v4 scheduler and independently reconstructs them from the complete
parent. Nominal macro dt remains 118.65234375 microseconds. It retains local
subdivision without a cumulative refinement count and does not broaden the old
multi-opening failure classifier to unrelated errors. Cached cost evidence is
SHA-pinned, and the complete immediate parent is reread using the old verifier.
Original 450-second forecast, 540/600-second child limits, 650-second launch
reserve, user pause and 0.22-second stationary target remain binding.

User continuation is interpreted as the previous one-hour campaign duration.
Actual launch must specify --wall-seconds 3600 and freeze one start/deadline
configuration. All children inherit that deadline; preparation does not create
a configuration or a process. No solver is launched by test/import/inspect.
The launch must run hidden and retain a detached process with kill-on-job-close
disabled. Complete pairs are checkpointed; partial or failed pairs are preserved
as evidence only. A startup or geometry-only PASS cannot be reported as new
physical progress or deposition. Real continuation remains NOT_RUN until launch.
