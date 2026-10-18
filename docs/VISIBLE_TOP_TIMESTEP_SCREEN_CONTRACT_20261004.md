# TOP event: larger physical time-step screening

Registered before these physical evaluations. This is a diagnostic temporal
screen of the unchanged HOT_QS_VISIBLE_SUPPORT_CHAIN_RULE_CONTACT_APPROXIMATE_V1
model, not a new accepted trajectory or a physical validation.

Only after the complete 42.875--44.875 ms / 32 x 62.5 us TOP window and its
independent readback PASS, replay the *same true 42.875 ms parent* with 250 us,
500 us and 1000 us respectively. Each is an independent bounded stage with its
own log/receipt; candidate branches remain NOT_RESTARTABLE. No old state or
failed candidate is promoted or overwritten. Trial outputs are compared against
the actual accepted 62.5 us trajectory at every shared endpoint, not merely the
last state. Sources, physical parameters, full parent energy/contact/material
history, passive IDs and existing model transition are unchanged.

Every trial increment retains fresh full-map acceptance, scalar/batch equality,
complete interval energy prefix, geometry/penetration/temperature/density gates,
finite contact-shell enclosure and retained-edge whole-path certificate.
Only SHOULDER and qualified TOP non-bore events are permitted. All reference
fine increments inside a trial increment contribute event/blocked/crossed IDs
and raw/filtered graph neighbors to the cumulative local screen.

The pre-existing global and local comparison limits apply unchanged: 5 um
position, 0.5 K temperature, 10% stress (local floor 1 MPa), 20% plastic strain
(local floor 1e-5). Strict crossing ID lists must match at every endpoint.
Save per-increment comparisons, reference IDs, hashes, raw diagnostics and
candidate metadata; no accepted checkpoint is written. TOP must be observed
in each completed diagnostic branch. This only demonstrates agreement with one
smaller-step reference on this event interval; it is not proof of asymptotic
convergence, later-event reliability, or contact-history convergence.

Each stage predicts total cost as 30 s + 1.5 * count * prior full fine-stage
cost-per-step * (dt / 62.5 us), also taking the maximum measured root cost.
Forecast must be <=450 s. Runtime limits remain 540 s internal/600 s external,
single-thread, CPU only. Larger dt changes discretization, not material time,
mass, modulus, temperature or prescribed tool/feed speeds. Abort on any failed
gate, preserve failure, and do not automatically change thresholds or retry
with additional physics. A successful study may inform a separately registered
subsequent same-parent coarse/fine continuation; study states cannot restart it.
