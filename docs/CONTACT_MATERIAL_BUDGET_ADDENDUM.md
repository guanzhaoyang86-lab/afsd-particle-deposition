# Evaluation cap reduction before the sole real root

2026-10-04, following the completed cost probe, before starting a root.
The fixed 32-evaluation plan did not satisfy its 450 s launch forecast:
1.5*32*9.0645269+30 = 465.0972912 s. Its root_allowed=false is retained.

The separate run_contact_material_root.py consumes the same frozen component,
contact and full-trial cost evidence. It reduces the cap to
min(32, floor((450-30)/(1.5*cost))) = 30 evaluations, forecast437.9037105 s.
Internal540/external600 s, dt500ns, beta.75/history6 and all acceptance unchanged.
There is one attempt only; no retry/dt shrink or starting from the rejected state.
The cap arithmetic is checked independently before launch. Original gate files,
source hashes and failed preflight remain unchanged. No process accuracy claim.
