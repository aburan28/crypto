# Next gate: isolated transport attribution and portable price

This is a pending protocol, not an executed run. It depends on the verified
[four-size local comparison](REPORT.md). It does not alter any V1–V3 record.

## Hypothesis and fixed controls

On the same four fields and eight V3 source targets per field, add a public
base-field isomorphism-only transport arm using
`ψ₂(x,y)=(4x,8y)` from `E:y²=x³+1` to `E_iso:y²=x³+64`.
Verify the inverse `ψ₂⁻¹`, target and generator images, curve order and
scalar checks with independent arithmetic. Keep the exact V3 degree-three
destination arm `E':y²=x³−27`, source relation arm, and negation rho.
Use a duplicate source-relation A/A arm as well as duplicate rho, with
balanced order across targets and rounds. Keep twelve factor-base points,
both relation streams, recovery, and the same collection/descent caps.
Use at least five measured rounds and one warmup. Freeze an exact executable
spec and source/input hashes before execution; retain failures and restarts.

The generic two-class step floor remains `sqrt(π/4)` in counted-group `S`;
matched negation rho on the same original source target is the reference.
The falsifiable mechanism question is whether degree-three destination
relations beat both direct source relations and the isomorphism-only model
after complete cold pricing. A source/destination difference that lies
inside relation A/A noise does not distinguish the mechanism.

## Pricing and release gate

Run one paired native session on an L2-or-better Linux isolab worker with
whole-core placement, calibrated native conversion for map field work,
roots, lookups, modular rows, rho canonicalization and hash operations,
and recorded memory/PMU vectors. Price each target cold, including every
arm's factor base, pair table, map, descent and verification. Report curve
selection/map construction separately from the public-input solve, and
charge any reusable setup by the exact number of targets that reuse it.
Keep a direct isolated wall/CPU table alongside the portable counted unit.

Require all original-source scalars verified, independent relation replay,
no unpriced unit in the portable `S` total, both A/A controls within their
frozen noise bands, a same-session exact replay certificate, and a second
host's independent full replay. Report paired target intervals per size;
do not pool the four sizes into a scaling claim. Stop and preserve an
attempt when a correctness, calibration, isolation, or independent-host
gate fails. The current environment has no configured isolab client, so
execution is pending an available worker and credentials; no worker job
has been submitted.
