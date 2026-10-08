# IC dashboard review, 5 October 2026

The [desktop readiness view](desktop-readiness-pipeline.png) shows the new
512-query F5/SAT outcome chart beside the pipeline admission table. The
[mobile overview](mobile-overview.png) shows the current decision and compact
navigation at 390px. These are screenshots of the generated
[scoreboard](../../index-calculus-scoreboard.html), not new measurement data.

The [browser receipt](browser-check.json) pins the page SHA-256 and records
passing desktop, 390px and 320px layout, dark mode, no-JavaScript rendering,
evidence search and preserved deep links. The renderer also passed
`node tools/render_ic_dashboard.mjs --check`, which verifies evidence hashes
and that the historical ledger remains unchanged. The 512-query counts come
from the source-pinned stage report; the selected two-instance SAT pilot is
explicitly excluded from natural-yield estimates.

The later [updated readiness view](desktop-readiness-updated.png) adds the
independently audited one-million-conflict SAT preparation. Its
[mobile view](mobile-readiness-updated.png) and [updated browser receipt](browser-check-updated.json)
check that the three distinct solver-budget rows remain readable at 390px.
The [updated overview](desktop-overview-updated.png) leads with the actual
decision: SAT has all 29 reusable logs and passing disclosed marked-stdin
transport, but no fresh target result or speed claim. The added counts and
transport status are pinned to their separate original audits; the
earlier 100,000-conflict SAT row remains visible as a separate experiment.
