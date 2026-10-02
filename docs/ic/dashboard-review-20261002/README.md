# IC dashboard readability and verification

The canonical scoreboard now opens with a bounded-autolab overview: current
outcome, paired tournament and rho graphs, six-step preparation/online diagram,
F5 cost ledger, next gates and searchable historical evidence. Every original
ledger-body byte is retained inside a native collapsed details element. Search
and old hash links open that ledger. With JavaScript disabled, the overview,
charts and full ledger disclosure still work; browser Find supplies search.

The graphs cite frozen round-three confirmation data. Per-cell dots are estimates,
not invented confidence intervals; the overall interval and stricter familywise
upper bound are distinct. The qualified `ic_online` denominator differs from the
incumbent cold role. The rho graph uses the actual `rho_online` role and paired
single-target costs; it is neither a new F5/SAT comparison nor batch amortization.
F5's new disclosed control is explicitly uncalibrated; its matched-rho speedup
stays unknown. This front page covers the bounded autolab, not all repository
research. The dated library panels retain their individual questions and scopes.

[Overview data](../dashboard-overview-data.json) includes the exact source files
and their content hashes. Regenerate the static overview with:

```sh
python3.12 tools/render_ic_dashboard.py
```

The renderer asserts the selected accepted decision/control states; do not point
it at a different outcome and silently keep these labels. When a new gate lands,
update the explicit source selection, state descriptions and date together.
The browser performs navigation only, not statistical recomputation. Existing
research updaters may continue adding their historical panel markers; place new
panels within the `legacy-evidence` container. Keep source reports immutable.

The retained [local browser checks](browser-checks.json) pass with HTTP(S)
requests blocked: collapsed landing, two data-bound graphs, positive/empty
search, individual report links, direct hash opening, desktop/mobile overflow
and no JavaScript exceptions. Screenshots show [desktop light](desktop-light.png),
[desktop dark](desktop-dark.png), [mobile labels](mobile-chart.png) and the
[pipeline and control ledger](pipeline.png). Visual inspection enlarged the
mobile SVG labels before the final capture.

[browser-check-source.mjs](browser-check-source.mjs) is the exact local CDP check
source. It connects only to the separate headless Chrome profile under private
tmp; adapt its explicit local paths for another checkout. It uses Node's built-in
WebSocket and installs no browser dependency. It performs no native scientific
solver or experiment. [preservation.json](preservation.json) records exact
historical-body preservation and idempotent rendering. Skill validation passes
with the existing Python 3.13/PyYAML environment; the initial Python 3.12 invocation
lacked PyYAML and was not counted as a pass. All scientific F5 replay is recorded
separately in the consumed control's result.
