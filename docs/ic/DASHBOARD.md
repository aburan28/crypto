# IC research dashboard

Open `docs/index-calculus-scoreboard.html` directly in a browser. The page is
self-contained: its overview styles, graphs and search are embedded, and it
requires no server, external fonts, package installation or solver execution.

The first screen answers the current decision. Subsequent sections show the
last challenger comparison and its uncertainty, the matched rho reference,
solver readiness, a complete pipeline diagram and the proposed tournament
process. Historical results by regime and the complete ledger are collapsed
under Evidence. Search and old fragment links open the relevant detail.

[Desktop preview](dashboard-preview-desktop.png) ·
[Phone comparison preview](dashboard-preview-mobile.png) ·
[Original local browser checks](dashboard-browser-check-20261003.json).

The primary question is the verified online cost of one unseen supplied point,
after reusable preparation, paired with rho on the same point and resources.
Preparation and cold costs remain separate. Disclosed controls and historical
multi-target results do not answer this question. Pending costs stay unknown.

## Update and verify

Use Node.js 22 or later; no dependencies are required:

```sh
node tools/render_ic_dashboard.mjs
node tools/render_ic_dashboard.mjs --check
node tools/check_ic_dashboard_browser.mjs
```

The browser check uses an installed Chrome or Chromium with a fresh temporary
profile. It checks chart visibility, the admitted SAT status, quiet evidence
search, old deep links, desktop and 390/320-pixel layouts, dark appearance,
JavaScript-disabled rendering and the absence of page network requests. It
does not use the user's browser profile. Set `IC_DASHBOARD_CHROME` to override
the executable and `IC_DASHBOARD_SCREENSHOTS` to save preview images and the
check receipt. Pull-request CI repeats the checks and retains screenshots.

The renderer reads the already accepted frozen evidence, checks its explicit
SHA-256 pins and evidence-level gates, and quotes recorded estimates and
confidence intervals. Chart coordinates and decimal formatting do not calculate
new statistics. It does not start research, dispatch solvers or modify run
artifacts. Updating a pin or claim requires review of the underlying evidence.

Edit `docs/ic/dashboard-overview.css` and
`docs/ic/dashboard-overview-interactions.js` for presentation. The renderer
embeds both into the standalone page. Edit the HTML template in the renderer
for overview content. Regenerate the page and its source-bound overview data
in the same PR.

The renderer retains the historical evidence and existing regime summary
verbatim. The legacy collapse handler is scoped to the ledger so it cannot
hide the overview's graphs. It records the ledger and summary hashes in the overview data and checks
preservation during rendering. Old anchors must continue to work. Do not use
the overview renderer to delete superseded evidence, turn an incomplete result
into a win, or compare diagnostic times from differently qualified controls.
