#set page(width: 297mm, height: 210mm, margin: (x: 12mm, y: 9mm))
#set text(font: "Helvetica Neue", size: 8pt)
#set heading(numbering: none)

= N83 dense-root 64-point panel
2026-10-13 10:18–10:50 UTC · same public points · Apple M4 Pro shared host

K1798 dense met the frozen lowest-mean-probe rule among completed arms: 169,069.8125 probes per target across 64 points, versus 300,021.90625 for K1500 dense and 371,728.25 for both K1271 layouts. All 256 target rows in the four completed arms passed independent scalar and four-point relation replay. The K1500 regular control again reached its 900 s cap.

#figure(image("panel.svg", width: 100%), caption: [Recorded mean probes, packed-slot allocation and exploratory online time. Run ID: run-20261013T1018Z. The capped K1500 regular arm has no target values.])

#table(
  columns: (auto, auto, auto, auto, auto, auto, auto),
  inset: 3pt,
  stroke: 0.4pt + rgb("b7c5d3"),
  [Arm], [Actual B], [Columns], [Slots], [Mean probes], [Online ms], [Status],
  [K1271 regular], [210,986], [1,271], [8 GiB], [371,728.25], [171.706], [64/64 replay],
  [K1271 dense], [210,986], [1,271], [4 GiB], [371,728.25], [151.812], [64/64 replay],
  [K1500 dense], [249,000], [1,500], [8 GiB], [300,021.91], [129.283], [64/64 replay],
  [K1500 regular], [missing], [missing], [missing measured], [missing], [missing], [timeout 124],
  [K1798 dense], [298,468], [1,798], [8 GiB], [169,069.81], [76.774], [64/64 replay],
)

The producer binary, rank seed 29, 12 rank threads, one target thread and 900 s cap were fixed across arms. Its point input SHA-256 was `e9ed7ddb0aeafdbfc7354d67279a8712572eea9d2aee3896f99633bdc2e937cc`. The separate known-scalar sidecar was used only by the native verifier. The K1798 extension was frozen before the panel started. Online times are exploratory on a shared host and are secondary to the one-target IC/rho run.

#pagebreak()

= Paired probe differences and replay

#figure(image("panel_paired.svg", width: 82%), caption: [Paired mean target-probe differences, left minus right. The intervals are 10,000-resample 95% percentile bootstrap intervals. Positive values favor the right-hand arm.])

K1500 dense minus K1798 dense: 130,952.09375 probes, interval [52,057.703125, 209,990.828125]. K1271 dense minus K1798 dense: 202,658.4375, interval [121,448.546875, 288,107.375]. K1271 dense minus K1500 dense: 71,706.34375, interval [−24,081.171875, 168,322.5]. K1271 regular and dense had exactly the same probe count and first witness on all 64 points. The K1500 regular/dense parity comparison remains censored by regular-arm timeout.

The native analyzer independently checked each completed arm's recovered scalar with the generic curve library; each archived base's representative subgroup membership, 166-point signed Frobenius orbit and label cycle, and point distinctness; and all 256 four-point target sums. All checks passed. Maximum online phase residual was 0.000501 ms. The three base-dump SHA-256 digests and raw run paths are listed in the Markdown report beside this PDF. Selection is by the frozen mean-probe rule among completed verified arms; no CPU isolation receipt is attached to these wall times.

Evidence: `../../analysis_panel.json`, `../../runs/panel/`, `../../PROTOCOL.md`, `../../CAPACITY_CEILING.md`, `../../PROOF_MAP.md`. The single-target same-point rho comparison is in `../run-20261013T0942Z/REPORT.pdf`.
