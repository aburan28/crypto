#set page(paper: "a4", margin: 18mm)
#set text(font: "Helvetica Neue", size: 9pt)
#set heading(numbering: "1.")

= Genus-128 GHS source-lock refresh

Run `run-20261015T0929Z` · 15 October 2026

The current-main curve and Jacobian-transfer registries were copied byte-for-byte into the branch before rebuilding the WeakCurves snapshot and native GHS audit. Their SHA-256 hashes are `85300a415cfaef50e561a09fe8e0be56fb2b6ecc93babf5906a0d0d91c90d5fc` and `ddc14e19741f017ef3880d004906f355053abd1499f23a8ad8f433642cb2f21d`. The prior clean Linux PR merge checkout had used these newer source bytes and rejected the older snapshot; its failed job log is retained.

#figure(
  image("ghs_cost_sweep.svg", width: 100%),
  caption: [Recorded conditional cost sweep after the source-lock refresh. Both named routes select degree bound 15.],
)

#table(
  columns: (1.4fr, 1fr, 1fr, 1fr), inset: 5pt,
  [Source curve], [Model log2 queries], [Model log2 total], [Nominal base],
  [`c2pnb176w1`], [51.872229], [71.336033], [48,921,052],
  [`c2pnb272w1`], [51.872229], [71.773568], [48,921,052],
)

The 84-bound sweep is numerically identical to the first run. Only the source hashes changed. These are cost-model outputs under the half-split-place, polynomial-smoothness, and operation-price assumptions; the current receipt retains null fields for actual usable places, ordinary-query relation yield, useful rank, and target descent.

Evidence: `research/weak_curve_leaderboard_20261012/ghs_genus128_cost_audit.json`; first-run receipt in `reports/run-20261015T0854Z/`; clean-runner failure in `validation/ghs128_ci_snapshot_mismatch_attempt1.log`. Reproduction: native `--check`, catalog snapshot `check`, graph and decomposition index `--check`, and the 81-file SHA-256 manifest. This refresh changes provenance only and adds no curve-level finding or explicit isogeny map.
