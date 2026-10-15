#set page(paper: "a4", margin: 18mm)
#set text(font: "Helvetica Neue", size: 9pt)
#set heading(numbering: "1.")

= Exact-integer genus-128 GHS model audit

Run `run-20261015T0854Z` · 15 October 2026

The native audit replays the source-locked conditional cost model for the GF(4) Type II, magic-eight routes of `c2pnb176w1` and `c2pnb272w1`. It counts monic irreducibles and degree-128 smooth monic polynomials exactly across 84 factor-degree bounds. Both routes select bound 15 and 48,921,052 nominal half-split places.

#figure(
  image("ghs_cost_sweep.svg", width: 100%),
  caption: [Recorded model cost by factor-degree bound. Two named curves; no measured relation yield.],
)

#table(
  columns: (1.4fr, 1fr, 1fr, 1fr),
  inset: 5pt,
  [Source curve], [Model log2 queries], [Model log2 total], [Nominal base],
  [`c2pnb176w1`], [51.872229], [71.336033], [48,921,052],
  [`c2pnb272w1`], [51.872229], [71.773568], [48,921,052],
)

The exact smooth numerator is `1374288545819365998587940548442445318355469588992205247403709040067002` over `4^128` monic polynomials. The model assumes half the base irreducibles yield usable divisor places, monic-polynomial smoothness predicts divisor smoothness, and the registry's operation prices and sparse linear algebra apply. Actual usable places, ordinary-query relation yield, useful rank, and target descent are unmeasured fields in the receipt.

Source: `research/weak_curve_leaderboard_20261012/ghs_genus128_cost_audit.json` (source SHA-256 fields); `docs/curves/jacobian-transfers.json`; `src/cryptanalysis/jacobian_transfer/cost.rs`. Producer: `research/weak_curve_leaderboard_20261012/src/bin/ghs_genus128_cost_audit.rs`. Reproduction: `cargo run --release --locked --bin ghs_genus128_cost_audit -- --check`. The previous Python model receipt is archived in `preliminary_python/`; the current receipt is native Rust.

The next gate is a smooth magic-eight cover with certified selected-subgroup image, actual divisor places, measured ordinary-query relations and rank, final linear algebra, and target descent. The conditional model does not add a curve-level finding or alter the isogeny graph.
