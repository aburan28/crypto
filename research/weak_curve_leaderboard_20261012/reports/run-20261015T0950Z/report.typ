#set page(paper: "a4", margin: 18mm)
#set text(font: "Helvetica Neue", size: 9pt)

= WeakCurves cover-source alignment

Run `run-20261015T0950Z` · 15 October 2026

The clean Linux catalog job passed its snapshot, native GHS audit and attack-map checks, then found a source-digest mismatch in the decomposition index. The `main` covers inventory changed from SHA-256 `93a0267f762e84037328636dd4d60f1b378462818500e711c7fb81e9e8c7452c` to `7d2f8c5c90b0b84aff79959223443285c5f6ce8e52cc5058be2e4699fdb6e884`. The current-main file was copied byte-for-byte and the index and site copy were regenerated.

#figure(image("dependency.svg", width: 100%), caption: [Source files and the derived export checked by the decomposition validator.])

The decomposition index still has 922 exact identities, 917 checked cover certificates, five formula replays and thirteen property filters. Only its `covers_sha256` field changed. The native genus-128 GHS cost sweep and its conditional prices remain identical to the previous run. Actual usable places, ordinary-query relation yield, useful rank and target descent remain unmeasured.

Evidence: `research/weak_curve_leaderboard_20261012/validation/ghs128_ci_decomposition_source_attempt2.log`, `docs/curves/covers.json`, `research/weak_curve_leaderboard_20261012/decomposition_properties.json`. The current-head Linux replay remains the final gate. This source-lock refresh adds no curve-level finding or isogeny edge.
