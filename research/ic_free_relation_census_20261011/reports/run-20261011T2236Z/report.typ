#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show raw: set text(size: 8.5pt)
= Free-relation census of the K = 600 pair-root tables — run report

*Study* `research/ic_free_relation_census_20261011` (aburan28/crypto). *Prediction* `PREDICTION.md`, commit `e51993fc` (before the fixture runs). *Run* 2026-10-11 ≈22:36Z UTC, commit `e51993fc`, binary sha256 `2e183196…b203`. Host: Claude Code cloud container, 4 vCPU, 15 GB, shared with a concurrent run; untimed study (wall seconds logged, not measurements).

== What ran
`examples/koblitz_free_relation_census.rs <base_header.jsonl>` on the three committed fixtures: n = 61 (`experiments/koblitz-single-target-n61-20260925/factor-base.jsonl`), n = 73 (`…n73-20261002/base_n73_K600.jsonl`), n = 83 (`…n83-20261006/base_n83_K600.jsonl`). 3 runs requested, 3 completed. Method: all unordered pair-states `P_i ± π^k P_j` keyed by least Frobenius rotation of x; sort; each collision verified as `S₁ = η·π^t S₂` by point comparison, row `Σ c·e` over `Z/r`, each kept row re-verified as `Σ [c]P = O`; rank mod r. Calibration: six synthetic bases (`runs/calibration.jsonl`).

== Recorded values
#figure(image("census.svg", width: 92%), caption: [Per rung: the note's prediction `n⁴K⁴/r`, the class-count prediction `K⁴n³/r` (both derived), the measured nontrivial rows and rank mod r. Dashed: rank needed, K − 1 = 599. Records `runs/census_n{61,73,83}.json`.])

#table(columns: 8, inset: 4pt,
  [*rung*], [*log₂ r*], [*states*], [*note n⁴K⁴/r*], [*classes K⁴n³/r*], [*rows (meas.)*], [*rank (meas.)*], [*needed*],
  [n = 61], [47.2], [21,960,000], [11,016], [181], [147], [49], [599],
  [n = 73], [56.3], [26,280,000], [43], [0.6], [3], [1], [599],
  [n = 83], [52.9], [29,880,000], [718], [8.6], [18], [6], [599],
)
Identity failures 0 and row-check failures 0 on every rung; 1,800 (= 3K) structural collisions per rung give rows ≡ 0 and are dropped. Rank is exactly one third of the rows at every rung: a 4-term zero-sum pairs up three ways, each a separate table collision. Calibration (synthetic, unsaturated): n = 37 → 2,698 rows vs 2,847 (`K⁴n³/r`) vs 105,329 (note); n = 43 → 8,585 vs 8,667 vs 372,674.

== Reading
The note's falsifier (≥ 500 verified relations at n = 61 and 83, full rank at n = 61) fails at every rung; its count is high by ≈ n. The rank stage is not free; the ledger's precompute charge stands. Making it free needs K ≳ (3r/n³)^(1/3) ≈ 1,300 / 8,600 / 3,600.

Raw records: `runs/{census_n61,census_n73,census_n83}.json`, `runs/calibration.jsonl`, `runs/{commit,sha256}`.
