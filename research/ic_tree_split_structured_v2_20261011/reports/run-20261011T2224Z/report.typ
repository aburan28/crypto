#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.5pt)
#show raw: set text(size: 8.5pt)
= Tree-split cut-out screen, round 2 — run report

*Study* `research/ic_tree_split_structured_v2_20261011` (aburan28/crypto). *Registration* commit `8a7109c4` (pushed before the run); cairn registration objective `sha256:49b79374…`. *Run* `runs/registered`, 2026-10-11 22:24:56Z – 23:01:18Z UTC, commit `6506a1e2`, harness binary sha256 `28fdd513…cfdb9`. Host: Claude Code cloud container, 4 vCPU, 15 GB, shared with a concurrent census run; untimed study (the metric is a degree; `secs` logged, not measured). Self-test passed.

== What ran
`run.sh` → `examples/tree_split_structured_v2.rs`, 102 cells requested, *102 completed*, 0 failed; 2 `Z` cells censored (|Z| = 111,628 and 132,870 at n = 15, ℓ = 5: degree 5 needs 174,437 monomials > the 90,000 cap; all arms censored, excluded as registered). Arms: `frob2` (Frobenius-stable kernels in no proper subfield), `nbw` (normal-basis windows), `poly` at n = 23, 25, the object `Z = Sym²(T2)` with controls R and ZR, and the `subfield` gate.

== Recorded values
#figure(image("ds_vs_ell.svg", width: 100%), caption: [Measured `D_s` against ℓ per series; solid = curve set, dotted = the smaller control. 2,000 samples per arm per degree. Source `runs/registered/cells.jsonl`; plotter `reports/plot_cells.py`.])

#table(columns: 5, inset: 4pt,
  [*group*], [*cells (uncensored)*], [*gaps*], [*slopes of D_s(curve) on ℓ*], [*registered verdict*],
  [G0 subfield gate], [2], [T2 ⊆ GF(2^7), deg-1 deficit ≥ 14], [—], [*PASS*],
  [G1 T2 / frob2, n = 21], [16 (16)], [all 0], [0.6 – 1.5], [*CLOSED*],
  [G2 T2 / nbw, n = 17–23], [40 (40)], [39 × 0, 1 × −1], [0.8 – 1.1], [*CLOSED*],
  [G3 T2 / poly, n = 23, 25], [20 (20)], [all 0], [0.6 – 0.8], [*CLOSED*],
  [G4 Z / poly, n = 13, 15], [8 (7)], [all 0], [no slope clause], [*CLOSED*],
  [G5 Z / nbw, n = 13, 15], [8 (7)], [all 0], [no slope clause], [*CLOSED*],
  [G6 Z / frob2, n = 15], [4 (4)], [all 0], [no slope clause], [*CLOSED*],
  [frob2, n = 35 (supplementary)], [4 (4)], [all 0], [—], [not in the verdict],
)
The one −1 (nbw, n = 21, a = 1, ℓ = 4) is a degenerate cell, |F| = |T2| = 4, where the random 4-point control drew one of its 2,000 samples from its 8-point affine closure (probability ≈ 0.4%): control noise.

== Checks
Lean `Verdict.lean` (kernel `decide`): the registered rule (`Rule.lean`, sha256 `ec513363…b1e3`) on the 102 recorded cells gives gate = true and CLOSED for G1–G6; mutating a verdict fails. cairn result objective `sha256:84dde9a6…` accepted by the kernel; `cairn audit` (epochs 5 s): chain intact. CI workflow `tree-split-structured-v2-lean.yml` regenerates `Data.lean` and the readout from the records.

Raw records: `runs/registered/{cells.jsonl, cells.log, readout.txt, selftest.txt, commit, sha256}`.
