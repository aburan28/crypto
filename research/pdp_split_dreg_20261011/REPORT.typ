#set page(paper: "a4", margin: 16mm)
#set text(size: 9.5pt)
#show raw: set text(size: 8.5pt)

= Split and hybrid point decomposition against matched rho

*Study* `research/pdp_split_dreg_20261011` · branch `claude/pdp-split-dreg` · 2026-10-11/12 · *stage diagnostic* (one PDP solve plus a relation-phase projection; no end-to-end IC, no speedup or `vs_rho` claim). Protocol: `PROTOCOL.md` (preregistered), `AMENDMENT_1.md`. Analysis: `ANALYSIS.json` (native `pdp_split_analyze`). Full text and tables: `README.md`.

== Question

Does cutting one decomposition `R = P₁+…+P_m` into summation-polynomial pieces joined by auxiliary abscissae lower the solving degree and the cost, does guessing `g` Boolean unknowns, and how does each compare with rho on the same curve?

== Method (measured)

Koblitz curves `icv1-f2m11-t67-05f5aa36`, `icv1-f2m13-t181-515ee569`, `icv1-f2m17-tm101-00378d4e`, `icv1-f2m19-tm797-9c54981b`; factor base `span(1, z, …, z^(ℓ−1))`. Plans: unsplit `S₄` (m = 3, 3ℓ unknowns), S₃ chain 2+1 (3ℓ + n), two-way 2+2 with an `S₄` link (m = 4, 4ℓ + n), three-way 2+1+1 (4ℓ + 2n). Hybrid `g ∈ {0,2,4}`, every branch solved. F4 over `F₂[v]/(v²−v)`; solving degree = highest F4 step degree at which the basis grew; cost = F4 word XORs. 8 cells × 4 targets → 192 runs; every verified decomposition set checked against exhaustive enumeration. Matched rho: ecbench session `ECBS1he23809badf5b` (8 public targets × 5 rounds per curve, audit OK with `--replay-all`, receipt `0b56a9b7…cc0c`).

== Correctness

First grid: one cell failed the exhaustive check — a solution read-out defect in the library extractor on a non-inter-reduced F4 basis, not F4 (Amendment 1). Fixed in `pdp_split::solutions_from_basis` with a regression test; amended grid: *192/192 runs complete and exact*. Across both grids every run's solving degree and word-XOR count is identical; only the three amended runs' verified counts differ. cairn replay objectives for one cell and the rho audit were accepted (same host).

#figure(image("reports/figures/fig1-solving-degree.svg", width: 100%), caption: [F4 solving degree per target (4 per cell) by plan and guess count `g`; amended grid `runs-v2/`. Recorded values from `ANALYSIS.json`.])

#figure(image("reports/figures/fig2-word-xors.svg", width: 100%), caption: [F4 word XORs per target, all `2^g` branches summed (log). The S₃ chain is dearer than unsplit at ℓ = 4 and cheaper at ℓ = 5 for `g = 0`; unsplit with `g = 4` is cheapest in every `m = 3` cell.])

#figure(image("reports/figures/fig3-projection-vs-rho.svg", width: 100%), caption: [Relation-phase *projection* in S, `(F+1)/λ · X̄ · ns_per_word_xor / √r` (pinned ratio from another host; linear algebra and descent omitted, so a lower bound on IC), against *measured* strong and signed-Frobenius rho from session `ECBS1he23809badf5b`. Two cells have `λ = 0` in 256 exhaustively counted targets.])

== Predictions

#table(columns: (auto, 1fr), align: left,
  [P1], [Chain degree < unsplit degree at m = 3: *held* (4 vs 7, all six cells).],
  [P2], [Chain dearer at g = 0, ratio growing in n at ℓ = 4: *held at ℓ = 4* (3.9×, 8.0×, 28×, 69×); *falsified at ℓ = 5* (chain 20.8× and 6.6× cheaper).],
  [P3], [Three-way degree < two-way at m = 4: *held* (4 vs 7); cost mixed (113× cheaper at n = 11, 1.3× dearer at n = 13).],
  [P4], [Guessing: unsplit cheaper at g ≥ 2 *held*; chain dearer at g = 4 *held only at ℓ = 5*.],
  [P5], [Exhaustive check: *failed* in the first grid (read-out defect), *held* after Amendment 1.],
  [P6], [No projected S below strong rho: *held*; closest 22.8× strong rho, 391× signed-Frobenius rho.],
)

== Interpretation (labelled)

Splitting lowered the solving degree every time and paid in `n` unknowns per cut; whether that was net cheaper depended on ℓ (a crossover between ℓ = 4 and ℓ = 5 at g = 0). Guessing — which removes unknowns instead of adding them — gave the cheapest variant in every m = 3 cell. No variant came within an order of magnitude of rho on these curves. Limits: `r ≤ 2^18`, four targets per cell, pinned pricing, descriptive wall time, no unsplit m = 4 (no binary `S₅`), standard rather than Frobenius-invariant factor base.
