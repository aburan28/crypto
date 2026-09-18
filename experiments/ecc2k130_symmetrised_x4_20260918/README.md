# X4 symmetrised end-to-end, 2026-09-18

X4 of `RESEARCH_ECC2K130_ROUTE_TARGETS.md`: unchained symmetrised
`S_{m+1}` through collection, every phase priced. Protocol was frozen
before the run (same file, earlier revision). This freeze is the
measurement, not a slope.

**Boundary, stated before the run.** Floor `Λ = m/n` for a
Frobenius-stable base. Unit `Λ = trials · C(|F_u|, m−1) / 2^n`. Ratio
`Λ · n / m`. Algebraic reductions are a stage diagnostic until T3.
Falsifier for advance: slope of `log₂(Λ · n / m)` vs `n` consistent with
zero is **engineering**; a decrease is **advance**. A fit requires ≥4
usable rungs in one frame. Mixing K0 with K1 is inadmissible.

**T7.** `DecompositionStrategy::Symmetrised` agrees with `Enumerate` on
every completed sampled target of every priced rung (12/12 at
`n = 7, 9, 15, 17`). Recovered logs satisfy `[k]P = Q` on the usable
rungs.

**Result.** No fit. Three usable K1 rungs recovered a verified log:

| n | ℓ | \|F_u\| | trials | Λ | Λ·n/m | Enumerate |
|---|---:|---:|---:|---:|---:|---|
| 7 | 4 | 29 | 4 | 12.6875 | 29.60 | verified, 4 trials |
| 15 | 5 | 61 | 18 | 1.0052 | 5.03 | verified, 18 trials |
| 17 | 9 | 409 | 11 | 7.0022 | 39.68 | verified, 14 trials |

`n = 9` agreed 12/12 (all refuted) and both oracles collected zero
relations in 20,000 trials: not a usable rung. `n = 5, 11, 19` skip
(`prepare_symmetrised_attack` is `None`). `n = 13` has no curve.
`n = 23` was started and killed after ten minutes with no cell.

The three ratios are not flat and do not decrease through the last
rung. That is not a least-squares fit and is not a classification.
The predicted class remains engineering; it is not confirmed.

Host: `ip-172-31-19-103`, rustc 1.98.1.
Runner: `examples/symmetrised_e2e.rs`.
Receipt: `summary.json`.
