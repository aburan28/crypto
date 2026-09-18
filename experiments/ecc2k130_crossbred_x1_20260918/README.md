# Crossbred X1 attempt, 2026-09-18

T4 selection is applied (`examples/crossbred_bench.rs`). X1 asks for
`Q / C(|F|, 2)` fitted over ≥4 agreeing `m = 3` rungs. This freeze is
that attempt, not a slope.

**Boundary, stated before the run.** `Q_enum = C(|F|, 2)` word-ops at
one word-op per pair (undercounts enumeration). Falsifier: slope of
`log₂(Q / Q_enum)` vs `ℓ` ≥ −0.1 (`α ≥ 1.9`). Success: slope ≤ −0.5
(`α ≤ 1.5`) with `agree = yes` on every rung.

**Result.** No fit. Only two usable agreeing rungs have `|F| ≥ 3`.
`n = 13` and `n = 23` have no determining space under T4. Degenerate
bases (`|F| = 1`) at `n = 7` and `n = 15` are excluded from the ratio.
`Q / Q_enum` is 83.2 at `ℓ = 4` and 319.2 at `ℓ = 6` — Crossbred is
costlier than the conservative pair count, and the two-point sketch of
`log₂(Q/Q_enum)` vs `ℓ` rises (not a slope, not a falsification).

Host: `ip-172-31-19-103`, rustc 1.98.1, source after T4
(`9664f9d0` plus skip diagnostics in this commit).
