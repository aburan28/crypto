# Crossbred X3, 2026-09-18

X3 of `RESEARCH_ECC2K130_ROUTE_TARGETS.md`: Crossbred on the
symmetrised `u`-frame. Protocol was frozen before the run (same file,
earlier revision). This freeze is the measurement, not a slope.

**Boundary, stated before the run.** `Q_enum = C(|F_u|, 2)` word-ops at
one word-op per pair (undercounts enumeration). Axis: `ℓ = dim V`.
Falsifier: slope of `log₂(Q / Q_enum)` vs `ℓ` ≥ −0.1 (`α ≥ 1.9`), or no
determining space on any `m = 3` rung with `|F_u| ≥ 3`. Success: slope
≤ −0.5 (`α ≤ 1.5`) with `agree = yes` on every usable rung. A fit
requires ≥4 agreeing `m = 3` rungs with `|F_u| ≥ 3`.

**Result.** No fit. Only two usable agreeing rungs have `|F_u| ≥ 3`,
both on `K_1`: `n = 7` (`ℓ = 4`, `Q/C = 4.653`) and `n = 15` (`ℓ = 5`,
`Q/C = 6.409`). `K_0` is degenerate (`|F_u| = 1`) on every small rung
under the frozen paired-oracle divisor. Larger rungs that do have a
real `F_u` (`n = 17`, `23`, `31`) have no determining space under T4.
The additional falsifier `kernel_dim = 0` throughout is **not** met —
a space exists at `n = 7` and `n = 15` — but there are not four rungs
to fit. `Q / Q_enum > 1` on both usable cells. Every priced cell has
`filters = 0`.

A two-point sketch of `log₂(Q/C)` vs `ℓ` between those rungs has slope
`+0.46` (`α ≈ 2.46`). That is not a least-squares fit and is not a
result.

Host: `ip-172-31-19-103`, rustc 1.98.1.
