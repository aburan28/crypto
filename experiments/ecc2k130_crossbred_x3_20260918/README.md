# Crossbred X3 protocol, 2026-09-18

X3 of `RESEARCH_ECC2K130_ROUTE_TARGETS.md`: Crossbred on the
symmetrised `u`-frame. This file freezes the protocol **before** the
run. Numbers go in `summary.json` after the ladder finishes.

**Boundary, stated before the run.** `Q_enum = C(|F_u|, 2)` word-ops at
one word-op per pair (undercounts enumeration). Axis: `ℓ = dim V`.
Falsifier: slope of `log₂(Q / Q_enum)` vs `ℓ` ≥ −0.1 (`α ≥ 1.9`), or no
determining space on any `m = 3` rung with `|F_u| ≥ 3`. Success: slope
≤ −0.5 (`α ≤ 1.5`) with `agree = yes` on every usable rung. A fit
requires ≥4 agreeing `m = 3` rungs with `|F_u| ≥ 3`.

**Selection.** T4: smallest `D` then smallest `k` with
`kernel_dim ≥ v − k`, `v − k ∈ 1..=64`, `solve_crossbred` on the probe
not exhausted, `k < min(v, 22)`.

**Instance.** `K_a` with `a` given (`--a`, default 0). Divisor
`divisor_for_dimension(n, (n+1).div_ceil(m))`, then
`build_symmetrised_factor_base` / `build_symmetrised_system`. This is
the paired-oracle convention, not chosen per rung.

```bash
cargo run --release --example crossbred_bench -- --sym --no-sweep \
  5:3 7:3 9:3 11:3 13:3 15:3 17:3 19:3 23:3
```
