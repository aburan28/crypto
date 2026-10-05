# Round 1: does the semi-regular law predict refutation degrees it has not seen?

Registered before any registered draw is measured. §2 lists what ran before registration.
Goal: [GOAL.md](GOAL.md).

## 1. The law under test

For an exported system with `v` Boolean unknowns and equation degrees `d_1 … d_k`, let
`d_reg` be the index of the first non-positive coefficient of
`(1 + z)^v / ∏ (1 + z^{d_i})`. [predict.py](predict.py) computes it per draw from the
exported file. **The law: the refutation degree lies in `[d_reg, d_reg + 1]`.** Every
draw's prediction is in [runs/registered/predictions.jsonl](runs/registered/predictions.jsonl),
written before measurement.

The rival is the fixed-`ℓ` fit of the earlier rounds: `x4` refutes at `ℓ + 5`, and `rr`
at its `n = 17, 19` value (6 at `ℓ = 4`, 8 at `ℓ = 5` and `ℓ = 6`), whatever `n` is.

## 2. What ran before registration (disclosed; not evidence)

- **The twelve measured cells** of [ic_gb_ladder_20261003](../ic_gb_ladder_20261003/RESULTS.md)
  and [ic_dense_ladder_20261004](../ic_dense_ladder_20261004/RESULTS.md) were compared with
  the series. Every reading is at `d_reg` (both arms at `ℓ = 2`; `rr` at `ℓ = 3, 4, 6`) or
  `d_reg + 1` (`x4` at `ℓ = 3…6`; `rr` at `ℓ = 5`). The law was found on these cells, so
  they are not evidence for it.
- **A smoke test of the runner on a cell later dropped from the panel.** Two draws of
  `K₁/2⁶¹ ℓ = 4` were measured: `rr` refuted at 5 and `x4` at 8, both inside their windows
  (`[5, 6]`, `[8, 9]`) and `x4` below the fit's 9. Because they were seen, that cell is
  excluded; `K₀/2⁶¹ ℓ = 4` replaces it. The two readings are reported, never counted.
- **An M4RI driver** (libm4ri PLUQ) gave the right verdict on nine calibration cases but
  was slower than `examples/macaulay_dense.rs` on these matrices. It is not used.

## 3. Panel

Exported by `examples/rr_degree_ladder.rs --dump-dir`, seed **20261005** (new, so no draw
was seen before), four rootless draws per cell. Hashes:
[runs/registered/dump.sha256](runs/registered/dump.sha256); the `.sing` files are not
committed and regenerate byte-identically.

| cell | arms | `x4` window | `rr` window | role |
|:--|:--|:--|:--|:--|
| `K₀/2⁶¹ ℓ = 6`, `K₁/2⁶¹ ℓ = 6` | `x4`, `rr` | `[9, 10]` | `[6, 7]` | **discriminating**: the fit says 11 and 8 |
| `K₁/2²⁹ ℓ = 5`, `K₁/2⁴⁵ ℓ = 5`, `K₁/2⁶¹ ℓ = 5` | `x4`, `rr` | `[9, 10]` | `[6, 7]` | **discriminating for `rr`** (the fit says 8); `x4` consistency |
| `K₁/2²³ ℓ = 5` | `x4`, `rr` | `[9, 10]` | `[7, 8]` (one trivial draw: `[0, 1]`) | boundary: the law moves `rr` between `n = 23` and 25 |
| `K₁/2²³ ℓ = 6` | `rr` | | `[7, 8]` | consistency |
| `K₁/2⁵⁹ ℓ = 6` | `x4`, `rr` | `[10, 11]` | `[7, 8]` | boundary: the law moves `x4` between 59 and 61 |
| `K₀/2⁶¹ ℓ = 4` | `x4`, `rr` | `[8, 9]` | `[5, 6]` | consistency |

`x4` on `K₁/2²³ ℓ = 6` is not run: both laws give the same window there and it costs
about 10 hours.

## 4. Measurement

`examples/macaulay_dense.rs`, unchanged since `ic_dense_ladder_20261004` apart from its
post-run clippy edit. [lanes.sh](lanes.sh) scans `D = 1, 2, …` per draw up to `d_reg + 2`,
so a reading one above the window is measured, not censored, and stops at the first
refuting `D`. Per process: 86,400 CPU-s and 13.5 GB, the engine refusing any basis above
12 GB before allocating. A killed or refused process censors the draw at `≥ D`; censoring
is never evidence, and there are no retries. Units run in parallel lanes, on this host and
on other 4-core hosts (each regenerates its systems and checks them against
`dump.sha256`); degrees are deterministic and wall times are advisory.

## 5. Predictions and decision rule

- **P1.** Every resolved draw reads inside its window.
- **P2 (the discriminating test).** On the five discriminating `x4`/`rr` cells, the median
  is inside the window and therefore below the fit: `x4` at 9 or 10 on `n = 61, ℓ = 6`;
  `rr` at 6 or 7 on `n = 29, 45, 61` at `ℓ = 5` and on `n = 61, ℓ = 6`.

Verdict, fixed now:

- **Generic law holds out of sample** if P2 holds on every discriminating cell with at
  least three resolved draws, and at least 90% of all resolved draws satisfy P1.
- **Fit holds, generic law refuted** if the median of a majority of the discriminating
  cells equals the fit's value.
- **Neither** otherwise; the readings are reported against both.

A reading below `d_reg` would be the more surprising failure (structure that helps the
attacker) and is reported first if it occurs.

## 6. Scope

Two Koblitz curves, `n` from 23 to 61, `ℓ` 4–6, `m = 3`, four draws per cell, one engine.
A stage diagnostic: it tests a predictor of one stage's degree. It does not claim an
end-to-end cost, it says nothing about `m ≥ 4`, and nothing at `n ≈ 83` or 131 beyond
[cost_n131.py](cost_n131.py), which is conditional on this round's verdict.
