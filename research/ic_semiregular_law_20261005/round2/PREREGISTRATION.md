# Round 2: a sharper rule for `x4`, and round 1's structural rule for `rr`, on fresh cells

Registered before any registered draw is measured. Goal: [../GOAL.md](../GOAL.md).
Round 1: [../RESULTS.md](../RESULTS.md) (registered verdict "neither").

## 1. What round 1 left

- `x4` sat inside `[d_reg, d_reg + 1]` on all 32 draws, but the window does not say which
  end. `rr` refuted below its window on two cells, and a post-hoc diagnostic traced that to
  duplicated equations plus extra top-form relations.
- After round 1, before this registration, one more rule was compared with every recorded
  reading (all rounds; disclosed, not evidence): **`d_aff`**, the first `k` at which the
  running sum of `(1 + z)^v / ∏ (1 + z^{d_i})` is non-positive, the series of the
  homogenized affine system. On the 56 recorded `x4` draws the reading equals `d_aff` on 49
  and is one below on 7 (the three `ℓ = 2` draws of `K₁/2¹⁷` and `K₁/2⁴⁵ ℓ = 5`). It is
  never above.

## 2. Predictions, computed per draw before measurement

[predict2.py](predict2.py) writes [runs/registered/predictions.jsonl](runs/registered/predictions.jsonl).

- **A (`x4`).** The refutation degree equals `d_aff`. **A′:** it is at most `d_aff`.
- **B (`rr`).** Reduce the equations by graded elimination (`../topforms.py`). Let `d_red`
  be the semi-regular `d_reg` of the reduced profile. The refutation degree equals
  `d_red − 1` when `n ≥ 29`, and `d_red` when `n < 29`. The rule was fitted on round 1's
  nine `rr` cells; this round tests it on new ones.

## 3. Panel

Seed **20261012** (new), four rootless draws per cell; hashes in
[runs/registered/dump.sha256](runs/registered/dump.sha256).

| cell | `x4`: `d_aff` (fit `ℓ + 5`) | `rr`: predicted | role |
|:--|:--|:--|:--|
| `K₀/2⁵³ ℓ = 3`, `K₁/2⁵⁷ ℓ = 3` | 7 (8) | not run | A against the fit |
| `K₀/2⁵⁷ ℓ = 4`, `K₁/2⁵⁹ ℓ = 4` | 8 (9) | 5 | A against the fit; B |
| `K₀/2⁵³ ℓ = 5`, `K₁/2⁵⁷ ℓ = 5` | 9 (10) | 5 | A against the fit; B |
| `K₀/2⁵⁷ ℓ = 6` | 10 (11) | 6 | A against the fit; B |
| `K₀/2⁴¹ ℓ = 5`, `K₀/2⁴⁵ ℓ = 5` | 10 (10) | 6 | A at the boundary where round 1 missed once; B |
| `K₀/2³¹ ℓ = 4` | 9 (9) | 5 | A consistency; B |
| `K₀/2²³ ℓ = 5` | 10 (10) | 7 | A consistency; B with `n < 29` |
| `K₀/2³⁷ ℓ = 6` | not run | 7 | B (`x4` there needs about 3 h per draw and both rules give 11) |

`rr` at `ℓ = 3` is not run: every draw reduces to a system containing a constant.
For the same reason one `rr` draw of `K₁/2⁵⁹ ℓ = 4` is excluded. `predict2.py` writes 0
for such draws, but the constant only appears at the equations' own degree, so 0 is not a
real prediction.

Smoke test before registration: two excluded `rr` draws of `K₁/2⁵⁷ ℓ = 3`, used only to
check the runner. Nothing registered was measured.

## 4. Measurement

As in round 1: `examples/macaulay_dense.rs` through [lanes.sh](lanes.sh), scanning
`D = 1 …` up to predicted + 2 and stopping at the first refuting `D`. Per process: 86,400
CPU-s and 13.5 GB. Censoring is never evidence, and there are no retries.

## 5. Decision rule, fixed now

- **A holds** if every `x4` cell's median equals `d_aff` and at least 90% of resolved `x4`
  draws equal it. **A′ holds** if no resolved `x4` draw is above `d_aff`.
- **B holds** if every `rr` cell's median equals its prediction and at least 90% of
  resolved `rr` draws equal it.
- Each is reported on its own. A failure of A at the two boundary cells, with A′ intact,
  would be reported as "`d_aff` is an upper bound, not exact".

## 6. Scope

Two Koblitz curves, `n` 23–59, `ℓ` 3–6, `m = 3`, four draws per cell, one engine. A
stage diagnostic about one decomposition step's degree; no end-to-end cost.
