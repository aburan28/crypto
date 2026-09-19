# Toy-suite rho parity for Koblitz pair-table index calculus

**Question.** On the five small Koblitz cells the candidate tournament
already solved, does the production pair-table pipeline cost no more
than matched signed-Frobenius rho in exclusive group operations, on
every useful completed verified case?

**Answer.** No. Iteration 1 reported the gate met (mean `α`
0.76–0.91) and that figure was an accounting error: the `|F| − 1`
point additions a windowed probe makes before its lookups were not
charged (§10). Charged, the five toy cells cost 16–23× rho and the
twelve-rung ladder (§9–§10, `2^{11}` to `2^{39}`) costs 8× to 600×
rho on every rung, 108/108 pairs verified on both arms. Iterations
0–2 are retained as before marks; their gates are withdrawn. Class of
the correction: **accounting**. Iterations 4–5 (§11–§12) are the
engineering that remains after the correction and are predicted, in
advance, not to reach the gate on any rung. Not degree-131, not a
refit of X1–X4; the product-law floor at `n = 131` is untouched.

Companion: [`RESEARCH_ECC2K130_ROUTE_TARGETS.md`](RESEARCH_ECC2K130_ROUTE_TARGETS.md)
T12. Scoreboard panel `#ecc2k130-rho-parity-20260918`.

## 1. Boundaries, written before the run

The instance family is `K_a : y² + xy = x³ + a x² + 1` over `F_{2^n}`,
prime-order subgroup of order `r`, automorphism group `⟨−1⟩ × ⟨π⟩` of
order `A = 2n`.

**Floor.** An `m = 3` factor base of `B` points yields at most
`B(B+1)/2` unordered pair sums and therefore uniform-target 3-summand
yield at most `min(1, B³/(6(r−1)))` up to the usual factorial and
repetition conventions. The table-row rule does not enlarge `B`. The
counting bound does not move.

**Reference.** Pollard rho with the same automorphisms, the shipped
`koblitz_signed_frobenius_rho_with_progress`, on the same curve,
subgroup, public target and algorithm seed. Expected walk length
`√(π r / (4 n))`. Setup (jump table and parallel-walk starts) is
charged.

**Unit.** Exclusive group operations

```
G = point additions + binary-method additions of every charged scalar multiplication
S = G / √r
```

Doublings count as additions. A scalar `s > 1` costs
`(⌊log₂ s⌋) + (hw(s) − 1)`. Rho's recorded scalar multiplications are
converted at the expected cost `1.5 · (bitlen(r) − 1)`. Linear algebra
is a stage diagnostic (column count is `O(1)` on these cells) and is
not dropped from the correctness budget: every column logarithm is
re-checked in the group and that check is charged. Wall-clock is a
footnote.

**Ratio.** `α = G_IC / G_rho`. The gate is `α ≤ 1` on **every useful
completed verified cell**. An aggregate cannot hide a losing cell.
This is the same all-cases shape as
[`research/rho_parity_20260915/`](research/rho_parity_20260915/), in
group operations rather than Valgrind instructions.

Class: **engineering**. `S` may fall against rho; the ratio
to the counting floor does not, because `B` and `m` are unchanged.
Iteration 1 confirmed that class: `S` fell, the floor ratio did not.

## 2. Frozen suite

| Cell | Curve | Recipe floor `B` | Role |
|---|---|---:|---|
| n13a0 | `K_0 / F_{2^13}` | 78 | useful |
| n17a1 | `K_1 / F_{2^17}` | 102 | useful |
| n19a0 | `K_0 / F_{2^19}` | 114 | useful |
| n19a1 | `K_1 / F_{2^19}` | 114 | useful |
| n23a0 | `K_0 / F_{2^23}` | 138 | useful |

These are the tournament's five cells. The recipe is
`SubgroupOrbits { seed: 43, points: 6n }`. Sampling stops at the first
Frobenius closure that reaches the floor; the built `|F|` may overshoot
by a batch of orbits. `K` is the number of signed Frobenius orbits of
the materialized base.

**Targets.** Three public hash-to-curve seeds per cell, unknown scalar,
domain `ic-workflow-public-target-v1`. Frozen seeds
`0x51EED001`, `0x51EED002`, `0x51EED003`. Three paired repetitions.
Order on pair `(cell, target, rep)` is IC-then-rho when
`(n + a + target_index + rep)` is even, else rho-then-IC.

**Algorithm seed.** `0x5EED_0007`. Rho's per-target seed is
`algorithm_seed XOR i · 0x9e3779b97f4a7c15`, the production worker's
rule.

## 3. Frozen candidate

Production path, not the tournament's specialised `tiny` worker:

1. `DecompositionStrategy::PairTable`, `m = 3`, `descent_m` unset.
2. `allow_direct_relation = false`.
3. `collapse_negation = true`.
4. Linear algebra: dense (column count is the orbit count, `O(1)`).
5. `max_trials = 4096`, collection batches of 64.
6. Folded pair table with
   `t = PairSumTable::optimal_folded_rows(K, |F|, r)` rows.
   The formula is `t · |F| + (K+3)/(λ · cov(t))` with
   `λ = |F|(|F|+1)/(2r)` and `cov(t) = 1 − ((K−t)/K)²`,
   `1 ≤ t ≤ K`. It is a function of the public base, not of a target.
7. Walked collection: `collection_window = |F| − 1`, the largest
   window that still takes the walked-probe path.
8. Descent: `IndividualLogSolver` walking probes (production).
9. Every column logarithm re-checked as `[x_o]G = R_o`. Recovered
   `d` re-checked as `[d]G = Q` in general arithmetic.

Cold single-target: factor-base construction, table build, relation
collection, certification, descent and final check are all charged to
that target. There is no warm amortisation across targets.

## 4. What is charged

**IC.** Table-build additions. Iteration 0 recomputed every stored sum
on the occupancy pass and again on the fill pass (`additions = 2 ×
stored`). Iteration 1 holds `(key, orbit)` from the first `add_many`
and scatters from that list (`additions = stored`). Charging one pass
while still adding twice is relabelling and is not this iteration.
Walked collection: one binary-method multiplication per 64-probe run
plus one addition per probe, and one stride multiplication per batch.
Column certification: binary-method cost of each `x_o`. Descent: three
binary-method multiplications and 63 additions to place 64 walks, plus
one addition per subsequent walked probe, plus the final `[d]G`.
Factor-base cofactor projections: `K` multiplications by the cofactor.

**Rho.** `rho_exclusive_group_ops`: `setup_group_additions +
walk_group_additions + expected_binary_method_group_ops(bitlen(r)) ·
(setup_scalar_multiplications + candidate_verification_scalar_multiplications)`.

Shared curve construction is charged to neither arm (both need it; it
cancels in the ratio only if it is equal, and it is not a group
operation of the solvers). Target hashing is charged to neither: the
target is an input.

## 5. Gate, success, abandonment

**Success.** Every useful cell has every repetition of every target
verified `[d]G = Q` on both arms, and the mean `α` on that cell is
`≤ 1`, and no individual verified pair on that cell has `α > 1`.

**Not a success.** A cell that solves but loses. A timeout. A skipped
verification. Pooling cells. Quoting wall-clock. Changing `B`, `m`,
the row formula, the window rule, the seed, or the conversion after
seeing a cell.

**Class, measured.** Engineering: the counting floor is unchanged.
Hitting the gate is a practicality result on five toy cells. It does
not move `α` at `n = 131` and is not claimed to.

**Abandon this line** if the first frozen measurement has a cell with
mean `α > 2` after the row rule and walked probes, with no remaining
production lever that changes `G` without changing `B` or `m`. Then
the leftover is the specialised tournament worker's instruction-count
win, which already exists and is a different unit.

**Inadmissible.** Degree-131 extrapolation. Mixing K0 with K1 in a
slope. Dropping table build, certification, or descent from `G`.
Using Valgrind instructions as this note's unit (that unit already
has a strict win in
[`research/ic_candidate_tournament_20260915/`](research/ic_candidate_tournament_20260915/);
this note asks a different question). Changing T4, `Q_enum`, or the
X4 divisor.

## 6. How to run

```bash
cargo run --release --example rho_parity_e2e
cargo run --release --example rho_parity_e2e -- --quick   # not a result
```

Receipt: `experiments/ecc2k130_rho_parity_20260918/`.
Iteration 0: `experiments/ecc2k130_rho_parity_20260918/iteration-0/`.
Iteration 1: `experiments/ecc2k130_rho_parity_20260918/iteration-1/`.
Iterations 2–5 (`--ladder`, `--batch`, `--recipe`): §9–§12.
Runner: `examples/rho_parity_e2e.rs`. From iteration 3 the runner
charges the summand additions (§10); a rerun of the iteration 0/1
command now reports the corrected `G_IC`, with the old subtotal kept
as `g_ic_before_summand_correction`.

```bash
cargo run --release --example rho_parity_e2e -- \
  --out experiments/ecc2k130_rho_parity_20260918/iteration-1
```

## 7. Iteration 0, two-pass table (measured 2026-09-18)

Host `ip-172-31-19-103`. Cited from
[`experiments/ecc2k130_rho_parity_20260918/iteration-0/summary.json`](experiments/ecc2k130_rho_parity_20260918/iteration-0/summary.json).
All 45 pairs recovered `[d]G = Q` on both arms. Class: **engineering**.
`B` and `m` are the frozen recipe. The counting floor did not move.

| Cell | `|F|` | `t` | `G_IC` | table adds | `G_rho` | mean `α` | max `α` | gate |
|---|---:|---:|---:|---:|---|---:|---:|---|
| n13a0 | 182 | 1 | 685 | 364 | 548–555 | 1.241 | 1.250 | unmet |
| n17a1 | 272 | 1 | 978 | 544 | 869–951 | 1.063 | 1.125 | unmet |
| n19a0 | 304 | 1 | 1069 | 608 | 897–903 | 1.187 | 1.192 | unmet |
| n19a1 | 304 | 1 | 1093 | 608 | 988–1106 | 1.054 | 1.106 | unmet |
| n23a0 | 368 | 2 | 1931 | 1380 | 1474–1670 | 1.238 | 1.310 | unmet |

`all_cases_gate = false`. IC cost is setup-dominated and constant per
cell. One pair on n19a1 already has `α = 0.988`; the cell still fails
because other targets on that cell are above 1. Collection is one
64-probe batch on every cell; descent is one walking probe after
64-walk placement. Neither is the dominant term. Table adds are
`2 × stored` (n13: 182, n17: 272, n19: 304, n23: 690).

Mean `α` is below 2 on every cell, so the abandonment clause in §5
does not fire. The remaining production lever that changes `G` without
changing `B` or `m` is to stop paying for the stored sums twice.

## 8. Iteration 1, one-pass table (measured 2026-09-18)

**Hypothesis, frozen before the run.** Occupancy already computes every
stored sum. Hold `(key, orbit)` per row and scatter into buckets from
that list. Do not call `add_many` a second time. Charge
`additions = stored`. The stored keys, the row rule
`t = optimal_folded_rows(K, |F|, r)`, `B`, `m`, the window, the seeds,
and the conversion are unchanged. Predicted class: **engineering**.

Host `ip-172-31-19-103`. Cited from
[`experiments/ecc2k130_rho_parity_20260918/iteration-1/summary.json`](experiments/ecc2k130_rho_parity_20260918/iteration-1/summary.json).
All 45 pairs recovered `[d]G = Q` on both arms.
`all_cases_gate = true`. Table adds equal stored entries (n13: 182,
n17: 272, n19: 304, n23: 690), not twice that. Collection and descent
counts are the iteration-0 values, so the whole drop is the second
`add_many` going away. The predicted `G_IC` in the hypothesis table
matched the measurement on every cell.

| Cell | `|F|` | `t` | `G_IC` | table adds | `G_rho` | mean `α` | max `α` | gate |
|---|---:|---:|---:|---:|---|---:|---:|---|
| n13a0 | 182 | 1 | 503 | 182 | 548–555 | 0.911 | 0.918 | met |
| n17a1 | 272 | 1 | 706 | 272 | 869–951 | 0.767 | 0.812 | met |
| n19a0 | 304 | 1 | 765 | 304 | 897–903 | 0.850 | 0.853 | met |
| n19a1 | 304 | 1 | 789 | 304 | 988–1106 | 0.761 | 0.799 | met |
| n23a0 | 368 | 2 | 1241 | 690 | 1474–1670 | 0.796 | 0.842 | met |

Worst pair is n13a0 target 1, `α = 0.918`. No verified pair has
`α > 1`. The §5 gate holds on every useful cell.

Class: **engineering**. `S` fell; the ratio to the counting floor did
not, because `B` and `m` did not. Iteration 0 remains the before mark
on the scoreboard. This is not n=131 parity and is not claimed to
move `α` at degree 131.

## 9. Iteration 2, the ladder (frozen before the run)

**Question.** The five cells span `log₂ r` = 11.0 to 21.0 and IC is
setup-dominated on all of them. Where does parity break, and what are
the fitted exponents of `G_IC` and `G_rho` in `r`? `AGENTS.md` §5 asks
for exponents over at least four sizes; the gate alone does not supply
them.

**Cells added.** Every Koblitz curve `KoblitzCurve::new(a, n)` returns
for `25 ≤ n ≤ 47`, with its prime subgroup order. The ladder is
indexed by `log₂ r`, not `n`, because the cofactor varies:

| Cell | `log₂ r` | recipe `|F|` | `K` | `t` | expected rho walk |
|---|---:|---:|---:|---:|---:|
| n29a1 | 15.37 | 464 | 8 | 1 | 34 |
| n31a0 | 20.46 | 496 | 8 | 1 | 191 |
| n39a0 | 26.03 | 624 | 8 | 5 | 1,176 |
| n37a0 | 27.78 | 592 | 8 | 7 | 2,212 |
| n43a1 | 32.11 | 688 | 8 | 8 | 9,210 |
| n47a1 | 36.64 | 752 | 8 | 8 | 42,242 |
| n41a0 | 39.00 | 656 | 8 | 8 | 102,621 |

`|F|`, `K`, `t` are from the public recipe (`SubgroupOrbits{43, 6n}`,
`optimal_folded_rows`) and are written here before any target is
solved. The five §2 cells are re-run in the same receipt so every
rung has one code hash.

**Unchanged.** Unit, gate shape, recipe, row rule, window rule,
`m = 3`, seeds, target domain, repetitions, order rule, conversion,
one-pass table.

**Changed, and why.** `max_trials` `4096 → 2²⁰` and rho
`max_iterations_per_restart` `2²⁰ → 2²⁴`. Both are caps, not costs:
every trial and iteration actually spent is charged. A cell that
exhausts either cap is a failure, retained as a failure.

**Diagnostics added, not charged.** Table lookups
(`collection_trials · window` and `descent_trials · |F|`). The unit
charges a walked probe one addition; the `|F|` memory probes it makes
are not group operations. Their count is written next to `G` so the
unit's blind spot is visible.

**Fit.** Least-squares slope of `log₂ mean(G)` against `log₂ r` over
every cell whose nine pairs all verified, for each arm. Rho's slope
should sit near `0.5` once the walk dominates its fixed setup.

**Prediction, written before the run.** `G_IC` is
`table + K·cofactor muls + collection + certification + descent`.
Collection probes scale as `2rK / |F|³` at fixed `|F| ≈ 6n`; so as
`r` grows with `|F|` nearly flat the IC slope tends to `1`, twice
rho's. On these cells `|F|³` is between `10⁸` and `4·10⁸`, so the
collection term stays below rho's walk through `log₂ r = 39`, and the
gate is predicted to **hold on every rung** while the fitted IC slope
is predicted to exceed rho's. The extrapolated crossing in this unit
is where `2rK/|F|³ ≈ √(πr/4n)`, about `r ≈ 2⁴⁸` at `|F| ≈ 650`,
past `MAX_N = 63`'s reach in one run. That crossing is an
extrapolation and will be labelled one.

**Class.** Engineering if the gate holds: the counting floor at these
`B` is unchanged. If the IC slope is at or below `0.5` on ≥ 4 rungs,
that contradicts the counting argument and must be re-examined as an
accounting error before it is called anything else.

**Inadmissible.** Choosing which of the twelve cells to fit after
seeing them. Growing `|F|` on a rung to rescue it. Quoting the
crossing as a measurement.

```bash
cargo run --release --example rho_parity_e2e -- --ladder \
  --out experiments/ecc2k130_rho_parity_20260918/iteration-2
```

**Measured 2026-09-19**, host `ip-172-31-19-103`, cited from
[`iteration-2/summary.json`](experiments/ecc2k130_rho_parity_20260918/iteration-2/summary.json).
108/108 pairs verified on both arms. Under the iteration 0–2
accounting the gate **fails on three of twelve rungs**: n29a1 (mean
`α` 1.245), n47a1 (2.18, max 4.44), n41a0 (1.89, max 2.57). Slopes
over twelve rungs: IC 0.273, rho 0.265. The prediction "gate holds on
every rung" is falsified, and the collection model it rested on
(`2rK/|F|³` probes) was wrong by the factor that §10 identifies. The
iteration-2 numbers are superseded by §10 and are kept as the before
mark; they are not the thread's result.

## 10. Iteration 3, the accounting correction (measured 2026-09-19)

**What was wrong.** A windowed `m = 3` probe `R` forms `R − P_k` for
every `P_k` in its window before it looks anything up
(`PairSumTable::witnesses_fast_window`, one `add_many` over the
window). Those are point additions: the same batched affine addition
the table build is charged one operation per output for. §4 charged
the probe one addition and the window nothing. The descent probe does
the same over the whole base (`decompose_fast(state, 3)`, blockwise,
stopping at its first witness). So iterations 0–2 dropped
`trials × window` additions from collection and about
`trials × |F|` from descent. `AGENTS.md` §6: a ratio improved by
dropping a cost from the budget does not count. The tournament's
Valgrind unit never had this hole, which is one reason its win and
this unit's must not be mixed.

**Correction.** `G_IC` now includes `collection_trials × window` and
`descent_trials × |F|`. The descent term overstates by at most one
scan (the successful probe stops early); on every rung that is under
`|F|` operations against a `G_IC` in the tens of thousands. The
uncorrected subtotal is kept in every receipt as
`g_ic_before_summand_correction`. Nothing else changed: same seeds,
recipe, rows, window, caps, code hash for the arithmetic.

**Class: accounting.** Numbers changed, the algorithm did not. No gain
is claimed for iterations 0–2; their gates are withdrawn. Because the
counts needed for the correction (`collection_trials`, `window`,
`descent_trials`, `factor_base_size`) are in the frozen receipts of
iterations 0, 1 and 2, the corrected `α` for those runs is derived
from them and is identical to the iteration-3 rerun, which is
deterministic on the same seeds.

Cited from
[`iteration-3/summary.json`](experiments/ecc2k130_rho_parity_20260918/iteration-3/summary.json).
108/108 pairs verified on both arms. `G_IC` is constant per cell;
`G_rho` is the mean over nine pairs.

| Cell | `log₂ r` | `|F|` | `t` | table | collection summands | descent summands | `G_IC` | `S_IC` | `G_rho` | `S_rho` | mean `α` | min–max `α` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| n13a0 | 10.97 | 182 | 1 | 182 | 11,584 | 182 | 12,269 | 274.1 | 552 | 12.33 | 22.2 | 22.1–22.4 |
| n29a1 | 15.37 | 464 | 1 | 464 | 29,632 | 464 | 31,132 | 151.1 | 833 | 4.04 | 37.4 | 36.5–37.9 |
| n17a1 | 16.00 | 272 | 1 | 272 | 17,344 | 272 | 18,322 | 71.5 | 922 | 3.60 | 19.9 | 19.3–21.1 |
| n19a0 | 17.00 | 304 | 1 | 304 | 19,392 | 304 | 20,461 | 56.6 | 900 | 2.49 | 22.7 | 22.7–22.8 |
| n19a1 | 18.00 | 304 | 1 | 304 | 19,392 | 304 | 20,485 | 40.0 | 1,039 | 2.03 | 19.8 | 18.5–20.7 |
| n31a0 | 20.46 | 496 | 1 | 496 | 31,680 | 496 | 33,310 | 27.8 | 1,230 | 1.03 | 27.1 | 25.6–28.6 |
| n23a0 | 21.00 | 368 | 2 | 690 | 23,488 | 368 | 25,097 | 17.3 | 1,563 | 1.08 | 16.1 | 15.0–17.0 |
| n39a0 | 26.03 | 624 | 5 | 2,340 | 39,872 | 1,456 | 44,453 | 5.37 | 3,992 | 0.48 | 11.3 | 9.8–12.7 |
| n37a0 | 27.78 | 592 | 7 | 2,590 | 37,824 | 5,723 | 46,921 | 3.09 | 5,681 | 0.37 | 8.4 | 6.7–9.6 |
| n43a1 | 32.11 | 688 | 8 | 3,096 | 439,680 | 69,259 | 514,405 | 7.55 | 16,700 | 0.25 | 31.2 | 25.1–35.6 |
| n47a1 | 36.64 | 752 | 8 | 3,384 | 21,052,032 | 825,195 | 21,956,488 | 67.2 | 58,516 | 0.18 | 611.7 | 225–1,261 |
| n41a0 | 39.00 | 656 | 8 | 2,952 | 21,253,440 | 6,886,469 | 28,241,615 | 38.1 | 57,538 | 0.08 | 532.1 | 301–759 |

Gate: **fails on every rung.** Best rung n37a0 at `α = 8.4`. Fitted
slopes over the twelve rungs: `G_IC ∝ r^{0.40}`, `G_rho ∝ r^{0.27}`
(rho's fixed setup still dominates its walk below `log₂ r ≈ 26`; on
the top four rungs alone the slopes are IC `0.88`, near the model's
`1` at fixed `|F|`, and rho `0.32`, low because n41a0's three targets
finished in 0.56× the expected walk; three targets per rung is a noisy
slope and it is reported as such). Corrected iteration
0 and 1 on the five toy cells: mean `α` 16.5–23.1 and 16.1–22.7. The
one-pass table was worth 1–3% of the corrected `G_IC`, not the 20–40%
the uncorrected figure showed.

**Corrected answer to the §1 question: no.** In exclusive group
operations the production pair-table pipeline costs 8× to 600× matched
rho on every Koblitz cell from `2^{11}` to `2^{39}`, and the
iteration-1 "parity" was the window's additions going unpriced.

## 11. Iteration 4, collection batch of 8 (frozen before the run)

**Why.** With the window charged, a collection batch is
`64 × (|F| − 1)` additions whatever the system needs. On the rungs
below `2^{27}` every probe decomposes and the solver needs about
`K + 2` relations; the first batch of 64 buys 50-odd relations it
never uses. Trying the system after every 8 probes stops there
instead. This is a runner granularity, not a change to `B`, `m`, the
window, the row rule, the seeds or the accounting; each stride
multiplication per batch is still charged.

**Prediction.** Cells whose 64-probe batch yielded `≥ 16` relations
(n13a0 through n39a0) fall to 16–24 trials, so their `G_IC` drops by
about 4× and `α` lands between 4 and 10. n37a0 (10 relations in 64)
and the three rungs above it are unchanged within one batch. Gate
still fails on every rung. Class: **engineering**.

```bash
cargo run --release --example rho_parity_e2e -- --ladder --batch 8 \
  --out experiments/ecc2k130_rho_parity_20260918/iteration-4
```

**Measured 2026-09-19**, host `ip-172-31-19-103`, cited from
[`iteration-4/summary.json`](experiments/ecc2k130_rho_parity_20260918/iteration-4/summary.json).
108/108 verified. As predicted: the eight rungs through n39a0 stop at
16–32 trials and land at mean `α` 4.8–7.8 (n29a1 15.2, its 24 trials
against a 464-point window); n37a0 stops at 56 trials (7.6); the top
three rungs move by less than one batch. Gate fails on every rung.
Class: **engineering**. Slopes over twelve rungs IC 0.48, rho 0.27.

| Cell | trials (was 64+) | `G_IC` | `G_rho` | mean `α` (iter. 3 → 4) |
|---|---:|---:|---:|---|
| n13a0 | 16 | 3,544 | 552 | 22.2 → 6.4 |
| n29a1 | 24 | 12,610 | 833 | 37.4 → 15.2 |
| n17a1 | 16 | 5,287 | 922 | 19.9 → 5.7 |
| n19a0 | 16 | 5,889 | 900 | 22.7 → 6.5 |
| n19a1 | 16 | 5,915 | 1,039 | 19.8 → 5.7 |
| n31a0 | 16 | 9,527 | 1,230 | 27.1 → 7.8 |
| n23a0 | 16 | 7,460 | 1,563 | 16.1 → 4.8 |
| n39a0 | 32 | 24,596 | 3,992 | 11.3 → 6.2 |
| n37a0 | 56 | 42,407 | 5,681 | 8.4 → 7.6 |
| n43a1 | 584 | 478,775 | 16,700 | 31.2 → 29.0 |
| n47a1 | 28,016 | 22,106,848 | 58,516 | 611.7 → 615.8 |
| n41a0 | 32,440 | 28,424,411 | 57,538 | 532.1 → 535.4 |

(n47a1 and n41a0 tick up by a few relations' worth of stride
multiplications, one per batch; that is the granularity's cost at the
top and it is charged.)

## 12. Iteration 5, `|F|` sized by `r` (frozen before the run)

**Model, from the corrected accounting.** With the full folded table
(`t = K`) and the window charged, one addition buys one lookup that
hits with probability `|F|²/(2r)`, so

```
G_IC ≈ |F|²/(4n)  +  (K+1) · 2r/|F|²  +  K · 1.5 · bits(r)  +  descent
       table         collection             certification
```

with `K ≈ |F|/(2n)`. The first two terms are minimised at
`|F|* = (2r)^{1/3}`, where `G_IC ≈ 3(2r)^{2/3}/(4n)`. The 6n recipe
puts `|F|` at 656 on n41a0 where `|F|*` is 10,322, and pays for it in
collection, which is 30× the table there.

**Recipe change.** `points = ⌈(2r)^{1/3}⌉` in place of `6n`, same
sampler, same seed 43. The sampler adds eight orbits per round, so
on rungs where `|F|*` is below eight orbits the base is the same
eight-orbit minimum as before and this iteration changes nothing.
Batch 8 from §11 is kept. `B` changes, so the counting floor per cell
changes with it by construction; the row states its `|F|`, and no
ratio to the old floor is quoted.

**Public parameters, derived before any target is solved**
(`build_subgroup_orbit_factor_base(kc, 43, ⌈(2r)^{1/3}⌉)`,
`optimal_folded_rows`, `build_folded_rows`; all target-independent):

| Cell | `⌈(2r)^{1/3}⌉` | `|F|` | `K` | `t` | table adds | model collection | model `G_IC` |
|---|---:|---:|---:|---:|---:|---:|---:|
| n13a0 … n31a0, n23a0, n39a0 | 16–516 | unchanged | 7–8 | as before | as before | as before | as before |
| n37a0 | 773 | 1,184 | 16 | 6 | 5,994 | 9,178 | ≈ 17,200 |
| n43a1 | 2,103 | 2,752 | 32 | 15 | 32,250 | 56,386 | ≈ 93,100 |
| n47a1 | 5,978 | 6,016 | 64 | 41 | 169,576 | 440,432 | ≈ 620,000 |
| n41a0 | 10,322 | 10,496 | 128 | 80 | 580,560 | 1,498,165 | ≈ 2,097,000 |

**Prediction.** Against the iteration-3 measured `G_rho`: n37a0
`α ≈ 3.0`, n43a1 `≈ 5.6`, n47a1 `≈ 10.6`, n41a0 `≈ 36`; the eight
lower rungs equal iteration 4. Gate fails on every rung. Class:
**engineering** against the rho reference (α falls, the method does
nothing a generic algorithm cannot).

**Extrapolation, marked as such.** At the model optimum the ratio to
rho's walk is `α ≈ 1.34 · r^{1/6} / √n`. It is below 1 only for
`r < (n/1.8)³`, which is `2^{8.5}` at `n = 13` and `2^{13.5}` at
`n = 41`, under every rung here. At `n = 131`, `r ≈ 2^{129}` it is
about `3.5 × 10^5`. This rests on the `r^{2/3}` cost of the pair
table, and on the measured slopes below, not on a measurement at
degree 131.

**Inadmissible.** Reading the eight unchanged rungs as evidence for
the recipe. Quoting the extrapolation as a measurement. Reporting
iteration 5 without the iteration-3 rows beside it.

```bash
cargo run --release --example rho_parity_e2e -- --ladder --batch 8 --recipe cbrt \
  --out experiments/ecc2k130_rho_parity_20260918/iteration-5
```

**Measured 2026-09-19**, host `ip-172-31-19-103`, cited from
[`iteration-5/summary.json`](experiments/ecc2k130_rho_parity_20260918/iteration-5/summary.json).
108/108 verified. The eight lower rungs are byte-identical to
iteration 4, as the sampler floor predicted. The four upper rungs got
the derived `|F|`, `K`, `t` and table adds exactly, and then needed
far more relations than `K + 1`:

| Cell | `|F|` | `K` | relations needed | `K ln K / 3` | trials | collection adds (model) | `G_IC` | `G_rho` | mean `α` (iter. 4 → 5, predicted) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| n37a0 | 1,184 | 16 | 30 | 15 | 48 | 56,784 (9,178) | 68,485 | 5,681 | 7.6 → **12.2** (3.0) |
| n43a1 | 2,752 | 32 | 43 | 37 | 80 | 220,080 (56,386) | 259,748 | 16,700 | 29.0 → **15.9** (5.6) |
| n47a1 | 6,016 | 64 | 96 | 89 | 344 | 2,069,160 (440,432) | 2,260,191 | 58,516 | 615.8 → **61.8** (10.6) |
| n41a0 | 10,496 | 128 | 285 | 207 | 1,072 | 11,250,640 (1,498,165) | 11,886,436 | 57,538 | 535.4 → **217.7** (36) |

**What the model missed.** A relation is one equation in three of
the `K` orbit unknowns, and the system cannot solve until every orbit
has appeared in at least one of them. With third summands landing
uniformly that is a coupon-collector count, about `K ln K / 3`, not
`K + 1`; the measured counts sit at or above it, and the solver also
had to wait for full rank past mere coverage. The relation count
enters the collection term linearly, so the collection came in 4–7×
above the model and the optimum `|F|` is smaller than `(2r)^{1/3}`
by the cube root of that factor. n37a0 is the case where the recipe
overshot: doubling `|F|` on a rung whose collection was already
cheap bought a larger table and twice the relations.

**Result.** Net of both levers the top three rungs improved 1.8×,
10× and 2.5× over iteration 3 and n37a0 worsened 1.6×. Gate fails on
every rung; best rung n23a0 at `α = 4.8` (unchanged from iteration
4). Class: **engineering**, with the quantitative prediction
falsified by the factor above and the sign of the prediction (α
falls on the rungs the recipe touches) holding on three of four.
Slopes over twelve rungs IC 0.40, rho 0.27; top four rungs IC 0.66,
rho 0.32 (rho's top-rung slope is depressed by n41a0's lucky targets,
§10).

## 13. Where this leaves the thread

In exclusive group operations, on every Koblitz cell from `2^{11}` to
`2^{39}`, the pair-table pipeline costs between 4.8× and 218× matched
signed-Frobenius rho after the accounting correction and the two
engineering rounds. The remaining structural cost is the collection
term `(relations) · 2r / |F|²`, in which each point addition buys one
lookup that hits with probability `|F|²/(2r)`; no choice of `|F|`
brings the sum with the `|F|²/(4n)` table under `√r` at these sizes,
and the model gap widens as `r^{1/6}`. The lever that is not
exhausted is the relation count: steering probes toward uncovered
orbits, or a table built so that every stored pair covers a
prescribed orbit, would bring `K ln K / 3` toward `K + 1` and recover
at most the factor between the model and measured columns above
(4–7×), which does not reach the gate on any rung. That is the
falsification: the §5 target is not reachable in this unit by any
lever this note has left, and the thread closes negative. The
tournament's Valgrind-instruction win is a different unit and is not
contradicted or confirmed by this.
