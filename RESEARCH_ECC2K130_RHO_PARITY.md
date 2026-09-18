# Toy-suite rho parity for Koblitz pair-table index calculus

**Question.** On the five small Koblitz cells the candidate tournament
already solved, does the production pair-table pipeline cost no more
than matched signed-Frobenius rho in exclusive group operations, on
every useful completed verified case?

**Answer.** Not yet measured. This note freezes the protocol before the
run. It is a toy-suite practicality gate, not degree-131 parity, not an
exponent claim, and not a refit of X1–X4.

Companion: [`RESEARCH_ECC2K130_ROUTE_TARGETS.md`](RESEARCH_ECC2K130_ROUTE_TARGETS.md)
T12. Scoreboard panel `#ecc2k130-rho-parity-20260918` is filled only
when the receipt exists.

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

Class, predicted: **engineering**. `S` may fall against rho; the ratio
to the counting floor does not, because `B` and `m` are unchanged.

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

**IC.** Table-build additions (both occupancy and fill passes).
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

**Predicted class.** Engineering: the counting floor is unchanged.
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
Runner: `examples/rho_parity_e2e.rs`.
