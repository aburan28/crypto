# Protocol: multiplier closure in the 4-sum pair-table index calculus

Frozen 2026-10-10, before the registered runs.  **Stage diagnostic on toy
Koblitz curves.**  `S`, end-to-end cost and speedup are **unset**.  Nothing
below is a claim about ECC2K-130 or any deployed curve.

## Question

The compact-orbit producer (`docs/ic/COMPACT_ORBIT_PRODUCER.md`, the landed
`n = 61 … 83` rungs) builds its summand set from `K` columns and the
multipliers `{±τ^k}`, `2n` per column, and finds 4-sum relations by probing
a pair table of `x(s_a + s_b)`.  `RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md`
§A1 proposes closing the multiplier set under small-norm elements of
`Z[τ]` ("loops of isogenies" at the crater).  The companion note
`docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md` derives that in this
pipeline any set of cheap known-log multiples plays the same role as the
Frobenius orbit in the *yield* model, so the `Z[τ]` framing should be
cosmetic: plain integer multiples `{1, …, c}` should do the same.  It also
derives that the multiplier closure cannot change the per-target probe
count at equal table memory; it can only trade columns for multipliers
(fewer relations, smaller linear algebra), and that this trade must stop
when the multiples of one column start to coincide additively.

## Instrument (`multipliers.rs`, standalone Rust, `rustc -O`)

Curve `y² + xy = x³ + 1` over `F_{2^n}`, `#E = cofactor · r`, `r` prime;
`G` of order `r`; `μ` the Frobenius eigenvalue on `⟨G⟩` (checked
`τ(G) = [μ]G`, `μ² + μ + 2 ≡ 0`, `μ^n ≡ 1`).  Each cell fixes `Kc`; every arm
has `K · c = Kc` columns times multipliers, so the summand set has
`L = 2n · Kc` points and the pair table the same nominal size:

| arm | multipliers per column | K |
|:--|:--|--:|
| `frob` | `{±τ^k}` | `Kc` |
| `tau` (c) | `{±α τ^k : α ∈ A_c}`, `A_c` the `c` smallest-norm elements `a + bτ` with `a` odd, `a > 0` | `Kc / c` |
| `int` (c) | `{±a τ^k : a ∈ {1, …, c}}` | `Kc / c` |

`α(P)` is computed from the point as `aP + b τ(P)` and self-checked against
`[a + bμ]P`.  The table stores `(x(s_a + s_b), a, b)` for `a < b`; a probe
draws a random pair `(c, d)`, forms `T = R − s_c − s_d`, and matches `±T`
against the table.  A hit is a 4-sum relation; its row over the `K`
columns is checked against the hidden column scalars before use.  Per arm:
`hits` targets `ρG` with fresh `ρ` (first hit each), echelon over `Z/r`
tracked after every row, columns solved and compared with the hidden
scalars, then one hidden target decomposed and its scalar verified.

Memory metric: the number of distinct Frobenius-canonical `x` among the
table entries, which is what a τ-canonicalised index (the landed design)
would store.  Model: `L² / (4n)`.

## Derivations stated before the run

- **Yield law.**  A probe hits with probability about `L² / (2r)` (each of
  the `L⁴/24` 4-subsets sums to `R` with probability `1/r`, and each is
  reached through 6 of the `L²/2` probed pairs).  Mean first-hit probes
  `≈ 2r / L²`, the same in every arm of a cell because `L` is the same.
- **Equal-memory trade.**  Relations needed equal `K`, so total probes to
  full rank scale as `K · 2r/L²`: the `tau`/`int` arms with `K = Kc/c` need
  `c` times fewer probes than `frob` for the same table.
- **Coincidence loss.**  Multiples of one column coincide additively
  (`aP + bP = (a+b)P`, and the `Z[τ]` analogue), so with few columns and
  many multipliers the set of distinct pair sums, and with it the hit
  rate, falls below the random-set model.  The loss should be visible in
  the canonical-count metric and in the probes of the `K = 1` arms.

## Disclosed smoke tests (not cited as results)

- `n = 19, Kc = 6, 16 hits, seed 1`: all arms recovered logs and target;
  `rows_to_full_rank = K` everywhere; mean probes 3.8 (`frob`) to 8.2
  (`tau K=1`) against `2r/L² = 5.0`.
- `n = 23, Kc = 12, 64 hits, seed 1`: `frob` 14.2 against `2r/L² = 13.8`;
  `tau`/`int` 13.2–16.7 at `K ≥ 4`, 16.6–18.9 at `K = 2`, 39.6/33.2 at
  `K = 1`; canonical count 3176 (`frob`) falling to 1834/1961 (`K = 1`).
  Both smoke runs were used to fix the band widths below.

## Registered cells

| cell | n | r | Kc | L | pairs | hits | seed |
|:--|--:|--:|--:|--:|--:|--:|--:|
| C23 | 23 | 2,095,853 | 12 | 552 | 152,076 | 64 | 20261010 |
| C37 | 37 | 230,603,167 (cofactor 596) | 36 | 2,664 | 3,547,116 | 64 | 20261010 |
| C41 | 41 | 549,756,390,943 | 72 | 5,904 | 17,425,656 | 48 | 20261010 |

Arms per cell are chosen by the program: `frob` at `K = Kc`, then `tau` and
`int` at four divisors of `Kc` (the largest proper divisor, two middle
ones, and `K = 1`).

```sh
rustc -O -o multipliers multipliers.rs
./multipliers 23 20261010 12 64 results/c23.json
./multipliers 37 20261010 36 64 results/c37.json
./multipliers 41 20261010 72 48 results/c41.json
python3 score.py
```

## Predictions (pass/fail)

- **M1 (yield law).**  In every cell the `frob` arm's mean first-hit
  probes lie in `[0.5, 2] × 2r/L²`.
- **M2 (multiplier equivalence at K ≥ 3).**  Every `tau` and `int` arm with
  `K ≥ 3` has mean first-hit probes within `[0.7, 1.5] ×` the `frob` arm's
  in the same cell.
- **M3 (Z[τ] is cosmetic).**  For every `(cell, K)` pair the ratio
  `tau / int` of mean first-hit probes lies in `[0.6, 1.6]`.
  *Decisive outcome that would revive the `Z[τ]` framing:* `tau ≤ 0.67 × int`
  at every `K ≥ 2` in every cell.  Not expected.
- **M4 (coincidence loss at K = 1).**  In every cell, both `K = 1` arms
  have mean first-hit probes `≥ 1.5 ×` the `frob` arm's, and canonical
  counts `≤ 0.75 ×` the `frob` arm's.
- **M5 (rows).**  `rows_to_full_rank ≤ K + 2` and no degenerate (all-zero)
  rows in every arm.
- **M6 (correctness).**  `logs_ok` and `target_ok` in every arm of every
  cell.
- **M7 (equal-memory trade).**  For every `tau`/`int` arm with `K ≥ 3`,
  `total_probes_to_full_rank(frob) / total_probes_to_full_rank(arm)` lies
  in `[0.4, 2.5] × (Kc / K)`.

## Decision rule

- M1–M7 pass and the decisive outcome of M3 does not occur: the
  multiplier closure is a constant-factor precompute trade (relations and
  linear algebra divided by `c`, per-target probes unchanged), available to
  any cheap known-log multiplier set and bounded below by the coincidence
  loss; `Z[τ]` buys only free evaluation.  Class: **stage diagnostic,
  closed as a yield lever; open only as a precompute-size engineering
  option at the landed rungs.**
- M3 decisive: `Z[τ]` multipliers have a structural advantage over
  integers; investigate the mechanism before any further claim.
- M2 or M4 fail: the coincidence loss sets in earlier or later than
  derived; report the measured `K` threshold and redo the derivation.
- M5 or M6 fail: a bug; halt.

## Inadmissible moves

Changing `Kc`, `hits`, the seed or the arm set after the run; reading a
smoke test as a result; any statement about `n = 131` beyond the
arithmetic of the companion note; any `S` or end-to-end figure.

## Amendment after the registered runs (2026-10-10, written after C23 and C37 completed and C41's first arm had reported)

Two registered choices were wrong and are kept on record as failures of
M5/M6 in `results/score.json`:

1. **C41 registered `hits = 48` for a `K = 72` arm.**  An arm cannot reach
   full rank with fewer rows than columns, so the `frob` arm of C41 fails
   M5 and M6 by construction.  This is a protocol error, not a property of
   the method.
2. **`hits = 64` at `K = 36` (C37) left one column untouched** in the
   registered seed: the replay with per-column diagnostics
   (`results/c37_replay_touches.json`, identical means to the registered
   run) shows column 5 with zero touches in 64 rows, rank 35.  With sparse
   4-term rows the probability that some column of 36 is missed in 64 rows
   is about 2 %; the seed-1 diagnostic (`results/diag37_seed1.json`, 200
   hits) reached rank 36 in exactly 36 rows.  M5's bound `K + 2` did not
   allow for this tail and also failed on C23's `frob` arm (18 rows for
   rank 12).

Follow-up cells, same `Kc` and arm set, new seed `20261011`, `hits ≥ 2K`
for every arm, scored with the same M1–M7 and reported beside the
registered cells, never in their place:

| cell | n | Kc | hits | seed |
|:--|--:|--:|--:|--:|
| C37b | 37 | 36 | 128 | 20261011 |
| C41b | 41 | 72 | 160 | 20261011 |

The instrument gained two output fields before the follow-ups
(`column_touches`, `rank_trace`); the algorithm and the random stream are
unchanged, verified by the C37 replay reproducing the registered means
(63.859 probes) to the last digit.
