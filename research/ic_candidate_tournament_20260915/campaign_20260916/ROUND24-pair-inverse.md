# Round 0024 — one inverse per ± pair, and the crossover cells at forty fixtures

Pre-registered before any round-0024 stage ran. Everything here that reads
measurements reads round 0023, which is closed, or the pre-round check in §4,
which draws no fixture from round 0024's seed.

## 1. Why this round exists

Round 0023 ([RESULTS.md](RESULTS.md#round-0023-the-crossover-cells-in-the-harness-and-the-base-rule-out-of-sample))
put the crossover cells in the harness and left two things open.

* **Past the crossing, the job is collection.** It is 65% of `scaled`'s
  instructions at `n37a0`, 81% at `n43a1`, and 78% at `n61a1`. Inside collection
  at `n43a1`, the batch inversion of the decomposition scan is 35% and orbit
  naming (`NormalBasis::canon`) is 31%. Anything that moves the crossover by a
  constant factor has to move one of those two.
* **Twelve fixtures do not pin rho at the crossover cells.** Prediction 5 missed
  at `n37a0` (1.085 against [1.15, 1.70]) because rho's median iteration count
  on those twelve fixtures ran 1.195× its model. Round 0019 made the same
  diagnosis at `n23a1` and fixed it by allocation. This round does the same.

**This round is not expected to beat rho, and says so here rather than after.**
The collector's cost over rho grows about `r^0.16` (round 0022). A saving in
collection alone buys a constant factor, and at `n43a1` that constant would
have to be 2.6.

## 2. The candidate

**`pairinv`** (`round24-pairinv.patch`, on top of `scaled`). The base stores
each point next to its negative, and `+P` and `-P` share an abscissa. The scan
computes the rest `T - Q` for every base point `Q`, and its denominator is
`T.x + Q.x`. That is therefore the same field element for `Q` and `-Q`.
`scaled` puts both into Montgomery's batch inversion. `pairinv` puts in one and
reads it twice.

The sharing is keyed on an **abscissa comparison with the previous block
position**, not on the base's layout. It is correct for any order, and falls
back to one inversion per point where neighbours differ. Scan order, rests,
witnesses, relations, factor base and logarithms are untouched.

Canonical orbit naming was considered and rejected as a lever for this round.
It is already a run-doubling search over nibble-table coordinates. The cheap
Frobenius invariants of its input (popcount, longest zero run) cannot filter,
because a pair table of about 25,000 orbit names fills every bin they define.

## 3. The round

| | |
|:--|:--|
| seed | `2026092424` |
| cells | `--cells 13a0,17a1,19a0,23a0,23a1,31a0,37a0,43a1` (round 0023's) |
| holdouts | `--holdout-cells 19a1,29a1,59a0,61a1` (round 0023's) |
| allocation | `--confirmation-cases n23a1=40,n37a0=40,n43a1=40`; 228 confirmation cases |
| profile | pilot, one target, objective `rho`, `--require-native-progress` |
| comparison | `--comparison-kind factor-base-policy` |
| limits | `--timeout 300` (round 0023's worst trial: 7.5 s), `--max-processes 6000` |
| arms | `incumbent`, `scaled`, `pairinv` (`round-0024-single-candidates.json`, built by `round24_candidates.py`) |
| trials | 4,824: aa 48, smoke 96, development 288, selection 288, confirmation 2,052, replay 2,052 |

**The baseline stays round 0023's incumbent.** Under objective `rho` a
challenger is promoted only if it beats rho at every cell, so `scaled` was not
promoted, although it passed every gate against the incumbent. Keeping the same
baseline keeps the incumbent ratios comparable across rounds 0023 and 0024.
`scaled` is carried as the reference that `pairinv` must reproduce exactly.

**The allocation raises `n37a0` and `n43a1` from twelve to forty**, as round 0019
did for `n23a1`, which stays at forty. Raising a count cannot buy a pass
(`parse_confirmation_allocation` refuses a lowering). It only resolves an
estimate. The holdouts stay at twelve: they are out-of-sample tests of the base
rule, which round 0023 already passed.

**Selection takes the two cheapest eligible development arms and confirms the
cheaper of them.** Only that arm enters confirmation and replay, so the
`pairinv`/`scaled` comparison is measured in development and selection and in
§4. Confirmation scores the winner of selection against the incumbent and rho.

## 4. Before the round: the arm's tests, and equivalence

`cargo test --release --offline --locked --lib cryptanalysis::koblitz_tiny_ic`
on the `pairinv` tree: **9 passed, 0 failed**, including the base-construction
test that caught round 0023's first patch.

[round24_pairinv_check.py](round24_pairinv_check.py) runs round 0023's frozen
`scaled` worker and a `pairinv` build on round 0023's own fixtures. It uses six
a cell from development and selection, or from confirmation for the two
holdouts, under the arm config round 0023 froze. It compares the whole parsed
report except `elapsed_seconds`, the one field that differs between two runs of
the same worker. Every `pairinv` report also goes through `oracle.py`.

**Identical reports on all 72 jobs** ([raw output](round24-pairinv-check.txt)), every one verified. Instructions:

| cell | `pairinv` / `scaled` | 95% band | min – max |
|:--|--:|:--|:--|
| `n13a0` | 0.9964 | [0.9961, 0.9968] | 0.9960 – 0.9971 |
| `n17a1` | 0.9956 | [0.9951, 0.9961] | 0.9947 – 0.9963 |
| `n19a0` | 0.9883 | [0.9846, 0.9920] | 0.9796 – 0.9922 |
| `n19a1` | 0.9856 | [0.9821, 0.9892] | 0.9800 – 0.9911 |
| `n23a0` | 0.9788 | [0.9756, 0.9821] | 0.9739 – 0.9836 |
| `n23a1` | 0.9689 | [0.9636, 0.9741] | 0.9594 – 0.9791 |
| `n29a1` | 0.9974 | [0.9972, 0.9977] | 0.9970 – 0.9979 |
| `n31a0` | 0.9788 | [0.9747, 0.9829] | 0.9731 – 0.9857 |
| `n37a0` | 0.9350 | [0.9306, 0.9393] | 0.9243 – 0.9388 |
| `n43a1` | 0.9170 | [0.9155, 0.9184] | 0.9155 – 0.9205 |
| `n59a0` | 0.9604 | [0.9553, 0.9655] | 0.9517 – 0.9699 |
| `n61a1` | 0.9124 | [0.9089, 0.9159] | 0.9088 – 0.9194 |

That is less than the profile suggested. If halving the batch's elements
halved `batch_inv`, the saving at `n43a1` would be about 0.81 × 0.35 / 2 ≈ 14%;
it is 8.3%. Where the gap goes is not measured here. The saving does follow
collection's share of the job: `n59a0`, where curve setup is half the job,
saves less than `n43a1` or `n61a1`.

## 5. Predictions

1. **`pairinv` is `scaled`, computed more cheaply.** On every development and
   selection fixture the two arms build byte-identical factor bases, and all
   arms recover the same logarithm on every stage.
2. **Below one everywhere, by the §4 amounts.** `pairinv`/`scaled` below one
   at every cell measured in development and selection, and in
   [0.90, 0.95] at `n43a1`, [0.91, 0.96] at `n37a0` and [0.95, 0.99] at
   `n23a1`. The bands are wider than §4's because the fixtures are new.
3. **Selection confirms `pairinv`.** It is the provisional challenger.
4. **Against the incumbent,** `pairinv`/incumbent over the whole panel on
   confirmation in [0.68, 0.87] in instructions, with the upper 95% limit below
   one. §4's ratios composed with round 0023's `scaled`/incumbent cells give
   0.772. The band is round 0023's interval width, because the crossover cells'
   `scaled`/incumbent ratio varies by fixture.
5. **Past the crossing, rho still wins, and now at forty fixtures.**
   `pairinv`/rho above one at `n37a0`, `n43a1`, `n59a0` and `n61a1`, with
   `n37a0` in [1.00, 1.50] and `n43a1` in [1.8, 2.6]. On the eight panel
   cells, below one in instructions, with `n23a1` below 0.95.
6. **The decision.** `pairinv` passes every gate against the incumbent on both
   final stages and fails `rho_gate` at the crossover cells: `retained`,
   `beats_rho_strict` false. **This is the expected outcome, not a failure of
   the round.**

## 6. Falsification and scope

Falsification: any bad, missing or unmatched certificate; any unfinished trial
in any arm or in rho; any logarithm that differs between arms on the same
fixture; any factor base that differs between `scaled` and `pairinv`.
`pairinv`/`scaled` above one at any cell voids prediction 2 and the patch's
stated mechanism.

Every conclusion is scoped to these twelve cells, one seed, the pilot profile,
this collector and this checker. §4 reads round 0023's fixtures, so its numbers
are in-sample for nothing in this round. Nothing here bears on sect163k1 or on
any curve outside the census, and no exponent claim follows from a constant
factor.
