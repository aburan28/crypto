# Round 0025 preparation: a pair-or-triple switch, and the checker amendment it needs

Nothing here is a tournament round. It prepares one arm and one checker
change, and checks both on round 0024's closed fixtures. No trial of any
round was run and no ledger record was written. Raw output:
[round25-switch-check.txt](round25-switch-check.txt).

## 1. The checker amendment

[round25-oracle-mixed-summands.patch](round25-oracle-mixed-summands.patch),
against round 0023's frozen `evaluator/oracle.py`. Apply it with `patch -p1`
inside a copy of that directory. [evaluator-r25/](evaluator-r25/) is the
result:

| file | sha256 |
|:--|:--|
| `evaluator-r25/oracle.py` | `ae8f9b2a21b54dd28c6eada6d84ebd27df77a0986cd96a13c24ffb3c9f99510b` |
| `evaluator-r25/tournament.py` | `1e5b27743bf3bf99a3834882ddde4992f658874438576e7532646feb655651cc` (round 0023's, unchanged) |

**The problem.** `verify` required `report['summands'] == summands`, and it
required every relation and every non-degenerate descent to have exactly
`summands` indices. The tournament passes each arm's `config['summands']`.
An arm at `summands: 4` therefore could not return 3-term relations on any
fixture.

**The change.** A report declares its own `summands`. It may equal the
configured count, exactly as before. It may also be an integer `length` with
`3 ≤ length < summands`. The declared count then binds the whole report:
every relation has exactly `length` indices, and every descent has 0 or
`length`. Nothing else in `verify` changes. `tournament.py` does not need to
change, because it already hands the arm's configured count to `verify`,
and that count is now the ceiling.

**Why it does not weaken soundness.** A relation's length was never part of
the soundness argument. Under the amendment, as before:

- every relation's points are re-added in the group and must equal `[a]G`;
- every row must hold against the certified column logarithms;
- the rows must have full rank;
- every column logarithm is checked by scalar multiplication;
- every recovered scalar must satisfy `[d]G = Q`, and must follow from its
  own descent relation, which is re-added in the group.

What the configured count carries is the method: which collector produced
the relations. The amendment keeps that visible. A report shorter than its
configuration adds `relation_summands` to the certificate the tournament
writes into each receipt.

**Mixing is refused.** One declared length binds every relation and descent
in a report, and a report with mixed lengths is refused. Mixing would be just
as sound. It is refused for these reasons:

- The switch never produces it, because it runs one collector per job.
- A single length keeps "the method this report used" one fact.
- It follows round 0017's precedent, where two certificate formats were
  admitted and mixing them in one report was refused.

**Fewer than three summands is refused** (`MIN_SUMMANDS = 3`, the pair
table). Shorter relations would be sound too, but no arm emits them, and
admitting only what a registered collector produces keeps the change as
small as possible. It also means an arm configured at three is checked
exactly as before.

**Additive, checked.** When the declared count equals the configured one,
`length` is the configured value itself. That value is not a copy, so even
`3.0 == 3` is read as before. The certificate is then byte-identical, and so
is every refusal message. Tested:

- `test_equal_counts_decide_identically_on_random_certificates` runs 40
  random honest and tampered certificates.
- `test_stored_receipts_verify_identically`: every VERIFIED report of rounds
  0014 through 0020 on disk goes through both checkers. That is 34,668 profile
  and native reports under `R25_FULL=1`. The old and new checker must agree,
  and on profile reports both must equal the certificate the receipt
  recorded. See §5 for the counts.

**Caveats.**

- **An arm at `summands: 4` can now return a 3-term report whatever its
  solver.** Correctness is unaffected. Attribution is affected: a
  `triple_counted` arm that silently ran the pair table would verify. The
  receipt's certificate shows `relation_summands: 3`, but the frozen
  tournament does not act on it. If a round wants the method pinned per arm,
  its analysis has to read that key.
- **The tournament compares only `solutions` between a trial's profile run
  and its native run.** A non-deterministic arm could run different lengths
  in the two runs and still pass. The switch is deterministic: it decides
  from `(n, r, points)` alone.

## 2. The switch collector

[round25-switch.patch](round25-switch.patch), `-p1` against
`round24-sources/counted`, which is the tree round 0024 froze as
`source_candidates/counted/source`. It adds the solver `pair_or_triple`,
configured at `summands: 4`.

- **`koblitz_tiny_ic::choose_collector(n, r, points)`** returns
  `Collector::Pair` or `Collector::TripleCounted` before any base point is
  drawn. It also returns the numbers it priced.
- **The worker** runs `TinyIc::new_with(…, chosen)` exactly as the
  `pair_table` and `triple_counted` solvers do. It then reports:
  - `summands` as the collector's own count (3 or 4);
  - `collector` (`"pair_table"` or `"triple_counted"`);
  - `collector_model` (orbits, rows and modelled work of both collectors).

  The report's `summands` used to be `cfg.summands`. It is now
  `ic.collector().summands()`, which is the same value for every other solver
  because the fast path admits each collector only at its own count. The two
  new fields are written only by the switch.
- **Neither collector changed.** The patch adds functions and tests to
  `koblitz_tiny_ic.rs`, and it does not edit any existing line of that file.

### The rule

Both collectors are priced in the same units, table entries built plus rests
scanned, by the model each already sizes itself with:

- **Pair.** `W₂ = t·|F| + (K + 3)/(λ₂·cov₂(t))`, with `λ₂ = |F|(|F|+1)/(2r)`
  and `cov₂ = 1 − ((K−t)/K)²`. This is the expected work the pair collector's
  own row rule minimises. It is taken at the pair plan:
  - `K = policy_orbits(r)` where that is more than one batch;
  - otherwise whole batches of eight orbits until `2n·K ≥ points`;
  - `t` minimises `W₂`.
- **Triple.** `W₃ = t·|F|(|F|+1)/2 + (K + 3)/(λ₃·cov₃(t))`, with
  `λ₃ = C(|F|+2, 3)/r` and `cov₃ = 1 − ((K−t)/K)³`. This is the triple
  table's own model (`triple_work`), taken at the orbits and rows the counted
  sizing chooses.

**Rule: run the triple table if and only if `W₃ < W₂`.** A tie keeps the
pair table. No constant is fitted, and the threshold is the two models'
break-even. With `points = 6n`, as the tournament's jobs ask:

- break-even is near `r ≈ 5·10⁷` at `n = 29`, `1.3·10⁸` at `n = 37`,
  `2.8·10⁸` at `n = 43` and `1.4·10⁹` at `n = 61`;
- for `n ≤ 23` it never switches.

(From a Python replica of the two functions, checked against the
`collector_model` numbers the worker printed on all 12 cells.)

**Provenance, stated plainly.**

- **The model is the code's own.** `W₂` and `W₃` are the formulas already in
  the source (round 0007's pair model; `research/ic_triple_table_20260923`'s
  `W₃`). The rule reads only `(n, r, points)`, and no round-0024 number is in
  the source.
- **Pre-round-24 data agrees with it.** `research/ic_triple_table_20260923`
  (`confirm.json`, 32 fixtures a cell, another seed) measured triple/`scaled`:
  - 1.493 at `n23a1`, where the model gives 2.01 and chooses pair;
  - 0.589 at `n37a0` (model 0.77, triple);
  - 0.409 at `n43a1` (model 0.46 at that arm's two orbits, triple).

  All three are on the side the model chooses. Round 0024's pre-round probe on
  round 0023's fixtures (`round24-amended-probe.txt`) is on the model's side
  at all twelve cells too.
- **The model was picked after I had seen round 0024's results.** The
  task gave the per-cell ratios, and I compared three pricings that are
  internally consistent:
  1. `W` for both collectors;
  2. the counted sizing's per-target count for both (a full `|F|` scan per
     success);
  3. a hybrid of the two.

  Pricings 1 and 3 order all twelve cells correctly. Pricing 2 chooses the
  triple table at `n13a0`, `n17a1` and `n23a1`, which is wrong. It charges a
  full scan to every successful probe, while at small `r` a probe stops at
  its first hit. I chose 1 because it is the pair collector's existing
  model. **So the per-cell agreement below is not an out-of-sample
  confirmation of the rule.** That needs a fresh seed, and preferably cells
  near the break-even (for example `n29`–`n31` at larger `r`, and `n37` just
  above `r ≈ 10⁸`), under a pre-registration.
- **The model is not calibrated in absolute terms.** Its `W₃/W₂` departs
  from the measured ratio by 0.27× to 1.46× across the twelve cells. It
  overstates the gap at small cells, where fixed costs pull the ratio toward
  one, and understates it at `n59a0`. The rule uses only the side of one, and
  the smallest margins are 2.0 on the pair side (`n23a1`) and 0.77 on the
  triple side (`n37a0`).

## 3. Equivalence on round 0024's closed fixtures

[round25_switch_check.py](round25_switch_check.py) ran on four fixtures a
cell, 48 in all. For the eight panel cells it used development and selection
fixtures, and for the holdouts (`n19a1`, `n29a1`, `n59a0`, `n61a1`) it used
confirmation fixtures. The switch worker was built from the patch, and the
references are round 0024's frozen workers: `worker`, the pair arm, and
`source_candidates/counted/worker`, the triple arm. Instructions are callgrind
`Collected:` under valgrind 3.22.0, with all three workers copied to
equal-length paths.

- **48/48** switch reports equal the chosen collector's report from its own
  worker, field for field. Only `elapsed_seconds`, `collector` and
  `collector_model` differ.
- **48/48** switch reports verify under `evaluator-r25/oracle.py` at the
  arm's `summands: 4`. The pair cells verify as 3-term reports with
  `relation_summands: 3` in the certificate, and the triple cells as 4-term.
- All three collectors recover the same logarithm on every fixture.

| cell | rule picks | `W₃/W₂` | r24 triple/pair (confirmation) | r24 better of two | match | switch/pair, instr | switch/triple, instr |
|:--|:--|--:|--:|:--|:--|--:|--:|
| `n13a0` | pair | 6.49 | 1.967 | pair | yes | 1.0103 | 0.504 |
| `n17a1` | pair | 6.63 | 2.205 | pair | yes | 1.0065 | 0.438 |
| `n19a0` | pair | 6.73 | 2.327 | pair | yes | 1.0067 | 0.431 |
| `n19a1` | pair | 5.17 | 2.115 | pair | yes | 1.0077 | 0.462 |
| `n23a0` | pair | 2.89 | 1.542 | pair | yes | 1.0087 | 0.655 |
| `n23a1` | pair | 2.01 | 1.323 | pair | yes | 1.0088 | 0.781 |
| `n29a1` | pair | 14.07 | 3.835 | pair | yes | 1.0036 | 0.264 |
| `n31a0` | pair | 7.45 | 3.260 | pair | yes | 1.0062 | 0.318 |
| `n37a0` | triple | 0.77 | 0.524 | triple | yes | 0.423 | 1.0005 |
| `n43a1` | triple | 0.48 | 0.291 | triple | yes | 0.242 | 1.0001 |
| `n59a0` | triple | 0.44 | 0.648 | triple | yes | 0.701 | 1.0000 |
| `n61a1` | triple | 0.44 | 0.314 | triple | yes | 0.381 | 1.0001 |

These are geometric means of four fixtures. The bands are in the raw output.
The ratio against the collector the rule did not pick is these four
fixtures' reading of the triple/pair gap, not a new measurement of it.

**What the switch costs over the collector it runs:**

- **Against round 0024's own arm:** 1.0036–1.0103 on the pair cells,
  geometric means. The largest is at `n13a0`, the smallest job, where the
  per-fixture maximum is 1.0109. On the triple cells it is 1.0000–1.0006.
- **Two parts, measured.** The raw output's `switch/same-tree` column
  separates them.
  - *The decision and its report fields.* Against the same tree's own
    `pair_table` run, the switch costs 1.0018–1.0093. That is about 5,100
    instructions a job at `n13a0`, all in `factor_base_and_tables`
    (per-phase callgrind dumps of one fixture). Part of it is formatting the
    `collector_model` numbers.
  - *The tree.* Where the switch runs the pair table, it runs the counted
    tree's `Collector::Pair`. The frozen triple worker, run as
    `pair_table`, costs 0.9998–1.0070 of the matched tree's pair binary.
    The largest difference is at `n23a1`. The reports are identical.
- **Instructions only.** Native time was not measured here.

## 4. Tests

- **`cargo test --release --offline --locked --lib cryptanalysis::koblitz_tiny_ic`**
  on the patched tree: **16 passed, 0 failed**. That is the 13 of `counted`
  plus three new ones:
  - `pair_work_matches_the_collectors_row_rule`: the switch's pair pricing
    predicts the rows the pair collector actually built, on every supported
    test curve at three base sizes, and the planned orbit count at
    `points = 6n`;
  - `the_switch_rule_at_the_tournament_cells`: the rule's output at the
    twelve cells, and its margins;
  - `the_switch_runs_the_chosen_collector`: the chosen collector builds the
    same base and table as when built alone, emits relations of its own
    length, and solves.
- **[test_round25_oracle.py](test_round25_oracle.py)**, under plain Python
  and pytest: **10 passed, 0 failed**. It covers:
  - hand-built certificates (see the next paragraph);
  - a 3-term report verifying under a 4-summand configuration, where the
    unamended checker refuses it only with `changed summand count`;
  - a 4-term report verifying identically under both checkers;
  - refusals of a declared count above the configured one, of one below
    three, of non-integer declarations, and of mixed relation or descent
    lengths;
  - seven tamperings at each length (relation scalar, index, shortened,
    lengthened, descent scalar, column log, logarithm);
  - real reports from round 0024's pair and triple workers at `n19a0` and
    `n37a0`, with the same tamperings;
  - the stored-receipt identity sweep.

  The hand-built certificates use the checker's own group arithmetic on
  `n13a0`, where every logarithm can be looked up. They therefore test the
  checker's acceptance logic, not a solver.

## 5. Scope and what remains

- Every number here is on round 0024's fixtures, four a cell, with
  instructions only.
- The switch's per-cell agreement with round 0024's better-of-two is
  in-sample (§2). Scoring it needs a pre-registered round. That round must
  name the amended checker (`evaluator-r25`, sha256 above) as its evaluator,
  and must state how `relation_summands` is read.
- Where a matched rho stands is unchanged. In round 0024 the better of the
  two collectors is above rho at `n23a1` and at every crossover cell in
  instructions. The switch takes that better of two, plus at most about 1%, so
  it cannot beat a matched rho where neither collector did.

## Reproduce

```sh
S=round24-sources/counted            # from round24_amended_candidates.py
cp -a $S /tmp/sw && (cd /tmp/sw && patch -p1 < $PWD/round25-switch.patch)
(cd /tmp/sw && cargo build --release --offline --locked --example ic_tournament_worker)
(cd /tmp/sw && cargo test --release --offline --locked --lib cryptanalysis::koblitz_tiny_ic)
cd .. && python3 campaign_20260916/round25_switch_check.py \
    /tmp/sw/target/x86_64-unknown-linux-musl/release/examples/ic_tournament_worker
python3 campaign_20260916/test_round25_oracle.py        # R25_FULL=1 for every stored receipt
```
