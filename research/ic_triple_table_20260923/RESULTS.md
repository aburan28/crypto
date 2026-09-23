# Results: a triple-sum table for the tournament's IC arm

**Class: engineering** candidate for the IC tournament (`AGENTS.md` §3).  It
changes the collector's cost law, not the verdict at 131 bits, and it is
offered to the tournament's harness rather than scored outside it.

Pre-registered in `PREREGISTRATION.md` (`2f26244`, before the collector
existed), with an addendum fixing the statistic and seeds (`ea5a1db`, after a
4-fixture probe, before this run).  Four-core Xeon container; instruction
counts are callgrind `Collected:` under valgrind 3.22.0, the tournament's own
unit, which does not depend on the host.  Both arms built by `build_arms.sh`
from committed inputs.

## The collector

`triple-table.patch`, against `scaled` (round-0020 `both` +
`round21-wide-pair-table.patch` + `round23-scaled-base.patch` — the tree round
0023 ran and round 0024's `pairinv` builds on):

- a `Collector` mode.  `Pair` is the incumbent, unchanged; `TinyIc::new` still
  builds it, so every existing caller and test is untouched;
- `Triple` stores `C_o + P_j + P_l` (`j ≤ l`) under the same Frobenius orbit
  names, scans one summand per probe exactly as the pair collector does, and
  carries a hit back to four base indices by the same rotation and sign;
- its own base size, the `K` minimising `W₃`, floored at two orbits because
  the checker needs more than one column, sampled as one batch of that size;
- the worker routes `solver: "triple_table"` with `summands: 4` to it, and
  only that pairing: a `pair_table` job at four summands still takes the
  general path, as it did before.

## Correctness

**96 of 96 fixtures verified by `oracle.py`** on both arms (the pair arm at
three summands, the triple arm at four), and **every pair recovered the same
logarithm**.  The solver's own tests pass 11 of 11 (`unit-tests.txt`),
including a new end-to-end test that re-adds every triple relation and every
descent witness in the library's general arithmetic on seven curves at three
seeds each.

## Confirmation: triple against `scaled`

32 paired fixtures a cell, seeds from `random.Random(20260924)`, disjoint from
the probe.  Geometric mean of the paired instruction ratios, 95% percentile
bootstrap interval (10,000 resamples) — the statistic fixed before the run.

| cell | `r` | `scaled` base | triple base | **triple/`scaled`** | median | range | median trials, `scaled` → triple |
|---|--:|--:|--:|--:|--:|--:|--:|
| `n23a1` | 4,196,903 | 8 orbits | 2 | **1.493** [1.411, 1.575] | 1.492 | 0.91–1.99 | 11 → 6 |
| `n37a0` | 230,603,167 | 16 | 2 | **0.589** [0.535, 0.643] | 0.585 | 0.29–1.05 | 32 → 32 |
| `n43a1` | 4,644,189,029 | 24 | 2 | **0.409** [0.337, 0.494] | 0.429 | 0.17–1.25 | 123 → 363 |

`confirm.json` holds every row.

## Scored against the pre-registration

| criterion | registered | measured | |
|---|---|---|---|
| correctness gate | every report verifies, same logarithm | 96/96, 96/96 | **met** |
| `n43a1` | in or below `[0.563, 0.649]` | 0.409, interval `[0.337, 0.494]` | **met — below** |
| the gain grows with `r` | `n23a1` > `n37a0` > `n43a1` | 1.493 > 0.589 > 0.409, intervals disjoint | **met** |
| refuted if | `n43a1` ≥ 0.85, or no growth with `r` | neither | not refuted |

**Verdict: supported.**

Two misses, stated:

- **`n23a1` misses its band, 0.938–0.951, and loses outright at 1.49.**  The
  cause was named before the probe and not quantified: the model's optimum
  there is one orbit, the checker needs two, and at two orbits the triple
  table costs more than the whole pair-table job.  Quantified after the probe
  (not a prediction), the model puts it at 1.65–1.81; it measured 1.49.
- **`n37a0` and `n43a1` beat their bands**, by 0.23 and 0.15.  The prediction
  held the non-collection share fixed and flagged that as conservative; with
  two orbits instead of 16 or 24, base construction and linear algebra shrink
  too.  The model is right about the direction and the growth and
  conservative about the size.

## What this does and does not show

- **It is a rate change, measured on three cells.**  The gain grows with `r`
  with disjoint intervals, which a constant-factor lever cannot do.  Three
  cells do not fit an exponent; the `r^{1/10}` is the derivation's, and this
  is consistent with it rather than a measurement of it.
- **It is not scored by the tournament.**  Its gates are native time as well
  as instructions, confirmation and replay stages, and the parity rule.  None
  of those was run here.  And `OPERATIONS.md` compares arms at matched summand
  counts, which a four-summand arm is not; admitting it may need an amendment,
  which is that lane's decision.
- **The obvious next arm is a switch, and it is not measured here.**  The pair
  table wins at `n23a1` and the triple table at `n37a0` and `n43a1`; a
  collector choosing between them by the model at run time would keep the
  first and take the other two.  That design came from this data, so it needs
  its own pre-registration before it is scored.

## Secondary, not pre-registered: against rho on the same fixtures

Added after the confirmation run.  **No decision above rests on it**, and the
tournament's own gates decide `beats_rho_strict`.  It exists because the
tempting alternative — multiplying the triple/`scaled` ratio by the
tournament's scaled/rho, measured on other fixtures — is an inference rather
than a measurement.

The worker's matched rho (`automorphism_order = 2n`, oracle-verified in its
own mode, 96 of 96 complete) on the same 96 fixtures, instruction counts only;
`rho-same-fixtures.json`.

| cell | `scaled`/rho | **triple/rho** | fixtures below rho, `scaled` → triple |
|---|--:|--:|--:|
| `n23a1` | 0.689 [0.621, 0.769] | 1.029 [0.940, 1.130] | 26 → 14 of 32 |
| `n37a0` | 1.381 [1.152, 1.649] | **0.813 [0.700, 0.941]** | 8 → 19 of 32 |
| `n43a1` | 2.281 [1.864, 2.754] | **0.933 [0.753, 1.159]** | 4 → 19 of 32 |

- **The harness reproduces the tournament's crossover.**  `scaled`/rho reads
  0.69, 1.38 and 2.28 against the tournament's 0.83, 1.40 and 2.21 (round
  0022), on different fixtures and a different build host.  That is the check
  that this measurement is measuring the same thing.
- **At `n37a0` the triple collector is below rho**, the whole interval under
  one, where `scaled` is 1.38× rho.  **At `n43a1` the point estimate is below
  rho** and the interval spans one, where `scaled` is 2.28×.
- **In instructions only.**  The tournament gates on native time too, and in
  fifteen rounds a candidate has repeatedly passed on instructions and failed
  on native time.  Nothing here says what that gate would do.
- **A switch would read 0.69, 0.81 and 0.93 — below rho at every cell measured.**
  That is the pair collector at `n23a1` and the triple at the other two, and it
  is arithmetic on this table, not a run.  Measuring it is the next step, under
  its own pre-registration.

## Reproducing

```sh
research/ic_triple_table_20260923/build_arms.sh /abs/work       # restores round 0020, applies the patches, builds both
python3 research/ic_triple_table_20260923/check.py \
    /abs/work/scaled/target/x86_64-unknown-linux-musl/release/examples/ic_tournament_worker \
    /abs/work/triple/target/x86_64-unknown-linux-musl/release/examples/ic_tournament_worker \
    --fixtures 32 --seed 20260924 --out confirm.json
python3 research/ic_triple_table_20260923/model.py              # the prediction
```

`check.py` draws seeds in the confirmation run's order, so it reproduces its
fixtures.  **What "reproduces" means here, measured:** the committed harness,
invoked the way the run was (relative worker paths from the same directory),
returns `scaled`'s count for the first fixture exactly and the triple arm's
within 6 instructions of 3.97 million.  Invoked with absolute paths, both arms
move by about 450 instructions — the worker's start-up sees a longer path —
which cancels in the ratio (1.38526 against 1.38532).  A count from this
directory is comparable to another count taken the same way, and a ratio is
comparable across invocations; a raw count is not comparable across
invocations to better than a few hundred instructions.
