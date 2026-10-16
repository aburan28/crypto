# R09: the scan keyed from slopes, pipelined

Declared on 2026-10-08, before any R09 timed run, on v3 (main's
`995ea207`). Its class is not the reference class: it is the AVX-512
class without VPCLMULQDQ, where the scan's subtraction is the scalar
product. This host, a Cascade Lake, is in it.

## Hypothesis

On this class the collection scan costs 38–41 ns a scanned summand at the
three target sizes, in four parts
([stage record](../../explorations/slope-keys-cl-20261008/README.md#the-scans-stages-on-this-host)):

| curve | `log₂ r` | subtract | key | filter | admitted |
|:--|--:|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 15.3 | 12.0 | 6.4 | 4.6 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 15.5 | 13.4 | 6.7 | 5.2 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 15.5 | 13.7 | 8.6 | 3.3 |

The candidate removes work from the first two parts, and overlaps the last
two with the arithmetic:
- **Keys from slopes.** A rest's abscissa is `x₃ = λ² + λ + x_R + x_k + a`.
  The basis change is linear and turns squaring into a rotation by one,
  so `coords(x₃) = c + rotl(c) + d_k + e` with `c = coords(λ)`,
  `d_k = coords(x_k)` kept once per base point, and `e = coords(x_R + a)`
  once per trial. The key needs one basis change, as before, and no
  abscissa: the squaring, the sum and its store are spared for every
  summand the filter rejects, and only the admitted ones are completed.
- **Eight chains in the slopes' batched inversion**, the same inverses: on
  this class a carry-less multiply takes seven cycles, three of them
  chained in each product, and four chains leave the multiplier idle.
- **The rotation step in fewer micro-ops:** `v + v` for `v << 1`, and a
  per-lane variable shift for the shift by a count in a register, which is
  two micro-ops on Skylake-SP and Cascade Lake. The same value.
- **The scans of a run's trials in three stages:** trial `i` is prepared
  while trial `i − 1`'s filter words come in, trial `i − 1` is filtered
  and its admitted keys' bucket entries sent for, and trial `i − 2` is
  resolved after its runs were sent for. Each trial's answer depends only
  on its own point and scan, and the answers come back in trial order.
- **The folded build keys its rows from slopes** the same way.

## The phase and its share

On this class the scan is 57–59% of cold time at `2^44.3` and `2^44.5` and
85% at `2^47.2`; the build is most of the rest (stage record).

## Class

Engineering, by AGENTS.md §3: every decision the pipeline makes is v3's.
The keys are the same keys, so the table and every relation, logarithm and
counter are the same.

## The candidate

[`candidate.patch`](candidate.patch): eight commits on `995ea207`, ending at
`8b3d24d7` (SHA-256 in `candidate.patch.sha256`); the last only runs rustfmt.
- `koblitz_fast.rs`: `FastCurve::slopes_many` and `sum_from_slope`,
  `has_lane_kernel`, `FrobeniusCanon::slope_keys_paced` (and
  `slope_keys_on`, which names the path for the tests), and the rotation
  step.
- `semaev_decomp.rs`: `Gf2::batch_inv8`.
- `koblitz_index_calculus.rs`: `PreparedScan`, `prepare_scan`,
  `filter_scan`, `prefetch_scan_runs`, `resolve_scan`, `finish_scan` and
  `scan_pipelined`; `collect_walked`'s runs through them; the folded
  build's rows from slopes; each base point's normal coordinates kept on
  the table.
- **Where it applies:** a folded table with the canonical key, on a curve
  whose subtraction is the scalar product. Where the eight-lane kernel
  runs (VPCLMULQDQ and a field it suits), the one-call scan and the keyed
  build keep it, unchanged.

The tests, on the candidate's tree, all pass:
- `pipelined_scans_answer_as_the_one_call_scan`: the stages against
  `decompose_fast_scan`, trial for trial and in order, with trials that
  do not scan and scans the stages refuse mixed in, for runs of every
  length from none to past the pipeline's depth;
- `slope_keys_equal_the_rests_keys`: at four degrees, with degenerate
  sums mixed in, the slopes are `add_many_lazy`'s, the slope keys are the
  rests' keys on the portable path and on the AVX-512 one, and a sum
  rebuilt from its slope is the sum;
- `the_partitioned_fold_build_stores_the_two_pass_table`: the slope-keyed
  build stores the scalar build's table;
- the existing key tests (`koblitz_fast::`) and the scan tests
  (`a_folded_table_s_scan_finds_what_the_full_table_s_scan_finds`,
  `an_index_scan_s_third_summand_is_always_one_it_was_given`,
  `an_aimed_unit_s_relations_still_verify_in_the_group`) unchanged.

## Hardware class

**R09's class:** x86-64 with AVX-512F and PCLMULQDQ, and without
VPCLMULQDQ (`icprog host-class r09`: `avx512f`, `pclmulqdq`,
`!vpclmulqdq`, the `!` marking a feature the class excludes).
- **The result holds for that class.** Skylake-SP and Cascade Lake are in
  it; the reference class is not.
- **The runner enforces it:** `icprog run r09` refuses every timed step on
  a host outside the class.
- **Elsewhere:** where VPCLMULQDQ is present the eight-lane kernel runs
  where it did, and the candidate's path runs only at fields that kernel
  does not take; no claim is made for those. On a CPU without AVX-512 the
  slope keys take the portable rotation, which the tests check against the
  rests' keys; no claim is made there either.

## Main has moved under v3 (disclosed)

The plan's standing rule (§8, 2026-10-05) asks for a drift round like R07
before the next declared round once main has moved under the baseline.
Main has moved from `995ea207` to `4b02be35` (335 commits) since R07.
R09 is declared on v3 without a timed drift round first, because none of
that drift is on the path a round times:
- **The scan and the build:** `koblitz_fast.rs` and `semaev_decomp.rs`
  are unchanged; in `koblitz_index_calculus.rs` main changed only the SAT
  decomposition's coordinate domains (`u128` codes), which `ic price`'s
  pipeline does not call.
- **The rest of `ic price`'s cold time:** unchanged. Main's other changes
  to `ic` add a `large-prime` subcommand and embed the build commit.
- **The strong rho** (`koblitz_strong_rho.rs`) changed. It is the rule
  comparison's reference, not a round's arm, and cold time excludes it.
- **The check:** before R09's timed steps, main's head `ic` is built and
  pinned on all 90 suite rows untimed, beside the candidate's pin. If any
  of its outputs differs from v0's, R09 stops and a drift round runs
  first.

When R09's results reach main, the results pull request re-bases the
candidate on main's head, and a drift round on this class re-bases the
class's baseline there.

## Arms

- **The base:** v3, main's `ic` at `995ea207`, the binary R07 accepted
  (`ic-main-995ea207`).
- **The candidate:** the base plus `candidate.patch`, built with
  `IC_BUILD_COMMIT` and v3's `Cargo.lock`.
- **The manifest** records both binaries' SHA-256 and the commits they
  were built from.

## Before any timed step

1. **The tests,** on the candidate's tree:
   `cargo test --release --lib -- koblitz_fast koblitz_index_calculus semaev_decomp`.
2. **The pin, untimed:** the candidate on all 90 suite rows. Every output
   must equal v0's.
3. **The A/A:** the base against a byte-identical copy on `M1`'s 22 rows,
   five rounds, on this host: the A/A bands are the round's own.

## Rows

- **The suite rows:** `M1`'s 22 rows, and the other six suite rows at
  each target size: 40 rows, five rounds, ABAB, isolated, 400 processes.
- **The fresh holdouts:** eight rows at each target size, by suite v1's
  own construction (`icprog holdouts`):
  - recipe seeds 226 to 229;
  - two targets a seed, `T143` to `T150`;
  - `public_hash_seed` 23143 to 23150, and rho seeds `0x230000 + 143`
    to `+ 150`.

  They are committed in [`holdouts/`](holdouts/) with their
  `SHA256SUMS`. No round and no exploration has run them. 24 rows, five
  rounds, ABAB, isolated: 240 processes.
- **The extension.** A set whose interval half-width exceeds 3% after five
  rounds gets rounds 6–10, pooled.

## Power

The final exploration's per-pair spread of the log ratio was 0.029, 0.057
and 0.035 at the three sizes. At forty pairs a set, that gives a 95%
half-width of about 0.9%, 1.8% and 1.1%. The interval clears 1.05 with
80% power when the true ratio is at least about 1.06, 1.08 and 1.07.

## Prediction

| curve | `log₂ r` | exploration, base over candidate | 95% interval | predicted |
|:--|--:|--:|:--|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 1.119 (8 pairs) | [1.093, 1.146] | **1.09–1.15×** |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 1.112 (8 pairs) | [1.060, 1.166] | **1.06–1.17×** |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 1.103 (8 pairs) | [1.071, 1.136] | **1.07–1.14×** |

Elsewhere the scan and the build are a smaller share of cold time, and
both are keyed the same way at every size, so no size is predicted to
regress.

## Success and stop

**Accepted** if all of the following hold:
1. The tests pass, and every pinned output is identical.
2. At each target size, the paired cold-time ratio's 95% interval lies
   above 1.05, on the suite rows and on the fresh holdouts separately.
3. No size regresses beyond its A/A band.

**Rejected** otherwise. **Stopped** if a test fails, an output differs,
or a verification fails.

**If accepted**, the candidate becomes the baseline of R09's class, with
its binary hash, host and `S` from this round's runs. The reference class
keeps v3 and its declared rounds (R06, then R08), which are measured on
v3's code; a later round brings the two lines together.

## Inadmissible

- Reusing an exploration's runs, pooling them with R09's, or choosing
  between them.
- Changing the candidate, the suite, the recipes or the targets, or drawing
  new holdout seeds after a run.
- Pooling contended runs.
- Quoting the scan, the build or collection as the speedup: the speedup is
  cold time.
- Crediting the gain to one part of the candidate: the round measures the
  five parts together; the exploration's micro-benchmark splits them, as a
  stage diagnostic.

## Cost

- **The tests, the pin and the A/A:** about 90 untimed and 220 timed
  processes, about an hour and a half.
- **The suite rows and holdouts:** 640 timed processes, about three
  hours, and up to 640 more if extended. They share the benchmark lock
  with Track B's campaign, which runs on this host at the same time; every
  process is isolated, and the lock keeps them from overlapping.

## Order

- **R09 runs on v3, now,** on this host. It does not wait for R06: the two
  are on different classes, and neither's code changes the other's
  measurement.
- **The commands** are R08's, with `r09` for the round:
  `icprog run r09 <step>` for `manifest`, `pin`, `aa`, `compare`,
  `holdout` and `extend`, then `icprog analyse r09`.

## Explorations before this declaration (disclosed)

All are in
[`../../explorations/slope-keys-cl-20261008/README.md`](../../explorations/slope-keys-cl-20261008/README.md),
on this host, on 2026-10-08:
- **The scan's stages on this class:** R04's probes on v3, three rounds of
  `M1`'s rows at the three target sizes.
- **A2 re-checked here:** huge pages on the pair tables, two rounds.
- **The pipelined scan alone** (`7a1d132e`), one round, stopped; and the
  slope-keyed scan without the build (`1e44a603`), one round, stopped.
- **The micro-benchmark** (`examples/pipe_scan_bench`), which split the
  scan's halves and priced each lever.
- **The final candidate's exploration** (`0b1afe38`, the candidate less
  the gate on the eight-lane kernel, which is inert on this class): four
  rounds of `M1`'s rows at the three target sizes. It gave the prediction
  above.
