# Compiled schedules show limited 2x gains; the universal objective is unmet

Every candidate fails the preregistered universal dramatic and incremental gates.
The strongest family is compiled full-word scanning: the 64-point kernel and its
fixed dispatcher each pass **3/18** >2x comparisons. These include two fresh n24
holdout groups. This is a bounded engineering result, not the required universal
gain or a cryptanalytic breakthrough. The primary is ineligible for confirmation.

The complete comparison covers **192 distinct systems, 49 methods and 75,264
observations**. Every result is verified: 147 SAT and 45 UNSAT fixtures per arm.
All 168 predecessor inputs and all 32 prior methods remain, with 24 new holdouts.
No prior figure, exception or frozen verdict is superseded by regrouping.

## Fixed acceptance results

The reference has 39 methods, including quiet, full-word, scalar-projection and
single-stage controls. A candidate removes only itself when it belongs to that
roster. All eighteen group comparisons and every completion are required. The
matched-control threshold is a separate 1.05 lower bound; it cannot replace the
strongest-reference gate. Every universal dramatic and incremental result is
REJECTED, even when some individual groups pass.

| Candidate | >2x strongest-reference groups | >1x strongest-reference groups | >1.05x matched-control groups |
|---|---:|---:|---:|
| `wide64_quiet` | 0/18 | 8/18 | 13/18 |
| `wide64_unrolled` | 3/18 | 13/18 | 13/18 |
| `word16_unrolled` | 0/18 | 12/18 | 16/18 |
| `word_dispatch` | 3/18 | 14/18 | 17/18 |
| `leaf16_word_unrolled` | 0/18 | 0/18 | 17/18 |
| `byte_simd` | 0/18 | 0/18 | 18/18 |
| `byte_planes` | 0/18 | 0/18 | 18/18 |
| `byte_unrolled` | 0/18 | 5/18 | 18/18 |
| `leaf16_byte_simd` | 0/18 | 0/18 | 18/18 |
| `leaf16_byte_planes` | 0/18 | 0/18 | 18/18 |
| `leaf16_byte_unrolled` | 0/18 | 0/18 | 18/18 |

The complete >2x subgroup results are retained below; none is presented as a pass
of the full protocol:

| Candidate | Split | Variables | Family | Paired median ratio | 95% interval |
|---|---|---:|---|---:|---|
| `wide64_unrolled` | regression | 24 | unplanted | 2.121658 | [2.111656, 2.135872] |
| `wide64_unrolled` | holdout | 24 | planted | 2.101325 | [2.083576, 2.131333] |
| `wide64_unrolled` | holdout | 24 | unplanted | 2.147599 | [2.072068, 2.231157] |
| `word_dispatch` | regression | 24 | unplanted | 2.126564 | [2.107128, 2.145454] |
| `word_dispatch` | holdout | 24 | planted | 2.109842 | [2.083912, 2.206426] |
| `word_dispatch` | holdout | 24 | unplanted | 2.146314 | [2.102378, 2.209661] |

The dispatcher also misses four incremental groups: regression n16 planted and
n20 cross-planted, and holdout n16 planted and cross-planted. Their lower bounds
are respectively 1.000000, 0.999973, 0.925926 and 0.779613. The previous n24
cross-planted seed2097153 exception remains at ratio **0.194210** for the dispatcher
and **0.093432** for direct SIMD byte filtering against the current reference.

## Projection work and measurement order

Across one deterministic solve per fixture, two-stage direct SIMD screening
processes **606,446,080 partial points**, makes **2,368,803 second-stage checks**,
and requires **9,516 complete syndrome checks**. Its single-stage control makes
2,368,803 complete checks on the same ordered points. This large work-count
reduction does not establish a complete-cost gain: setup, coefficient maintenance,
mask extraction and all failed filtering remain charged. These are single-grid
counts, not all eight timing repetitions and not collected independent relations.

The randomized/reverse protocol addresses timing differences noticed between
equivalent wrappers during cyclic discovery. Under the primary protocol, the
median component/dispatcher ratios across fixture medians are **1.002699** at
n16, **0.998049** at n20 and **1.001462** at n24. Their mathematical work and first
models match exactly. These results are consistent with measurement-context or
code-layout effects in the discovery differences; they do not prove a specific
hardware cause. No instruction-cache-cold claim is made.

## Cold complete standalone solve costs

Units are milliseconds including encoding, every image and schedule, filtering,
complete checks, recovery, destruction and result validation. The table pools
two fresh n24 holdout seeds and eight repetitions. Its ratio is a descriptive
ratio of pooled planted costs relative to retained direct SIMD. Acceptance uses
paired pointwise fastest-reference costs, a different statistic. Every row is
class **engineering** and has correctness PASS.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Direct SIMD / method, planted | Correctness |
|---|---:|---:|---:|---:|---|
| Search-only | 82.267875 | 59.109792 | 133.809917 | 0.025 | PASS |
| Degree-3 flat kernel | 389.806354 | 236.008417 | 544.432396 | 0.005 | PASS |
| Degree-3 sparse bucket | 1434.410437 | 766.203542 | 1832.271167 | 0.001 | PASS |
| Degree-3 hybrid kernel | 676.465416 | 361.093209 | 753.735458 | 0.003 | PASS |
| Selective degree-2 flat | 93.188875 | 67.080521 | 143.119875 | 0.022 | PASS |
| Selective degree-2 one-word | 86.461167 | 61.242230 | 135.998333 | 0.024 | PASS |
| Ordered-specialization search | 70.492688 | 50.582021 | 114.221813 | 0.029 | PASS |
| Fixed-quadratic state | 29.291563 | 20.800375 | 47.691500 | 0.070 | PASS |
| Retained packed state | 17.407000 | 12.600105 | 28.807396 | 0.117 | PASS |
| Full RREF, list | 204.686375 | 87.435541 | 435.574125 | 0.010 | PASS |
| Full RREF, wide | 42.473230 | 19.855313 | 95.332063 | 0.048 | PASS |
| Echelon/affine tail, list | 148.274250 | 75.484792 | 243.890874 | 0.014 | PASS |
| Echelon/affine tail, wide | 45.624187 | 24.110333 | 76.576750 | 0.045 | PASS |
| Packed without diagnostic hashing | 12.495354 | 8.820979 | 20.444416 | 0.163 | PASS |
| Recursive affine/products, list | 516.991584 | 206.844083 | 623.828562 | 0.004 | PASS |
| Recursive affine/products, compact | 43.627188 | 21.257562 | 64.053354 | 0.047 | PASS |
| Direct Gray, scalar | 17.462333 | 13.337146 | 20.751812 | 0.117 | PASS |
| Direct Gray, SIMD | 2.038874 | 2.003750 | 3.220125 | 1.000 | PASS |
| 12-variable leaves, scalar | 18.364604 | 12.039438 | 28.516708 | 0.111 | PASS |
| 12-variable leaves, SIMD | 9.562000 | 5.834291 | 12.479625 | 0.213 | PASS |
| 16-variable leaves, scalar | 16.459479 | 14.677500 | 21.557334 | 0.124 | PASS |
| 16-variable leaves, SIMD | 3.197521 | 2.606396 | 3.860292 | 0.638 | PASS |
| Transported blocks, scalar | 15.147854 | 12.674250 | 20.539084 | 0.135 | PASS |
| Transported blocks, SIMD | 2.149083 | 1.806938 | 2.875146 | 0.949 | PASS |
| Initial restriction, list/scalar | 15.908042 | 6.974771 | 20.867188 | 0.128 | PASS |
| Initial restriction, packed/SIMD | 2.117771 | 0.926500 | 2.877042 | 0.963 | PASS |
| Transported 16-variable leaves, scalar | 14.241354 | 14.636875 | 21.519687 | 0.143 | PASS |
| Transported 16-variable leaves, SIMD | 3.247459 | 2.570146 | 3.755355 | 0.628 | PASS |
| Fibers with row elimination | 5.986562 | 5.545438 | 9.374667 | 0.341 | PASS |
| Fibers with scalar column membership | 5.037688 | 4.532063 | 8.839833 | 0.405 | PASS |
| Fibers with SIMD screening | 1.966979 | 1.670146 | 3.652000 | 1.037 | PASS |
| Fibers with zero-only SIMD screen | 2.962292 | 2.578375 | 5.166291 | 0.688 | PASS |
| Direct 16-point, quiet | 1.937105 | 1.748355 | 2.667021 | 1.053 | PASS |
| Full-word 64-point, quiet | 1.461688 | 1.135646 | 1.789583 | 1.395 | PASS |
| Projected bytes, scalar | 5.921750 | 4.971187 | 7.801771 | 0.344 | PASS |
| Projected bytes, SIMD | 1.238875 | 1.104625 | 1.951437 | 1.646 | PASS |
| Leaf16 full words, quiet | 2.179438 | 2.294042 | 3.499896 | 0.936 | PASS |
| Leaf16 projected bytes, scalar | 5.501542 | 5.617355 | 8.762771 | 0.371 | PASS |
| Leaf16 projected bytes, SIMD | 1.786313 | 1.448959 | 2.682500 | 1.141 | PASS |
| Single-stage bytes, quiet | 1.779834 | 1.332938 | 2.341958 | 1.146 | PASS |
| Leaf16 single-stage bytes, quiet | 2.332917 | 1.737416 | 3.146209 | 0.874 | PASS |
| Projected bit planes | 1.536125 | 1.161938 | 2.018687 | 1.327 | PASS |
| Leaf16 projected bit planes | 2.117979 | 1.531125 | 2.864583 | 0.963 | PASS |
| Compiled projected bytes | 1.140313 | 0.900521 | 1.573042 | 1.788 | PASS |
| Compiled full-word 64-point | 0.698354 | 0.560230 | 0.839333 | 2.920 | PASS |
| Leaf16 compiled projected bytes | 1.467542 | 1.215521 | 2.348187 | 1.389 | PASS |
| Compiled full-word 16-point | 0.984708 | 0.969562 | 1.537875 | 2.071 | PASS |
| Fixed compiled dispatcher | 0.578354 | 0.565666 | 0.838083 | 3.525 | PASS |
| Leaf16 compiled full words | 1.792749 | 1.588083 | 2.321146 | 1.137 | PASS |

[run_01/results.json](run_01/results.json) retains every interval, model, work count
and exclusive phase cost. [SUMMARY.json](SUMMARY.json) adds individual exceptions,
projection totals and same-kernel alias ratios. Confidence intervals concern
these fixed fixtures, not a population-wide or worst-case guarantee.

## Validation and custody

The producer passes **60 Rust tests**. Thirty-two Python evidence tests pass,
including exact frozen replay, original-equation checks, retained inputs and
methods, compiled-schedule equivalence, random/reverse ordering, null traces,
filter conservation, source/counter corruption and censored costs. The mathematical
contract is in [README.md](README.md), with the producer argument in
[CORRECTNESS.md](CORRECTNESS.md). These are not independent external review.

The first local evidence-test attempt could not copy fixtures because the system
temporary volume was full. The test harness now creates scratch copies on the
checkout volume and registers cleanup before copying. The rerun passed; no timed
source or frozen measurement changed. [VALIDATION_ATTEMPTS.json](VALIDATION_ATTEMPTS.json)
retains both observed outcomes and explicitly leaves unavailable execution
timestamps null.

The complete campaign took **1,556.231 seconds**. Whole-worker peak RSS was
**29,851,648 bytes**, including all methods and reference preparation. Candidate-
specific allocation and calibrated operation counts remain unmeasured. No measured
cell is censored, and no confirmation was launched. [RUN_LEDGER.json](RUN_LEDGER.json)
binds the seven discovery bundles, completed primary, source lineage, report and
validation attempts. Executed artifacts remain immutable. Production, full
index-calculus and rho costs remain **null**. The broad dramatic-gain goal is open.

## Next bounded question

The most useful next question is whether direct construction of the compiled
full-word schedule can materially lower complete cost. The current timing fields
do not separate schedule construction from scanning. A discovery-only cost split
must first establish whether a construction-only change could possibly meet the
strongest-reference target; a counting argument alone is insufficient.
[NEXT_EXPERIMENT.md](NEXT_EXPERIMENT.md) records that prospective cost-bound gate
and the coefficient-cancellation contract. No implementation or gain is claimed.
