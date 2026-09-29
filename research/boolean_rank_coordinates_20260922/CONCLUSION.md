# Exponential lookup removed; larger-system gains remain workload-dependent

The coordinate index now uses combinatorial prefix counts rather than a table
with 2^n entries. It is verified through **36 Boolean variables**. At n=36,
degree 3, the rank table has **148 entries** (1,184 bytes of usize payload on
this host). The ambient monomial basis still has 7,807 columns and is charged.
This is a structural memory improvement, not evidence of a full-solver speedup.

Two bounded experiments were completed. Neither satisfies its universal
dramatic-construction gate. The broader complete-algorithm objective remains
open; the earlier 6–12-variable gains cannot simply be extrapolated.

## Initial scaling test

`run_01` compares sorted construction, packed construction with binary search,
and packed construction with combinatorial rank. A dense lookup runs only at
n=12. Ranked indexing beats binary search in all nine larger-size holdout
comparisons, by paired medians **1.20–1.28x**. But only **11/18** total gates
pass: just two of the nine sorted/ranked comparisons establish the required
2x lower bound. At n=36, paired ratios against sorted are 1.606x for quadratic,
1.490x for linear-drop and 1.219x for restricted-cycle fixtures.

This run has 144 cells, 5,616 batch-arm samples and **76,752** independently
oracle-checked outputs. The worker, protocol and all raw data are retained.

## Sparse-intermediate successor

`run_02` tests a new mechanism on new holdouts: one scratch bitmap per
application, with a touched-word list. Each intermediate row retains only
nonzero word blocks. Scratch/marker state resets even after cancellation, and
the final dense output is allocated only at its actual compact width. All
scratch work, sparse allocations, output construction and validation are charged.

The successor keeps fresh sorted, binary and dense-intermediate ranked arms.
It passes **22/27** gates, so universal promotion is still **REJECTED**. The
five failures are all three comparisons against ranked at n=20, the sorted
comparison for n=20 restricted-cycle (lower bound 1.992x), and the sorted
comparison for n=36 restricted-cycle (lower bound 1.767x).

At n=36 the sparse arm is a useful finite improvement over ranked in every
family. Paired median ratios and 95% bootstrap intervals are:

| Family | Sorted / sparse [interval] | Ranked / sparse [interval] |
|---|---:|---:|
| Quadratic | 2.159 [2.135, 2.182] | 1.347 [1.330, 1.363] |
| Linear drop | 2.206 [2.184, 2.219] | 1.448 [1.439, 1.464] |
| Restricted cycle | 1.807 [1.767, 1.840] | 1.463 [1.416, 1.491] |

Cold batch32 timings below are copied from the corrected successor report.
They are milliseconds, including setup and exact output validation; medians
pool two holdout seeds and **twenty** balanced repetitions. The ratio column is
a descriptive ratio of pooled medians for the quadratic family. Acceptance
uses paired per-size/family intervals, not that pooled ratio. Every measured
row has oracle equality PASS and is classified as engineering.

| Variant at n=36 | Quadratic (ms) | Linear drop (ms) | Restricted cycle (ms) | Sorted / arm, quadratic | Correctness |
|---|---:|---:|---:|---:|---|
| Sorted | 3.581979 | 4.658271 | 0.849187 | 1.000 | PASS |
| Binary lookup, dense intermediate | 2.868584 | 3.937604 | 0.856833 | 1.249 | PASS |
| Ranked lookup, dense intermediate | 2.232979 | 3.075438 | 0.688917 | 1.604 | PASS |
| Ranked lookup, sparse intermediate | 1.666313 | 2.124417 | 0.469021 | 2.150 | PASS |
| Exponential dense lookup | NOT_EXECUTED | NOT_EXECUTED | NOT_EXECUTED | null | not measured |

Sparse intermediates regress at n=20 compared with ranked dense intermediates:
paired ranked/sparse medians are 0.888x, 0.935x and 0.957x. No post-holdout
hybrid or selected-family success claim is made. A future representation policy
would require new holdouts and its own charged comparison.

The successor has 144 cells, 12,240 batch-arm samples and **167,280** checked
outputs. Together both runs verify **244,032** matrix outputs. No measured
cell failed or was censored. Dense-lookup cells above n=12 are deliberately
NOT_EXECUTED with null costs, not measured failures or hypothetical timings.

## Validation, memory and accounting correction

The initial seven Rust tests and successor eight tests pass. They include all
6,144 small coefficient/degree/active-mask cases across the applicable arms,
independent rank enumeration through n=36, high variable bits, exact row order,
context changes, current caps, and sparse scratch reset after word cancellation.
A full-width constant generator at n=36 exceeds the shared 4,096-row cap and
has an explicit resource-refusal regression. It is not mistaken for a zero
matrix or a mathematical refutation.

Twelve evidence-integrity tests verify and replay both runs, including all 30
manifest file hashes and the predecessor linkage. At n=36 the ranked/sparse
contexts retain 141,160 bytes with full multiplier masks, and 68,584 bytes with
the restricted mask. Persistent context size does not measure the intermediate
allocation saving: those rows and scratch buffers are temporary. Whole-worker
peak RSS reached 18,120,704 bytes in both runs and includes the common oracle
corpus and every variant. It is not candidate-specific peak memory.

**Additive legend correction:** the frozen `run_02/RESULT.md` inherited the
phrase “twelve balanced repetitions.” Its actual protocol, worker commands,
sample count and analysis used **twenty**. `REPORT_CORRECTION.json` binds the
original, corrected report, authoritative protocol and unchanged numerical
results by SHA-256. `run_02_REPORT_CORRECTED.md` changes only that legend.
The original frozen files remain intact. The current analyzer formats the
repetition count from the protocol; replay tests require identical numerical
results and exactly the corrected report. No measurement or gate outcome changes.

## Scope and next unresolved work

All work remains standalone generic Boolean matrix construction. There is no
curve input, key recovery, production solver integration, full-solving cost or
rho comparison. The WDSat/full IC suite is inapplicable; missing full-pipeline
quantities remain null. Bootstrap intervals describe fixed seeds on this host,
not independent reproduction or a population claim.

The exponential lookup obstacle is resolved in this bounded constructor.
The next unresolved costs are the remaining ambient-basis/context work on
restricted inputs, the density-dependent intermediate representation, and the
fraction of complete algebraic solving spent on construction. Those require
measurement before claiming a dramatic complete-algorithm improvement.
