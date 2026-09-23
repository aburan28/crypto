# Conditional linear fibers are correct, but the complete-cost gain is rejected

The SIMD fiber candidate passes **0/18** dramatic strongest-reference comparisons
and **1/18** incremental comparisons. Its matched scalar-column backend gate passes
17/18 groups and is also rejected as a universal claim. All three complete gates
remain **REJECTED**, and the primary is ineligible for confirmation.

| Comparison | Groups passed | Required | Universal verdict |
|---|---:|---:|---|
| >2x versus the pointwise fastest fixed reference | 0 | 18 | REJECTED |
| >1x versus the same reference | 1 | 18 | REJECTED |
| >1.05x versus the matching scalar-column method | 17 | 18 | REJECTED |

The current reference contains all 28 predecessor methods and the row, scalar-column
and zero-only SIMD controls. There are 32 total methods and 31 reference methods.
The complete run contains **168 distinct systems and 37,632 observations**, with
128 SAT and 40 UNSAT fixtures per arm. All complete with verified results. It retains
all 144 prior inputs and adds 24 unused holdouts. Every retained method preserves
its models, logical counters and traces on the shared inputs.

## What changed, and why the counting reduction was insufficient

An exact maximum independent variable set turns every outside assignment into a
small changing linear system. Scalar and SIMD screens reject only explicit zero-row
contradictions, including fixed XOR combinations of original equation rows. Every
survivor receives exact column-span membership or direct row elimination, followed
by original-variable recovery and full model verification.

Across one deterministic solve per fixture, the SIMD method screens **25,304,736
outside assignments** and makes **105,861 full linear queries**, about **0.4183%**
of screened prefixes. This is a count over one grid, not all seven repeated timing
samples. The screen work and coefficient updates still have to be paid. Rank sums
are ranks of repeated conditional matrices, not independent collected relations.

Paired median ratios against the fastest reference range from **0.403962 to 1.127409**.
The only incremental pass is the n24/unplanted regression group, with a lower bound
of **1.000156**; that marginal fixed-group result does not establish a robust gain.
The old n24/cross-planted/seed2097153 exception remains, with ratio **0.121522** against
the current reference.

The matched scalar/SIMD median ratios range from 1.780589 to 2.942116, but the n24
cross-planted holdout group's interval is **[1.024697, 3.404951]**, below its 1.05
lower-bound threshold. Three systems select six variables and use the documented
scalar-screen fallback, including one fresh cross-planted fixture. Their results
are retained; the implementation was not changed after seeing holdout timings.

## Cold complete standalone solve costs

Units are milliseconds, including encoding, exact selection, all setup and screens,
linear solving, recovery, destruction and result validation. This table pools two
fresh n24 holdout seeds and seven repetitions. Its descriptive last numeric column
uses pooled planted costs relative to retained direct SIMD. Acceptance instead
uses paired pointwise fastest-reference costs on every designated group. Every row
is class **engineering**, with correctness PASS.

| Method | Planted (ms) | Cross-planted (ms) | Unplanted (ms) | Direct SIMD / method, planted | Correctness |
|---|---:|---:|---:|---:|---|
| Search-only | 20.276459 | 26.065313 | 126.155937 | 0.128 | PASS |
| Degree-3 flat kernel | 68.772166 | 50.899687 | 513.035354 | 0.038 | PASS |
| Degree-3 sparse bucket | 235.042229 | 170.269083 | 1760.644084 | 0.011 | PASS |
| Degree-3 hybrid kernel | 129.278416 | 89.139417 | 808.516667 | 0.020 | PASS |
| Selective degree-2 flat | 23.387625 | 29.699292 | 148.162521 | 0.111 | PASS |
| Selective degree-2 one-word | 20.576042 | 26.291375 | 131.962458 | 0.126 | PASS |
| Ordered-specialization search | 17.889646 | 22.693646 | 111.392792 | 0.145 | PASS |
| Fixed-quadratic state | 7.578250 | 9.555167 | 45.521791 | 0.343 | PASS |
| Retained packed state | 4.477958 | 5.729937 | 27.890125 | 0.581 | PASS |
| Full RREF, list | 233.075979 | 69.649104 | 402.975000 | 0.011 | PASS |
| Full RREF, wide | 49.221750 | 15.239688 | 90.753042 | 0.053 | PASS |
| Echelon/affine tail, list | 34.739730 | 23.479625 | 210.481458 | 0.075 | PASS |
| Echelon/affine tail, wide | 10.178250 | 7.162834 | 65.169875 | 0.256 | PASS |
| Packed without diagnostic hashing | 3.074937 | 4.047271 | 19.832270 | 0.846 | PASS |
| Recursive affine/products, list | 172.130958 | 159.895625 | 421.451542 | 0.015 | PASS |
| Recursive affine/products, compact | 18.144458 | 15.882125 | 45.553271 | 0.143 | PASS |
| Direct Gray, scalar | 17.895333 | 16.017187 | 20.737167 | 0.145 | PASS |
| Direct Gray, SIMD | 2.601812 | 2.535916 | 3.105292 | 1.000 | PASS |
| 12-variable leaves, scalar | 3.599833 | 5.628312 | 26.941083 | 0.723 | PASS |
| 12-variable leaves, SIMD | 1.604062 | 2.678062 | 12.073583 | 1.622 | PASS |
| 16-variable leaves, scalar | 2.904375 | 5.419834 | 21.291812 | 0.896 | PASS |
| 16-variable leaves, SIMD | 0.529583 | 0.968480 | 3.802375 | 4.913 | PASS |
| Transported blocks, scalar | 15.881062 | 15.826667 | 20.688625 | 0.164 | PASS |
| Transported blocks, SIMD | 2.081021 | 2.191229 | 2.849396 | 1.250 | PASS |
| Initial restriction, list/scalar | 15.594229 | 7.857625 | 20.597230 | 0.167 | PASS |
| Initial restriction, packed/SIMD | 2.456291 | 1.116708 | 2.844771 | 1.059 | PASS |
| Transported 16-variable leaves, scalar | 2.921396 | 5.488542 | 21.114354 | 0.891 | PASS |
| Transported 16-variable leaves, SIMD | 0.509312 | 0.938187 | 3.716312 | 5.108 | PASS |
| Fibers with row elimination | 2.786709 | 2.241522 | 7.505500 | 0.934 | PASS |
| Fibers with scalar column membership | 2.525792 | 2.043166 | 7.146187 | 1.030 | PASS |
| Fibers with SIMD screening | 0.953291 | 0.743417 | 2.458729 | 2.729 | PASS |
| Fibers with zero-only SIMD screen | 1.453480 | 1.221750 | 3.667479 | 1.790 | PASS |

[run_01/results.json](run_01/results.json) retains every cell, exclusive phase cost,
semantic work count and interval. [SUMMARY.json](SUMMARY.json) adds individual-case
ratios, dimension counts and the aggregate single-grid work. These intervals concern
fixed fixtures, not a population-wide guarantee. Earlier figures and verdicts remain
in their frozen bundles and the canonical scoreboard.

## Validation, accounting and custody

The producer passes **49 Rust tests**, including all **1,157,359** linear systems
through four rows and columns with all right-hand sides. Tests also cover exact
selection, every coefficient and SIMD lane, redundant-row rejection soundness,
changing rank, smallest witnesses, full equation words, complete solves and caps.
Twenty-six Python evidence/census tests pass, including exact replay, original
models, source and counter corruption, retained reference methods, all shared
predecessor behavior, optimal selection and complete UNSAT domain accounting.
These are producer checks, not an independent external audit.

The [contract](README.md) and [correctness argument](CORRECTNESS.md) distinguish
outside prefixes, screen lane-rounds, linear queries and rejected full extensions.
No rejected extension is relabelled as an evaluated full assignment. UNKNOWN remains
censored. Each method appears nine or ten times at each position in each size's
42-fixture cohort; no within-fixture position-balance claim is made.

The campaign took **1,071.082 seconds**. Whole-worker peak RSS was **29,818,880 bytes**
and includes all methods and reference preparation. Candidate-specific allocation
and calibrated operation counts remain unmeasured. No cell was censored, and no
confirmation run was launched. [RUN_LEDGER.json](RUN_LEDGER.json) binds the three
discovery probes, completed primary, support census, source lineage and report.
Executed sources and manifests remain unchanged. Full index-calculus, production
solver and rho costs remain **null**. The dramatic-gain objective remains open.

## Structural follow-up, without a speed claim

The exact discovery support census finds no equation independent of the currently
selected variables on all six n24 systems and on five of six n16 systems. Thus a
simple static untouched-row screen would provide no rejection on the current n24
choices. Alternative sets can trade a smaller eliminated dimension for more
untouched equations; construction and the larger outside domain must be charged.
All twelve censuses finish within the cap, after 32–232 independent subsets.

A separate representation hypothesis is to evaluate a fixed subset of up to eight
equations with byte syndromes and perform full-word checks only for survivors.
That could increase SIMD lane count, but requires exact full-equation recovery,
complete cost accounting and controls for changed batching and diagnostic hashing.
[NEXT_EXPERIMENT.md](NEXT_EXPERIMENT.md) records the prospective contract. It is
unimplemented and unmeasured; no gain or asymptotic conclusion is assigned to it.
