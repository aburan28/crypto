# Combined construction and reduction reject both new candidates

Charging complete canonical linear reduction changes the construction-only
conclusion. At 36 variables the strongest retained staged path is about
**1.38–1.50x** faster than sorted construction plus the same reducer on the
successor holdouts. Neither immediate streaming nor cached nonzero-word and
pivot-incidence reduction improves on that staged reference. Both strict
dramatic-combined-workload gates pass **0/9** comparisons.

These are complete measurements of the bounded **construction-plus-RREF** task.
They are not measurements of complete polynomial solving or index calculus.
The broader thread objective remains open.

## Experiment 1: streaming in ambient coordinates

`run_01` retains sorted, ranked, sparse-intermediate and small-system dense
construction controls, each followed by the same pivot-indexed reducer. The
streaming candidate generates and eliminates rows immediately in ambient
coordinates, then compacts the canonical result. It counts every nonzero source
row before elimination, including dependent rows.

All 53,040 reduced outputs match an independent construction and Gauss-Jordan
oracle. Streaming loses to the pointwise fastest staged path in every larger
size/family comparison: paired reference/candidate medians are 0.533–0.772x,
or approximately 1.30–1.88x slower. At n=36 restricted-cycle, it executes about
12 times as many word XORs while preserving the same logical row eliminations.
Removing the intermediate matrix did not compensate for wider elimination.

The original sorted pipeline's n=36 linear-drop construction fraction is
about 45.7–45.9%. Setting just that timed phase to zero, with every other phase
held fixed, projects only 1.84–1.85x total improvement. This is a conditional
diagnostic, not a bound on coupled implementations or complete solving.

## Experiment 2: nonzero words and backward incidence

`run_02` uses fresh holdouts and keeps every previous pipeline as a fresh
same-binary control. The new arm retains compact coordinates, caches nonzero
word positions in immutable forward pivot rows, and computes backward target
lists from the echelon pivot incidence. Those lists remain valid because
higher-pivot elimination cannot change lower-pivot columns. Each pivot word
list is refreshed after all higher pivots have been removed.

All 98,280 outputs and logical row-XOR counts match. Nevertheless, the new
candidate loses all nine comparisons: paired fastest-reference/candidate
medians are 0.585–0.869x, or approximately 1.15–1.71x slower. Fewer word XORs
are not fewer total operations or lower wall time; cache/list preparation,
indirection and cleanup are also charged. The run does not isolate the
individual overhead contributions.

At n=36, the new reducer performs about 40% fewer word XORs on quadratic
fixtures, yet reduction time rises from 0.533188 to 1.172268 ms per batch.
For linear-drop fixtures it performs about 57% fewer word XORs, while reduction
time rises from 1.185915 to 1.456624 ms. Neither counter reduction is promoted
as a speedup.

## Matched cold costs at 36 Boolean variables

The following milliseconds are copied from `run_02/RESULT.md`. Each cold batch
contains eight inputs and includes setup, construction, complete forward and
backward elimination, compaction, exact output validation and destruction.
Values are medians over two holdout seeds and thirty balanced repetitions.
Ratios in the table use pooled medians for the quadratic column; acceptance
uses paired per-size/family ratios against the fastest prior pipeline.
Every measured row has exact oracle equality PASS and class engineering.

| Variant | Quadratic (ms) | Linear drop (ms) | Restricted cycle (ms) | Sorted / arm, quadratic | Correctness |
|---|---:|---:|---:|---:|---|
| Sorted construction + common reducer | 1.489562 | 2.441646 | 0.272854 | 1.000 | PASS |
| Ranked construction + common reducer | 1.168688 | 2.074896 | 0.257354 | 1.275 | PASS |
| Sparse construction + common reducer | 0.996125 | 1.771021 | 0.183812 | 1.495 | PASS |
| Ambient-coordinate streaming | 1.912479 | 2.882667 | 0.343083 | 0.779 | PASS |
| Sparse construction + incidence reducer | 1.634313 | 2.035062 | 0.256021 | 0.911 | PASS |
| Dense lookup construction + common reducer | NOT_EXECUTED | NOT_EXECUTED | NOT_EXECUTED | null | not measured |

Dense lookup controls execute only at the n=12 bridge. Their larger-size costs
remain null. Streaming's separate construction/reduction times also remain
null because only its fused phase is observable; the fused cost is charged.

## Verification and memory accounting

The two complete grids contain 288 cells, 34,920 batch-arm samples and
**151,320 oracle-verified RREF outputs**. No measured cell failed or was censored.
Both frozen manifests retain source/protocol and executable hashes, fixtures,
raw samples, process receipts and test output. Previous experiment files are
unchanged. The successor's source hash binds its predecessor without implying
that the two implementations are the same.

Eight Rust tests cover all 5,050 small binary matrices with explicit row-span
enumeration, all 6,144 small Boolean coefficient/degree/mask combinations,
larger family fixtures, dependent source-row caps, current column caps,
word boundaries, context changes and 400 systems over all eight Boolean points.
Fifteen evidence-integrity tests replay both runs and reject false phase zeros,
missing reductions, changed logical XOR counts/ranks, hidden incidence storage,
altered sources, incomplete runs and fabricated dense controls.

`basis_max_bytes` is maximum end-of-insertion basis storage over the batch,
not full peak memory. The incidence candidate separately reports its prepared
pivot-mask/list capacities, excluding later cache growth and other temporaries.
The run_02 n=36 quadratic incidence bases use 170,740–186,716 bytes across the
two seeds, plus 26,032–26,760 bytes of incidence storage. The common reducer's
bases use 115,984–123,988 bytes. Whole-worker peak RSS reached 10,108,928 bytes
in run_01 and 10,600,448 bytes in run_02, including all arms and references.

## Source-contract audit and the next comparison

`SOURCE_CONTRACT.json` records a read-only audit at the experiment's base commit.
The repository has **two different output contracts**: the public counted
matrix routine returns full reduced rows, while the specialized solver consumer
needs contradiction/forced-variable consequences. At that base, the latter
already selects a linear/constant row-space intersection and flat matrix
storage by default at 24 or more variables, subject to its explicit overrides.
The full-RREF experiment must not be substituted for that specialized path.

The next solver-relevant mathematical comparison must use the actual requested
output contract and retain the existing linear-tail/flat-storage method as a
strong control. Claiming a gain merely by omitting work that consumer already
omits would be invalid. No production solver, curve input or scalar workflow
was run or changed here. Complete polynomial-solving cost, full IC cost and
rho ratio remain null.
