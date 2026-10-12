# Coefficient-aware reuse works; cold performance promotion is rejected

The parameterized support-envelope contract is implemented and validated on
bounded generic Boolean matrices. It produces new matrices from changed
coefficient vectors while correctly handling product cancellations, generator
degree drops, newly required multipliers, constant and zero generators, and
support escapes. This answers the correctness question left open by the
preceding exact-support experiment.

The eager envelope representation is **not promoted for performance**. It
passes **0 of 48** preregistered holdout gates. Every gate's median control /
envelope cold-time ratio is below one (range 0.5801–0.9322). No production
solver integration follows from these measurements.

## Retained evidence

- `run_01/protocol.json`: frozen contract, new holdout seeds and 48 gates.
- `run_01/worker.rs`: actual compiled worker, identical to the top-level source.
- `run_01/results.json` and `RESULT.md`: all 256 cells and derived comparisons.
- `run_01/raw-n*.jsonl`: fixture assignments and exact raw output for every arm.
- `run_01/receipts.json`: fresh-worker exit status, time and peak RSS.
- `run_01/manifest.json`: 15 source/result/receipt file hashes.
- `run_01/test_receipt.json`: seven passing Rust tests, including 8,192 exhaustive
  coefficient/degree/active-mask cases. Ten Python integrity tests replay the
  complete evidence and reject missing data, altered sources, false hit/fallback
  claims, suppressed degree-drop multipliers and silently expanded envelopes.

The fixed run contains **12,800 batch-arm samples** and **272,000 oracle-verified
matrix outputs**, with **35,520 changed-input envelope hits** and **3,360 exact
support-escape fallbacks**. These hit counts include repetitions. Before arm and
repetition multiplication, the complete fixture grid requires **60,816 additional
multipliers** due to generator degree drops. Exact-support controls have zero
changed-input hits. The grid completed without censored or failed workers.

## Charged construction comparison

The following numbers are copied from the frozen `run_01/RESULT.md`. Times are
cold batch milliseconds, including validation, for batch size 64, pooled across
four variable sizes, two holdout seeds and ten repetitions. Ratios are descriptive
ratios of pooled medians; acceptance uses paired per-size intervals in
`results.json`. Every row is a construction-stage engineering diagnostic with
correctness PASS, not a full-method speed claim.

| Variant | Repeat (ms) | Changing coefficients (ms) | Degree cycle (ms) | Support escape (ms) | Direct / arm, coefficients | Matrix cache / arm, coefficients |
|---|---:|---:|---:|---:|---:|---:|
| Direct construction | 1.657479 | 1.594291 | 1.825187 | 1.567375 | 1.000 | 0.993 |
| Verified layout reuse | 1.399625 | 1.659000 | 2.159729 | 1.562709 | 0.961 | 0.954 |
| Exact product schedule | 0.215458 | 1.631145 | 1.816083 | 1.550729 | 0.977 | 0.970 |
| Exact packed-matrix cache | 0.165813 | 1.582396 | 1.739000 | 1.615521 | 1.008 | 1.000 |
| Parameterized support envelope | 1.610500 | 2.081271 | 2.184021 | 2.195896 | 0.766 | 0.760 |

The envelope's median retained storage over this holdout set is **100,715 bytes**,
compared with **4,878 bytes** for the packed-matrix cache. At n=12 the average of
the two seed medians is **319,016 bytes** versus **11,328 bytes**. Storage excludes
allocator metadata. Whole-worker peak RSS ranged from 2,293,760 to 4,505,600 bytes;
it includes all variants and the common reference matrices, so it cannot be
attributed to the candidate alone.

Compilation cost is not the sole problem: on the coefficient-changing family
at batch 64, the envelope's median application cost exceeds direct construction
on **each of the eight held-out size/seed cells**. Multiplying the batch size to
amortize compilation therefore has no measured route to a win for this
implementation on those cells. This observation is finite, not an impossibility
theorem about support envelopes.

## Interpretation and remaining degree of freedom

The predecessor's exact key fixed every coefficient, so it could only reuse an
identical matrix. This candidate instead retains linear maps from coefficient
bits to product coefficients, and genuinely reuses them across changing systems.
That resolves the mathematical contract issue. It does not make the current
allocation and compaction strategy efficient.

One concrete code-level hypothesis remains untested: this eager implementation
visits every planned multiplier and filters by degree, allocates surviving rows,
then sorts and compacts output support again. Bucketing the plans by multiplier
degree would avoid visiting ineligible plans on quadratic inputs. Reusing a
predeclared column index with an occupancy map could avoid repeated sorting,
while still compacting the actual support and honoring current caps. These are
**unmeasured follow-up hypotheses**, and would need a new fixed experiment with
new holdouts and the same direct and exact-matrix controls. The current run is
closed as a negative performance result; its holdouts must not become tuning data
for another claim against those same holdouts.

Full-pipeline cost, operation calibration and rho ratio remain null. The
standalone worker does not invoke WDSat or Gröbner solving, so the full IC suite
does not apply. No asymptotic or cryptanalytic conclusion is drawn.
