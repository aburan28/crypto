# Results: each summation-polynomial path pursued

**Retained:** borrowed symbolic expansions plus once-per-curve fixed-slot grouping
and Horner specialization, with corrected operation accounting. The independent
fixed-accumulation candidate is validated and preserved, but its slower per-target
path is superseded by grouped Horner. The prior row/cache allocation optimization
remains in place.

The three local results relative to the corrected reference are:

* Borrowing cached expansions lowers symbolic setup time by about 1.0%.
* Fixed-slot accumulation lowers specialization time by about 49.4% without
  changing the field-operation count.
* Grouped Horner lowers specialization time by about 55.4%, and lowers actual
  specialization multiplications. The combined setup/Horner candidate retains
  those evaluation gains and lowers setup time by about 1.2%.

These are **local CPU measurements**, not a new end-to-end ECDLP speedup.
The independent full-DLP confirmation below establishes neither an improvement
nor a regression. Its CPU-time control explains why the apparent initial
full-DLP gains should not be promoted to an overall speedup claim.

## Original frozen six-way ablation

One unit in every numerical cell: **candidate time / corrected-reference time**,
with input-cluster bootstrap 95% intervals in brackets. All variants perform the
same verified workload within a column. Micro timings separate symbolic setup
from specialization; residual and full-DLP columns include entire child processes.
The reference is explicitly 1.0. Absolute times and every pair remain in
[the comparison](run-01/comparison.json).

| Variant | Symbolic setup | Specialization | Residual corpus | Cold full-DLP corpus |
|---|---:|---:|---:|---:|
| Legacy accounting | 0.999876 [0.997804, 1.001606] | 1.001221 [0.994232, 1.007177] | 0.996691 [0.990552, 1.000328] | 0.962602 [0.928205, 0.993459] |
| Corrected reference | 1.000000 [1.000000, 1.000000] | 1.000000 [1.000000, 1.000000] | 1.000000 [1.000000, 1.000000] | 1.000000 [1.000000, 1.000000] |
| Borrowed setup | 0.989942 [0.987686, 0.991640] | 1.013252 [1.006896, 1.018728] | 0.982021 [0.949233, 0.998981] | 0.937438 [0.890032, 0.987369] |
| Fixed slots | 1.003492 [1.001240, 1.005884] | 0.506250 [0.504315, 0.508242] | 1.037073 [0.997631, 1.085220] | 0.935409 [0.890692, 0.982371] |
| Grouped Horner | 1.002204 [0.999964, 1.004330] | 0.446125 [0.442805, 0.449634] | 1.021846 [0.999810, 1.066536] | 0.936030 [0.889020, 0.985526] |
| Combined (retained) | 0.988336 [0.986704, 0.990133] | 0.446360 [0.443583, 0.448725] | 0.998004 [0.997311, 0.998958] | 0.963040 [0.924996, 0.998961] |

The original 342 runs completed with no failures. This comprises 144 micro runs
(589,824 timed specializations), 108 stage runs (6,480 verified residuals), and
90 verified cold full-DLP runs (16,200 independent oracle cross-checks).
Every specialized-polynomial digest agrees, as do complete decomposition sets,
relation statistics, rho results and all counters apart from exactly reconciled
specialization charges. All 21 prior completed stage/full-DLP outputs replay
exactly in the legacy control, ignoring only timing fields. No outliers are
removed. See [raw records](run-01/raw.jsonl) and [provenance](run-01/provenance.json).

## Why the initial full-DLP timing is not the verdict

The accounting-only legacy control also appeared faster than the corrected
reference, even though it executes the same polynomial arithmetic. Individual
paired runs showed descheduling outliers: for example, seed 11 repetition 0 took
0.6083 seconds for reference and 0.4395 for combined, whereas the other two pairs
were about 0.427–0.429 seconds each. This is an unresolved timing confound, not
proof that algebra became substantially cheaper.

Before selecting a default, [a second contract](confirmation-contract.json)
froze five repetitions on the five previous curves plus two fresh curves,
seeds 907 and 1009. It reran legacy, reference, Horner and combined and retained
whole-process wall time **and child CPU time separately**. All 140 runs recovered
the correct scalar and passed 23,940 independent oracle checks; exact outputs
and specialization-cost deltas agree. The original evidence is retained unchanged.

One unit below: time / corrected-reference time, with input-cluster 95% intervals.

| Variant | Whole-process wall time | Child user+system CPU time |
|---|---:|---:|
| Legacy accounting | 0.985827 [0.936350, 1.048988] | 0.999646 [0.997372, 1.002353] |
| Corrected reference | 1.000000 [1.000000, 1.000000] | 1.000000 [1.000000, 1.000000] |
| Grouped Horner | 0.985358 [0.930106, 1.041011] | 1.001691 [0.998953, 1.004212] |
| Combined (retained) | 0.962640 [0.905238, 1.021824] | 1.000787 [0.998405, 1.003021] |

All non-reference confirmation intervals include 1. The retained combined
candidate's CPU-time ratio is 1.000787; this supports
neither a full-DLP speedup nor a regression. The large, repeatable local
specialization benefit is retained, with no claim that it materially accelerates
complete scalar recovery. See [confirmation raw records](confirmation-01/raw.jsonl),
[comparison](confirmation-01/comparison.json), and [provenance](confirmation-01/provenance.json).

## Accounting and limits

At p=67, seed 11, specialization actually costs 1342 Fp multiplications per
reference call and 1199 with Horner: 10.66% fewer products in this
phase. Symbolic setup uses 41,360 products in both variants. The historical
estimate was 1,770 per specialization, which overcharged the reference by 428.
That correction is separate from the additional 143 products saved by Horner.
The runner checks the exact per-curve delta using the measured microbenchmark
charges and number of solver calls, for every stage and full-DLP report.

This small phase's reduction is diluted by Macaulay elimination, normal forms,
root extraction and relation collection. Full common-operation totals, normalized
S, cost/rho and cost/floor remain unmeasured. Existing reported solver totals are
partial metrics; they include the corrected measured specialization charge, but
do not newly calibrate cold setup and runtime bookkeeping. Classification is
**engineering**, plus an explicitly separated **accounting correction**. The
summation-polynomial degree, relation-yield ceiling and asymptotic boundary do
not change. No extrapolated crossover or GPU/RDMA improvement is claimed.

The validation corpus is small and the host is shared. Confidence intervals are
clustered by curve/seed, with repetitions kept together. The confirmation has
seven input clusters, not 35 independent curves. Full-DLP times include independent
MITM and rho auditing. These limits preclude a general performance guarantee.

## Implementation and reproducibility

All 15 Gaudry unit tests pass for the retained implementation, including direct
termwise specialization checks over every element of F_(7³), sampled larger-field
inputs, zero/base-field targets, symbolic S4 identities, border retries, independent
MITM equivalence and dense/sparse/large-prime scalar recovery. [Test log](gaudry-tests.log).
All source variants and their hashes are preserved under [sources/](sources).
[selection.json](selection.json) records the default and its source hash.

The current shared checkout also passes `cargo check --lib`. Build commands,
workload definitions and timing boundaries are in [README.md](README.md).
`reference.patch` records the accounting change; the remaining patches compare
each optimization to that corrected reference. The standalone fixed-slot path
remains reproducible even though grouped Horner is the selected implementation.
