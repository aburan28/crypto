# Point counting: measured changes and formal scope

2026-10-16; PR #1923. This report supplements the historical
[implementation study](../sea_point_count_20261015/README.md), which made no
timing claim. The present study measures the original source cf274521 against
the optimized source 44f8d17f, using identical native workers and verified
orders. All 140 requested measurement records completed; none was excluded
by the visible-process contention policy.

## Observed performance

The three BSGS point-count workloads have paired median baseline/candidate
ratios 1.98652, 2.79019 and 3.24310. Each exceeds its workload's observed A/A
deviation envelope. SEA and the complete P-192 API do not establish an
improvement under that predeclared noise policy. Cached SEA at 40 bits was
slower in every A/B pair; its median ratio is 0.81228, within a particularly
wide 49.308% A/A envelope. That observation is retained rather than hidden.

The following table is copied from the [native summary](evidence/run-20261016T0205Z/summary.md).
Times are measured microseconds per call, with five A/B pairs per workload.
Ratios are medians of paired ratios, not ratios of the two displayed medians.

| Workload | A median us/call | B median us/call | A/B median ratio | Paired range | A/A max deviation | A/B pairs |
| --- | ---: | ---: | ---: | --- | ---: | ---: |
| bsgs24 | 61.747 | 31.495 | 1.98652 | 1.88259..2.14328 | 26.575% | 5 / 5 |
| bsgs32 | 327.251 | 117.994 | 2.79019 | 2.73356..3.14247 | 17.594% | 5 / 5 |
| bsgs40 | 1936.899 | 591.519 | 3.24310 | 3.07922..3.42764 | 8.710% | 5 / 5 |
| sea40_cached | 264.650 | 314.525 | 0.81228 | 0.80232..0.93396 | 49.308% | 5 / 5 |
| sea192_cached | 1904284.547 | 1825125.454 | 1.02393 | 0.90254..1.10098 | 6.766% | 5 / 5 |
| count192_cold | 98449066.339 | 97670607.289 | 1.00469 | 0.99858..1.01812 | 2.869% | 5 / 5 |
| cm256_control | 19048.878 | 17948.105 | 1.03405 | 0.98103..1.17123 | 19.957% | 5 / 5 |

The host is an AMD EPYC 9V74 VM with five visible CPUs. Children were pinned
to CPU 4, movable threads were evicted, and visible-process CPU use was
monitored under the repository benchmark lock. PSI is unavailable; some idle
system threads could not be moved, and cloud neighbors and frequency were
uncontrolled. This is not an exclusive-host measurement. The A/A envelope is
a descriptive repeatability control, not a confidence interval.

The cached rows exclude canonical Phi-cache import and one warm-up call;
their setup costs are retained separately. The original cache preparation
recorded 82.296 seconds of setup. The cold P-192 row includes fresh polynomial
construction, counting and curve/twist certification. Warm-stage timing is
not substituted for complete API timing. The secp256k1 row is an unchanged
CM point-count control; it does not time CM curve generation.

## Algorithm and correctness changes

The implementation retains every baby index for a repeated group point,
including repeated infinity, uses exact integer sizing and bounds the final
giant block. Later sampled points can test only surviving candidates when
the scalar-test cost estimate is lower than a complete search. Candidate
intersections use ordered results. The small-field count also hoists the
fixed giant-step multiple. These are engineering changes to the existing
algorithms; no asymptotic complexity improvement is asserted.

The TLA+ sample-only negative control prompted a concrete constructor
regression. For CM data p=11, D=-7, v=2, t=-4, the old constructor returned
an order-8 curve when the requested order was 16. For p=59, D=-23, v=2,
t=-12, it returned order 48 instead of 72. Both failed for all eight tested
seeds. It also accepted invalid norm data t=-36 at p=59. The retained
[pre-fix run](evidence/pre-fix-cm-regression.jsonl) has 17 failures among
32 signed-trace/seed cases plus one invalid-input check.

The constructor now checks probable primality and the ordinary CM norm
t^2 - D*v^2 = 4p before constructing a curve. For fields below 2^22 it
compares an exhaustive count; for larger fields it requires successful
curve and twist Hasse certificates. If neither twist is verified, it returns
an error, including when factorization is insufficient. It exposes the
verification evidence on CmCurve. The [post-fix run](evidence/post-fix-cm-regression.jsonl)
passes all 33 checks. A 40-bit integration test exercises the certificate
path for both trace signs, as well as the existing CM isogeny-loop check.

CM generation was already implemented for ordinary, nonexceptional j
invariants through Hilbert class polynomials and twist selection. This patch
repairs its order guarantee; it does not add general support for j=0 or 1728,
prove the analytic polynomial rounding, or certify the full endomorphism ring.
The constructor's GLS metadata h and r are not certified by its base-field
order check. These scope limits are also stated in the API usage guide.
The constructor fix happened after timing and is not included in either
frozen benchmark arm. No generation speedup is claimed.

## Checked mathematical argument

Let the lawful additive group have step P and target Q. For positive baby
width B, every index k has the unique decomposition k=iB+j with 0<=j<B,
by quotient and remainder. The scalar equation [k]P=Q is equivalent to
[j]P=Q-[iB]P by addition and cancellation. A complete baby multimap therefore
finds every solution in 0<=k<C, including repeated group points. Testing k<C
excludes excess entries in the last block. Conversely every emitted match
satisfies the scalar equation. This is Lean's bsgs_complete theorem, with
smul_add and subtract_left proved from the stated group laws.

Within one block, insertion order increases j. Every valid index of an
earlier block is smaller than every index of a later block, since
iB+j < (i+1)B <= i'B+j'. A positive signed trace stride preserves this
ordering. These are giant_blocks_ordered and signed_progression_ordered;
the concrete table's insertion and traversal remain implementation obligations.

If the true candidate starts in C and every sampled test T accepts it,
intersection C with T retains it. Induction proves this for any number of
samples. Directly filtering C by T and intersecting a complete search for T
are the same predicate. A singleton surviving set must consequently contain
the true candidate. These are refinement_preserves, refinements_preserve,
direct_equals_full_intersection and unique_sound. The premise that real SEA
constraints and point tests accept the true trace is explicit, not proved here.

Finally, if two numbers lie in the same interval of width w and both are
multiples of a genuine divisor d>w, they are equal: unequal multiples differ
by at least d, contradicting the width. This is interval_unique, the
uniqueness step of the Hasse certificate. Establishing genuine group-order
divisibility and the arithmetic/primality bridge remains a separate obligation.

## Formal execution and limits

TLC checked Fixed (165,750 distinct states), Wide (310,500), and Certification
(58,956), including invariants, ordering and bounded termination. The legacy
single-index table and sample-only certificate controls both failed with
the expected nontrivial counterexamples and exit status 12. Lean 4.19.0
checked all ten theorems without admitted proofs. Axiom inspection reports
only standard propext and Quot.sound where used. All source, pinned tool
hashes, finite bounds and [final logs](formal/logs/final/) are retained.

These are finite-state model checks and unbounded abstract lemmas, not a
machine-checked refinement of compiled Rust. The [proof map](VERIFICATION.md)
names the exact assumptions and [open implementation obligation](../../problems/sea-rust-group-refinement-proof.md):
prove the actual field/point operations form the group used by Lean, then
prove SEA residue generation and the divisor/primality connection. Analytic
CM construction, integer bounds, compiler behavior and timing resistance are
also outside the checked claims. An independent lean4checker executable is
unavailable; ordinary Lean kernel checking and axiom inspection completed.

## Requirement status and reproducibility

| Requested work | Status | Evidence or exact remaining gap |
| --- | --- | --- |
| Improve point counting and SEA implementation | Implemented; executable checks pass | Exhaustive helper tests, actual curve regressions, release API/SEA/kernel tests; speed claims limited to measured BSGS cases |
| Determine CM generation support and implement if absent | Existing support audited and corrected | Signed-trace regression; exhaustive or Hasse order evidence; exceptional j and analytic certification remain unsupported/unproved |
| Measure speedup | Complete for all seven frozen workloads | 140 records, five A/A and five A/B pairs each, matched orders/work, no excluded pairs |
| Model in TLA+ and check proofs | Checked within documented scope | Three positive TLC configurations, two negative controls, ten Lean lemmas, pinned CI replay |
| Full cryptographic implementation verification | Partial; open | Concrete Rust/group refinement, SEA congruences, certificate divisibility/primality and CM analytic proof are missing |

Release validation passes 41 tests, with one pre-existing ignored test.
The native analyzer also tests rejection of bad/mismatched records. Native
and library Clippy checks, touched-file formatting, and the two CM example
builds are recorded in evidence. The existing crate-wide formatting failure
in lib.rs predates this branch; it is not hidden by editing unrelated code.

Replay commands are in [README](README.md), [measurement accounting](MEASUREMENT.md)
and the [run report](reports/run-20261016T0205Z/REPORT.md). The chart source is
render_reports.py, plotting only the native summary; both SVG and 150 dpi PNG
copies accompany the PDFs under benchmarks/20261016-sea-point-count.
The historical SEA report, performance-gains data, IC scoreboard/timeline and
curve-transfer figures are unchanged: this point-count result does not alter
their DLP workloads or recorded claims. No DLP or curve-security claim is made.
