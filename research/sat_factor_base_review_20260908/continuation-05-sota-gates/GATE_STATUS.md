# Koblitz index-calculus SOTA gate status

Current through Stage 135, 2026-09-21. The historical optimization chain is
`stage-99-optimization-chain-20260913/verification.json`. Stage 108 adds the
machine-replayable four-shard and direct-routing chain, five host-identified
routing comparisons, the selected five-pair `n=53` panel, and the refreshed
public unknown-scalar run. Stage 109 composed the current seven-gate audit on
2026-09-13 without changing any prior frozen artifact. Stage 124 re-sealed that
audit over the boundary ledger as updated on 2026-09-21 by the index-calculus
boundary ledger (operation-counted rows for the prime and binary regimes,
oracle pricing and the whole-process count on the Koblitz rows;
`research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md`), chaining to
the Stage 109 seal. Stage 125 re-seals it again over that ledger's Round 2
of the same day (folded pair tables, walk targets, the exact counting
ceiling and a balanced Koblitz base; the note's §10, with the first-round
records kept in the ledger's history), chaining to the Stage 124 and Stage
109 seals. Stage 126 re-seals it once more over that ledger's Round 3, also of
the same day (a walk restart that costs one group operation instead of
thirty-two scalar multiplications, bases sized at the derived family optimum
on the prime and binary ladders, and the family shape law reported per row;
the note's §11), chaining to the Stage 125, Stage 124 and Stage 109 seals. No
Koblitz gate fact moved at any of these re-seals and no prior frozen artifact
changed.

Stage 127 supersedes the current audit without mutating those snapshots. It
hash-pins the frozen 160-input Phase-B matrix and the current scalar-blind
`n=53` workflow evidence, and keeps them explicitly separate: Phase B is the
common `n, ell, m` Semaev solver matrix, while the optimized `n=53` run uses a
subgroup-orbit factor base that has no equivalent WDSat/Magma `ell` encoding.
It also incorporates the current one-core refresh, native projected-predicate
operation counts, every accepted and rejected optimization process, and the
unchanged Magma and unaffiliated-review gaps.

Stage 128 adds a current unified `n=41, ell=6, m=3` single-target series over
the fixed public two-torsion-saturated Frobenius-union base.  It charges the
same process from algebraic materialisation through projected predicate,
relations, verification, sparse linear algebra, descent and automorphism rho,
under five default-thread and five one-worker repetitions.  It replaces the
older online-only reading with a current full-cost loss in both thread classes.

Stage 129 materialises the exact non-Frobenius-closed standard
`n=59, ell=9, m=3` Phase-B factor base in the end-to-end workflow.  A bounded
1.2-million-probe attempt reaches zero natural relations, as the public counting
ceiling predicts, while same-target signed-Frobenius rho completes in five
metered repetitions.  The result is censored at its cap: it is larger-regime
execution evidence, not a completed IC run or an UNSAT certificate.

Stage 130 counts the public cofactor classes of every standard `n=59` width
from `ell=9` through `16`, charges those counters and the width constructions,
and binds the result to the exact compact-table memory formula.  It establishes
a finite no-go frontier for materialized standard bases under the current 4 GiB
budget; other algebraic families, implicit bases and larger storage remain open.

Stage 131 completes one public `n=59` unknown-scalar workflow with the
cofactor-projected image of the standard `ell=14` base.  The recipe uses only
public cofactor multiplication and point identity deduplication; it retains no
scalar preimages or discrete-log labels and does not enumerate the target
subgroup.  Both the default-thread and one-worker runs recover and verify the
same scalar from the same 26,796 relations.  Fully charged IC remains much
slower than automorphism rho in both thread classes.

Stage 132 widens that same algebraic recipe to `ell=15` on the same target,
solver, collection window, sparse filter and rho seed.  It halves probes and
summand scans, reducing default-thread IC wall by 4.81 percent and one-worker
wall by 17.23 percent.  Peak RSS rises by about 4.8 times, and full IC still
loses to rho by three orders of magnitude on one worker.

Stage 133 holds the selected `n=59, ell=15` instance fixed while tuning its
collector.  Equal 102.4-million-summand arms reject windows 512 and 2048 in
favour of 1024.  Two source-pinned buffer-reuse pilots preserve every relation
but regress relation-unit wall, so that implementation is rejected and archived.

Stage 134 adds an explicit witnessed-compact pair table on that same instance.
Each four-byte compact rest carries a four-byte packed pair witness; every
candidate is re-added in the exact fast group, then exact candidates are sorted
before selection.  Both complete modes reproduce Stage 132's relation hash and
reduce full wall and core cost while raising peak memory by about 1.5 times.

Stage 135 caches each first-pass packed pair-sum key in row-major order and
derives its witness from the row and offset.  Scatter then avoids the second
quadratic group-addition pass.  Full time and CPU fall again in both thread
classes, while the temporary 4.34 GB key cache raises construction peak RSS to
about 10 GB.

The campaign has target-independent algebraic factor bases; matched native-XOR,
WDSat, CryptoMiniSat, direct-MITM, GGMP, and signed-Frobenius-rho controls; a
balanced 160-instance PDP panel through `n=59`; public unknown-scalar end-to-end
runs at degrees 23, 31, 41, 53, and 59; and same-target known-answer and
scalar-blind construction, rank, solve, and rho comparisons at `n=53`. It has not passed all seven gates and
does not establish a new state of the art.

| Gate | Status | Current evidence | Remaining requirement |
|:--|:--|:--|:--|
| 1. Charge every stage and resource | **Partial overall** | The corrected Phase-B matrix charges 21,038.596137 core-seconds. The current `n=53` campaign retains 318 measured processes / 1,395.640868 sequential wall-seconds / 3,192.181904 core-seconds. Stages 132--134 retain every selected and rejected `n=59` process. Stage 135 adds three cached-build processes: 743.666397 wall-seconds, 1,576.843696 core-seconds and 10,073,538,560 B maximum RSS. | Licensed Magma process resources are absent. Preinstalled OS/toolchain acquisition remains an explicit exclusion. |
| 2. WDSat, CryptoMiniSat, Magma F4, MITM, GGMP | **Partial** | Native XOR SAT, WDSat, CryptoMiniSat, and direct MITM ran on the exact 160-input packet. Standard `n=31`/`n=41` and GGMP `n=31` are represented. The Stage-32 successor removes the two original WDSat buffer errors without rewriting Stage 26. | Execute all 160 Stage-22 Magma inputs on a licensed host under the frozen one-thread/no-retry contract, seal the return before truth scoring, and report F4 resources. |
| 3. Single-core, core-seconds, memory, conflicts, wall | **Partial because Magma is absent** | Every executed Phase-B arm retains wall, core-seconds, peak RSS, conflicts or operations, tree memory, and workflow wall. The selected cached `n=59` run reports default-thread IC 76.342293 s / 774.241342 core-seconds / 10,072,309,760 B RSS and one-worker IC 648.780325 s / 635.617997 core-seconds / 10,044,899,328 B RSS. | Supply the same fields for licensed Magma F4. SAT conflicts remain inapplicable to exact pair-table arms and are reported as null. |
| 4. Scale through `n=31`, `n=41`, and a larger PDP regime | **Satisfied for finite execution coverage including one completed larger IC run** | Phase B covers `n=31`, GGMP `n=31`, `n=41`, and `n=59`. Stages 129--130 retain the standard `n=59` cap and exact width frontier. Stages 131--132 complete and optimize the cofactor-projected `n=59, m=3` public unknown-scalar workflow at `ell=14` and `15`. | The evidence is finite and toy-sized; it is not an asymptotic scaling law or a literature-scale speed record. |
| 5. Unknown scalar with no constructed factor-base logs | **Satisfied for finite degrees 23, 31, 41, 53, and 59** | Stage 108 archives the `n=53` public hash-seed-53001 run. Stage 132 retains public hash-seed-59001 at `n=59`: both thread modes derive all 16,344 logs from the same 54,749 relations, recover `d=17861472351607`, and verify `[d]G=Q`. Neither target scalar nor factor-base logs are supplied. | Repeat on independent public seeds and obtain unaffiliated replay; these strengthen rather than replace the finite gate-5 execution. |
| 6. Full cost against automorphism-optimized Pollard rho | **Failed for current `n=59`; `n=41` loss; `n=53` default-thread wall pass only** | Stage 135 retains the same 2.6M probes / 2.6624B scans / 54,749 relations and exact verification while caching construction keys. Full IC/rho wall improves to 140.153x default and 1,158.501x one-worker; construction peak RSS rises to 10,073,538,560 B and 10,044,899,328 B. | Reduce the same `n=59` full cost without hiding its 10 GB construction peak, preserve the `n=53` result, and obtain independent replay. |
| 7. Independent external reproduction and novelty review | **Missing** | Issue [#97](https://github.com/aburan28/crypto/issues/97) and [mtrimoska/EC-Index-Calculus-Benchmarks#1](https://github.com/mtrimoska/EC-Index-Calculus-Benchmarks/issues/1) now include the five-run selected panel, current-head source pin, unknown-scalar result, exact verifier boundary, and the `CONCUR` / `QUALIFIED` / `BREAKS` format. Stage 11 archives GitHub Actions run [34428320022](https://github.com/aburan28/crypto/actions/runs/34428320022) as a project-authored Linux degree-23 reproduction (`STAGE11_RESULTS.md`); it is not independent review. Reviewers bind each gate to per-instance records with [the evidence guide](EXTERNAL_REVIEW_EVIDENCE.md) and [review template](EXTERNAL_NOVELTY_REVIEW_TEMPLATE.json); empty fields remain requests, not completed review. | An unaffiliated reviewer must return a sealed reproduction and source-pinned novelty/correctness assessment. Project-authored CI and replays do not satisfy independence. |

The local measurements above are limited to their stated fields. The meter's
`single_core_seconds` field equals `total_core_seconds`, both user plus system
CPU; it is not a separate measurement of single-core elapsed time. Requested
thread counts and maximum per-process RSS do not by themselves establish
observed single-core execution or aggregate parallel peak memory. Inclusive
outer receipts must not be added to their nested process charges.

Field-size coverage also does not establish literature-instance equivalence.
The n41 ell5 m3 cell differs from WDSat's prominent n41 ell20 m2 experiment.
The n31 standard/GGMP comparison uses different curve coefficients and actual
base sizes (a=1, 31 points versus a=0, 63 points), so its timing ratio cannot
isolate the construction's effect. See [the parameter and sampling comparison](LITERATURE_MATCHING_20260910.md).

External reviewers can bind each gate assessment to individual backend and
instance observations using [the evidence guide](EXTERNAL_REVIEW_EVIDENCE.md)
and [review template](EXTERNAL_NOVELTY_REVIEW_TEMPLATE.json). The packet includes
WDSat, GGMP and 2025 SATIC prior-art questions; empty fields and invitations
remain requests for evidence, not completed external review.

## Current Phase-B matrix

The corrected 480-row SAT view has 141 true positives, 120 true negatives, 219
inconclusive outcomes, no false classifications, and no solver errors. Standard
`n=31` and `n=41` resolve all 40 targets under all three SAT backends. The
`n=31` GGMP cell yields 2 native, 1 WDSat, and 18 CryptoMiniSat true positives;
the rest are inconclusive. All three SAT backends are inconclusive on the 40
`n=59` targets under their caps. Direct MITM classifies all 160 inputs: 80 true
positives and 80 true negatives.

| Backend | Core-seconds | Summed process wall | Peak RSS | Conflicts or operations |
|:--|--:|--:|--:|--:|
| Native XOR SAT | 958.931674 | 1,278.114043 s | 140,333,056 B | 9,233,732 conflicts |
| WDSat, corrected view | 7,404.749336 | 9,524.938771 s | 24,252,416 B | 7,599,249 conflicts |
| CryptoMiniSat | 6,237.273520 | 8,137.581155 s | 122,703,872 B | 19,480,743 conflicts |
| Direct MITM | 759.710741 inclusive outer | 766.538655 s one-CPU elapsed | 85,479,424 B tree | 4,804,252 additions; 4,768,800 pair entries |

The original Stage-26 cells remain immutable. Stage 32 adds 426.718011 charged
core-seconds to turn the two `n=59` WDSat buffer assertions into clean
timeout-inconclusive terminals. The classifications do not change.

## Current n=41 and n=53 end-to-end results

Stage 39 fixes the `n=41` factor base algebraically as the two-torsion
saturation of the Frobenius union generated by masks `[1,2,4,8,16,32]`.
Construction uses zero target samples or discrete-log labels: 2,380 abscissae,
4,759 rational points, 60 signed orbits, and 29 projected columns.

It derives and certifies all 29 logs from 35 relations over 114,688 probes.
Five domain-separated public hash targets solve in 0.259641 seconds versus
0.909281 seconds for signed-Frobenius rho. Precomputation is 4.310245 seconds,
so the online ratio crosses while the amortized ratio remains 5.025825. The
four-core ABBA control reduces median precomputation to 2.010618 seconds and
complete five-target IC wall to 2.479844 seconds, while increasing median
process CPU by 28.7 percent.

Stage 68 retains the current point-defined `n=53` base with 9,964 points and 94
orbit columns, selected without scalar labels or target-subgroup enumeration.
Packed witnesses reduce the support table to 738,197,504 bytes. Rank-aware
collection reaches augmented rank 95 with exactly 95 relations, and ordered,
pipelined support expansion plus fixed/fused PCLMUL arithmetic preserves every
relation hash and group check. On the same public point
`Q=(2565091273463387,5885236316843894)`, fused direct takes 6.039837 wall
seconds / 11.269810 core-seconds / 1,051,987,968 B RSS versus 5.106080 seconds
for signed-Frobenius rho: a 1.182872 loss. Fresh build plus direct is 19.474746
times rho.

Stage 72 then retains the optimized unrelated public hash target. It constructs
no scalar, derives all factor-base logs from 95 relations, and has direct IC and
rho independently recover `7892094459170` with `[d]G=Q`. Unknown-scalar direct
takes 6.571510 seconds and remains 1.467611 times rho; fresh build plus direct
is 25.359925 times rho. The corrected Linux single-core receipt pins parent and
child to CPU 0: direct takes 10.945707 wall / 8.576129 core-seconds versus rho
5.490656 wall / 4.334277 core-seconds, a 1.993515 wall and 1.978676 CPU loss.

Stage 89 retains the post-optimization selected stack. It keeps generic initial
pair sums and generic Frobenius support expansion while selecting split Bloom
indices, insert-hash reuse, shared dual-sign denominators, 16-byte compact query
scratch, combined n=53 dispatch, and the exact direct-PCLMUL pair-query batch.
The specialized query batch improves collection by 17.66--19.31 percent across
three same-binary A/B runs while preserving every relation hash and counter.
Five selected-stack process ratios against the same signed-Frobenius rho target
are 0.946483, 1.019315, 1.094974, 0.932363, and 1.018776: two wins, three losses,
and a 1.018776 median loss. The current-head process reports 5.227982 wall
seconds / 9.193539 core-seconds / 1,053,990,912 B RSS versus 5.131633 seconds
rho. Fresh build plus direct remains 20.262764 times rho. The optimized public
unknown-scalar process takes 5.629163 seconds and remains 1.279112 times rho.

Stage 92 refreshes the forced-single-core comparison on that selected stack.
Linux affinity restricts both parent and child to CPU 0, and the direct arm uses
one query thread with the exact specialized pair batch. Direct takes 9.405716
wall seconds / 6.868316 core-seconds / 863,047,680 B RSS versus 6.945058 wall
seconds for rho: a 1.354303 wall and 1.340180 core loss. Fresh build plus direct
remains 14.691007 times rho.

Stage 99 archives the subsequent optimization chain. Direct low/high
coordinate windows replace the mixed Bloom hash; candidates that already pass
the filter skip a redundant second filter read; and an Itoh--Tsujii addition
chain reduces each degree-53 inversion from 105 field products to 59. The
matched Stage 94 arm reduces wall by 5.09 percent and core cost by 6.96 percent
with identical relation hashes. Stage 95 retains width 4096. Stage 96 rejects a
one-word blocked filter despite reducing exact misses by 69.27 percent because
its query and process wall regress. Stage 97 then records five online direct
wins, but its 0.969444 median misses the predeclared 0.95 gate. Stage 98 reduces
repeated evidence output by 98.13 percent and wall by 1.20 percent while
retaining the base hash, all relation hashes, and every mathematical and
reference validation. Host-to-host rho variation remains material.

Stage 108 archives the four-shard support-table and direct-routing chain.
Sharding reduces matched process wall by 17.32 percent and support setup by
32.55 percent. Replacing the per-query mixed shard hash with xor routing then
wins all five host-identified comparisons, with 0.904540 median wall and
0.912199 median core ratios against mixed routing. The selected EPYC 9V74
direct/rho panel ratios are 0.847095, 0.858937, 0.839190, 0.829081, and
0.838930. The public hash-derived unknown-scalar run takes 4.307977 seconds
versus 4.362556 rho and reconstructs the target without a supplied scalar or
constructed base logs. The exact support still occupies 738,197,504 bytes;
selected median direct core is 2.261565 times rho and fresh build plus direct
is 17.636692 times rho.

Stage 127 adds the distinct current `ic workflow` route on public hash seed
53001.  Its factor base is the target-independent subgroup-orbit recipe with
36,464 points and 344 projected columns; no factor-base logs or target scalar
are constructed.  The selected one-unit stream has 588 relations with canonical
SHA-256 `d8e18605b4a27f307101ecbb5ba9973b50fc01ba272d6f39d216c4fd6a615d37`
and recovers `7892094459170`.  Predicate construction is charged once and emits
344 cofactor multiplications, 1,896,128 derived-coordinate squarings, 18,232
negations and 72,928 canonical-orbit coordinate squarings.  Five matched
one-worker runs have median IC 5.195076 seconds versus 1.580426 rho, a 3.29
times loss; five default-thread runs have median IC 0.891978 versus 1.580089
rho, a 1.777988 rho/IC wall advantage.  Whole-process CPU remains dearer, so
this does not close the full-cost gate.

Stage 128 refreshes one target from the fixed-algebraic `n=41` holdout under
the current workflow.  The public recipe saturates the Frobenius union generated
by masks `[1,2,4,8,16,32]` with two-torsion; it materialises 4,759 points and
29 projected columns from zero target samples or scalar labels.  Fourteen
8,192-probe units yield the same 35 relations in all runs, canonical SHA-256
`8e51f7adab62f4bc351bcb35225fd6522d7cd56eeeac4f162bc00ef6c9b8058a`.
Public hash seed 41301 has no constructed scalar; relation-derived logs recover
`281099696942` and verify `[d]G=Q`.  Five default-thread processes have median
full IC 0.503511 seconds versus 0.118236 rho (4.258367 times slower), 3.336178
whole-process core-seconds and 13,221,888 B RSS.  Five one-worker processes have
median IC 2.829142 seconds versus 0.118650 rho (23.817013 times slower), 2.823121
whole-process core-seconds and 10,944,512 B RSS.  The twelve-process retained
series charges 21.784403 sequential wall-seconds, 37.073066 core-seconds and
14,352,384 B maximum RSS.  This current full-cost result fails the rho gate even
though the older five-target experiment had a fast online descent after shared
precomputation.

Stage 129 adds the exact standard Phase-B `K_1/F_(2^59)` base to the workflow:
the polynomial span `x = sum(c_i z^i), 0 <= i < 9` has 512 abscissae, 483
rational points, 242 negation columns and 231 public projected columns.  It is
not Frobenius-closed, so the implementation refuses the folded table and uses
an exact 116,886-entry compact table.  Public hash seed 59001 constructs no
scalar.  Eight 150,000-probe units scan 307,200,000 summands and return zero
relations; logs records zero solve attempts and a censored rank failure.  The
full capped process costs 18.155840 wall seconds / 78.395618 core-seconds /
14,696,448 B RSS.  There are at most `C(485,3)=18,896,570` unordered triples
against subgroup order 25,179,555,920,633, so even ideal uniform coverage is at
most `7.5047e-7`; 1.2M probes have expected yield at most 0.9006 and 231
relations need at least 307.8M probes before window, collision, dependency and
rank losses.  Five metered same-target rho runs all recover `17861472351607`;
median rho is 3.728387 wall / 3.719220 core-seconds and maximum RSS is
88,522,752 B.  The spent IC cap is already 4.87 times median rho wall and 21.08
times median rho core without reaching one relation.

Stage 130 makes the cofactor obstruction exact across standard widths.  In the
small cyclic quotient of order 22,894, only `4.354e-5` to `4.593e-5` of
unordered triples cancel for every measured `ell=9..16`, close to the exact
`1/22894` penalty.  The best class-aware probe lower bounds among bases fitting
the 4 GiB compact-table budget are 25.73B at `ell=13`, 6.45B at `ell=14`, and
1.582B at `ell=15`; the last table already occupies 2,709,116,566 bytes.  At
`ell=16`, the lower bound falls to 399.24M probes but 2,157,161,086 pair entries
need 9,975,660,343 bytes, over budget.  A charged `ell=14` pilot builds the
132,999,895-entry table and scans 102.4M summands in 85.949930 wall /
146.766002 core-seconds / 896,630,784 B RSS with zero relations, consistent with
its at-most 0.125 expected yield.  The full frontier charges 16 processes,
119.367630 sequential wall-seconds, 185.617302 core-seconds and 896,630,784 B
maximum RSS.  This is a finite no-go result for materialized standard bases
under this budget, not for implicit bases or index calculus in general.

Stage 131 removes the standard base's cofactor-class obstruction without using
secret labels: mapping every `ell=14` parent point through public `[h]` produces
16,308 prime-subgroup points and 8,094 projected columns.  Each complete run
collects 26,796 verified relations over 5.3 million probes and 5.4272 billion
summand scans.  Sparse filtering reduces this to an 802-column core; block
Wiedemann uses 611 products, all 8,094 logs are group-certified, and both modes
recover `17861472351607` for public hash seed 59001.  The default-thread workflow
uses 132.782432 IC wall seconds / 1,387.731459 whole-process core-seconds /
895,205,376 B peak RSS versus 0.553665 seconds rho.  One worker uses
1,219.842388 IC wall seconds / 1,190.875026 whole-process core-seconds /
850,935,808 B peak RSS versus 0.546318 seconds rho.  Online descent takes
0.024512 and 0.023714 seconds respectively, but those online ratios exclude the
large precomputation and do not satisfy the full-cost gate.

Stage 132 selects `ell=15` for this one target.  The public cofactor image has
32,934 points, 16,467 negation orbits and 16,344 projected columns; its exact
compact table stores 542,340,645 pairs.  Both complete modes collect the same
54,749 relations, canonical SHA-256
`e2004ac6e81979caf984e4fa745dc1a1ee99d13892a3f60615a5662179e3ad99`,
over 2.6 million probes and 2.6624 billion scans, then recover and verify the
same scalar as Stage 131.  Against `ell=14`, default IC wall falls from
132.782432 to 126.398334 seconds and whole CPU from 1,387.731459 to
1,194.982194 seconds.  One-worker IC wall falls from 1,219.842388 to
1,009.627483 seconds and whole CPU from 1,190.875026 to 1,000.180872 seconds.
Peak RSS grows by 4.888 times default and 4.798 times one-worker.  The selected
full-cost ratios remain 229.710 and 1,811.815 times rho.

Stage 133 runs two equal-scan neighboring windows and two repetitions of a
source-pinned buffer-reuse candidate.  Window 512 reaches 95.52 percent and
window 2048 reaches 91.47 percent of the selected window-1024 relation
throughput.  Reusing window scratch buffers preserves the 2,157-relation pilot
hash but takes 1.1235 and 1.1725 times the selected relation-unit wall.  The four
fresh rejected processes charge 115.333875 sequential wall-seconds,
1,078.171690 core-seconds and 4,378,951,680 B peak RSS.  The selected Stage 132
source remains unchanged.

Stage 134 stores one packed pair witness beside each of the 542,340,645
compact rests under an 8 GiB ceiling.  Exact group re-addition rejects truncated
rest collisions, and sorting exact candidates makes the parallel build
deterministic.  Default IC falls from 126.398334 to 85.747666 seconds and whole
CPU from 1,194.982194 to 853.002930 seconds, with RSS increasing from
4,375,724,032 to 6,565,560,320 B.  One-worker IC falls from 1,009.627483 to
724.574469 seconds and CPU from 1,000.180872 to 720.258242 seconds, with RSS
increasing from 4,082,794,496 to 6,253,510,656 B.  Both modes reproduce canonical
relation SHA-256
`e2004ac6e81979caf984e4fa745dc1a1ee99d13892a3f60615a5662179e3ad99`
and recover the same verified unknown scalar.  Full cost remains 152.615 and
1,329.575 times rho.

Stage 135 caches 542,340,645 packed keys (4,338,725,160 temporary bytes) from
the counting pass.  Default build wall falls from 21.389824 to 11.712773 seconds
and full IC from 85.747666 to 76.342293 seconds; one-worker build falls from
198.245643 to 101.933057 seconds and full IC from 724.574469 to 648.780325
seconds.  Both full modes preserve the 54,749-relation hash and verified scalar.
Peak RSS rises to 10,072,309,760 B default and 10,044,899,328 B one-worker.
Full cost remains 140.153 and 1,158.501 times rho.

The narrow supported conclusion is unchanged: this is strong internal
engineering and finite public toy-research evidence. The known SAT-based
point-decomposition and Frobenius-invariant-factor-base literature remain the
prior-art baseline. The finite `n=41` online descent is faster after
precomputation, and the selected `n=53` panel passes the predeclared online wall
threshold on the recorded host. Its memory, core, and fresh-build costs remain
above rho. Licensed Magma and unaffiliated
reproduction/novelty review remain open. This is not a new Koblitz
index-calculus SOTA result.
