# Koblitz index-calculus SOTA gate status

Current through Stage 128, 2026-09-21. The historical optimization chain is
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

The campaign has target-independent algebraic factor bases; matched native-XOR,
WDSat, CryptoMiniSat, direct-MITM, GGMP, and signed-Frobenius-rho controls; a
balanced 160-instance PDP panel through `n=59`; public unknown-scalar end-to-end
runs at degrees 23, 31, 41, and 53; and same-target known-answer and
scalar-blind construction, rank, solve, and rho comparisons at `n=53`. It has not passed all seven gates and
does not establish a new state of the art.

| Gate | Status | Current evidence | Remaining requirement |
|:--|:--|:--|:--|
| 1. Charge every stage and resource | **Partial overall** | The corrected Phase-B matrix charges 21,038.596137 core-seconds. The current `n=53` campaign retains 318 measured processes / 1,395.640868 sequential wall-seconds / 3,192.181904 core-seconds / 245,972,992 B maximum RSS, including rejected variants. Selection emits native projected-predicate counts and every workflow stage emits wall, CPU and cumulative RSS. | Licensed Magma process resources are absent. Preinstalled OS/toolchain acquisition remains an explicit exclusion. |
| 2. WDSat, CryptoMiniSat, Magma F4, MITM, GGMP | **Partial** | Native XOR SAT, WDSat, CryptoMiniSat, and direct MITM ran on the exact 160-input packet. Standard `n=31`/`n=41` and GGMP `n=31` are represented. The Stage-32 successor removes the two original WDSat buffer errors without rewriting Stage 26. | Execute all 160 Stage-22 Magma inputs on a licensed host under the frozen one-thread/no-retry contract, seal the return before truth scoring, and report F4 resources. |
| 3. Single-core, core-seconds, memory, conflicts, wall | **Partial because Magma is absent** | Every executed Phase-B arm retains wall, core-seconds, peak RSS, conflicts or operations, tree memory, and workflow wall. The current scalar-blind `n=53` workflow has a five-pair one-worker panel: median IC 5.195076 wall, rho 1.580426 wall, whole-process 6.757385 core-seconds and 104,644,608 B RSS. The default-thread panel reports IC 0.891978 versus rho 1.580089. | Supply the same fields for licensed Magma F4. SAT conflicts remain inapplicable to exact pair-table arms and are reported as null. |
| 4. Scale through `n=31`, `n=41`, and a larger PDP regime | **Satisfied for finite execution coverage** | Phase B covers `n=31`, GGMP `n=31`, `n=41`, and PDP-only `n=59`. Stage 128 refreshes a unified current `n=41` unknown-scalar full workflow from algebraic base construction through rho. Stage 42 and the current workflow cover exact same-target `n=53`. | The evidence is finite and toy-sized; it is not an asymptotic scaling law. The `n=59` arm remains PDP-only rather than end-to-end IC. |
| 5. Unknown scalar with no constructed factor-base logs | **Satisfied for finite degrees 23, 31, 41, and 53** | Stage 108 archives public hash-seed-53001 without constructing or supplying its scalar, derives all 94 factor-base logs from 95 verified relations, and has direct IC and rho independently recover `d=7892094459170` with `[d]G=Q`. Selected direct is 4.307977 seconds versus 4.362556 rho, ratio 0.987489. | Repeat on independent n=53 public seeds and obtain unaffiliated replay; these strengthen rather than replace the finite gate-5 execution. |
| 6. Full cost against automorphism-optimized Pollard rho | **Current `n=41` full-cost loss; `n=53` default-thread wall pass only** | Stage 128 charges the current `n=41` workflow in one process: default-thread IC/rho is 4.258367 and one-worker IC/rho is 23.817013. The current scalar-blind `n=53` workflow has default-thread rho/IC 1.777988 but one-worker IC remains about 3.29 times rho and total core cost remains higher. | Run the larger PDP regime end to end, reduce the current one-core and total-core ratios to at most 1 while retaining the `n=53` default-thread wall result, and obtain independent replay. |
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

The narrow supported conclusion is unchanged: this is strong internal
engineering and finite public toy-research evidence. The known SAT-based
point-decomposition and Frobenius-invariant-factor-base literature remain the
prior-art baseline. The finite `n=41` online descent is faster after
precomputation, and the selected `n=53` panel passes the predeclared online wall
threshold on the recorded host. Its memory, core, and fresh-build costs remain
above rho. Licensed Magma and unaffiliated
reproduction/novelty review remain open. This is not a new Koblitz
index-calculus SOTA result.
