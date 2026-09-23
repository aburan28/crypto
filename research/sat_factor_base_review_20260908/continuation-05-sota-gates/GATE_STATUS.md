# Koblitz index-calculus SOTA gate status

Current through Stage 128, plus the additive Stage 159 native-F4
single-target supplement, Stage 160 fixed-X1 construction specialization, and
Stages 161 and 162 trusted-mask hashing and grouped critical-pair selection,
plus Stage 163 shape-selected block-4 M4RI elimination, and Stage 164 dense
pair/column indexing, batch UPDATE, and trimmed block-8 M4RI, plus Stage 165
target-independent sparse fixed-X1 ordering, plus Stage 166 generation-tagged
dense cancellation before F4 monomial-product sorting, plus Stage 167
deterministic parallel fixed-X1 batches with an exact single-thread control,
plus Stage 168 one-batch scheduling of the selected 39-system prefix,
2026-09-23. The historical optimization chain is
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
the note's §11), chaining to the Stage 125, Stage 124 and Stage 109 seals.
Stage 127 re-seals it a fourth time over that ledger's Round 4, also of the
same day: the unit's conversion from native counters to group additions now
uses ratios pinned in `docs/ic/calibration.json` rather than factors measured
on the host at the start of each run (the note's §12). That round reprices
every row by up to 9.5 per cent without moving a single native counter on any
of the 537 rows compared, and chains to the Stage 126, 125, 124 and 109 seals.
Stage 128 re-seals it a fifth time over that ledger's Round 5, on 2026-09-22:
the walk's sixteen restart offsets are drawn on first use instead of at setup,
so a row is charged for the offsets it took rather than for a pool it may never
reach (the note's §13). Because the offsets come off their own random stream and
are taken in a fixed cycle, the before and after walk bit-identical
trajectories, and all 537 rows are identical in every counter but the pool's
with zero trajectories moved; the saving is confined to the small rungs. It
chains to the Stage 127, 126, 125, 124 and 109 seals. No Koblitz gate fact moved
at any of these re-seals and no prior frozen artifact changed.

The campaign has target-independent algebraic factor bases; matched native-XOR,
WDSat, CryptoMiniSat, direct-MITM, GGMP, and signed-Frobenius-rho controls; a
balanced 160-instance PDP panel through `n=59`; public unknown-scalar end-to-end
runs at degrees 23, 31, 41, and 53; and same-target known-answer and
scalar-blind construction, rank, solve, and rho comparisons at `n=53`. It has not passed all seven gates and
does not establish a new state of the art.

| Gate | Status | Current evidence | Remaining requirement |
|:--|:--|:--|:--|
| 1. Charge every stage and resource | **Partial overall** | The corrected Phase-B matrix charges 21,038.596137 core-seconds. Stages 99 and 108 retain separately metered build, base/support construction, direct, rho, wall, core-second, RSS, collection, query, validation, LA, and output charges. The selected Stage 107 panel has 318.511978 available build+science core-seconds / 112.437872 sequential wall-seconds / 1,720,258,560 B maximum sampled tree RSS. Through Stage 168 the native-F4 campaign's measured lower bound is 8,381.199829 core-seconds / 7,988.745489 sequential wall-seconds / 6,310,576,128 B maximum process RSS across 223 metered components. | Licensed Magma process resources are absent. The complete native-F4 campaign total is `null` because some worktree setup, focused checks and profiling, score composition, documentation, and Git operations lack complete receipts. The Stage-168 build is reused from and already charged in Stage 167 rather than double-counted. |
| 2. WDSat, CryptoMiniSat, Magma F4, MITM, GGMP | **Partial** | Native XOR SAT, WDSat, CryptoMiniSat, and direct MITM ran on the exact 160-input packet. Standard `n=31`/`n=41` and GGMP `n=31` are represented. The Stage-32 successor removes the two original WDSat buffer errors without rewriting Stage 26. Stages 159 through 168 retain one exact-target `n=59` true-positive from the repository's native F4, with the same-target SAT and MITM receipts, exact source validation, a target-independent but post-hoc-selected outer schedule, exact construction controls, and matched parallel/single-thread F4 execution. | Execute all 160 Stage-22 Magma inputs on a licensed host under the frozen one-thread/no-retry contract. The one-target native-F4 supplement is a different formulation and does not replace licensed Magma or a full native-F4 panel. |
| 3. Single-core, core-seconds, memory, conflicts, wall | **Partial because Magma is absent and selected CPU-0 is stale** | Every executed Phase-B arm retains wall, core-seconds, peak RSS, conflicts or operations, tree memory, and workflow wall. Stage 168's 12-thread native F4 records 4.987010 wall / 28.179809 total core-seconds / 2,763,997,184 B RSS; the same source's Stage-167 one-thread control records 19.828459 wall / 19.805719 single-core seconds / 904,085,504 B RSS. Both charge 16,821,055,616 word XORs and report conflicts as `null`. The selected Stage 107 panel reports median direct 3.630413 wall / 9.739701 core-seconds / 1,055,776,768 B RSS versus rho 4.308264 wall / 4.307199 core-seconds. Stage 92 pins CPU 0, but predates the four-shard route. | Supply the same fields for licensed Magma F4 and refresh forced-single-core on the selected stack. Requested one- and twelve-thread execution and CPU time are recorded; observed affinity was not pinned for Stage 168. |
| 4. Scale through `n=31`, `n=41`, and a larger PDP regime | **Satisfied for finite execution coverage** | Phase B covers `n=31`, GGMP `n=31`, `n=41`, and PDP-only `n=59`. Unknown-scalar end-to-end controls cover `n=31` and `n=41`. Stage 42 constructs, ranks, solves, and verifies an exact same-target known-answer `n=53` instance. | The evidence is finite and toy-sized; it is not an asymptotic scaling law. The `n=59` arm remains PDP-only. |
| 5. Unknown scalar with no constructed factor-base logs | **Satisfied for finite degrees 23, 31, 41, and 53** | Stage 108 archives public hash-seed-53001 without constructing or supplying its scalar, derives all 94 factor-base logs from 95 verified relations, and has direct IC and rho independently recover `d=7892094459170` with `[d]G=Q`. Selected direct is 4.307977 seconds versus 4.362556 rho, ratio 0.987489. | Repeat on independent n=53 public seeds and obtain unaffiliated replay; these strengthen rather than replace the finite gate-5 execution. |
| 6. Full cost against automorphism-optimized Pollard rho | **Online wall gate passed; full-cost/core gate false** | Stage 108 selects the four-shard direct route after five direct/mixed wins with 0.904540 median wall and 0.912199 median core ratios. Its identified EPYC 9V74 panel has five direct/rho wins and 0.839190 median wall ratio. Median direct core remains 2.261565 times rho, retained support is 738,197,504 B, and fresh build plus direct is 17.636692 times rho. | Reduce memory below 512 MiB, core ratio to at most 1, and full-build cost while retaining the online wall result; refresh selected CPU-0 and obtain independent replay. All improvements are finite constants, not exponent changes. |
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

The narrow supported conclusion is unchanged: this is strong internal
engineering and finite public toy-research evidence. The known SAT-based
point-decomposition and Frobenius-invariant-factor-base literature remain the
prior-art baseline. The finite `n=41` online descent is faster after
precomputation, and the selected `n=53` panel passes the predeclared online wall
threshold on the recorded host. Its memory, core, and fresh-build costs remain
above rho. Licensed Magma and unaffiliated
reproduction/novelty review remain open. This is not a new Koblitz
index-calculus SOTA result.

## Stage 159 native-F4 single-target supplement

Stage 159 runs the repository's full Boolean `f4-f2` implementation on one
authenticated `n=59, ell=9, m=3` Phase-B target. The solver receives no truth
class, witness, discrete-log label, or target-subgroup enumeration. A direct
three-summand S4 expansion exceeded the watchdog at 5.25--6.31 GB RSS in two
censored attempts. The retained formulation fixes `X1`, runs full F4 over the
remaining 18 Boolean variables, and checks every returned root by exact curve
lifting and group addition.

The first complete arm took 205.702959 wall seconds, 204.982191 core-seconds,
and 950,878,208 bytes RSS. Skipping only `X1` masks with no rational curve lift
reduced it to 100.179380 wall seconds, 99.920352 core-seconds, 921,681,920
bytes RSS, 65 F4 calls, and 87,513,949,370 elimination word XORs while
returning the same exact witness. This is a 2.05x solver-stage engineering
gain. A trace-homomorphism equation was retained as a rejection: it raised
wall by 18.2 percent and word XORs by 58.8 percent despite lowering RSS.

On the same local host and exact target, direct MITM completed in 2.999661
seconds / 2.991745 core-seconds / 45,842,432 bytes, so selected F4 remains
33.40x slower by wall and CPU and 20.11x larger by RSS. Native XOR SAT reached
its 100,000-conflict cap in 8.906472 seconds. The authenticated Linux receipts
for this target leave WDSat and CryptoMiniSat censored at their 120-second
watchdogs. Native F4 is a distinct implementation and formulation; it does
not satisfy the still-missing licensed Magma F4 gate. Its conflict field is
`null` because F4 does not expose SAT conflicts.

The Stage 159 measured research lower bound, including five distinct clean
builds, six F4 attempts, local controls, failed builds, and dependency-fetch
receipts, is 1,596.064123 sequential wall-seconds, 1,565.839519 core-seconds,
and 6,310,576,128 bytes maximum process RSS. Complete campaign cost remains
`null` because some early incremental checks lacked an outer meter. The full
160-input native-F4 panel, licensed Magma, end-to-end IC/rho cost, and
unaffiliated review remain open. The canonical result is
[`stage-159-native-f4-single-target-20260922/result.json`](stage-159-native-f4-single-target-20260922/result.json).

## Stage 160 fixed-X1 constant-linear construction

Stage 160 preserves the exact Stage 159 factor base, blind target, X1 schedule,
65 F4 calls, 3,864,601 Boolean terms, equation fingerprint, 87,513,949,370
elimination word XORs, algebraic roots, and group-verified witness. It changes
only how each fixed-X1 S4 system is constructed. The shared symbolic `X2*X3`
product is formed once, while fixed `X1` and target powers act through exact
linear multiplication maps in the polynomial basis. A unit gate compares the
complete specialized and generic polynomial vectors across 60 field/value
cases.

The clean blind run reduces fixed-X1 construction from 26.323585 seconds to
0.973558 seconds, a 27.04x construction improvement. Whole-target F4 wall
falls from 100.179380 to 73.757990 seconds and core cost from 99.920352 to
73.585054 seconds, a 1.36x wall improvement. Peak RSS falls slightly from
921,681,920 to 910,934,016 bytes. The clean one-job build plus custody run is
241.262254 wall seconds / 236.554363 core-seconds / 1,683,472,384 bytes peak
RSS.

Direct MITM remains 24.59x faster by same-host wall, 24.60x cheaper by CPU,
and 19.87x smaller by RSS. Four other metered ideas are retained as rejected:
global matrix ordering, trailing-zero row trimming, a hashed reducer index,
and a dense reducer index. Through this stage the measured campaign lower
bound is 2,319.762066 sequential wall-seconds, 2,281.647854 core-seconds, and
6,310,576,128 bytes maximum process RSS. Complete campaign cost remains
`null` because the incremental compiles and focused tests were not outer
metered.

This remains one finite toy-PDP target. Licensed Magma, the full native-F4
panel, end-to-end IC/rho cost, and unaffiliated reproduction remain open. The
canonical additive result is
[`stage-160-fixed-x1-constant-specialisation-20260922/result.json`](stage-160-fixed-x1-constant-specialisation-20260922/result.json).

## Stage 161 trusted-mask hashing

A charged 20-second sample of the Stage 160 binary collected 3,334
top-of-stack samples: 1,143 in dense echelon elimination, 533 in critical-pair
installation, 238 in column packing, 225 directly in SipHash writes, and 118
in queued-pair retention. Stage 161 replaces SipHash only for internal maps
and sets keyed by solver-constructed `u64` monomial masks. The deterministic
SplitMix64 hash changes bucket placement only; hashbrown still resolves every
collision by exact key equality. Tuple keys and other data retain their prior
hashers.

The clean blind run preserves the equation fingerprint, 3,864,601 terms, all
65 F4 calls, matrix dimensions, pair counters, 87,513,949,370 elimination word
XORs, roots, and group-verified witness. Symbolic/matrix build falls from
30.523303 to 17.204156 seconds, a 1.77x phase improvement. Whole-target F4
wall falls from 73.757990 to 60.645599 seconds and core cost from 73.585054 to
60.498754 core-seconds, a 1.22x improvement. Peak RSS rises 5.85 percent to
964,263,936 bytes.

Same-host direct MITM remains 20.22x faster by wall and CPU and 21.03x smaller
by RSS. The clean one-job build plus custody run costs 228.347900 wall seconds,
224.535510 core-seconds, and 1,692,057,600 bytes peak RSS. Through this stage
the measured campaign lower bound is 2,684.534288 sequential wall-seconds,
2,641.558033 core-seconds, and 6,310,576,128 bytes maximum process RSS.
Complete cost remains `null` because one sandbox-denied profiler attempt was
interrupted before its receipt and incremental compiles and focused tests were
not all outer metered.

The result remains one finite toy-PDP target. Licensed Magma, the full
native-F4 panel, end-to-end IC/rho cost, and unaffiliated reproduction remain
open. The canonical additive result is
[`stage-161-fast-mask-hash-20260923/result.json`](stage-161-fast-mask-hash-20260923/result.json).

## Stage 162 grouped critical-pair selection

A charged post-hash profile put 695 of 3,327 top-of-stack samples in
`State::insert`. Binary disassembly placed its dominant samples in the
quadratic scan that tests whether any other new-pair LCM divides the current
LCM. Stage 162 groups equal LCMs and finds proper divisor LCMs by exact submask
lookup. It then emits the same lowest-index duplicate representative in the
same pair order as the previous UPDATE implementation. A 3,600-case randomized
test compares selected pairs, order, and chain/product counters directly with
that quadratic reference.

The clean blind run preserves every equation, term, F4 call, matrix dimension,
pair counter, 87,513,949,370 elimination word XORs, root, and exact witness.
Whole-target wall falls from 60.645599 to 57.862894 seconds and core cost from
60.498754 to 57.736423 core-seconds, a 1.048x improvement. Peak RSS changes by
0.40 percent to 968,146,944 bytes.

Same-host direct MITM remains 19.29x faster by wall, 19.30x cheaper by CPU,
and 21.12x smaller by RSS. The clean build plus custody run costs 225.353725
wall seconds, 221.602144 core-seconds, and 1,692,352,512 bytes peak RSS. The
measured campaign lower bound through Stage 162 is 3,029.348391 sequential
wall-seconds, 2,981.632960 core-seconds, and 6,310,576,128 bytes maximum process
RSS. Complete campaign cost remains `null`.

This remains a finite one-target solver-stage result. Licensed Magma, the full
native-F4 panel, end-to-end IC/rho cost, and unaffiliated reproduction remain
open. The canonical additive result is
[`stage-162-grouped-pair-selection-20260923/result.json`](stage-162-grouped-pair-selection-20260923/result.json).

## Stage 163 shape-selected block-4 M4RI

Stage 163 routes matrices with at least 128 rows and 256 columns and at most
four times as many columns as rows through block-4 Method of Four Russians
elimination. Other shapes retain streaming elimination, and
`PQ_F4_DISABLE_M4RI=1` provides a same-binary control. Four large synthetic
matrix shapes across word boundaries return identical pivot-column sets and
canonical row spaces under the two kernels.

The same-binary control takes 58.370052 wall seconds / 57.840795 core-seconds;
M4RI takes 54.538915 / 54.413164, a 1.070x wall and 1.063x CPU improvement.
The clean selected run takes 54.773779 wall seconds, 54.662604 core-seconds,
and 915,668,992 bytes peak RSS. It routes 196 matrices through M4RI and reduces
counted elimination-and-table work from 87,513,949,370 to 45,879,309,338 word
XORs. The equation fingerprint, target classification, roots, and exact
group-verified witness are unchanged; a few intermediate basis-path row and
reducer counters differ slightly and are retained rather than called
identical.

Direct MITM remains 18.26x faster by wall, 18.27x cheaper by CPU, and 19.97x
smaller by RSS. The selected clean build plus run costs 223.248497 wall
seconds, 218.415925 core-seconds, and 1,689,485,312 bytes peak RSS. The measured
campaign lower bound through Stage 163 is 3,589.954839 sequential wall-seconds,
3,532.297699 core-seconds, and 6,310,576,128 bytes maximum process RSS.
Complete campaign cost remains `null`.

The result remains a finite one-target solver-stage improvement. Licensed
Magma, a full native-F4 panel, end-to-end IC/rho cost, and unaffiliated
reproduction remain open. The canonical additive result is
[`stage-163-m4ri-echelons-20260923/result.json`](stage-163-m4ri-echelons-20260923/result.json).

## Stage 164 dense native-F4 pipeline

Stage 164 keeps the same algebraic factor base, blind `n=59, ell=9, m=3`
target, 65 fixed-X1 F4 calls, 3,864,601 Boolean terms, equation fingerprint,
two algebraic roots, and exact three-point curve witness. It adds dense
epoch-stamped LCM grouping on the 18-variable systems, exact restricted
submask covers, reusable pair-selection storage, indexed symbolic reducers,
dense monomial-to-column maps, order-preserving batch UPDATE, cache-local
M4RI tables, consecutive-pivot extraction, and trailing-zero table trimming.
Eleven focused tests include 3,600 randomized UPDATE comparisons, 230,400
indexed-reducer comparisons, randomized batch-vs-sequential state equality,
100,000 DegRevLex comparator comparisons, and exact large-matrix row-space
equality against streaming elimination.

The clean blind F4 process takes 36.450481 wall seconds, 36.319425
core-seconds, and 946,388,992 bytes peak RSS. It charges 27,264,366,281
elimination-and-table word XORs. The same clean binary takes 48.447073 seconds
with M4RI disabled and 38.698749 seconds with batch insertion disabled; all
arms retain the same equations, pair-criterion totals, roots, and witness.
Relative to Stage 163, selected wall improves by 1.503x and core cost by
1.505x. A separately metered lockfile dependency fetch costs 0.316796 wall
seconds / 0.102986 core-seconds, while dependency acquisition plus the clean
one-job build and custody run cost 214.044317 wall seconds / 208.873448
core-seconds / 1,689,780,224 bytes peak RSS.

Direct MITM remains 12.15x faster by wall, 12.14x cheaper by CPU, and 20.64x
smaller by RSS. The measured native-F4 campaign lower bound through Stage 164
is 6,464.947902 sequential wall-seconds, 6,383.990738 core-seconds, and
6,310,576,128 bytes maximum process RSS across 132 components. Complete cost
remains `null` because some development checks were not outer-metered.

This remains one finite toy-PDP solver-stage improvement. Licensed Magma, the
full native-F4 panel, end-to-end IC/rho cost, and unaffiliated reproduction and
novelty review remain open. The canonical additive result is
[`stage-164-native-f4-dense-20260923/result.json`](stage-164-native-f4-dense-20260923/result.json).

## Stage 165 target-independent sparse fixed-X1 order

Stage 165 changes only the order in which the symmetric direct-S4 solver fixes
the first factor-base coordinate. Coefficient masks are sorted by Hamming
weight and then numeric value. This permutation is defined solely by the
published polynomial-basis coefficients: it does not inspect the target,
truth class, witness, subgroup, discrete-log labels, or a known scalar.
`PQ_F4_X1_ORDER=ascending` retains the previous same-binary order, and a unit
test checks exact permutation and ordering for dimensions zero through twelve.

On the clean blind target, the selected process visits 85 masks, constructs
and completes 39 rational fixed-X1 F4 systems, and returns the same exact three
curve points as the ascending schedule. It takes 21.512089 wall seconds,
21.484441 core-seconds, 948,502,528 bytes peak RSS, and 16,821,055,616 charged
word XORs. The same clean binary in ascending order visits 138 masks, completes
65 systems, and takes 34.446599 wall / 34.412474 core-seconds, for a 1.601x
wall improvement. Three interleaved development pairs have 21.368530 versus
34.230726 wall-second medians, also 1.602x.

The schedule itself is target-independent, but it was selected post hoc after
inspecting this target. The result is therefore target-specific engineering,
not evidence that the order improves expected time over fresh or adversarial
targets. Direct MITM remains 7.17x faster by wall, 7.18x cheaper by CPU, and
20.69x smaller by RSS. The clean build plus custody run costs 200.076092 wall
seconds / 195.904177 core-seconds / 1,723,154,432 bytes peak RSS.

The measured native-F4 campaign lower bound through Stage 165 is 6,902.647655 sequential wall-seconds, 6,817.038349 core-seconds, and 6,310,576,128 bytes
maximum process RSS across 150 components. Complete cost remains `null`.

Licensed Magma, the full native-F4 panel, end-to-end IC/rho cost, and
unaffiliated reproduction and novelty review remain open. The canonical
additive result is
[`stage-165-x1-weight-order-20260923/result.json`](stage-165-x1-weight-order-20260923/result.json).

## Stage 166 dense F4 product cancellation

Stage 166 changes the construction of F4 rows after monomial multiplication.
For Boolean systems with at most twenty variables, one generation-tagged
`u32` array now records the parity of every mapped output mask. The first
occurrence stores a compact key whose natural order is descending degree and
ascending mask. The solver sorts only odd-parity survivors instead of sorting
all mapped terms and cancelling adjacent duplicates afterward.
`PQ_F4_DISABLE_DENSE_MUL=1` restores the exact previous path. A randomized
test compares the two products through twenty variables, including generation
wraparound, and the existing F4/Buchberger and M4RI/streaming certificates
continue to pass.

The clean blind process takes 19.887753 wall seconds, 19.858236 core-seconds,
and 900,169,728 bytes peak RSS. Across 818,171 monomial products it maps
651,502,431 input terms to 549,909,017 surviving terms, cancelling 101,593,414
terms before sorting. The exact clean disabled control takes 21.582929 wall /
21.552617 core-seconds / 939,032,576 bytes, for a 1.085x wall improvement.
Three interleaved development pairs have 20.107567 versus 21.526022
wall-second medians, a 1.071x improvement. Equation fingerprint, F4 calls,
matrix and algebraic counters outside the new product counters, roots, and the
exact curve-verified witness agree.

Direct MITM remains 6.63x faster by wall, 6.64x cheaper by CPU, and 19.64x
smaller by RSS. The clean build plus custody run costs 193.162520 wall seconds
/ 187.303839 core-seconds / 1,703,051,264 bytes peak RSS. The failed default
offline-cache build is retained and charged; it stopped at vendoring because
`redis 0.32.7` was absent and did not fall back to the network.

The measured native-F4 campaign lower bound through Stage 166 is 7,400.453979 sequential wall-seconds, 7,307.628306 core-seconds, and 6,310,576,128 bytes
maximum process RSS across 180 components. Complete cost remains `null`.

Licensed Magma, the full native-F4 panel, end-to-end IC/rho cost, and
unaffiliated reproduction and novelty review remain open. The canonical
additive result is
[`stage-166-dense-product-cancellation-20260923/result.json`](stage-166-dense-product-cancellation-20260923/result.json).

## Stage 167 deterministic parallel fixed-X1 batches

Stage 167 exploits only the independence already present between rational
fixed-X1 systems. The backend gathers the next deterministic batch in the
algebraic schedule, solves every launched system through a bounded Rayon
pool, waits for the complete batch, charges every result, and then inspects
results in schedule order. It does not receive target labels, a known scalar,
the witness, subgroup enumeration, or factor-base logarithms.

The clean thirteen-thread F4 process takes 5.589951 wall seconds, 36.302718
total core-seconds, and 3,033,563,136 bytes peak RSS. The same clean binary at
batch one takes 19.828459 wall seconds, 19.805719 single-core seconds, and
904,085,504 bytes. Parallel wall improves 3.547x while total CPU rises 1.833x
and RSS rises 3.355x. Three interleaved development pairs give 5.534149 versus
19.708587 wall-second medians, a 3.561x speedup. Both clean arms visit 85
masks, complete the same 39 F4 systems, charge 16,821,055,616 word XORs, have
the same equation fingerprint, roots and exact witness, and execute no
speculative system after the successful one.

Batch width thirteen was selected post hoc: this target's first valid
relation is in rational system 39, so three batches end exactly there. The
mechanism is label-blind, but this width result is target-specific rather than
fresh-target or expected-time evidence. Direct MITM remains 1.86x faster by
wall, 12.13x cheaper by CPU, and 66.17x smaller by RSS. Clean build plus the
parallel custody run costs 173.190335 wall seconds / 199.814249 core-seconds /
3,033,563,136 bytes peak RSS. The full-cost gate remains false.

The measured native-F4 campaign lower bound through Stage 167 is 7,948.233270 sequential wall-seconds, 8,144.184517 core-seconds, and 6,310,576,128 bytes
maximum process RSS across 215 components. Complete cost remains `null`.

Licensed Magma, the full native-F4 panel, fresh-target validation, end-to-end
IC/rho cost, and unaffiliated reproduction and novelty review remain open. The
canonical additive result is
[`stage-167-parallel-x1-batches-20260923/result.json`](stage-167-parallel-x1-batches-20260923/result.json).

## Stage 168 one deterministic batch of 39 systems

Stage 168 puts the same selected 39 rational fixed-X1 systems into one
deterministic batch, removing the two barriers in Stage 167. Twelve Rayon
workers execute the full batch dynamically; every system is charged before
the completed results are inspected in algebraic schedule order. The
factor base, equations, roots, work counters, and exact witness are unchanged.

The clean F4 process takes 4.987010 wall seconds, 28.179809 total core-seconds,
and 2,763,997,184 bytes peak RSS; `single_core_seconds` is `null`. Relative to
Stage 167's clean batch-13 arm, wall improves 1.121x, CPU falls to 0.776x, and
RSS falls to 0.911x. Three twelve-thread development runs have a 4.744985
wall-second median. The selected source and clean build are unchanged, and
the build cost is reused rather than charged twice.

Batch size 39 is exactly the already known successful prefix and was selected
post hoc. It is not a generic stopping rule: an unknown fresh target could
require additional batches or no relation within the cap. Direct MITM remains
1.66x faster by wall, 9.42x cheaper by CPU, and 60.29x smaller by RSS. Clean
build plus the selected custody run is 172.146798 wall seconds / 191.677868
core-seconds / 2,763,997,184 bytes peak RSS. The full-cost gate remains false.

The measured native-F4 campaign lower bound through Stage 168 is 7,988.745489
sequential wall-seconds, 8,381.199829 core-seconds, and 6,310,576,128 bytes
maximum process RSS across 223 components. Complete cost remains `null`.

Licensed Magma, the full native-F4 panel, fresh-target validation, end-to-end
IC/rho cost, and unaffiliated reproduction and novelty review remain open. The
canonical additive result is
[`stage-168-single-batch-39-20260923/result.json`](stage-168-single-batch-39-20260923/result.json).
