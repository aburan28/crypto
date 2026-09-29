# Further 2× matrix-F5 complete-call experiment

## Boundary, hypothesis, and scope

The accepted opt-in one-thread reference is the combined matrix-F5 mode
archived in `research/gf2_table_reuse_20260929/RESULT.md`: 125.89 ms marginal
median for `f5_n24_m24_d4` on its pinned AMD EPYC 7763 runner. The unmodified
reduced-output default was 254.87 ms on that same run. A **further 2×** means
the same complete call, including criterion, row build, reduction, and full
`Vec<F2BoolPoly>` output, at no more than half the paired reference time.
The output form remains selective echelon. Rank, row space, and pruning must
match the reference; a pivot-choice arm may change raw echelon rows, term
count, and reduction operation count, which must all be reported.
No IC online-time or DLP speedup follows from a matrix-F5 stage result.

The frozen source reference is `aa38ca38eff9fcce5b51fcdb70910aea8343616e`.
The SHA-256 digests of its relevant files are `20d2c0933daea0673f2b642226c056f8bdef8ad933614d994737e96923a1d951`
(`gf2_elim.rs`), `e8c66b4d9163160ed32b848820beb5013eb384ddde09a8bc7b4919565ae93959`
(`matrix_f5_f2.rs`), `934dba95f3ec85bcb4108708d1d0f1c16e334cae65962ae6deb5fcf9d565071e`
(`pq_groebner_f2.rs`), and `35cafca3dc64e8cf25a1ec63b1338eddb8da4725d0ca1c0f9dedc98313cb16aa`
(`f4_f2_bench.rs`).

The first bounded hypotheses are: (1) choosing a lighter pivot from a small
fixed candidate window may limit fill, reduce the 13.7-million-term output,
and shorten both reduction and unpacking; (2) wider pivot blocks assembled
from more eight-bit tables may reduce repeated matrix sweeps. The first
hypothesis can change raw echelon rows while preserving rank and canonical
row space, so it requires its own arm and an explicit output contract check.
Neither hypothesis is assumed beneficial before paired measurement.
The third bounded hypothesis is that six- or seven-bit Gray-code tables
reduce table-construction work enough to offset additional matrix passes;
compare widths 8 (reference), 7, and 6 at four tables per pass.
The fourth hypothesis is that rows remain sparse across the leading
degree-four monomial band. Convert packed rows to sorted column indices,
echelonize that band with a lowest-weight pivot per leading-column bucket,
then use the existing packed kernel on the lower-degree suffix. Return the
complete echelon basis as canonical `F2BoolPoly` rows and exact rank. Charge
the sparse phase's term visits separately from dense word XORs; its output
and operation counts may differ, while canonical row space must match.
The fifth hypothesis is that a stable lightest-first ordering of packed
input rows gives the dense kernel lighter pivots and less fill, without the
billions of list merges seen in a full sparse pass. Charge the sort within
the reduction phase and keep output, rank, and row-space checks unchanged.

Before reserving an eligible x86-64 runner, run a nonpromoting Apple ARM64
screen on the frozen seed and holdout A. Use one thread, one binary, one warmup
per arm, and three alternating reference/candidate pairs per seed. Compare
four, six, and eight tables, then pivot-candidate windows 1 and 4, then table
widths 8, 7, and 6, then dense reference versus hybrid sparse leading band,
then unsorted versus stable lightest-first packed rows as separate arms. This
screen can reject a regression or a correctness failure, but its ratios
cannot establish the requested gain. Preserve the raw screen receipt.

## Frozen workload and accounting

Use the seven F5 cases in `examples/f4_f2_bench.rs` with the four seeds
`0`, `badc0de1`, `5eed2026`, `f5c02a28`, each case run in a separate process
under one release binary. The primary is `f5_n24_m24_d4`. Reference and
candidate share `KIC_F5_ECHELON=2`, `KIC_F5_FUSED_BUILD=1`,
`KIC_F5_DIRECT_PACK=1`, `KIC_GF2_FORCE_AVX2=1`,
`KIC_GF2_REUSE_TABLE=1`, `KIC_GF2_SIMD=1`,
`KIC_GF2_DEFER_ABOVE=0`, and `KIC_GF2_WORD_BATCH=0`.
The reference disables each new option. Use one pinned CPU on a Linux
x86-64 host with AVX2 and BMI2 and `RAYON_NUM_THREADS=1`; do not claim a
Linux gain from Apple ARM64 screening. Record unsupported hosts explicitly.

Warm each arm once, then run five reference/reference A/A pairs and five
rotating reference/candidate pairs per seed. Preserve every process status,
full JSON line, source and binary SHA-256, CPU/affinity, wall and exclusive
phase times, row-space and raw signatures, rank, row/column/pruned counts,
word operations, and output terms. Include failures, timeouts, and OOMs.
The `wall_ms` interval runs from the start of the matrix-F5 call through
unpacking; fixture construction, fingerprints, and process startup are
outside it. Report the median of paired ratios and exact five-pair bootstrap
95% intervals. Report A/A ranges, all holdouts, and smaller-case ratios.

## Decisions

Promote a combined change as the requested further 2× only if the frozen
primary paired median and lower 95% bound are at least 2.00, every holdout
primary median exceeds its A/A maximum, no smaller-case median is below its
own A/A minimum, and all required correctness fields match. If an arm changes
raw echelon rows, require equal canonical row space and explicitly report
the changed representation; do not call it an exact-output gain.
An incremental candidate may be retained opt-in only if its primary paired
median and lower bound exceed 1.03 with the same guards. Archive negative
results and regressions. Stop an arm after a correctness failure or a primary
paired ratio below 0.95 in a complete predeclared screening panel; choose the
next structural hypothesis without a timing-driven rerun of that arm.
