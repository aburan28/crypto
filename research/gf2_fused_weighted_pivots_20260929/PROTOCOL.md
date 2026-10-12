# Fuse weighted pivot discovery into the GF(2) strip update

## Hypothesis and frozen reference

The current-weight pivot experiment in
`research/gf2_minweight_pivots_20260929/RESULT.md` cut the n24 output term
count roughly in half while preserving rank and canonical row space, but
slowed the complete call. Its implementation scans every remaining strip
row once to find each pivot and again to clear that pivot in the strip.
Select the next pivot during the strip-clearing pass that the packed kernel
already needs. This removes one full scan per pivot without changing the
weighted choice, its output rows, or counted row XORs. It does not remove
the row popcounts after table updates; those remain charged to reduction.

The frozen reference is main commit
`c9d3b9525355face849c59d954c1864923cdbc43`, with source SHA-256
`20d2c0933daea0673f2b642226c056f8bdef8ad933614d994737e96923a1d951`
for `src/cryptanalysis/gf2_elim.rs`,
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`
for `src/cryptanalysis/matrix_f5_f2.rs`, and
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`
for `examples/f4_f2_bench.rs`. The weighted arm applies the archived
`rejected_candidate.patch` from the prior experiment unchanged. The fused
arm adds only the discovery fusion. The accepted one-thread fast reference
is the full-output mode measured at 125.89 ms on its eligible x86 host in
`research/gf2_table_reuse_20260929/RESULT.md`; do not compare that absolute
time directly to another host.

## Workloads and cost

Use the seven matrix-F5 cases in `examples/f4_f2_bench.rs`, seeds `0`,
`badc0de1`, `5eed2026`, and `f5c02a28`, primary `f5_n24_m24_d4`. All arms
use selective echelon output, fused build, direct packed rows, AVX2 XOR
where supported, table reuse, direct scalar unpack, no deferred-above or
word-batch option, and `RAYON_NUM_THREADS=1`. Three arms in one binary are:
`fast` (ordinary pivot choice), `weighted` (prior complete current-weight
choice), and `fused` (same current-weight choice with next-pivot discovery
inside strip clearing). Use explicit environment flags to select arms.

First run an Apple ARM64 nonpromoting screen on frozen and holdout A, one
warmup per arm and three rotating triads per seed. Preserve every process
status and full output, source/binary hashes, host data, phases, term counts,
row signatures and counted operations. Stop if the fused arm differs from
the weighted arm in rank, raw rows, canonical row space, term count or
counted work, or if its paired complete-call median against `fast` remains
below 0.95 on either seed. A local positive screen is not an x86 claim.

If it survives, run on a Linux x86-64 host with AVX2 and BMI2, one pinned
allowed CPU and all four seeds. Warm all arms once, take five `fast`/`fast`
A/A pairs and five rotating three-arm sets per seed. Report paired ratio
medians with exact five-pair bootstrap 95% intervals and A/A ranges. Check
all smaller cells. The measured `wall_ms` covers criterion, row building,
reduction and full `Vec<F2BoolPoly>` unpacking. Exclude fixture construction,
process startup and fingerprinting. Preserve failures, timeouts and OOMs.

The requested further 2× passes only if the frozen fast/fused complete-call
paired median and lower interval bound both reach 2.00, each holdout primary
median exceeds its A/A maximum, no smaller case regresses below its own
A/A minimum, and all correctness checks pass. An incremental opt-in may
remain only if weighted/fused median and lower bound exceed 1.03 and the
fast/fused primary and smaller-case guards pass. Otherwise archive the
negative result and remove the experimental runtime option. Any result is
a matrix-F5 solver-stage diagnostic, not an IC online-time or DLP speedup.
