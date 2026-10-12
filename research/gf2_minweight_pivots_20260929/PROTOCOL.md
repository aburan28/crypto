# Current-weight pivots in the packed GF(2) F5 kernel

## Hypothesis and reference

The accepted one-thread fast F5 call was 125.89 ms on the frozen
`f5_n24_m24_d4` case in `research/gf2_table_reuse_20260929/RESULT.md`.
Direct scalar unpack subsequently reached 1.048× against that fast mode in
`research/f5_further_2x_20260929/RESULT.md`, with unchanged rank, raw rows,
and row space. The remaining cost on the latter x86 runner was about 60 ms
reduction and 52 ms unpack. A sparse leading-band experiment produced only
5.10 million output terms versus 13.73 million in the packed kernel, but
spent 4.25 billion term visits and was about 110× slower. Selecting the
lightest current pivot row inside the existing packed block kernel may
capture part of that fill reduction while retaining table updates.

The frozen source reference is main commit
`45c3114d5d5da81fb3daae7ea47e50b749543fa0`. Its SHA-256 digests are
`20d2c0933daea0673f2b642226c056f8bdef8ad933614d994737e96923a1d951`
for `src/cryptanalysis/gf2_elim.rs`,
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`
for `src/cryptanalysis/matrix_f5_f2.rs`, and
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`
for `examples/f4_f2_bench.rs`.

Maintain one weight per unpivoted packed row. Initialize it with the full
row popcount; after each table application, update the row's weight from the
modified row, and swap weights with rows during pivot swaps. At a candidate
column, inspect all eligible strip rows and choose the lowest-weight row
among ties. The current word may have been reduced inside the block while
the cached weight is from the previous block; this is a deterministic
approximation. Charge weight maintenance inside reduction timing. The
output remains full `Vec<F2BoolPoly>` echelon rows and an exact rank. Raw
rows, output terms, and reduction operation counts may change; canonical
row space, rank, and the criterion/build counts must match.

## Fixed local screen and eligible-host run

Use the seven F5 cases in `examples/f4_f2_bench.rs` with seeds `0`,
`badc0de1`, `5eed2026`, and `f5c02a28`; primary `f5_n24_m24_d4`. Each
benchmark process emits all seven cases. All arms use selective echelon,
fused build, direct packed rows, AVX2 XOR where supported, table reuse,
`KIC_GF2_DEFER_ABOVE=0`, `KIC_GF2_WORD_BATCH=0`, and one Rayon thread.
The three arms are (0) accepted fast mode with ordinary scalar unpack,
(1) the same mode with direct scalar unpack, and (2) direct unpack plus
`KIC_GF2_MIN_WEIGHT_PIVOT=1`. Keep all other settings equal. The A/A
reference is arm 0.

First run a nonpromoting Apple ARM64 screen on frozen and holdout A with one
warmup and three rotating triads per seed. Preserve all process outputs,
status, source and binary hashes, host data, phases, row signatures, term
counts and operations. Stop this hypothesis after any correctness failure,
or after a complete screen whose arm-1/arm-2 primary paired median is below
0.95 on either seed. A positive ARM64 screen is not an x86 speedup claim.

If the screen survives, run the same binary as separate arm processes on a
Linux x86-64 host with AVX2 and BMI2, one pinned allowed CPU, all four seeds,
one warmup per arm, five arm-0/arm-0 A/A pairs, and five rotating triads per
seed. Record every failure, timeout and OOM. `wall_ms` is the full call from
criterion through unpack; fixture construction, process startup and
fingerprinting are outside it. Report medians of paired ratios with exact
five-pair bootstrap 95% intervals, A/A ranges, all holdouts and all smaller
cases. Preserve raw JSON and source/binary hashes in the PR.

## Gate

The requested **further 2×** passes only if frozen arm-0/arm-2 complete-call
median and lower interval bound both reach 2.00, each holdout primary median
exceeds its A/A maximum, no smaller-case arm-0/arm-2 median falls below its
own A/A minimum, and every row-space/rank/count check passes. An incremental
candidate may remain opt-in only if frozen arm-1/arm-2 median and lower bound
exceed 1.03 with the same holdout and smaller-case guards. Otherwise archive
the negative result and remove the experimental option. Report any changed
raw rows and terms explicitly. This is a matrix-F5 solver-stage diagnostic,
not an IC online-time or DLP speedup.
