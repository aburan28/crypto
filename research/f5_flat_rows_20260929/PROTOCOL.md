# Contiguous GF(2) rows for the one-thread matrix-F5 call

## Hypothesis and frozen reference

The accepted fast F5 kernel stores each 203-word n24 matrix row in a
separate `Vec<u64>`. On the frozen degree-4 case it builds 6,924 rows and
12,951 columns. The accepted one-thread x86 complete-call median with
direct scalar unpack was 122.96 ms in
`research/f5_further_2x_20260929/RESULT.md`; reduction and unpack are
the major phases. A contiguous row arena may reduce allocator and TLB
costs during the repeated full-row table XORs. Preserve the same pivot
ordering, table algorithm, full `Vec<F2BoolPoly>` output, row fingerprints,
rank and counted word XORs. Copying/flattening the built rows is charged
inside the complete call and reduction phase; it is not hidden setup.

The frozen source reference is main commit
`c369f42d06a4b505cbf2df5dad873e39858b11a7`. SHA-256 hashes are
`20d2c0933daea0673f2b642226c056f8bdef8ad933614d994737e96923a1d951`
for `src/cryptanalysis/gf2_elim.rs`,
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`
for `src/cryptanalysis/matrix_f5_f2.rs`, and
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`
for `examples/f4_f2_bench.rs`.

## Frozen workloads and costs

Run the seven matrix-F5 cases in `examples/f4_f2_bench.rs` with seeds `0`,
`badc0de1`, `5eed2026`, and `f5c02a28`; primary
`f5_n24_m24_d4`. The arms use one binary and separate processes. Both
enable selective echelon output, fused row counting, direct packed rows,
forced AVX2 XOR where available, table reuse, direct scalar unpack, one
Rayon thread, and disable deferred-above and word batching. The only
difference is `KIC_F5_FLAT_ROWS=0|1`. A flat arm must return exact raw rows,
canonical row space, rank, output terms, F5 build/criterion counts and
reduction word XORs. Other matrix solvers keep their existing layout.

First run a nonpromoting Apple ARM64 screen on frozen and holdout A, one
warmup per arm and five alternating pairs per seed. Record every process
status, full output, source/binary digests, host and phase timings. Stop on
any correctness mismatch, or if the local frozen or holdout complete-call
paired median does not reach 1.05. The local screen is not an x86 claim.

If the local screen survives, run the same arms on a Linux x86-64 host with
AVX2 and BMI2 and one pinned allowed CPU. Warm both arms, take five
prior/prior A/A pairs and five alternating prior/flat pairs on each seed.
Report medians and exact five-pair bootstrap 95% intervals for the full
call and phases, along with A/A ranges, all smaller cases, failures,
timeouts and OOMs. `wall_ms` includes criterion, build, conversion,
reduction and full row unpack; it excludes fixture construction, process
launch and fingerprinting.

The requested further 2× passes only if the frozen full-call prior/flat
paired median and lower interval bound both reach 2.00, every holdout
primary exceeds its A/A maximum, no smaller case falls below its A/A
minimum, and all exactness checks pass. An incremental opt-in may remain
only if the frozen median and lower bound exceed 1.05 with the same
correctness, holdout and smaller-case gates. Otherwise archive the
negative result and remove the option. This is a matrix-F5 solver-stage
diagnostic, not an IC online-time or DLP speedup.
