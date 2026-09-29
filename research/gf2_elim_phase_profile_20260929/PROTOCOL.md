# Locate the one-thread GF(2) F5 elimination cost

## Frozen reference and question

The accepted matrix-F5 n24 degree-4 call on x86 spends about 60 ms
in elimination and 57 ms unpacking its full polynomial rows. Several
row-order choices saved output terms but increased elimination time;
AVX2 nibble compaction slowed unpacking. Before another redesign,
measure how elimination time splits among strip loading, pivot
selection and intra-block pivot reduction, Gray-code table build,
and clearing nonpivot rows. The instrumentation must leave the
algorithm, output rows, rank and counted work unchanged. It is a
diagnostic, not a speed candidate.

The frozen reference is main commit
`5e05d9d94a6116627142bd58650944510647d618`, with SHA-256
`20d2c0933daea0673f2b642226c056f8bdef8ad933614d994737e96923a1d951`
for `src/cryptanalysis/gf2_elim.rs`,
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`
for `src/cryptanalysis/matrix_f5_f2.rs`, and
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`
for `examples/f4_f2_bench.rs`. Prior accepted x86 timing is in
`research/gf2_avx2_table_build_20260929/RESULT.md`; its EPYC 7763
host is not paired with any new host in this experiment.

## Frozen workload, checks and decision

Use the seven matrix-F5 cases in `examples/f4_f2_bench.rs`, including
primary `f5_n24_m24_d4`, with seed XORs `0`, `badc0de1`,
`5eed2026`, `f5c02a28`. One release binary runs separate processes
for `KIC_GF2_PROFILE=0|1`. Both modes use selective echelon output,
fused F5 row counting, direct packed rows, direct scalar unpack,
forced AVX2 row XOR, Gray-code table reuse, one Rayon thread, and
disable deferred-above, word-batch, AVX-512 unpack and experimental
table builders. On Linux x86-64 with AVX2/BMI2, pin one allowed CPU.
One warmup per mode, five profile-off/profile-off A/A pairs and five
alternating off/on pairs per seed retain all outputs, failures,
timeouts and OOMs. The profile-on arm emits one `GF2_PHASES` JSON
line for each shared-eliminator call. Keep the complete call and
phase intervals in the receipt, along with source/binary hashes,
host, load, affinity, row fingerprints, rank, output terms and
counted work. Report five-pair A/A and off/on ratios to show the
instrumentation's own perturbation; do not treat them as a speedup.

Check that all seven F5 cases match bit-for-bit raw rows, canonical
row space, rank, criterion/build counts, output terms and reduction
word operations between modes. For the primary case, report median
profiled elapsed time and exclusive strip, pivot, table-build and
row-clear times; report their sum and residual. If the timer overhead
or host contention prevents a stable ordering, label the breakdown
inconclusive. Preserve the instrumentation patch and receipt, then
remove the runtime option from production. Use the largest measured
phase as the next optimization target. No 2× complete-call or IC
online-time claim follows from profiling alone.
