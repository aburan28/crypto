# Split matrix-F5 pivot work on one x86 thread

## Frozen reference and question

The pinned EPYC 7763 profile in
`research/gf2_elim_phase_profile_20260929/RESULT.md` found that
the n24 degree-4 shared GF(2) eliminator spends about 27.4 ms of
59.3 ms in its pivot loop, more than in row clearing (22.0 ms) or
Gray-code table construction (8.0 ms). That pivot bucket includes
searching the strip, XORing the selected pivot row with earlier
pivots in the same block, XORing earlier pivots by the new row, and
clearing the pivot bit from the strip. Measure these four exclusive
subphases before changing the kernel. Instrumentation must leave
all rows, rank and counted work unchanged; it is a diagnostic,
not a speed candidate.

The frozen reference is main commit
`ce79b74645c8e564ec361918dfb72bd990cc56d6`, with SHA-256
`20d2c0933daea0673f2b642226c056f8bdef8ad933614d994737e96923a1d951`
for `src/cryptanalysis/gf2_elim.rs`,
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`
for `src/cryptanalysis/matrix_f5_f2.rs`, and
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`
for `examples/f4_f2_bench.rs`.

## Frozen workload, checks and decision

Use the seven F5 cases in `examples/f4_f2_bench.rs`, including
primary `f5_n24_m24_d4`, with seed XORs `0`, `badc0de1`,
`5eed2026`, `f5c02a28`. One release binary runs separate processes
for broad-only (`KIC_GF2_PROFILE=1`,
`KIC_GF2_PIVOT_PROFILE=0`) and fine-pivot (`1`, `1`) modes.
Both modes use selective echelon output, fused F5 row counting,
direct packed rows, direct scalar unpack, forced AVX2 row XOR,
Gray-code table reuse, one Rayon thread, and disable deferred-above,
word-batch, AVX-512 unpack and experimental table builders. On
Linux x86-64 with AVX2/BMI2, pin one allowed CPU. Run one warmup
per mode, five broad/broad A/A pairs and five alternating broad/fine
pairs per seed. Preserve every output, failure, timeout and OOM,
source/binary hashes, host, load and affinity, row fingerprints,
rank, terms and counted work. Report the A/A and broad/fine ratios
only to characterize timer perturbation.

All seven F5 cases must match raw and canonical row fingerprints,
rank, criterion/build counts, output terms and reduction word ops
between modes. For the primary, report median fine-profiled pivot
time and exclusive scan, selected-row XOR, earlier-pivot XOR, and
strip-clear times; report their sum and residual. If the timer
overhead or host contention prevents a stable ordering, label the
breakdown inconclusive. Preserve the patch and receipt, then remove
profiling from production. Choose the largest measured subphase as
the next implementation target. No 2× complete-call or IC
online-time claim follows from profiling alone.
