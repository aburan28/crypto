# Probe 32 spread-out rows for a lighter F5 pivot

## Frozen hypothesis and reference

Full-strip sampled-weight pivot search nearly halved n24 matrix-F5
output terms but slowed reduction; scoring the first eight matching
rows kept reduction nearer the accepted path but saved only 1–2% of
terms (`research/gf2_sampled_pivot_20260929/RESULT.md`,
`research/gf2_bounded_pivot_20260929/RESULT.md`). This candidate
first uses the accepted search to find the lowest available pivot
column. When that column is the current lowest undecided bit, it
then probes 32 positions spread evenly across the remaining row
range: position `pivot_row + floor(j * remaining_rows / probes)` for
`j = 0..probes-1`, where `probes = min(32, remaining_rows)`.
Among probed rows carrying that pivot bit, choose the one with the
lowest cached score. The score counts eight evenly spaced packed
suffix words and is refreshed once after each pivot block. These
fixed probes should consider a wider set of rows than the first
eight while bounding selection cost. If the lowest undecided bit
has no row, retain the accepted full search and its selected pivot.
The candidate may change echelon rows, term counts and counted
reduction work. Rank, canonical row space, criterion and row-build
counts must remain exact.

The frozen reference is main commit
`4ccaf6e56edf78656d4af1dbed9c5fa15200c395`, with SHA-256
`20d2c0933daea0673f2b642226c056f8bdef8ad933614d994737e96923a1d951`
for `src/cryptanalysis/gf2_elim.rs`,
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`
for `src/cryptanalysis/matrix_f5_f2.rs`, and
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`
for `examples/f4_f2_bench.rs`.

## Frozen work and gates

The new arm sets `KIC_GF2_STRATIFIED_WEIGHT_PIVOT=1`; prior sets
`0`. Both use selective echelon output, fused F5 row counting,
direct packed rows, direct scalar unpack, forced AVX2 row XOR where
available, Gray-code table reuse and one Rayon thread. Disable
deferred-above, word-batch, AVX-512 unpack and experimental table
builders. Charge score refresh, pivot search, reduction, build and
unpack to the complete `wall_ms`. Fingerprints and fixture
construction are outside the timed interval.

First test rank and canonical row space on random sparse, dense and
partial-word matrices, then all seven F5 cases. The local ARM64
rejection screen uses seeds `0` and `badc0de1`, one release binary,
separate arm processes, one warmup per arm, five prior/prior A/A
pairs and five alternating prior/new pairs per seed. Preserve every
call, failure, timeout, OOM, source/binary hashes, host and load,
output signatures, phase costs, terms and counted work. A local
rejection is sufficient if either primary complete-call paired
median is below 1.05, any smaller case below 0.95, or correctness
fails. Wide A/A ranges make a small local gain inconclusive.

If the local gate passes, repeat the seven cases on Linux x86-64
with AVX2/BMI2 and one pinned CPU on four seeds (`0`, `badc0de1`,
`5eed2026`, `f5c02a28`), five A/A and five alternating pairs.
Report exact five-pair bootstrap 95% intervals, A/A ranges and all
smaller cases. The requested further 2× passes only if the frozen
complete-call median and lower bound both reach 2.00, every holdout
beats its A/A maximum, no smaller case falls below its A/A minimum,
and exactness passes. An incremental opt-in remains only if the
frozen median and lower bound exceed 1.05 with the same guards.
Otherwise archive the tested patch and receipt and remove the
runtime option. This is a matrix-F5 solver-stage diagnostic,
not an IC online-time or DLP speedup.
