# Bound the sampled-weight pivot search to eight ties

## Frozen hypothesis and reference

The sampled current-weight pivot experiment cut n24 matrix-F5 output
from 13.7 million to about 7.5 million terms, but its full-strip
search increased elimination time and slowed the complete call
(`research/gf2_sampled_pivot_20260929/RESULT.md`). This experiment
scores at most the first **eight** rows whose current strip has the
lowest available pivot bit. It then selects the lightest sampled row
among those eight. Each score counts eight evenly spaced packed
suffix words, refreshed once after each pivot block. If the current
lowest bit has no candidate, search the full strip to find the next
available pivot bit. The bounded search should preserve much of the
output sparsity while retaining the early exit that made the accepted
pivot search fast. The candidate may change echelon rows, output term
counts and counted reduction work. Rank, canonical row space,
criterion and row-build counts must remain exact.

The frozen source reference is main commit
`f75503bb1b97415ff5de2d3698e43bfb046dac5b`, with SHA-256
`20d2c0933daea0673f2b642226c056f8bdef8ad933614d994737e96923a1d951`
for `src/cryptanalysis/gf2_elim.rs`,
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`
for `src/cryptanalysis/matrix_f5_f2.rs`, and
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`
for `examples/f4_f2_bench.rs`.

## Frozen work and gates

The new arm sets `KIC_GF2_BOUNDED_WEIGHT_PIVOT=1`; prior sets `0`.
Both use selective echelon output, fused F5 row counting, direct
packed rows, direct scalar unpack, forced AVX2 row XOR where
available, Gray-code table reuse and one Rayon thread. Disable
deferred-above, word-batch, AVX-512 unpack and experimental table
builders. Charge score refresh, pivot search, reduction, build and
unpack to the complete `wall_ms`. Fingerprints and fixture
construction are outside the timed interval.

First test rank and canonical row space on random sparse, dense and
partial-word matrices, and verify the seven F5 cases. For the local
ARM64 rejection screen, use seeds `0` and `badc0de1`, one release
binary, separate arm processes, one warmup per arm, five prior/prior
A/A pairs and five alternating prior/new pairs per seed. Preserve
every call, failure, timeout and OOM, source/binary hashes, host and
load, output signatures, phase costs, terms and counted work.
A local rejection is sufficient if either primary complete-call
paired median is below 1.05, any smaller case falls below 0.95, or
correctness fails. Wide A/A ranges make a small local gain
inconclusive, not an eligible x86 speed claim.

If the local gate passes, repeat the seven cases on Linux x86-64
with AVX2/BMI2 and one pinned allowed CPU on four seeds (`0`,
`badc0de1`, `5eed2026`, `f5c02a28`), five A/A and five alternating
pairs. Report exact five-pair bootstrap 95% intervals, A/A ranges
and all smaller cases. The requested further 2× passes only if the
frozen complete-call median and lower bound both reach 2.00,
every holdout beats its A/A maximum, no smaller case falls below its
A/A minimum, and exactness passes. An incremental opt-in remains
only if the frozen median and lower bound exceed 1.05 with the same
guards. Otherwise archive the tested patch and receipt and remove
the runtime option. This is a matrix-F5 solver-stage diagnostic,
not an IC online-time or DLP speedup.
