# Sampled current-weight pivots for one-thread matrix-F5

## Frozen hypothesis and reference

Selecting the lightest *current* row among equal leading-column
candidates roughly halved returned n24 degree-4 terms, but the
full-weight implementation slowed reduction because it recounted
entire rows after block updates. Fusing pivot discovery into strip
clearing recovered some of that cost, yet the complete call remained
slower (`research/gf2_fused_weighted_pivots_20260929/RESULT.md`).
Sorting rows once by initial density produced only 3.5–3.9% fewer
terms (`research/f5_initial_density_order_20260929/RESULT.md`).

This experiment samples eight evenly spaced packed words in each
remaining row, from the current pivot word through the row suffix,
once after each pivot block. The sample score is used only to break
ties among rows with the same leading column. The fused strip pass
discovers the next pivot. Sampling should retain much of the output
sparsity of live weights while avoiding their full-row recounts.
The candidate may change echelon rows, output terms and reduction
word operations. Rank, canonical row space, criterion and row-build
counts must remain exact.

The frozen source reference is main commit
`809fd355c553a77e21bfdae22dcd1099cc1fd6e9`, with SHA-256
`20d2c0933daea0673f2b642226c056f8bdef8ad933614d994737e96923a1d951`
for `src/cryptanalysis/gf2_elim.rs`,
`bf35680f431ab7ac13a7615d334ec9b060217b84bdde36f03ce00319491d42da`
for `src/cryptanalysis/matrix_f5_f2.rs`, and
`d40f9f0c76e5c56c9f4a52596b9edef41f321aa3d70cfb023c4ba21fdde70336`
for `examples/f4_f2_bench.rs`.

## Frozen work and gates

The new arm enables `KIC_GF2_SAMPLED_WEIGHT_PIVOT=1`; the prior arm
sets it to `0`. Both use selective echelon output, fused F5 row
counting, direct packed rows, direct scalar unpack, forced AVX2 row
XOR where available, Gray-code table reuse, and one Rayon thread.
Disable deferred-above, word-batch, AVX-512 unpack and experimental
table builders. Charge sample maintenance, pivot search, reduction,
build and unpack to the complete `wall_ms`. Fingerprints and fixture
construction are outside the timed interval.

First test exact rank and canonical row space on random sparse, dense
and partial-word matrices, along with the seven frozen F5 cases.
For the local ARM64 rejection screen, use seed XORs `0` and
`badc0de1`, one release binary, separate arm processes, one warmup
per arm, five prior/prior A/A pairs and five alternating prior/new
pairs per seed. Preserve every call, failure, timeout, OOM, source
and binary hashes, host and load, output signatures, phase costs,
terms and counted work. A local rejection is sufficient if either
primary seed has a complete-call paired median below 1.05, any
smaller case is below 0.95, or correctness fails. A high-load or
wide-A/A local result is inconclusive for a small speed claim but may
still reject a large regression. This screen cannot establish an
x86 performance gain.

If the local gate passes, repeat with the same seven cases and four
seeds (`0`, `badc0de1`, `5eed2026`, `f5c02a28`) on Linux x86-64
with AVX2/BMI2 and one pinned allowed CPU, five A/A and five
alternating pairs. Report exact five-pair bootstrap 95% intervals,
A/A ranges and every smaller case. The requested further 2× passes
only if the frozen complete-call median and lower bound both reach
2.00, every holdout beats its A/A maximum, no smaller case falls
below its A/A minimum, and exactness passes. An incremental opt-in
remains only if the frozen median and lower bound exceed 1.05 with
the same guards. Otherwise archive the tested patch and receipt and
remove the runtime option. This is a matrix-F5 solver-stage
diagnostic, not an IC online-time or DLP speedup.
