# EXP7: Fourier flatness of cheap prime-field point predicates

**Frozen before execution:** 2026-10-07. **Baseline source:**
`aburan28/crypto@30da958ed337fb3d7522d7fd691aa15e6c733a97`.
This is a diagnostic on toy prime-order elliptic curves, not a discrete-log
solver or a claim about cryptographic-size curves.

## Question and hypothesis

For each predicate on the affine coordinates of `P = [k]G`, enumerate the
entire cyclic group, form its membership indicator `f(k)`, and inspect every
nontrivial Fourier coefficient on `Z/n`. The null hypothesis is that the
maximum coefficient is distributed like the maximum for a uniformly sampled
set of the same size with the same `k <-> -k` symmetry. The prediction is no
adjusted excess for a cheap coordinate predicate; a log interval must give a
clear positive control. This sweep also computes the full pair-sum count
vector from the Fourier coefficients, so it retains the information needed
for the proposed pair-sum experiment.

## Frozen panel and membership rules

Use `ic_boundary::find_prime_order_curve(bits, 20261007)` at `bits =
8, 10, 12, 14`, in that order. This repository constructor searches for a
nonsingular short-Weierstrass curve of prime order and a generator using
`StdRng`; the code revision and seed fix its output. Record `p, a, b, n,
G, cofactor, ICV1`. Enumerate `k = 0,...,n-1` by repeated point addition,
verify that every point is distinct and `[n]G = O`, and exclude `O` from all
sets. The predicates, chosen without looking at spectra, are:

| Case | Membership for affine point with integer `0 <= x < p` |
| --- | --- |
| smooth | `x > 0` and every prime factor of `x` is at most `B = max(2, floor(p^(1/4)))` |
| cf16 | `x > 0` and **every** nonzero partial quotient of the finite simple continued fraction of `x/p` is at most 16 |
| cantor3 | every base-3 digit of `x` is 0 or 2 |
| hamming | `popcount(x) <= floor((bit_length(p)-1)/2)` |
| farey | there exist `1 <= b <= H`, `-H <= a <= H` with `x*b = a (mod p)`, where `H = max(2, floor(sqrt(p)/3))` |
| legendre+++ | each of `x, x+1, x+2` has Legendre symbol `+1` modulo `p` |
| sha-x | first two bits of SHA-256 of `"exp7-sha-x-v1" || p_le64 || x_le64` are zero |

`random-pairs` is one seeded uniform set of `floor((n-1)/8)` inversion
pairs, independent of coordinates. `log-interval` is `1 <= k <=
floor((n-1)/4)` and is the deliberately log-aware positive control. Both
are frozen controls; neither is a cheap coordinate predicate. Cases with
membership 0 or `n` are retained as degenerate, with no normalized statistic.

## Spectrum, nulls, and decision

Compute `F(j) = sum_{k=0}^{n-1} f(k) exp(-2*pi*i*j*k/n)` for every `j`.
At `j != 0`, centering `f` has no effect. Record the largest `|F(j)|`
over `1 <= j <= floor(n/2)`, its index, and the normalized peak
`Z = |F(j)| / sqrt(m*(1-m/n))`, where `m = |S|`. Use an arbitrary-length
DFT (Bluestein over radix-2 FFT), with direct DFT checks on small vectors.
The inverse DFT of `F(j)^2` gives the ordered-pair count at each log index;
record its maximum absolute deviation from `m^2/n`, and check a direct
convolution on small test cases.

For each case on each curve, draw **512** independent null sets of exactly
`m` members. For coordinate cases, `sha-x`, and `random-pairs`, sample
`m/2` of the `(n-1)/2` inversion pairs uniformly without replacement.
For `log-interval`, sample `m` of `1,...,n-1` without replacement. Use
`StdRng` seeded by SHA-256 of `"exp7-null-v1" || p_le64 || case_name ||
replicate_le64`; the one `random-pairs` case uses replicate 512, so it is
disjoint from its calibration ensemble. Nulls preserve `m` and the
inversion symmetry forced by x-coordinate membership. Record null median,
95th percentile, maximum, and the conservative empirical peak p-value
`(1 + #{null peak >= observed peak})/513` for every case.

The predeclared candidate family contains the six non-hash coordinate cases
on four curves, hence **24** cells. A lead requires the candidate's empirical
p-value to be at most `0.05/24` (Bonferroni). Keep any lead open for an
independent seed and larger-null replay; a toy peak alone is no DLP result.
Check that `log-interval` exceeds its own null 95th percentile on all four
curves; otherwise diagnose the spectral pipeline. `sha-x` and
`random-pairs` must be reported even if they are outliers. Do not discard
failed, sparse, or outlying cases.

## Execution, evidence, and limits

Implement in Rust as `examples/exp7_fourier_flatness.rs`, using the existing
prime-curve arithmetic and SHA-256 implementation. Use `cargo test --locked
--release --example exp7_fourier_flatness` for DFT, group-enumeration, and
predicate checks. Run `cargo run --locked --release --example
exp7_fourier_flatness -- research/prime_fourier_flatness_exp7_20261007/results.json`.
Freeze the code in a commit before the measured sweep. Commit the exact
JSON output, source revision, command, curve records, analysis, diagram,
quantitative visual, and PDF. CPU wall time is only operational; this
experiment makes no speed claim and needs no timed-ratio comparison.

The work is `O(n log n)` per spectrum plus the complete `n`-point
enumeration, and it begins with **known** `k` labels. Charge this enumeration
and the 512-null calibration to the diagnostic. The Fourier test is not a
cheap unknown-log oracle: learning the memberships in log order already
requires the toy enumeration. Character-sum bounds and hardcore-bit
reductions have hypotheses not established by this screen; do not infer
universal flatness or a DLP break from its outcomes.
