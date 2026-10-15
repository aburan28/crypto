# Item 4 protocol: full Fourier spectrum of decompositions on the Kummer tower

Frozen 2026-10-14, before any data of this item was collected.
Design source: `reports/Index calculus leads next steps.md`, "### 4.".
Builds on `research/ic_exponent_falsifiers_20261013/` (C4, plain exponent sum,
falsified) and requirements-document candidate C2 (sparse-spectrum zero finding).

## Prior-work check (5 minutes)

- `research/ic_exponent_falsifiers_20261013/README.md`: C4 ran one character
  family only (`(1,1,1)·N/2^k`, the plain sum mod `2^k`); "C2 ... Not run yet";
  weighted forms untested.
- `research/prime_fourier_flatness_exp7_20261007/`: a Fourier flatness sweep
  of coordinate predicates over `Z/n` indexed by discrete logs. Different object
  (spectrum over the log group, not over the Kummer exponent group of a
  decomposition), so not a duplicate.
- `git log` (last 300 commits) and `research/`, `benchmarks/`: no other
  Fourier/spectrum test on Kummer-tower decompositions.

## Hypothesis

H1 (C2 alive): for some character `χ_c(e) = ω_N^{c·e}` of the exponent group
`(Z/N)^m`, `N = 2^t`, the value `c·e mod N` of the exponent tuples of the
decompositions of a target is concentrated (or anti-concentrated) given the
target, i.e. the target-conditioned spectrum of the decomposition indicator
has a coefficient far above what an unstructured base gives.
H0 (C2 closed on these toys): the spectrum is flat, indistinguishable from
the same statistic on relabelled and on random factor bases.

## Object and statistic

Base: PKM's Kummer tower `V = μ_N ⊂ F_p`, abscissa `ζ^e` labelled `e ∈ Z/N`,
both signs (ecbench `prime-kummer:t=T`). For an `m`-term decomposition
`R = ±P_1 ± … ± P_m` the x-only indicator of target `x(R)` is
`f_R(e_1..e_m) = 1` iff `S_{m+1}(ζ^{e_1},…,ζ^{e_m}, x(R)) = 0`; it is
`S_m`-symmetric, so it is a union of orbits `O` (one per unordered exponent
multiset).

Enumeration is exhaustive: every multiset of `m` signed base points is
summed once, so every target's full decomposition set is known (no target
sampling). Consequence, stated before the data: the *unconditional* multiset
of all decomposition tuples is then exactly the set of all base `m`-multisets,
whose spectrum is fixed by the base's marginal (which labels are on the
curve) and carries no target information. The informative statistic is the
**target-conditioned cross-orbit spectrum**

  `Â(c) = Σ_R Σ_{O≠O' ⊂ Z_R} f̂_O(c) · conj f̂_O'(c)`,  `f̂_O(c) = Σ_{σ∈S_m} ω^{-c·σ(e_O)}`,

computed for all `N^m` characters at once by one FFT of the symmetrised
cross-pair difference histogram. Under H0 (independent uniform exponents)
`E Â(c) = 0` and `Var Â(c) = 2·|Stab_{S_m}(c)|²·(m!)²·P` (`P` = number of
unordered cross-orbit pairs; ×2 for real characters `2c = 0`), so the
normalised coefficient is `Z(c) = Â(c)/sqrt(Var)`. A planted relation
`c·e ≡ s(R)` makes `Z(c)` large and positive. Statistics, over the
canonical classes of nonzero characters (up to `S_m` and `c ↦ −c`):

- `maxZ = max |Z(c)|` (two-sided) and the arg-max class;
- `Z_C4 = max over k of |Z((1,…,1)·N/2^k)|` (the C4 family, for continuity);
- top-k mass `M_k = Σ top-k Z² / Σ Z²`, `k = 1, 10, 100`;
- marginal: max normalised coefficient of the on-curve label indicator on `Z/N`.

## Arms

1. Tower (observed).
2. Permutation null, 200 per cell: the same decompositions with the labels
   relabelled by a uniform random bijection of `Z/N` (preserves every target,
   orbit count and pair count; destroys multiplicative meaning of labels).
3. Random-factor-base null, 100 per cell: `N` distinct uniform abscissas of
   `F_p^*`, labelled `0..N−1` in random order, same curve, same exhaustive
   enumeration, same statistic.
4. Positive control (planted bias in one character), on curve 0 of each
   (size, m): the tower's targets and orbit counts kept, every orbit's tuple
   replaced by a uniform tuple except that with probability `ρ` it is drawn on
   `c*·e ≡ s_R (mod N)`, `s_R` uniform per target; `c* = (1,3)` (m = 2),
   `(1,3,5)` (m = 3), i.e. a weighted linear form; `ρ ∈ {1, 0.3, 0.1, 0.03}`.
   Read against the cell's permutation-null threshold.

## Sizes

Primes: smallest `p > 2^b`, `p ≡ 1 mod 4096`: `p20 = 1073153`, `p22 = 4206593`,
`p24 = 16801793`. Five curves per prime, `y² = x³ + ax + b`, `(a, b)` uniform
from `numpy.random.default_rng(20261014000 + 100·size_index + curve_index)`
(`size_index` 0, 1, 2 for 20, 22, 24 bits), redrawn until `4a³ + 27b² ≠ 0`.

- `m = 2`, `t = 11` (`N² = 2^22` characters), all three primes.
- `m = 3`, `t = 7` (`N³ = 2^21` characters), all three primes.
- Continuity cell: the C4 curve `p = 946177, a = 304973, b = 530845`,
  `m = 3`, `t = 8` (`2^24` characters), same arms (200 perm, 100 random base).

The report's suggestion `t = 2,3,4` is not used: with `N ≤ 16` almost no
target has two decompositions at `p ≥ 2^20`, so the conditioned spectrum is
empty; `t` is chosen so that targets with ≥ 2 orbits are plentiful and
`N^m ≤ 2^24` fits the 700 MB cap. Targets: all of them (exhaustive), which is
≥ the design's 10^4–10^5 informative targets in every cell.

## Kill criterion (from the report) and decision

- **NEGATIVE (C2 closed on these toys)**: in every cell the tower's `maxZ` is
  ≤ the permutation null's 99.9th percentile (empirical quantile of the 200
  null maxima, linear interpolation), and `M_1, M_10, M_100` are inside the
  null's [0.1, 99.9] percentile band; and the positive control at `ρ = 1`
  exceeds the threshold in every (size, m).
- **POSITIVE**: some cell exceeds, and the exceedance survives a replay with
  2000 fresh permutation nulls at Bonferroni level `p < 0.05/31` and ranks
  above the random-base nulls.
- An exceedance that does not survive that replay is reported as a
  multiplicity false alarm (31 cells × several statistics are tested).
- **INCONCLUSIVE** if the positive control fails at `ρ = 1` anywhere.
- Sensitivity is reported as the smallest planted `ρ` detected per (size, m).

## Instruments

`examples/ic_next10_04_kummer_fourier.rs` (exhaustive enumeration, binary
output of targets with ≥ 2 orbits); `analysis/spectrum.py` (histogram, FFT,
nulls, controls). Single-threaded, `nice`, `prlimit --as=700MB`.
