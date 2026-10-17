# Pre-registration: new features for the NIST seed skip filter (2026-10-17)

Written and committed **before** any feature below was computed. Known at this
point: the 21 accepted candidates and their labels (10 published P, 11
skipped S), their curve orders `n` (in
`../nist_seed_provenance_20261016/census/orders_*.txt`), and the feature results
already reported in `../nist_seed_provenance_20261016/README.md` (H1–H6). None of
the features below appears there.

## Population and labels

The 21 candidates that `classify.gp` accepts, labelled P/S as in that README.
Every feature is a function of `(p, n)` only (`t = p + 1 − n`,
`D = t² − 4p`, twist order `n' = 2p + 2 − n`), so `b` and its sign play no
part.

## Hypotheses

- **H7 endomorphism ring.** `Z[π]` maximal at small primes: no `l²` divides
  `D` for odd primes `l < 10^6` (`D ≡ 5 mod 8`, so 2 never divides it).
  Also the count of such `l`.
- **H8 twist security.** `L(n')` = bit length of the largest prime factor of
  the twist order. Binary: `L(n') ≥ 160`; `L(n') ≥ bitlen(p) − 32`;
  `n'` = cofactor < 2^32 times a prime. Continuous: `L(n') / bitlen(p)`.
- **H9 `n ± 1`.** `L(n − 1)`, `L(n + 1)` (same binary thresholds:
  `≥ 160`, `≥ bitlen(p) − 32`) and `/ bitlen(p)` as continuous.
- **H10 a shortcut MOV test.** A 1997 shortcut such as testing
  `p^((n−1)/2) ≠ 1 (mod n)` would reject exactly the curves where `p` is a
  square mod `n`. Binary: `kronecker(p, n) = 1`; `kronecker(−p, n) = 1`;
  `p` a primitive root mod `n`; the embedding degree `k = ord_n(p)` equals
  `n − 1`; `(n − 1)/k ≤ 2`.
- **H11 characters of n and t in the other modulus.** `kronecker(x, n) = 1` for
  `x ∈ {−1, 2, 3, t, n'}` and `kronecker(y, p) = 1` for `y ∈ {n, t, D, n'}`.

Factorisations use PARI/GP `factor` with a per-number time limit of 600 s.
A number that is not fully factored is recorded as unfactored. Its largest
factor is then reported as an interval (the largest proven prime factor, up to
the bit length of the composite cofactor), and a binary threshold on it counts
only when the interval decides it.

## Decision rule (unchanged from the 2026-10-16 protocol)

- A binary feature explains the filter only if it separates P from S
  **perfectly**. By chance this happens with probability ≈ 2·2^−21 ≈ 1e−6 per
  balanced feature.
- Continuous features: exact Mann–Whitney U with a Bonferroni correction over
  the continuous features in this file. A perfectly splitting threshold is
  reported only as a candidate.
- **Negative control.** 100,000 random relabellings (10 P out of 21,
  uniformly). Report the fraction in which at least one binary feature in this
  file splits perfectly. That fraction is the family-wise false-positive rate.
- **Positive control.** A planted label (P := the 10 candidates with the
  smallest `n mod 2^20`) must be detected as a perfect split of a feature
  included for the purpose, with the same code.
- Anything short of a perfect split is reported as "not explained". Partial
  separations are listed, not claimed.
