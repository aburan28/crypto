# Toward a Gram–Schmidt Analysis of the secp256k1 LLL Degeneracy

A focused theoretical analysis attempting to derive — from first
principles — why the Boneh–Venkatesan lattice for secp256k1 fails
LLL where P-256, secp192k1, and secp224k1 succeed. The cross-
Koblitz empirical refutation of the near-power-of-2 hypothesis
(see [`RESEARCH_GLV_HNP.md`](RESEARCH_GLV_HNP.md) Phase 1.5
status) has narrowed the possibilities but not closed them; this
note derives what we can derive and identifies what's still open.

> **Status**: theoretical analysis sketch. The derivation goes as
> far as the structure of the Gram matrix and identifies the
> candidate source of degeneracy, but the precise eigenvalue
> distribution of the Gram matrix at secp256k1 vs. generic curves
> is left as an open computation.

## 1. The Boneh–Venkatesan basis

For an HNP attack on `m` biased ECDSA signatures with parameters
`(a_i, t_i)` and bias bound `B = 2^{k_{\text{bits}}}`, the basis
matrix is

```
       ┌─────────────────────────┬────────┬────────┐
       │ n² · I_m                │  0     │  0     │
B  =   ├─────────────────────────┼────────┼────────┤
       │ n · a_1, ..., n · a_m   │  B     │  0     │
       ├─────────────────────────┼────────┼────────┤
       │ n · t_1, ..., n · t_m   │  0     │ n · B  │
       └─────────────────────────┴────────┴────────┘
```

(dimensions `(m+2) × (m+2)`).

This is the implementation in
`crypto/src/cryptanalysis/hnp_ecdsa.rs::hnp_recover_key_with_reduction`.

The basis is integer-valued with entries of bit-length:
- Diagonal of upper block: `2 · n_{\text{bits}}` bits (the `n²` terms)
- Off-diagonal entries: `n_{\text{bits}} + \log_2(a_i)` bits (≈ `2 n_{\text{bits}}` since `a_i \in [0, n)`)
- `B` slot: `k_{\text{bits}}` bits ($\approx 0.75 n_{\text{bits}}$)
- `n · B` slot: `n_{\text{bits}} + k_{\text{bits}}` bits

So all entries have bit-length between `k_{\text{bits}}$` and `2 n_{\text{bits}}`. For secp256k1 with `n_{\text{bits}} = 256`, this is between 192 and 512 bits.

## 2. The Gram matrix

LLL's convergence depends on the Gram matrix `G = B B^T` (entries
are inner products of basis rows). For our `B`:

- `G_{ii}` for `i ≤ m`: `n^4` (= `(n²)²`)
- `G_{ij}` for `i, j ≤ m`, `i ≠ j`: `0` (diagonal block has no off-diagonal)
- `G_{i, m+1}` for `i ≤ m`: `n^2 · n · a_i = n^3 · a_i`
- `G_{i, m+2}` for `i ≤ m`: `n^3 · t_i`
- `G_{m+1, m+1}`: `(n · a_1)^2 + (n · a_2)^2 + ... + (n · a_m)^2 + B^2`
  `≈ m · n^2 · ⟨a²⟩ + B²` where ⟨a²⟩ is the avg `a_i^2`
- `G_{m+1, m+2}`: `n^2 · ⟨a · t⟩`
- `G_{m+2, m+2}`: `m · n^2 · ⟨t²⟩ + n^2 · B^2`

For generic `a_i, t_i` (uniform in `[0, n)`), we expect
`⟨a^2⟩ ≈ n²/3`, `⟨t²⟩ ≈ n²/3`, `⟨a · t⟩ ≈ n²/4`. So:

- `G_{m+1, m+1} ≈ m · n^4 / 3 + B^2 ≈ m · n^4 / 3` (since `B² ≪ n^4`)
- `G_{m+1, m+2} ≈ m · n^4 / 4`
- `G_{m+2, m+2} ≈ m · n^4 / 3 + n^2 · B^2`

The Gram matrix has roughly the structure:

```
       ┌─────────────────┬──────────────┬──────────────┐
       │ n^4 · I_m       │  ~n^3 · a    │  ~n^3 · t    │
G  ≈   ├─────────────────┼──────────────┼──────────────┤
       │ ~n^3 · a^T      │  m n^4/3     │  m n^4/4     │
       ├─────────────────┼──────────────┼──────────────┤
       │ ~n^3 · t^T      │  m n^4/4     │  m n^4/3     │
       └─────────────────┴──────────────┴──────────────┘
```

## 3. Eigenvalue gap

LLL converges in `O(\text{dim}^3 \log B)` iterations *when* the
Gram matrix's eigenvalues are well-separated. The Lovász condition
`δ · μ_{i-1}^2 ≤ μ_i^2 / μ_i^2` requires successive ratios of
orthogonalisation lengths to stay bounded.

For our `G`, the eigenvalues split into roughly three groups:

1. The `n^4 · I_m` block contributes `m` near-degenerate
   eigenvalues at `n^4`.
2. The 2×2 block in the lower-right has two eigenvalues governed
   by `m · n^4 / 3 ± m · n^4 / 4 = 7 m n^4 / 12` and `m n^4 / 12`.
3. Coupling terms (`n^3 · a^T`, `n^3 · t^T`) perturb these by
   `O(n^6 · m)` corrections.

The ratio between the largest and smallest eigenvalue at generic
`(a, t)` is `O(m)` — small, well-conditioned.

## 4. The candidate degeneracy mechanism

LLL's per-iteration cost is dominated by **size reduction**: for
each pair `(i, j)` with `i > j`, compute the Gram–Schmidt
coefficient `μ_{ij} = ⟨b_i^*, b_j⟩ / ‖b_j^*‖²` and reduce `b_i`
by `⌊μ_{ij}⌉ · b_j`.

If a `μ_{ij}` is "just barely" larger than `1/2` (the size-
reduction threshold), `⌊μ_{ij}⌉` rounds to 1, the reduction happens,
then in a subsequent iteration `μ_{ij}` becomes "just barely"
smaller than `−1/2`, leading to another reduction. This
oscillation can take many iterations before terminating, even
when the basis is provably LLL-reducible.

**Hypothesis**: for secp256k1's specific `(a, t)` distribution
(which is curve-specific because the `a_i, t_i` are derived from
ECDSA signing with that curve's $n$, $G$, and signing equation),
the Gram–Schmidt coefficients land in this oscillation regime
disproportionately often, while for P-256, secp192k1, secp224k1
they don't.

## 5. Why is secp256k1 different?

The `(a, t)` values come from:

```
   a_i = s_i^{-1} · r_i  (mod n)
   t_i = s_i^{-1} · z_i  (mod n)
```

where `r_i = (k_i · G).x mod n` and `s_i = k_i^{-1} (z_i + d · r_i) mod n`.

For different curves, the `(r_i, s_i, z_i)` distributions interact
differently with `n`. In particular, `r_i = (k_i · G).x mod n`
depends on `n` via the final modular reduction — and if `n` is
near `p` (as it is for secp256k1, where `n_{\text{secp}}` is
strikingly close to `2^{256} - 2^{32}` minus small terms), the
reduction `mod n` rarely changes the value of `r_i`, leaving
`r_i ≈ (k_i · G).x` directly. For P-256, `n` is further from
`p` (different bit-pattern), so the reduction more often
non-trivially adjusts `r_i`.

**This** might be the source: secp256k1's `r_i = (kG).x` values
land in a narrow window relative to `n`, making the `a_i, t_i`
distribution non-uniform in a way that produces clustering in the
Gram matrix. The cluster is the source of `μ`-oscillation.

## 6. The decisive empirical test

If hypothesis (§5) is correct, then **any curve where `n` is
strikingly close to `p` should exhibit the LLL-degeneracy**, while
curves with `n` further from `p` should not.

Let's compute `(p − n) / p` for the deployed curves:

| Curve         | `p` (bits) | `n` (bits) | `(p − n)/p`            |
|---------------|------------|------------|------------------------|
| secp192k1     | 192        | 192        | ≈ `1.5 × 10^{−29}`     |
| secp224k1     | 224        | 224        | ≈ `5.9 × 10^{−35}`     |
| secp256k1     | 256        | 256        | ≈ `3.7 × 10^{−39}` ←   |
| P-256         | 256        | 256        | ≈ `2.7 × 10^{−39}`     |

(approximate values — to be confirmed numerically.)

If secp256k1's `(p−n)/p` is markedly different from secp192k1's
and secp224k1's BUT close to P-256's, the hypothesis is REFUTED
(since P-256 doesn't fail). If all four `(p−n)/p` values are
comparable, then the closeness-of-n-to-p hypothesis is NOT the
discriminator either.

I don't have rigorous numerical comparisons here; computing them
exactly is the next step.

## 7. Empirical update — the bit-length sweep

The probe was extended to a wider curve set; results in
[`tests/lll_degeneracy_probe.rs::probe_lll_sweep_by_bit_length`](tests/lll_degeneracy_probe.rs):

| Curve            | n bits | LLL outcome     | Time             |
|------------------|--------|-----------------|------------------|
| secp192k1        |    192 | ✓ RECOVERED     | 389 ms           |
| secp224k1        |    224 | ✓ RECOVERED     | 617 ms           |
| **secp256k1**    |    256 | ✗ iteration cap | 10019 ms         |
| P-256            |    256 | ✓ RECOVERED     | 710 ms           |
| brainpoolP256r1  |    256 | ✓ RECOVERED     | 714 ms           |
| **P-384**        |    384 | ✗ iteration cap | 20392 ms         |
| **brainpoolP384r1** | 384 | ✗ iteration cap | 20131 ms         |
| **P-521**        |    521 | ✗ iteration cap | 20662 ms         |

The data identify **two distinct failure modes**:

- **Failure A (curve-specific)**: secp256k1 at 256 bits, while
  P-256 and brainpoolP256r1 succeed.  The §5 hypothesis
  ($r_i = (kG)_x \bmod n$ landing in a narrow window for
  secp256k1) remains viable but unproven.
- **Failure B (uniform at ≥ 384 bits)**: all three 384/521-bit
  curves fail.  None has the secp256k1 arithmetic; the failure
  is uniform.  **Most likely cause: the LLL iteration cap of
  $500 d^2 \cdot 8 + 10^4 \approx 4.1 \times 10^5$ is too small
  relative to the expected $O(d^4 \log B) \approx 7.7 \times 10^6$
  iterations needed for 384-bit-entry bases.**  Raising the cap
  to $\approx 10^7$ should patch Failure B (but not A).

Confirming Failure B's cap-hypothesis is a one-line change to
`crypto/src/cryptanalysis/lattice.rs`:

```rust
// Change:
let max_iter = 500 * n * n * 8 + 10_000;
// To something like:
let max_iter = 500 * n * n * (B_bits + 8) + 10_000;
```

where `B_bits` is the max-entry bit-length of the basis.  This
would scale the cap linearly with $\log B$, matching the LLL
worst-case complexity.

## 8. Reframed open questions

1. **Failure A (secp256k1 only at 256 bits)**: the candidate
   mechanism is $\mu$-oscillation (§4), and the candidate
   curve-specific source is the narrow-$r$ distribution (§5).
   The Gram-matrix eigenvalue spectrum comparison between
   secp256k1 and P-256 would empirically test this.

2. **Failure B (all $\geq$ 384-bit)**: most likely the iteration
   cap is undersized.  One-line fix above.  If LLL still fails
   after raising the cap, the genuine basis-degeneracy claim
   strengthens.

3. **Cross-cutting**: do additional $j = 0$ curves at 256 bits
   exhibit Failure A?  Test on custom-constructed Koblitz curves
   with $j = 0$ over primes other than secp256k1's specific $p$.

## 8. References

- Lenstra, Lenstra, Lovász, "Factoring polynomials with rational
  coefficients", Math. Annalen 1982.
- Boneh, Venkatesan, "Hardness of computing the most significant
  bits of secret keys in Diffie--Hellman and related schemes",
  CRYPTO 1996.
- Nguyen, Stehlé, "LLL on the average", ANTS 2006. (Average-case
  convergence analysis.)
- Schnorr, "A hierarchy of polynomial time lattice basis
  reduction algorithms", Theor. Comp. Sci. 1987. (BKZ background.)

## 9. Conclusion

The Gram–Schmidt analysis identifies **`μ`-oscillation in
size-reduction** as the likely failure mode. We've narrowed the
search to: clustering in `(a_i, t_i)` distribution that produces
near-`1/2` `μ_{ij}` values for secp256k1's specific arithmetic.
The next concrete step is to compute Gram-matrix eigenvalue
spectra empirically and confirm or refute the clustering claim.

This note documents the theoretical groundwork; the empirical
follow-up is left as future work (see §7).

---

## 10. RESOLUTION (2026-05-22): Both failures were numerical artifacts

### 10.1 Failure A → Overflow (resolved 2026-05-21)

The secp256k1 vs P-256 LLL discrepancy (§5, §7) was **not** a
genuine lattice-theoretic phenomenon.  It was a floating-point
overflow bug: `secp256k1.n ≈ 2^256`, so `n^4 ≈ 2^1024 > f64::MAX
≈ 2^1023`.  The GS norm computation produced `Inf → NaN`,
causing LLL to hit the iteration cap.  P-256 escaped only because
its `n` is smaller by `~2^192`, putting `P-256.n^4` just under
`f64::MAX`.  Fix: global scaling by `2^(max_bits − 500)`.  All
256-bit and 384-bit curves now recover 3/3 seeds.  The §5
μ-oscillation hypothesis, §4 Gram-matrix eigenvalue analysis, and
§7 r-distribution analysis are **refuted as root causes**.

### 10.2 Failure B → Catastrophic cancellation (resolved 2026-05-22)

P-521 LLL still failed after the overflow fix.  Root cause:
**catastrophic cancellation** in the f64 GS orthogonalization step.

The P-521 HNP basis has entries spanning `2^384` to `2^1042`
(658-bit dynamic range).  For the initial basis:

```
Row i:   n² × e_i   (i = 0..m-1)        [entries ~ 2^1042]
Row m:   (n·a₀, ..., n·a_{m-1}, 2^384, 0)
```

The GS orthogonalization must compute:
```
b*_m[l] = n·a_l - (a_l/n) × n²  =  n·a_l - n·a_l  =  0
```
Both terms are `~2^499` in scaled f64.  IEEE 754 subtraction of
two nearly-equal 53-bit values leaves a residual `~2^446`.  The
true `b*_m[m] = 2^384 ≈ 2^{-158}` (scaled) is `2^604` times
smaller than the phantom residual.  LLL's Lovász condition sees
`||b*_m||² ≈ m × 2^892` instead of `2^{-316}`, causing
pathological swap behavior.

### 10.3 Fix: 2048-bit fixed-point GS (`lll_reduce_hp`)

Implemented in `src/cryptanalysis/lattice.rs` as `lll_reduce_hp`.
Uses `BigInt`-based fixed-point arithmetic with `HP_PREC = 2048`
bits.  Each GS value stored as `BigInt × 2^2048`.

Key operations (`hp_mul`, `hp_div`, `hp_round`) and
`gram_schmidt_hp` function give `~2^{-1006}` residual for the
P-521 cancellation step vs `2^{446}` in f64 — a `2^{1452}`
improvement in precision.

**Empirical results** (`probe_p521_lll_hp`, m=8, k_bits=384):

| Reduction | Seed 0xC0FFEE | Seed 0xDEADBEEF | Seed 0x12345678 | Time/probe |
|-----------|---------------|-----------------|-----------------|------------|
| f64 LLL   | ✗ fail        | ✗ fail          | ✗ fail          | ~147s      |
| HP LLL    | ✓ RECOVERED   | ✓ RECOVERED     | ✓ RECOVERED     | ~79–82s    |

Note: HP LLL is **~1.8× faster** than f64 LLL because correct GS
enables proper LLL convergence (f64 LLL was looping uselessly to
the iteration cap).

### 10.4 Revised status of open questions (§8)

- **§8.1** Failure A: CLOSED.  Was overflow, not μ-oscillation.
- **§8.2** Failure B (P-384, P-521): CLOSED.  Was catastrophic cancellation.
  - P-384: overflow fix sufficient (resolved 2026-05-21).
  - P-521: required HP GS (resolved 2026-05-22).
- **§8.3** Koblitz cross-curve: CLOSED.  secp256k1 now recovers 3/3; no
  special Koblitz behavior remains.

### 10.5 Efficiency at higher m — CLOSED (2026-06-06)

Empirical results with incremental GS swap (implemented 2026-05-28):

| m  | dim   | per-probe time | 3-seed total | outcome          |
|----|-------|----------------|--------------|------------------|
| 8  | 10×10 | ~14s           | ~42s         | 3/3 recovered    |
| 16 | 18×18 | ~23s           | ~69.7s       | 3/3 recovered ✓  |
| 32 | 34×34 | ~57s           | —            | 1/1 recovered ✓  |

Actual scaling 8→32: 4.1× measured vs 11.6× theoretical O((m+2)²).
The incremental GS swap changes the dominant cost; the real complexity is
sub-quadratic in practice.  m=32 at ~57s/probe is fully practical.

**§10.5 CLOSED**: HP LLL scales to at least m=32 on P-521 with
sub-quadratic empirical growth.  No further optimisation needed for
signature counts ≤ 32.  For m > 32 the BKZ-style block reduction
is the natural next step (out of scope for this analysis).

