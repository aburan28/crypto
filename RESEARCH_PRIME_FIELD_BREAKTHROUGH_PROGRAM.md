# Research Program: Where a Prime-Field ECDLP Breakthrough Would Have To Come From

**Date.** 2026-10-08.  **Status.** Analysis plus two measured probes.
**No breakthrough is claimed.** This note turns the question "pursue a
breakthrough for ECDLP over prime fields" into (i) a quantitative map of what
is *provably or heuristically closed*, (ii) the exact condition any new
method must satisfy, and (iii) the cheapest falsifiers, two of which were run
today with the native engine from `RESEARCH_PRIME_FAST_ECDLP.md`.
**Code.** `examples/prime_ecdlp_breakthrough_probe.rs`; raw output in
`research/prime_fast_bench_20261008/probeA.md|json` (preprocessing rho,
corrected online walker) and `probe.md|json` (lattice probe; its probe-A
section is the pre-fix run).

---

## 0. The one-paragraph map

Everything that treats `E(F_p)` as a black-box group plus hashing of
coordinates — rho, kangaroo, BSGS, *and* every index-calculus variant whose
decomposition step enumerates or table-looks-up factor-base elements —
obeys the preprocessing lower bound `S·T² ≳ n` (Corrigan-Gibbs–Kogan 2018)
and the uniform bound `T ≳ √n` (Shoup 1997). The only non-generic handle a
prime field offers is that coordinates are *integers mod p* with no subfield,
no Frobenius and no nontrivial automorphism (for generic `j`). The single
known way to exploit integers-mod-`p` structure algebraically is small-root
lattice reduction (Coppersmith), and for Semaev's summation polynomials its
reach falls short of what index calculus needs by a factor in the exponent
that grows with the decomposition size. A breakthrough therefore needs one
of: a small-root method far beyond Coppersmith for a specific structured
family of polynomials; a factor base with an *algebraic* membership
condition of degree ≪ |F| (which no subset of `F_p` has); or a transfer of
the problem to a setting with more structure (isogeny / cover) that the
public probes in this repository have not found. Each is listed in §4 with a
falsifier.

## 1. Generic ceiling, in the units this repository already uses

Let `n ≈ p` be the prime group order. For any generic algorithm with `S`
words of target-independent advice and `T` online group operations per
target, `S·T² ≥ Ω̃(n)`.

| method | advice `S` | per-target `T` | `S·T²` |
|---|---:|---:|---:|
| rho (negation map) | 0 | `√(πn/4)` | `≈ n` (uniform bound) |
| preprocessing rho (Bernstein–Lange) | `n^{1/3}` | `≈ 2n^{1/3}` | `≈ 4n` — **tight** |
| IC-2, direct subtraction (PR #1560) | `B` logs | `p/(2B)` | `p²/(4B)` ≥ `p/4` |
| IC-3, pair table (PR #1560) | `B²` table | `0.75p/B²` | `0.56p²/B²` ≥ `0.56p` |
| IC-`m`, `(m−1)`-sum table | `B^{m−1}` | `(m!/2^m)·p/B^{m−1}` | `≥ p` |

The index-calculus rows are on or above the line for every `B`; the pair
table is just the generic time/memory trade-off in different clothing. This
is why the ledger's regime-C `vs_rho` target ("charged IC < rho") cannot be
met by any variant whose decomposition is generic, however fast the
arithmetic — PR #1560 made the constants ~10³× better and the gap to rho
from 32 bits upward is still > 40× and widening as `p/B` vs `√p`.

**What non-genericity must look like.** With `B` precomputed logs (`S = B`)
the per-target cost must be `T ≪ √(p/B)`; with `B ≈ p^{1/3}` that means
`T ≪ p^{1/3}`, i.e. decomposing a target in time *sub-linear in the factor
base* (the enumeration costs `B = p^{1/3}` by itself). Only an algebraic
decomposition can do that.

## 2. The algebraic lever and its reach

### 2.1 Small-x factor bases → small roots of summation polynomials

With `F = {P : x(P) ≤ B}`, `B = p^δ`, an `m`-decomposition of `R` is a
root `(x_1, …, x_m) ∈ [1,B]^m` of `S_{m+1}(x_1, …, x_m, x_R) ≡ 0 (mod p)`,
a polynomial of degree `2^{m−1}` in each variable. Trials per relation are
`≈ m!·p/(2B)^m`, so to get the *total* relation cost `B·(trials)` below `√p`
one needs `B^{m−1}·2^m/m! ≳ √p`, i.e. roughly

    δ ≥ 1 / (2(m−1))        (δ ≥ 1/2 for m=2, 1/4 for m=3, 1/6 for m=4),

and even then each trial must cost `O(1)` lattice reductions.

### 2.2 What Coppersmith-type lattices reach for this shape

For a modular polynomial whose monomials fill the hypercube `[0,d]^m` (here
`d = 2^{m−1}`), the Jochemsz–May basic lattice at level `ℓ` has dimension
`(dℓ+1)^m`, and the determinant condition gives, as `ℓ → ∞`,

    δ_max = 2 / (m(m+1)·d) = 1 / (m(m+1)·2^{m−2}).

(Derivation: Σ over monomials of `(ℓ − k)`, with `k = min_i ⌊α_i/d⌋`, is
`(dℓ)^m·ℓ·m/(m+1)`; Σ of each exponent is `(dℓ)^m·dℓ/2`; Howgrave-Graham
requires `p^{Σ(ℓ−k)}·X^{mΣα} < p^{dim·ℓ}`.) Numerically:

| m | polynomial | `δ_max` (lattice) | `δ` needed (§2.1) | shortfall |
|---|---|---:|---:|---:|
| 2 | `S₃`, degree 2 per variable | **1/6** (level 1: 1/18, level 2: 1/10, level 3: 0.119) | 1/2 | 3× |
| 3 | `S₄`, degree 4 | 1/24 | 1/4 | 6× |
| 4 | `S₅`, degree 8 | 1/80 | 1/6 | 13× |
| 5 | `S₆`, degree 16 | 1/240 | 1/8 | 30× |

Extended Jochemsz–May strategies move these by tens of percent, not by
factors of 3–30, and the symmetrised (elementary-symmetric) rewriting of
`S_{m+1}` lowers the monomial count, not the degree reach. **Probe B** (§3.2)
measures the `m = 2` row empirically so the derivation is not taken on faith.

### 2.3 Why algebraic factor bases do not rescue it

The decomposition is cheap over extension fields `F_{q^n}` because
"`x ∈ F_q`" is the *low-degree* condition `x^q = x` whose Weil restriction
turns `S_{m+1}` into `n` equations over `F_q` in `m·…` unknowns with a small
solution set. Over `F_p` any subset `F ⊂ F_p` of size `B` has membership
polynomial `L(X) = Π_{x∈F}(X − x)` of degree exactly `B`, with no
restriction-of-scalars trick to lower it; the system `S_{m+1} = 0, L(x_i) = 0`
has `≈ (2B)^m/m!·…` solutions over the algebraic closure, so any solver
that enumerates solutions costs ≥ `B^{m−1}`, which is enumeration again.
Multiplicative-subgroup bases (`x^B = 1`) make `L` sparse but not
lower-degree; `j = 0` orbits divide `B` by 3; neither changes the exponent.

## 3. Probes run on 2026-10-08

### 3.1 Probe A — preprocessing rho on prime-field curves (control)

The one known sub-`√n`-per-target method, never before run on prime curves
in this repository. `S = 2^s ≈ n^{1/3}` table entries from `G`-only walks,
`2^dp ≈ √(n/S)`, online walk from `[a₀]G + Q`.

| bits | S | dp | precompute steps | precompute wall s (14 thr) | online steps, mean of 16 | `2√(n/S)` | plain rho `√(πn/4)` (measured) | per-target step ratio | online wall ms | plain rho wall s | break-even targets |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 32 | 2^11 | 10 | 2.225e6 | 0.0 | 5.372e3 (16/16) | 2.896e3 | 5.808e4 (measured 6.029e4) | 11× | 7.80 | 0.061 | 1 |
| 36 | 2^12 | 12 | 2.005e7 | 0.3 | 6.140e3 (16/16) | 8.192e3 | 2.323e5 (measured 1.702e5) | 38× | 8.58 | 0.160 | 2 |
| 40 | 2^13 | 14 | 3.046e8 | 10.3 | 2.195e4 (16/16) | 2.317e4 | 9.293e5 (measured 7.567e5) | 42× | 38.71 | 0.617 | 18 |
| 44 | 2^15 | 15 | 3.188e9 | 64.5 | 6.188e4 (16/16) | 4.634e4 | 3.717e6 (measured 5.994e6) | 60× | 161.56 | 0.719 | 116 |

All 64 targets recovered and checked against the planted scalar. (A first
run whose online walker restarted at fruitless 2-cycles instead of taking
the precomputation's doubling escape sat 2–7× above `2√(n/S)`; the fix is
in the example and is a useful reminder that the online and precomputation
walks must be *identical* maps, escapes included. The first run is kept as
`probe.md` in the evidence directory.)

Reading: online steps track `2√(n/S)`, i.e. `n^{1/3}`-scaling per target,
against `√(πn/4)` for plain rho — the per-target factor grows like
`n^{1/6}`. The break-even column is how many targets on the same curve must
be solved before the precomputation pays for itself. This is the honest
"sub-rho per target" number for the ledger: non-uniform, generic, on the
`S·T² ≈ n` line.

### 3.2 Probe B — empirical small-root reach for `S₃`

Planted `R = F_i ± F_j` with `x_i, x_j ≤ p^δ`; level-`m` lattice of dimension
`(2m+1)²`; LLL (δ = 0.99); success means at least two reduced polynomials
vanish at the planted root over `Z` (the condition under which the
resultant step recovers it).

| bits | m | dim | δ | B = ⌊p^δ⌋ | ≥2 vanishing / trials | Howgrave-Graham holds (first two rows) | mean vanishing rows | ms per LLL |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 40 | 1 | 9 | 0.04 | 3 | 5/5 | 4/5 | 5.0 | 0 |
| 40 | 1 | 9 | 0.06 | 5 | 7/7 | 5/7 | 3.0 | 0 |
| 40 | 1 | 9 | 0.08 | 9 | 8/8 | 8/8 | 2.1 | 0 |
| 40 | 1 | 9 | 0.10 | 15 | 8/8 | 6/8 | 2.1 | 0 |
| 40 | 1 | 9 | 0.12 | 27 | 7/7 | 0/7 | 2.0 | 0 |
| 40 | 1 | 9 | 0.14 | 48 | 5/8 | 0/8 | 1.5 | 0 |
| 40 | 1 | 9 | 0.16 | 84 | 1/8 | 0/8 | 0.5 | 0 |
| 40 | 1 | 9 | 0.20 | 255 | 0/8 | 0/8 | 0.0 | 0 |
| 40 | 2 | 25 | 0.04 | 3 | 3/3 | 3/3 | 24.0 | 139 |
| 40 | 2 | 25 | 0.06 | 5 | 7/7 | 7/7 | 22.9 | 112 |
| 40 | 2 | 25 | 0.08 | 9 | 8/8 | 8/8 | 21.9 | 108 |
| 40 | 2 | 25 | 0.10 | 15 | 8/8 | 8/8 | 21.4 | 370 |
| 40 | 2 | 25 | 0.12 | 27 | 8/8 | 0/8 | 7.5 | 509 |
| 40 | 2 | 25 | 0.14 | 48 | 6/8 | 1/8 | 3.9 | 291 |
| 40 | 2 | 25 | 0.16 | 84 | 2/8 | 0/8 | 0.9 | 233 |
| 40 | 2 | 25 | 0.20 | 255 | 0/8 | 0/8 | 0.4 | 194 |
| 40 | 3 | 49 | 0.04 | 3 | 8/8 | 8/8 | 47.0 | 6476 |
| 40 | 3 | 49 | 0.06 | 5 | 5/5 | 5/5 | 47.0 | 7214 |
| 40 | 3 | 49 | 0.08 | 9 | 8/8 | 7/8 | 41.0 | 7883 |
| 40 | 3 | 49 | 0.10 | 15 | 7/7 | 7/7 | 42.6 | 4600 |
| 40 | 3 | 49 | 0.12 | 27 | 7/7 | 0/7 | 43.3 | 4231 |
| 40 | 3 | 49 | 0.14 | 48 | 7/8 | 0/8 | 17.0 | 2283 |
| 40 | 3 | 49 | 0.16 | 84 | 5/8 | 0/8 | 1.5 | 1747 |
| 40 | 3 | 49 | 0.20 | 255 | 0/8 | 0/8 | 0.0 | 1320 |
| 56 | 1 | 9 | 0.04 | 4 | 6/6 | 6/6 | 4.8 | 0 |
| 56 | 1 | 9 | 0.06 | 10 | 6/6 | 5/6 | 2.3 | 0 |
| 56 | 1 | 9 | 0.08 | 22 | 8/8 | 8/8 | 2.0 | 0 |
| 56 | 1 | 9 | 0.10 | 48 | 8/8 | 4/8 | 2.0 | 0 |
| 56 | 1 | 9 | 0.12 | 105 | 3/8 | 0/8 | 1.1 | 0 |
| 56 | 1 | 9 | 0.14 | 229 | 0/8 | 0/8 | 0.1 | 0 |
| 56 | 1 | 9 | 0.16 | 497 | 0/8 | 0/8 | 0.1 | 0 |
| 56 | 2 | 25 | 0.04 | 4 | 7/7 | 7/7 | 23.7 | 67 |
| 56 | 2 | 25 | 0.06 | 10 | 5/5 | 5/5 | 22.4 | 43 |
| 56 | 2 | 25 | 0.08 | 22 | 7/7 | 7/7 | 22.0 | 76 |
| 56 | 2 | 25 | 0.10 | 48 | 7/7 | 7/7 | 18.4 | 103 |
| 56 | 2 | 25 | 0.12 | 105 | 8/8 | 1/8 | 8.8 | 50 |
| 56 | 2 | 25 | 0.14 | 229 | 4/8 | 0/8 | 1.5 | 57 |
| 56 | 2 | 25 | 0.16 | 497 | 0/8 | 0/8 | 0.0 | 61 |
| 56 | 3 | 49 | 0.04 | 4 | 7/7 | 7/7 | 47.1 | 2130 |
| 56 | 3 | 49 | 0.06 | 10 | 5/5 | 5/5 | 38.6 | 956 |
| 56 | 3 | 49 | 0.08 | 22 | 7/7 | 7/7 | 47.1 | 4159 |
| 56 | 3 | 49 | 0.10 | 48 | 8/8 | 8/8 | 47.0 | 6063 |
| 56 | 3 | 49 | 0.12 | 105 | 7/7 | 0/7 | 30.1 | 1786 |
| 56 | 3 | 49 | 0.14 | 229 | 3/8 | 0/8 | 1.5 | 1005 |
| 56 | 3 | 49 | 0.16 | 497 | 0/8 | 0/8 | 0.1 | 597 |

Reading: the rigorous (Howgrave-Graham) criterion holds up to δ = 0.10 at every
level on both curves and never beyond; the heuristic criterion (two reduced
polynomials vanishing at the root over Z) fades between δ = 0.12 and 0.16 and
is zero by 0.16–0.20, sharpening toward the predicted 1/6 as `p` grows (the
56-bit curve fails earlier than the 40-bit one at every level). Raising the
level from 1 to 3 buys at most 0.02–0.04 in δ at 25–75× the cost per
reduction. So the lattice lever reaches `B ≈ p^{0.10–0.14}` for the
2-decomposition, against the `p^{1/2}` it would need to beat rho, and each
reduction costs 0.05–8 s against ~1 µs per enumeration step. The exponent
gap is 3–5×, not a constant, exactly as §2.2 predicts; higher `m` makes it
worse. Both the theory and this measurement say: **no breakthrough via
Coppersmith on the summation polynomials, at any size.**

## 4. Candidate levers, each with a falsifier

1. **Small roots beyond Coppersmith for the Semaev shape.** The `S_{m+1}`
   polynomials are not generic: they are symmetric, of a fixed resultant
   form, and their small-root instances come with the *curve structure*
   (the roots are x-coordinates of points whose sum is fixed). A method that
   reaches `δ ≥ 1/(2(m−1))` for them would be a breakthrough.
   *Falsifier:* Probe B's framework with any proposed lattice/embedding; the
   target row is `δ ≥ 1/4` at `m = 3` on a 40-bit curve. Nothing in the
   literature (Coppersmith 1996, Jochemsz–May 2006, Petit–Kosters–Messeng
   2016) reaches it.
2. **A factor base with sub-`B`-degree membership.** Needs a subset of
   `F_p` of size `B` cut out by equations of degree `≪ B` in a way that
   survives substitution into `S_{m+1}`. Candidates people have tried: images
   of low-degree maps (`x = g(t)`, `t` small — that is just small roots of
   `S_{m+1}(g(t_1), …)` with higher degree, worse by §2.2), subgroups,
   arithmetic progressions. *Falsifier:* count solutions of the combined
   system over `F̄_p` for a toy `p`; if it is `≥ B^{m−1}`, the lever is dead.
3. **Transfer to a structured setting.** Isogenies to curves with
   exploitable endomorphisms, covers by higher-genus curves over `F_p`,
   lifts to characteristic 0 (anomalous only). This repository's isogeny
   walker, CM audit and cover-route notes are the current evidence that
   nothing within reachable degree exists for the deployed curves.
   *Falsifier:* a weak target reachable at degree `ℓ` with
   `cost(ℓ-isogeny walk) + cost(attack on target) < √n`.
4. **Non-uniform / multi-target accounting.** Not a breakthrough but the
   only place sub-`√n` per target is real: preprocessing rho (Probe A),
   Kuhn–Struik multi-target rho. *Use:* make the regime-C ledger charge
   precomputation explicitly and record `S·T²` per row.
5. **Quantum.** Shor's algorithm is the actual breakthrough for every
   group; this repository's `quantum_estimator` already prices it. Not a
   classical lever.

### 4.6 Also examined on 2026-10-08 and closed on paper

- **Symmetric rewriting of `S₃` in `(e₁, e₂) = (x₁+x₂, x₁x₂)`.** Six
  monomials instead of nine, bounds `(2B, B²)`. Level-1 reach improves from
  `1/18` to `1/12`, but the asymptotic reach for the total-degree-2 triangle
  with these bounds is `1/9 < 1/6`: the smaller monomial set is outweighed
  by the `B²` bound on `e₂`. Not worth a run.
- **Multiplicative-subgroup factor bases** (`x^B = 1`, `B | p−1`). The
  membership polynomial becomes sparse, and `Π_{h∈H} S₃(h, X, z)` can be
  built by cyclic-norm folding in `O(log B)` polynomial products — but the
  output has degree `2B`, so every route (fold modulo `X^B − 1`, gcd, batched
  evaluation of `Ψ(Z) = Π S₃(x_i, x_j, Z)` over many targets) costs `Õ(B)`
  per target or returns to the `B²`-advice generic trade-off. No exponent
  change; it needs `p − 1` to have a divisor near `p^{1/3}`, which the
  deployed primes do not offer anyway.
- **Global lifting (xedni-style).** Take `E/Q` with small coefficients, use
  reductions of small-height rational points as the factor base so that
  relations in `E(Q)` reduce to free relations mod `p`. All those logs are
  spanned by `rank(E(Q))` unknowns, and the target still has to be written
  as a bounded combination of the generators — an `r`-dimensional
  bounded-coefficient DLP that is generic again. This is Silverman's xedni
  calculus; Jacobson–Koblitz–Silverman–Stein–Teske (2000) showed it does not
  beat rho.
- **Embedding into an extension to buy a subfield.** `E(F_p) ⊂ E(F_{p^n})`
  makes the subfield factor base `x ∈ F_p` algebraic, but Gaudry's cost
  there is `Õ(p^{2−2/n})` against the original `√p`: `p` for `n = 2`,
  `p^{4/3}` for `n = 3`. The embedding enlarges the problem faster than the
  structure helps.

## 5. Decision

Pursue (1) only with a concrete new lattice/embedding idea that changes the
exponent, never by tuning; keep (2)–(3) as closed-unless-new-structure; use
(4) to make the ledger's accounting honest. Everything else in regime C is
constant-factor engineering, which PR #1560 has largely done.
