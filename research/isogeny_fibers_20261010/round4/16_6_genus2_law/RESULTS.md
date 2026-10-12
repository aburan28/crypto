# Lead 16.6: the affine saturation law on a genus-2 Jacobian (round 4, 2026-10-11)

Code: `examples/r4_16_6_genus2_census.rs`. Raw outputs: `raw/`. Commands:
`run_all.sh`. Graphs: `law_measured_vs_predicted.svg`, `A_pair_correction.svg`
(script `plot_law.py`, data `plot_data.json`, `finite_size_A.json`). PDF:
`report.pdf` (source `report.typ`).

## Question, hypothesis, falsifier (§16.6)

*Question.* Does the affine saturation law of Prop. 15.2 (cost factor
|M|^{m−2}/(|A|·i) for the mitm pipeline) hold on a genus-2 Jacobian, where
the automorphism group is larger (order 10 on y² = x⁵ + b) and rational
(2,2)-kernels exist?
*Hypothesis.* The law is dimension-free: an A-closed base gains |A| (here 5
beyond the negation that every base already folds), a base saturated by a
rational kernel M of order μ is neutral at m = 2 and loses μ at m = 3, an
image base of index i gains i.
*Falsifier.* A measured factor outside the law beyond sampling error.

## Verdict: confirmed, with one structural difference

1. **The numbers follow the law in genus 2** for every arm we could build:
   A-closed with |A| = 10 (gain 5.0 beyond ±), kernel saturation with
   μ = 2, 4 (a Richelot (2,2)-kernel and a non-isotropic order-4 group, which
   behave the same), 5, the combined A + M₅ base (cost μ^{m−2}/|A|: 1/5 at
   m = 2, 1 at m = 3), and image bases with i = 5 and i = 16. At the
   largest field (p = 4111, r ≈ 10⁶ to 3·10⁶), the m = 3 cost ratios are
   0.198–0.200 (predicted 0.2), 1.89–1.98 (2), 3.82–3.92 (4), 4.85 (5) and
   0.977 (1). The class factors are 1.00 (1), 0.255 (0.25), 0.064 (0.0625)
   and 0.041 (0.04). All are measured; the intervals are in the tables below.
2. **One leading-order deviation, explained and not genus-specific.** On the
   curve whose cofactor group is fixed by ζ (C1, H = ⟨T₅⟩), the A-closed base
   gains extra m = 2 hits. The measured c₂ ratio is 1.32, 1.12–1.21 and 1.09
   at |F| = 50, 100 and 200, outside the law's leading term. The derived
   correction 1 + (h−1)(|A|−2)/(2|F|) = 1 + 16/|F| predicts 1.32, 1.16 and
   1.08. It comes from the pairs {ζᵃP, ζᵇ(−P)}, whose class (ζᵃ − ζᵇ)·cls(P)
   is 0 whenever A acts trivially on H. It vanishes as 1/|F|, is absent on
   C2 (ζ fixes no nonzero 2-torsion class; measured c₂ ratio 0.95–1.04), and
   has the same form for the j = 0 elliptic census of §15.2 (derived +5.6%
   there; §15.2 measured +2% in one census per size). This is a finite-size
   correction that grows with |A|, not a new saturation regime.
3. **Genus 2 differs in realizability, not in the law.** A base of degree-1
   divisors [P − ∞] can never be closed under translation by a nonzero t.
   C − ∞ and C − ∞ + t share at most 2 points (proof below). So every
   kernel-saturated base contains (μ−1)/μ of its elements as degree-2
   divisors (recorded: `degree1_points` = |F|/μ in every M-arm census). The
   census treats them as abstract group elements, so the law still holds
   for mitm. A Gaudry/Nagao-type decomposition oracle factors over degree-1
   divisors with x in a set, and cannot use such a base at all. In genus 2
   the §14.2 linearity barrier becomes total for kernel saturation, while
   automorphism folding stays free: any ζ-stable x-set gives an A-closed
   degree-1 base.

So the law is dimension-free in the census. Higher genus brings no extra
saturation regime. It brings a larger |A|, a larger finite-size correction
proportional to |A|, and the loss of kernel saturation as a degree-1 base
option.

## What ran

All runs are single-threaded on the shared 4-core host, about 6.5 min of
contended wall time in total (wall time is not used as evidence).

```sh
export CARGO_TARGET_DIR=…/scratchpad/r4target
cargo build --release --example r4_16_6_genus2_census
B=$CARGO_TARGET_DIR/release/examples/r4_16_6_genus2_census
$B search --lo 4000 --hi 4300 --maxh 16          # also 1000–1200, 2000–2300
$B census --p 4111 --b 2172 --size 80  --reps 40 --seed 1 --arms plain,A,M2,M4iso,M4non,G --json raw/census_c2_p4111_s80.json
$B census --p 4111 --b 2172 --size 160 --reps 40 --seed 3 --arms plain,A,M2,M4iso,M4non --json raw/census_c2_p4111_s160.json
$B census --p 4111 --b 3984 --size 100 --reps 40 --seed 2 --arms plain,A,M5,AM5,G --json raw/census_c1_p4111_s100.json
$B census --p 4111 --b 3984 --size 50|200 --reps 40 --seed 4 --arms plain,A --json raw/census_c1_p4111_s{50,200}_Aonly.json
# and the same census at p = 1051 and 2131 (run_all.sh has every line)
```

### Curves (toy instances; named by equation, no ledger slug exists for genus 2)

All curves are y² = x⁵ + b over F_p with p ≡ 1 (mod 10), so ζ₅ ∈ F_p and
A = ⟨ζ, ι⟩ ≅ Z/10 acts on J(F_p). They were chosen by `search`, which
computes all ten twists per prime (`raw/search_*.jsonl`).

| family | p | b | #J(F_p) | h | r | rational kernels |
|:--|--:|--:|--:|--:|--:|:--|
| C2: −b a 5th power, b non-square, J[2] ⊂ J(F_p) | 1051 | 1042 | 988 016 = 2⁴·61 751 | 16 | 61 751 | M2, M4iso, M4non |
| C2 | 2131 | 32 | 4 382 416 = 2⁴·273 901 | 16 | 273 901 | same |
| C2 | 4111 | 2172 | 16 434 896 = 2⁴·1 027 181 | 16 | 1 027 181 | same |
| C1: b square, −b not a 5th power | 1051 | 49 | 1 116 205 = 5·223 241 | 5 | 223 241 | M5 = ⟨(0,√b) − ∞⟩ |
| C1 | 2131 | 4 | 4 675 355 = 5·935 071 | 5 | 935 071 | M5 |
| C1 | 4111 | 3984 | 16 582 255 = 5·3 316 451 | 5 | 3 316 451 | M5 |

On C2, J(F_p) ≅ (Z/2)⁴ × Z/r, so 2J(F_p) = G and the image arm "inside
[2]J" coincides with the G arm (i = 16). It is reported once, as G. On C1,
T₅ = (0, √b) − ∞ has order 5 (5T₅ = div(y − √b)) and is ζ-fixed.

### Validation (every census file, field `validation`; all passed)

* #J by an **exhaustive enumeration of all reduced Mumford pairs**. Every
  u of degree ≤ 2 is visited and every candidate v is verified by u | f − v²,
  with 0 solver failures. It equals L(1) from #C(F_p) and #C(F_{p²}) for all
  six curves (`order_exhaustive` = `order_lpoly`).
* [#J]D = 0 for 40 random classes. Associativity, commutativity and
  D + (−D) = 0 hold on 300 random triples.
* Cantor addition agrees with the library's
  `crypto_lib::cryptanalysis::jv_cover::Hyp` (over F_{p²} ⊃ F_p) on 300
  random pairs per curve.
* ζ is a homomorphism (200 pairs), ζ⁵ = id, and ζ acts on G as a scalar λ
  with λ⁵ ≡ 1 (mod r) (`zeta_eigenvalue_on_g`).
* The kernels are closed subgroups. The Richelot kernel
  {0, e₅ − ∞, e₁ + e₂ − 2∞, e₃ + e₄ − 2∞} sums correctly. The 2-Sylow
  subgroup has order 16 and [J : 2J] = 16 on C2.

### Census (exactly §15.2's, ported)

Every base F is closed under negation. Columns U are the orbits of [h]F
under A = ⟨ζ, −1⟩. D_m is the exact number of distinct nonzero classes of G
that are sums of m base elements, with repetition, as in
`round3/automorphism_vs_translation.py`. The class factor is
c = D_m·m!/|F|^m. The cost model is U·(r/D_m)·L_m with L₂ = 1 and L₃ = |F|.
The census also records R_m (multisets hitting G) and D_m restricted to
pairwise distinct summands ("strict").

Arms (base elements are Mumford pairs):

* plain: degree-1 divisors [P − ∞], both signs, at random x (and the §15.2
  small-x base as `smallx`).
* A: S ∪ ζS ∪ … ∪ ζ⁴S, with S of size |F|/5 drawn from distinct ζ-orbits.
* M2, M4iso, M4non, M5: S ∪ (S + t) for t ∈ M, with |S| = |F|/μ.
* AM5: closed under both ζ and M5, with |S| = |F|/25.
* G: degree-1 divisors lying in G (i = h).

Each (curve, |F|) has 40 random bases plus the small-x base. Ratios to plain
are ratios of means over the 40 bases. The intervals are 95% percentile
bootstraps over bases (2000 resamples, computed in Rust).

**Normalization.** Each table also gives the cost with the plain base's
columns set to their nominal |F|/2 ("cost nominal"). The reason is that the
x-set over which y² = x⁵ + b has points is ζ-stable. A random plain base
therefore contains some ζ-related pairs by accident: about 8·C(|F|/2, 2)/p
of them. Measured plain columns are 0.83, 0.92 and 0.95 of |F|/2 at
p = 1051, 2131 and 4111. The raw cost ratios are in the JSON files.

## Measured (ratio to the plain base on the same curve, same |F|; prediction in parentheses)

p = 4111 (the largest field), 40 bases each:

| curve | \|F\| | arm | U | c, m = 2 (strict) | c, m = 3 | cost nominal m = 2 | cost nominal m = 3 | hit rate m = 3 |
|:--|--:|:--|:--|:--|:--|:--|:--|--:|
| C2 | 160 | A | 0.217 [0.215, 0.219] (0.2) | 0.980 [0.937, 1.025] (1) | 0.999 [0.987, 1.010] (1) | 0.203 [0.196, 0.211] (0.2) | 0.200 [0.198, 0.203] (0.2) | 0.041 |
| C2 | 160 | M2 | 0.523 [0.516, 0.530] (0.5) | 0.503 [0.482, 0.525] (0.5) | 0.255 [0.251, 0.259] (0.25) | 0.958 [0.925, 0.992] (1) | 1.888 [1.850, 1.922] (2) | 0.010 |
| C2 | 160 | M4iso | 0.266 [0.263, 0.270] (0.25) | 0.232 [0.225, 0.240] (0.25) | 0.064 [0.063, 0.065] (0.0625) | 1.043 [1.015, 1.069] (1) | 3.825 [3.761, 3.885] (4) | 0.003 |
| C2 | 160 | M4non | 0.266 [0.263, 0.270] (0.25) | 0.241 [0.232, 0.251] (0.25) | 0.064 [0.063, 0.066] (0.0625) | 1.013 [0.978, 1.047] (1) | 3.819 [3.736, 3.907] (4) | 0.003 |
| C2 | 80 | A | 0.208 [0.206, 0.211] (0.2) | 1.017 [0.912, 1.133] (1) | 1.007 [0.965, 1.052] (1) | 0.198 [0.183, 0.213] (0.2) | 0.199 [0.190, 0.207] (0.2) | 0.005 |
| C2 | 80 | M2 | 0.509 [0.501, 0.516] (0.5) | 0.529 [0.493, 0.567] (0.5) | 0.247 [0.234, 0.260] (0.25) | 0.937 [0.891, 0.987] (1) | 1.975 [1.878, 2.091] (2) | 0.001 |
| C2 | 80 | M4iso | 0.259 [0.255, 0.262] (0.25) | 0.225 [0.208, 0.242] (0.25) | 0.065 [0.062, 0.069] (0.0625) | 1.069 [1.011, 1.128] (1) | 3.821 [3.599, 4.033] (4) | 0.000 |
| C2 | 80 | M4non | 0.259 [0.256, 0.262] (0.25) | 0.241 [0.221, 0.264] (0.25) | 0.063 [0.058, 0.068] (0.0625) | 1.023 [0.960, 1.084] (1) | 3.922 [3.634, 4.287] (4) | 0.000 |
| C2 | 80 | G (i = 16) | 0.611 [0.596, 0.626] (1) | 15.98 [15.41, 16.56] (16) | 15.27 [15.00, 15.53] (16) | 0.050 [0.049, 0.052] (0.0625) | 0.038 [0.038, 0.039] (0.0625) | 0.080 |
| C1 | 100 | A | 0.211 [0.209, 0.213] (0.2) | 1.209 [1.147, 1.274] (1; corrected 1.16) | 1.010 [0.979, 1.047] (1) | 0.165 [0.157, 0.174] (0.2) | 0.198 [0.191, 0.204] (0.2) | 0.010 |
| C1 | 100 | M5 | 0.209 [0.207, 0.212] (0.2) | 0.205 [0.203, 0.207] (0.2) | 0.041 [0.041, 0.041] (0.04) | 0.987 [0.976, 1.000] (1) | 4.851 [4.803, 4.889] (5) | 0.000 |
| C1 | 100 | AM5 | 0.042 [0.042, 0.043] (0.04) | 0.205 [0.203, 0.207] (0.2) | 0.041 [0.041, 0.041] (0.04) | 0.199 [0.197, 0.201] (0.2) | 0.977 [0.976, 0.978] (1) | 0.000 |
| C1 | 100 | G (i = 5) | 0.824 [0.809, 0.839] (1) | 5.022 [4.972, 5.069] (5) | 4.906 [4.900, 4.912] (5) | 0.156 [0.153, 0.159] (0.2) | 0.159 [0.157, 0.162] (0.2) | 0.049 |

The m = 2 class factor for C2 is the strict one. On C2 the cofactor group
(Z/2)⁴ has exponent 2, so every doubling 2P of the base lands in G. That
makes 29% (|F| = 80) and 17% (|F| = 160) of the plain base's m = 2 hits
doublings (`plain_diag_frac2`). The non-strict G-arm c₂ ratio is therefore
11.6 rather than 16, while the strict one is 16.0.

Smaller fields, same arms (full rows in `plot_data.json`):

| curve, \|F\| | arm | p = 1051: c₃; cost nom. m = 3 | p = 2131: c₃; cost nom. m = 3 | p = 4111: c₃; cost nom. m = 3 | predicted |
|:--|:--|:--|:--|:--|:--|
| C2, 80 | A | 1.021; 0.196 | 1.005; 0.199 | 1.007; 0.199 | 1; 0.2 |
| C2, 80 | M2 | 0.253; 1.83 | 0.255; 1.90 | 0.247; 1.98 | 0.25; 2 |
| C2, 80 | M4iso | 0.062; 3.90 | 0.063; 3.92 | 0.065; 3.82 | 0.0625; 4 |
| C2, 80 | M4non | 0.063; 3.80 | 0.066; 3.71 | 0.063; 3.92 | 0.0625; 4 |
| C2, 80 | G | (not enough degree-1 points in G at p = 1051; arm skipped) | 14.26; 0.021 (hit rate 0.27) | 15.27; 0.038 (hit rate 0.08) | 16; 0.0625 |
| C1, 100 | A | 0.994; 0.201 | 0.982; 0.204 | 1.010; 0.198 | 1; 0.2 |
| C1, 100 | M5 | 0.044; 4.38 | 0.041; 4.77 | 0.041; 4.85 | 0.04; 5 |
| C1, 100 | AM5 | 0.044; 0.917 | 0.041; 0.966 | 0.041; 0.977 | 0.04; 1 |
| C1, 100 | G | 3.86; 0.092 (hit rate 0.54) | 4.69; 0.138 (0.17) | 4.91; 0.159 (0.05) | 5; 0.2 |

As in §15.2, the saturated m = 3 costs (M5: 4.38 → 4.77 → 4.85; AM5: 0.917
→ 0.966 → 0.977) climb to the prediction as r grows. The G arm's c₃ climbs
to i as the m = 3 hit rate falls (saturation of small r).

### Deviations from the leading-order law, each measured and accounted for

| effect | where | size measured | account (derived) | genus-specific? |
|:--|:--|:--|:--|:--|
| A-structured pairs {ζᵃP, ζᵇ(−P)} land in G when A acts trivially on H | C1, A arm, m = 2 | c₂ ratio 1.316 [1.227, 1.407], 1.209 [1.147, 1.274], 1.092 [1.059, 1.131] at \|F\| = 50, 100, 200 (p = 4111) | 1 + (h−1)(\|A\|−2)/(2\|F\|) = 1.32, 1.16, 1.08 | no: j = 0 elliptic has the same term with \|A\| = 6 |
| plain base partially ζ-folded (ζ-stable x-set) | every curve | plain U/(\|F\|/2) = 0.83 / 0.92 / 0.95 at p = 1051 / 2131 / 4111 | expected ζ-related pairs ≈ 8·C(\|F\|/2, 2)/p | no (j = 0 elliptic, ω-stable x-set) |
| G-filtered base partially A-closed | G arms | U ratio 0.43–0.82 (C1), 0.32–0.61 (C2) | G is A-stable; its few degree-1 points cluster in ζ-orbits | no |
| doublings 2P ∈ G for an exponent-2 cofactor | C2, m = 2 | 17–29% of plain hits | \|F\| doublings vs \|F\|²/(2h) generic hits | no (any h with exponent 2) |
| multiset repetition at m = 3 for M arms | M arms | c₃ 2–10% above 1/μ² (e.g. M4: 0.064 vs 0.0625; M5: 0.041–0.044 vs 0.04) | μ·C(\|F\|/μ+2, 3)/C(\|F\|+2, 3) ≈ μ⁻²(1 + 3(μ−1)/\|F\|) | no |

None of these is a regime change. Each is O(·/|F|) or O(·/r), and each has
the same form for the elliptic census.

### Proof of the structural statement (point 3 of the verdict)

Let t ∈ J(F_p), t ≠ 0, and suppose [P − ∞] + t = [Q − ∞] for points P, Q
(∞ allowed for Q, giving the class 0). Then t = [Q − P] = [Q + ιP − 2∞].
The reduced Mumford representative of t is unique. If it has degree 2, say
R₁ + R₂ − 2∞, then Q + ιP ∼ R₁ + R₂. Since R₁ + R₂ is not a fibre of x,
|R₁ + R₂| = {R₁ + R₂}, so {Q, ιP} = {R₁, R₂}: at most two choices of P. If it
has degree 1, R − ∞, then Q + ιP ∼ R + ∞. That divisor is not a fibre of x
either, so {Q, ιP} = {R, ∞}. The only solution with P a finite point is
P = ιR, Q = ∞. Hence |(C − ∞) ∩ (C − ∞ − t)| ≤ 2. A base of degree-1 divisors
closed under translation by t therefore has at most 2 elements. ∎ (The
elliptic case differs because there C − ∞ is all of E.)

## Limits (what was not done)

* Toy sizes only: p ≤ 4111, r ≤ 3.3·10⁶, |F| ≤ 200, one curve per (family, p),
  40 bases per cell. No statement about any deployed curve and no sub-rho
  claim. The cost is the §15.2 mitm *model* (columns × r/D_m × L_m), not a
  measured solve.
* No algebraic (Gaudry/Nagao) oracle was run. The realizability statement
  for such oracles is a proof about the base, not a measurement.
* The image arm via an explicit Richelot isogeny ψ: J′ → J was **not**
  built (no codomain or map was constructed). The only image arm is
  [2]J(F_p) = G on C2 (i = 16), plus G on C1 (i = 5).
* The M arms use rational subgroups (Richelot kernel, a non-isotropic
  order-4 group, ⟨T₅⟩) as translation groups. Nothing was computed on the
  quotient surfaces.
* `cargo clippy --release --example r4_16_6_genus2_census -- -D warnings`
  fails in this worktree because of pre-existing warnings in `src/`
  (`koblitz_factor_base_search.rs`, `koblitz_index_calculus.rs`,
  `koblitz_groebner.rs`), not because of this example. Plain
  `cargo clippy` reports zero warnings for the example, and
  `rustfmt --edition 2021 --check` passes.

## Next step

1. Run the same census through an algebraic decomposition oracle for
   genus 2 (Nagao-style, degree-1 base with x in an F_p-subspace of
   F_{p^n}) to test whether the A-fold of 5 survives there. The census
   predicts it does, because ζ preserves any ζ-stable x-subspace.
2. Build an explicit Richelot isogeny to get an index-4 image arm
   different from [2]J.
3. Feed the structured-pair correction 1 + (h−1)(|A|−2)/(2|F|) back into
   §15.2: it is a small, free m = 2 gain whenever A fixes the cofactor
   classes.

## What should move into `src/`

The example carries a self-contained genus-2 toolkit over F_p (u64): Mumford
pairs, Cantor addition with reduction, ζ-action, exhaustive Mumford
enumeration with per-candidate verification, and #J from L(1) using
#C(F_p) and #C(F_{p²}). It is validated against `jv_cover::Hyp`. A
`cryptanalysis::genus2_fp` module with these pieces would let later leads
reuse them instead of copying the example.
