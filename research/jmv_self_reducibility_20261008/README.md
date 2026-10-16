# JMV self-reducibility of ECDLP across an isogeny class: what is proved, what is open, and one obstruction

Status: **literature verdict, one proved obstruction, pilot toy evidence;
protocols J1 to J5 preregistered in [`PROTOCOL.md`](PROTOCOL.md), none run**

Date: 2026-10-08

Result class: **structural and asymptotic; no ECDLP speedup.**  `S` and
every ECDLP cost figure are **unset**.  Nothing here changes the cost of any
discrete logarithm.  The subject is whether DLP instances can be moved
between curves of the same order, and at what cost.

## The question

Two ordinary curves over `F_q` with the same number of points are
isogenous (Tate), so "same order" means "same isogeny class".  The class
splits into levels by endomorphism ring, arranged as ℓ-volcanoes for the
primes `ℓ` dividing the conductor `f_π` of `Z[π]` in `O_K`.  Jao, Miller and
Venkatesan asked whether DLP is equally hard across the class [JMV05].

## Verdict

| statement | status | source |
|:--|:--|:--|
| Random self-reducibility within one level | THEOREM under GRH | [JMV05] |
| Pairwise explicit isogeny within one level | `Õ(q^{1/4})` classical; polynomial given the connecting ideal | [Gal99, GHS02]; [PR23] |
| Across levels, conductor gap with all primes polynomial in `log q` | THEOREM under GRH | [JMV05] |
| Across a whole class, flat volcano or tall with small crater | rigorous `Õ(q^{1/4})` isogeny between any two curves | [Gal24] |
| Across a whole class, worst case | heuristic `Õ(q^{0.4})` | [Gal24] |
| Polynomial-time vertical step at a large prime `ℓ \| f_π` | OPEN | this note |

Consequences.  (1) A disproof in the strong sense, two same-order curves
whose DLP costs differ by more than `Õ(q^{1/4})`, is excluded for random
curves and for pairing-style classes with small crater.  (2) A disproof of
the polynomial claim would need a lower bound on computing a vertical
isogeny of large prime degree; none exists, and output-size bounds do not
apply because higher-dimensional representations never write the kernel
polynomial.  (3) A proof of the polynomial claim needs exactly one new
algorithm: the vertical step below.

Corollary for registered P-256 (from PR #1565's factorisation `−D_π = 3 ·
5 · q₁q₂q₃`, squarefree): `f_π = 1`, the class is a single level, and the
within-level theorem covers the entire class under GRH.  No vertical step
exists there.  Pairwise explicit isogenies inside that class cost
`Õ(h^{1/2}) ≈ 2^{64}` class-group operations by [GHS02], far below the
`2^{128}` DLP cost, which is a statement about equal hardness and not an
attack.

## The vertical step, made precise

Fix `ℓ` with `v_ℓ(f_π) = 1`, a surface curve `E₀` with `End = O_K`, and a
floor curve `E₁` one ℓ-step below, `End(E₁) = O′ = Z + ℓO_K` locally at `ℓ`.
Write `λ ∈ Z` with `λ ≡ t/2 (mod ℓ)`.

- THEOREM (Kohel).  Every `F_p`-isogeny `E₀ → E₁` has degree divisible by
  `ℓ`.  The gap between the levels is one vertical ℓ-isogeny.
- THEOREM.  On the floor, `π` is not scalar on `E₁[ℓ]`; its single
  eigenline is the unique `F_p`-rational ℓ-subgroup `G`, and `E₁/G ≅ E₀`.
  In ideal terms the ascending isogeny is `φ_I` for the kernel ideal
  `I = I(G) = {α ∈ O′ : α(G) = 0}`.
- LEMMA (ℓ² obstruction).  `I(G) = ℓ·O_K`.  Hence every endomorphism of
  `E₁` that kills `G` has degree divisible by `ℓ²`.

Proof of the lemma.  `π − λ` has trace `t − 2λ ≡ 0 (mod ℓ)` and norm
`((2λ − t)² − D_π)/4 ≡ 0 (mod ℓ²)`, so `γ = (π − λ)/ℓ` has integral trace
and norm and lies in `O_K`.  If `γ ∈ Z + ℓO_K` then `π ∈ Z + ℓ²O_K`, which
contradicts `v_ℓ(f_π) = 1`; so `Z + Zγ + ℓO_K` strictly contains `Z + ℓO_K`,
which has index `ℓ` in `O_K`, and therefore equals `O_K`.  Thus `(ℓ, π −
λ)O′ = ℓ(O′ + γO′) = ℓO_K`.  The kernel ideal `I(G)` contains `ℓ` and `π −
λ`, so it contains `ℓO_K`, which has index `ℓ` in `O′`; and `O′/I(G)` embeds
in `End(G) = F_ℓ`, so `I(G)` has index at least `ℓ`.  Hence `I(G) = ℓO_K`
and every `α ∈ I(G)` is `ℓβ` with `β ∈ O_K`, of norm `ℓ²N(β)`.  ∎

What the lemma rules out.  The polynomial-time route that closed the
analogous problems elsewhere, in QFESTA, SQIsign2D and Clapoti, starts
from an endomorphism `α` with `N(α) = ℓ · d`, `gcd(ℓ, d) = 1`, whose kernel
contains the target subgroup, and extracts the degree-`ℓ` factor by Kani's
lemma in dimension 2.  On the floor no such `α` exists.  The same holds
for maps from `E₁` to other floor curves, which are the only other maps
Clapoti can evaluate, under the proof obligation recorded as J5.3.  Any
higher-dimensional route to a vertical step must therefore start from
torsion images, which live over `F_{p^r}` with `r = ord(t/2 mod ℓ)`, up to
`ℓ − 1`.  This is a negative result for one method class under the stated
model, not for the problem.

Where the vertical step reappears after Clapoti (derivation, to be tested
as J5).  For two floor curves `E₁, E₁′` under the same surface curve `E₀`,
the composite `ψ = φ̂₁′ ∘ φ₁ : E₁ → E₁′` has cyclic kernel of order `ℓ²`
containing `G`, and its kernel ideal `(ℓ², π − μ)` is an invertible
`O′`-ideal exactly when `v_ℓ(N(π − μ)) = 2`.  Counting shows `ℓ − h` of the
`ℓ` cyclic order-`ℓ²` groups over `G` give invertible ideals, where `h = 1 +
(D_K/ℓ)` is the number of horizontal edges at `E₀`; the other `h` are the
composites that land on surface neighbours of `E₀`.  So floor-to-floor
passage through the surface is a horizontal degree-`ℓ²` action, evaluable
in polynomial time by [PR23], and the vertical step becomes: *given
efficient representations of two such composites `ψ, ψ″` with `ker ψ ∩ ker
ψ″ = G`, compute an efficient representation of `φ_G`.*  The naive route
through the ideal sum `I(ψ) + I(ψ″) = ℓO_K` hits the lemma again.

## Pilot toy evidence (run before this note was written; not preregistered)

[`jmv_floor_toy.sage`](jmv_floor_toy.sage), Sage 10.10, seed 20261008.
Surface curves are built from Hilbert class polynomials, so `End = O_K` is
known by construction; one descending step gives the floor.  Raw data:
[`results/pilot_results.jsonl`](results/pilot_results.jsonl).

| check | instances | result |
|:--|--:|:--|
| DLP on the floor transfers through the ascending ℓ-isogeny and pulls back | 13 | 13 pass |
| `E₀[ℓ]` rational exactly over `F_{p^r}`, `r = ord(t/2 mod ℓ)` | 13 | 13 pass |
| ascending kernel equals `E₁[ℓ] ∩ ker(π − λ)` | 13 | 13 pass |

Degrees `ℓ ∈ {3, 5, 7, 11, 13, 17, 19}`, `p ≤ 1013`, `D_K ∈ {−8, −11, −20}`.
One `ℓ = 17` instance at `p = 659` aborted in the torsion-basis sampler
(60 draws were too few when `E(F_{p^8})[17^∞]` is not elementary); the
shipped instrument draws 600 and J1 reruns it.  The time to classify the
`ℓ + 1` edges through `F_{p^r}` torsion ran from 0.12 s at `r = 1` to
45.3 s at `ℓ = 17, r = 16`, the poly(`r`) cost in miniature.  TOY-EVIDENCE.

## Next steps

- Conservative (J1, J3): extend the pilot to depth-2 volcanoes and measure
  the torsion route's cost in `r` and `ℓ` at CM-constructed classes with a
  large prime conductor.
- Representation-changing (J4, J5): replicate [Gal24]'s `q^{1/4}` vertical
  step at toy sizes, and test the invertible-ideal derivation above, which
  turns the vertical step into a gcd-of-isogenies problem.
- Speculative (J5.4): attack that gcd problem from torsion images in
  dimension 2 with the curve's actual small-extension `2^e`-torsion,
  proof obligation first, code after.

## References

- [JMV05] D. Jao, S. Miller, R. Venkatesan, *Do all elliptic curves of the
  same order have the same difficulty of discrete log?*, ASIACRYPT 2005;
  arXiv math/0411378.
- [JMV09] *Expander graphs based on GRH with an application to elliptic
  curve cryptography*, J. Number Theory 2009.
- [Gal99] S. Galbraith, *Constructing isogenies between elliptic curves
  over finite fields*, LMS JCM 2 (1999).
- [GHS02] S. Galbraith, F. Hess, N. Smart, *Extending the GHS Weil descent
  attack*, EUROCRYPT 2002.
- [Gal24] S. Galbraith, *Climbing and descending tall isogeny volcanos*,
  ePrint 2024/924; Res. Number Theory 2025.
- [PR23] A. Page, D. Robert, *Introducing Clapoti(s): evaluating the
  isogeny class group action in polynomial time*, ePrint 2023/1766.
- [JS10] D. Jao, V. Soukharev, *A subexponential algorithm for evaluating
  large degree isogenies*, ANTS IX (2010).
- [Rob22] D. Robert, *Evaluating isogenies in polylogarithmic time*,
  ePrint 2022/1068; and *On the efficient representation of isogenies*,
  ePrint 2024/1071.
- [Sut13] A. Sutherland, *Isogeny volcanoes*, ANTS X (2013).
- [Koh96] D. Kohel, *Endomorphism rings of elliptic curves over finite
  fields*, PhD thesis, Berkeley 1996.
