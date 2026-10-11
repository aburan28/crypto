#set page(paper: "a4", margin: (x: 2cm, y: 2cm), numbering: "1")
#set text(size: 10pt)
#set par(justify: true)
#set heading(numbering: "1.")
#show heading: set block(above: 1.2em, below: 0.6em)

#let thm(name, body) = block(
  width: 100%, inset: 8pt, radius: 3pt, stroke: 0.5pt + rgb("#2f5d9b"), fill: rgb("#f4f7fc"),
  [*#name* #body],
)
#let proof(body) = block(width: 100%, inset: (left: 8pt), [_Proof._ #body #h(1fr) $square$])

#align(center)[
  #text(size: 15pt, weight: "bold")[Lead 16.7: isogeny oracles in the algebraic group model]
  #v(2pt)
  Round 4 of the isogeny-fibers study · 2026-10-11 · theory round with a native toy sanity check \
  Source note: `RESULTS.md` (full statements and proofs) · raw output: `out_seed7.json`
]

#block(inset: 8pt, fill: rgb("#fff8e6"), radius: 3pt, width: 100%)[
  *Status labels.* Results called Theorem, Proposition, Lemma or Corollary are proved here or in `RESULTS.md`.
  _Heuristic_ means the class-factor uniformity heuristic of §14.1. _Conjecture_ means not proved.
  Numbers in §6 are *measured* (exact counts on one toy curve pair) or *derived* (ratios of measured counts).
  Toy curves only; no sub-rho claim; no claim about a deployed curve.
]

= Verdict

- *Model.* AGM-Iso(Ψ) is defined in §2. The adversary reads coordinates and may compute anything. It also has oracles for a fixed set Ψ of rational isogenies and their duals across one class, and every point it _submits to an oracle_ must come with a representation over the points it has received.
- *Theorem 1.* Any such adversary runs without oracles, with exactly the same success probability. The additive cost is O(Q(ℓ + log r)), and *no isogeny is ever evaluated*. The simulator needs only the class-joint instance (P#sub[j], xP#sub[j]). *Theorem 2* gives the converse, so the oracles are worth exactly the transported challenges.
- *Corollary 3.* If transport along a spanning tree is computable, oracle-DLP on E#sub[1] reduces to plain DLP on E#sub[1];. The loss is |J|−1 tree evaluations plus O(Q(ℓ + log r)). Non-tree isogenies, even of huge degree, enter only through a scalar.
- *Limits.* For DLP the AGM is vacuous without oracles: every algorithm outputs a scalar, so every algorithm is algebraic. The model has no lower-bound content and *cannot separate index calculus from rho*. Without the algebraic restriction on queries, the oracle becomes a transport oracle for factor-base points (Prop. 5). Simulating it then costs Q × (degree-dependent evaluation).
- *Summation polynomials.* Isogenies push decompositions forward exactly (Thm 6). *The kernel-saturation law of §14.1(b) is the image law of §14.1(a) on the image curve* (Thm 7). Image bases need no oracle (Prop. 8).
- *Constants vs exponents.* The oracle can only choose the representative: oracle-DLP ≡ CJ-DLP ≤ min#sub[j] DLP(E#sub[j]). "Constants only" is *false in general*, because GHS classes over F#sub[2#super[n]] with composite n contain both weak and non-weak members. It is *Conjecture 10* over prime fields.
- *Sanity check (measured).* The simulator answered 5000/5000 algebraic queries correctly after one evaluation of φ. All fiber-sum, injective and saturated identities held, with 0 mismatches. The exact equality D#sub[m];(G) = D#sub[m];(G′) held (6195 = 6195; 360 = 360). The heuristic ratios came out at 1.97–2.05 against a prediction of 2, and at 0.509 and 0.255 against predictions of 1/2 and 1/4.

#figure(image("reduction_flow.svg", width: 100%), caption: [Diagram of the reduction of Theorem 1, the class-joint instance (Thm 2), the computable-transport branch (Cor. 3) and the unrestricted query that breaks the bookkeeping (Prop. 5). This is a diagram of the proof, not data. For this theory round it is the required visual.])

= Setting and model

F#sub[q] is a finite field. E#sub[1], …, E#sub[k] are ordinary curves in one F#sub[q]-isogeny class with N = h·r, r prime, r ∤ h. G#sub[j] = E#sub[j];(F#sub[q])[r] has generator P#sub[j];. τ#sub[ij]: G#sub[i] → G#sub[j] is the abstract isomorphism P#sub[i] ↦ P#sub[j];. Ψ is a finite set of F#sub[q]-rational isogenies, closed under duals and *fixed before the instance is drawn*. The challenge is X = xP#sub[1];. CJ-DLP#sub[J] is the problem: given (P#sub[j], X#sub[j] = xP#sub[j]) for j ∈ J, find x.

#thm("Definition 1 (AGM-Iso(Ψ)).")[
The adversary A is any algorithm. It receives the curves, the P#sub[j], the names of the isogenies in Ψ and X, all as coordinates. The game keeps a list L#sub[j] of the points of E#sub[j] that A has received; initially L#sub[1] = (P#sub[1], X) and L#sub[j] = (P#sub[j]). A query O(ψ, Z, z) returns ψ(Z) and appends it to L#sub[t(ψ)];. A is *oracle-algebraic* if every query has Z = Σ#sub[t] z#sub[t] L#sub[s(ψ)];[t]. A outputs x\* ∈ Z/r. Its internal computation is unconstrained: rho, index calculus, Semaev polynomials, Weil descent, Gröbner bases and its own isogeny evaluations are all allowed. *COM(Ψ)* is the same model without the algebraicity requirement.
]

This is the Fuchsbauer–Kiltz–Loss AGM, applied to a family of groups with one representation list per curve. Oracle queries count as outputs, as in AGM-plus-ROM analyses. With Ψ = ∅ every DLP algorithm is algebraic, so in the AGM DLP hardness is an assumption and never a conclusion.

= Elliptic lemmas

#thm("Lemma 1.")[For F#sub[q]-rational ψ: E#sub[i] → E#sub[j], ψ(G#sub[i]) ⊆ G#sub[j] and ψ = λ#sub[ψ]·τ#sub[ij] on G#sub[i] for a unique λ#sub[ψ] ∈ Z/r. If r ∤ deg ψ, then λ#sub[ψ] ≠ 0.]
#proof[ψ commutes with π, so it maps rational points to rational points and r-torsion to r-torsion. A homomorphism between cyclic groups of order r is a scalar. If r ∤ deg ψ, then ker ψ meets G#sub[i] trivially.]

#thm("Lemma 2 (scalar calculus).")[(a) λ#sub[ψ₂ψ₁] = λ#sub[ψ₂]λ#sub[ψ₁];. (b) λ#sub[ψ̂]·λ#sub[ψ] = deg ψ mod r. (c) If α ∈ End(E) for ordinary E, with α = (a + bπ)/f and r ∤ f, then λ#sub[α] = (a + b)/f mod r. (d) Normalise P#sub[j] := γ#sub[j];(P#sub[1]) along a spanning tree with r ∤ deg on its edges. Then tree edges have λ = 1, their duals have λ = deg, and any other ψ has λ#sub[ψ] = λ#sub[α]/deg γ#sub[j] with α = γ̂#sub[j]ψγ#sub[i] ∈ End(E#sub[1]).]
#proof[(a) Evaluate both sides at P. (b) ψ̂ψ = [deg ψ]. (c) fα = a + bπ, and π is the identity on G, which consists of rational points. (d) Combine (a), (b) and (c).]

#thm("Lemma 3 (cofactor split).")[Let s·h ≡ 1 mod r and e#sub[G] = sh mod N. Then [e#sub[G]]Z ∈ G, Z − [e#sub[G]]Z ∈ E(F#sub[q])[h], and ψ(Z) = λ#sub[ψ]τ([e#sub[G]]Z) + ψ(Z − [e#sub[G]]Z).]

_Remark (hypothesis implicit in Thm 15.1)._ [x] is a rational isogeny of degree x² with λ = x. The oracle family must therefore be independent of the instance. Thm 15.1 meets this because the adversary must name its isogenies.

= Main results

#thm("Theorem 1 (oracle elimination).")[Let A be oracle-algebraic in AGM-Iso(Ψ), with Q queries, representation length ≤ ℓ, time T and success ε. Let J be the set of curves touched by Ψ, together with E#sub[1];. Then there is an oracle-free B, given the instance-independent scalars λ#sub[ψ], that solves CJ-DLP#sub[J] with success exactly ε. B runs in time T + O(Qℓ) operations in Z/r + O(Q log r) group operations, and never evaluates any ψ.]
#proof[B stores (u#sub[t], v#sub[t]) with L#sub[j];[t] = u#sub[t]P#sub[j] + v#sub[t]X#sub[j];. For a query (ψ: E#sub[i] → E#sub[j], Z, z) it computes (u, v) = Σ z#sub[t];(u#sub[t], v#sub[t]). By algebraicity Z = uP#sub[i] + vX#sub[i];. By Lemma 1 and τ(X#sub[i]) = xτ(P#sub[i]) = X#sub[j], ψ(Z) = λ#sub[ψ];(uP#sub[j] + vX#sub[j]). B returns this point and stores (λu, λv). Every answer equals the true oracle's, so A's view is identical to the real game.]

#thm("Theorem 2 (converse).")[If Ψ contains a spanning tree of J with r ∤ deg on its edges, and the generators are normalised as in Lemma 2(d), then CJ-DLP#sub[J] reduces to oracle-DLP using |J|−1 queries. Hence oracle-DLP on E#sub[1] and CJ-DLP#sub[J] are equivalent.]
#proof[Query the tree edges outward from X. This returns γ#sub[j];(X) = X#sub[j];.]

#thm("Corollary 3 (computable transport).")[If the tree kernels are known, oracle-DLP on E#sub[1] reduces to DLP on E#sub[1] with success preserved. The additive loss is Σ#sub[tree] C#sub[eval] + O(Q(ℓ + log r)), where C#sub[eval];(ℓ) = O(ℓ) with Vélu or Õ(√ℓ) with √élu per prime step. Kernel finding is instance-independent preprocessing.]

#thm("Corollary 4 (horizontal classes; heuristic path finding).")[If every curve of J has the same endomorphism ring, a path E#sub[1] → E#sub[j] is instance-independent. It is found in heuristic Õ(q#super[1/4]) time (Galbraith 1999; Galbraith–Hess–Smart 2002), after which the per-instance loss is polynomial. Jao–Miller–Venkatesan reduces to a _random_ curve of the level, not to a given E#sub[j], so it is not what is used here. Vertical moves across a conductor gap with a large prime factor are where Theorem 1 ends at CJ-DLP.]

#thm("Proposition 5 (unrestricted queries).")[COM(Ψ) is equivalent to oracle-free coordinate access plus transport oracles τ#sub[ij] on _arbitrary_ points of G#sub[i], plus constant-size instance-independent tables on E(F#sub[q])[h] (Lemma 3). With computable transport, COM-DLP reduces to DLP on E#sub[1] with loss Q·max#sub[j] C(γ#sub[j]) + O(Q log N), which grows with Q × degrees. Unlike AGM-Iso, COM moves non-algebraic points such as factor-base points between curves. Whether this ever helps beyond CJ-DLP is *open*; §14.4 gives a heuristic argument that it is neutral at best.]

#table(
  columns: (auto, auto, 1fr, auto),
  stroke: 0.4pt, inset: 4pt,
  [*model*], [*reduces to*], [*additive loss*], [*where*],
  [AGM-Iso(Ψ)], [CJ-DLP#sub[J]], [O(Q(ℓ + log r)); no degree dependence], [Thm 1, tight by Thm 2],
  [AGM-Iso, tree known], [DLP on E#sub[1]], [Σ#sub[tree] C#sub[eval] + O(Q(ℓ + log r))], [Cor 3],
  [AGM-Iso, horizontal], [DLP on E#sub[1]], [Õ(q#super[1/4]) once per curve (heuristic), then poly per instance], [Cor 4],
  [COM(Ψ), tree known], [DLP on E#sub[1]], [Q·max C(γ#sub[j]) + O(Q log N)], [Prop 5],
  [COM(Ψ), no transport], [none known], [—], [open],
  [GGM + isogeny oracles], [lower bound], [Ω(√r)], [Thm 15.1],
)

*Where the elliptic structure is used.* The homomorphism property linearises the answers. Prime order and rationality of G make ψ|G a scalar times a transport. The dual identity ψ̂ψ = [deg] makes duals free. Because π fixes G pointwise, the eigenvalues of endomorphisms in the class are computable. The idempotent split E(F#sub[q]) = G ⊕ E[h] separates a tabulable part from the transport.

= Summation polynomials in the model

Let φ: E → E′ be rational, M = ker φ ∩ E(F#sub[q]), μ = |M|. For a base F, d#sub[F]#super[(m)];(R) counts ordered m-tuples of F summing to R, and w#sub[F];(Q′) = \#(φ#super[−1];(Q′) ∩ F).

#thm("Theorem 6 (pushforward).")[For every Q′: Σ#sub[φ(R)=Q′] d#sub[F]#super[(m)];(R) = (w#sub[F])#super[\*m];(Q′). (a) If φ is injective on F: Σ#sub[t∈M] d#sub[F];(R + t) = d#sub[φ(F)];(φ(R)). (b) If F + M = F: d#sub[F];(R) = μ#super[m−1]·d#sub[φ(F)];(φ(R)) for every R.]
#proof[d#sub[F] = 1#sub[F]#super[\*m], and pushforward along a homomorphism commutes with convolution. In (b), translating one summand by t ∈ M is a bijection between decompositions of R and of R + t.]

#thm("Theorem 7 (saturation law = image law transported).")[Assume r ∤ deg φ and F + M = F with |F| = B. (i) Exactly, D#sub[m]#super[G];(F) = D#sub[m]#super[G′];(φ(F)). (ii) Heuristically, φ(F) sits in an image of index μ on E′, so c(φ(F)) ≈ μ·c#sub[plain];. (iii) Hence c(F) ≈ c#sub[plain]/μ#super[m−1], and the mitm cost factor is (1/μ)·μ#super[m−1] = μ#super[m−2];. These are Prop. 14.1(b) and the |M|#super[m−2] term of Prop. 15.2.]

#thm("Proposition 8 (image bases need no oracle).")[Let H = ψ(E′(F#sub[q])) with r ∤ deg ψ. Then R ∈ H iff (1 − e#sub[G])R ∈ H ∩ E(F#sub[q])[h], a fixed subgroup of size h/i. It equals i·E(F#sub[q])[h] when that group is cyclic. Membership therefore costs one scalar multiplication and a lookup, so the image gain i, the automorphism fold |A| and the saturation penalty are all coordinate-level effects that need no oracle.]

#thm("Proposition 9 and Example (constants vs exponents).")[In AGM-Iso(Ψ), cost(oracle-DLP) ≤ min#sub[j∈J] cost(DLP on E#sub[j]) + (|J|−1) queries. The oracle therefore moves an exponent iff the class has members with different best exponents. *Example:* over F#sub[2#super[n]] with composite n, GHS-weak and non-weak members coexist in one class (Gaudry–Hess–Smart 2002; Galbraith–Hess–Smart 2002), so "only constants" is false in general. In GHS02 the transport was computable by an isogeny walk, so there the oracle only removes the walk cost.]

#thm("Conjecture 10 (prime fields).")[For ordinary classes over F#sub[p], the best known exponent is a class invariant. Oracles then change the cost of known attacks by at most √3 for rho, and by the Prop. 15.2 factor for summation-polynomial methods. Evidence: Prop. 3.2 and Thm 15.1. Proving it would need lower bounds for non-generic algorithms.]

= Computational sanity check (measured)

`cargo run --release --example r4_16_7_agm_pushforward -- --pmin 1048576 --queries 5000 --seed 7`, single-threaded; the contended wall time of about 4 s is not a result. The curves are E: y² = x³ + 2x² + 2x over F#sub[p], p = 1048583, with N = 4·262321 (`ICV1:fp-1048583:-700:1049284:128:unk:unk:r:d62ba3a6bad8`), and E′ = E/⟨(0,0)⟩: y² = x³ − 4x² − 4x (`ICV1:fp-1048583:-700:1049284:10976:unk:unk:r:07d957b2581c`). φ and φ̂ are the explicit 2-isogeny formulas of Silverman, AEC III.4.5. Their slugs are not yet registered.

#table(
  columns: (1fr, auto),
  stroke: 0.4pt, inset: 4pt,
  [*check*], [*result (measured)*],
  [maps on curve, homomorphism, φ̂φ = [2], φφ̂ = [2] (2000 random points each)], [2000/2000 each],
  [λ#sub[φ̂] = deg φ/λ#sub[φ] = 2 on G, with P′ = φ(P)], [holds],
  [AGM simulation: random oracle-algebraic queries (2495 φ, 2505 φ̂)], [*5000/5000* agree],
  [simulator isogeny evaluations: per instance / instance-independent / naive], [1 / 1 / 5000],
  [cofactor split φ(Z) = φ(Z#sub[G]) + T(Z#sub[h]), 2000 random Z], [2000/2000],
  [fiber-sum identity, 6 base/m combinations, every Q′ in the support], [0 mismatches],
  [injective law (m = 2: 95894 R; m = 3: 15540 R)], [holds at all],
  [saturated law (m = 2: 24726 R; m = 3: 1408 R)], [holds at all],
  [D#sub[m];(G) saturated on E vs its image on E′ (exact prediction: equal)], [6195 = 6195; 360 = 360],
)

#figure(image("class_factor_ratios.svg", width: 92%), caption: [Ratios of class factors c = D#sub[m];(G)·m!/|F|#super[m], derived from exact census counts recorded in `out_seed7.json` (check 5), against the heuristic predictions (dashed). One curve pair and one seed; no error bars are available. Base sizes: 316/158 at m = 2 and 36/32/16 at m = 3.])

The plain base happened to be φ-injective at this x-range, so the "plain" and "injective" rows coincide. This is recorded as it fell, not engineered.

= What this does not say

- It is not a hardness result. The AGM gives no lower bound for DLP, and Thm 15.1 gives Ω(√r) only generically.
- It does not separate index calculus from rho; both are admissible adversaries.
- It does not cover non-algebraic oracle queries unless the transport is computable (Prop. 5).
- It does not settle CJ-DLP against DLP on E#sub[1] when transport is not computable (large-prime conductor gaps).
- "Only constants" is not a theorem. It is false for GHS classes in characteristic 2 and a conjecture over prime fields.
- The class-factor ratios are heuristic and come from one toy curve pair. Theorem 6, Theorem 7(i) and Proposition 8 are exact.

= Related work

*Confident:* Shoup (EUROCRYPT 1997); Nechaev (1994); Maurer (IMA C&C 2005); Fuchsbauer–Kiltz–Loss, _The algebraic group model and its applications_ (CRYPTO 2018); Fuchsbauer–Plouviez–Seurin (EUROCRYPT 2020, AGM + ROM); Bauer–Fuchsbauer–Loss (CRYPTO 2020); Semaev (ePrint 2004/031); Gaudry (J. Symb. Comput. 2009); Diem (Compositio 2011); Galbraith–Gaudry, _Recent progress on the ECDLP_ (DCC 2016); Gaudry–Hess–Smart (J. Cryptology 2002); Galbraith–Hess–Smart (EUROCRYPT 2002); Galbraith (LMS JCM 1999); Jao–Miller–Venkatesan (ASIACRYPT 2005); Vélu (1971); Bernstein–De Feo–Leroux–Smith (ANTS 2020); Silverman, AEC.
*Recalled, details unverified:* Menezes–Teske–Weng, _Weak fields for ECC_ (CT-RSA 2004); Agrikola–Hofheinz–Kastner (EUROCRYPT 2020); Zhandry, _To label, or not to label_ (CRYPTO 2022); Duman–Hartmann–Kiltz–Kunzweiler–Lehmann–Riepel, _Generic models for group actions_ (PKC 2023), which concerns group actions on curves, not point-evaluation oracles.
*Not searched:* no web search was run. We know of no prior AGM formulation with point-evaluation isogeny oracles, but that is a statement of the author's knowledge only.
