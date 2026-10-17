# Isogeny-computation algorithms: what exists, what this crate implements, what is missing

Scope: algorithms that *compute* an isogeny (or the data needed to), in the sense of the request
"all the approaches to generating isogenies between curves (Couveignes, Kohel, Galbraith, ...)".
This is a working survey for planning, not a literature review.

**Provenance of the entries.** Marked ✔ in the "checked" column: the statement (algorithm, input,
output, complexity) was read in this session from the source — the Bostan–Morain–Salvy–Schost paper
(arXiv cs/0609020, full text) for everything in §P2, and the abstracts/bibliographic records
retrieved by web search for Couveignes (ANTS-II 1996, LNCS 1122, pp. 59–65), Couveignes–Morain
(ANTS-I 1994), Galbraith–Stolbunov (arXiv 1105.6331) and Costello–Hisil (ASIACRYPT 2017). All other
references are cited from memory; their details have **not** been re-verified. Where a formula was
taken from memory, the test named in the last column is what establishes that the code is right,
not the citation.

**Status vocabulary** (following the repository's reporting rules): *implemented* = code + a test that
checks the output independently; *partial* = implemented with a stated restriction; *not implemented
– implementation gap* = nothing prevents it but effort; *not implemented – resource limit* = needs
hardware or a model of computation that is not available here (e.g. a quantum computer). No entry
is withheld for any other reason.

**Deliveries.** `V1` = first delivery, `V2` = kernel-only/BMSS extension (PR #1526, first commits),
`V3` = the speed and coverage round (arithmetic, big fields, char 2, quaternions, genus 2, SEA,
radical isogenies, relation lattice, curve models), `V4` = the fixed-degree round (Benčina–Kutas–Merz–
Petit–Stopar–Weitkämper's algorithms end to end, the remaining √élu paper algorithms, CSURF; study:
`research/fixed_degree_isogenies_20261016/`). Module paths are under `src/`.

## The problems, and which algorithms solve which

| | Problem | Typical algorithms |
|---|---|---|
| P1 | kernel → isogeny (codomain and map) | Vélu, Kohel, √élu, Montgomery/Edwards/Huff x-only formulas, chains with strategies, radical isogenies |
| P2 | two ℓ-isogenous curves (and maybe σ) → isogeny | the BMSS family: Stark, Elkies 1992/1998, Atkin, fastElkies, … ; small characteristic: Couveignes 1996, Lercier |
| P3 | curve + ℓ → its ℓ-isogenies (kernels / neighbours) | Φ_ℓ roots + Elkies, division-polynomial factoring, isogeny cycles (Atkin primes); application: SEA |
| P4 | two curves → *some* isogeny between them | Galbraith, GHS, Galbraith–Stolbunov, Kohel (volcanoes), Couveignes (hard homogeneous space), Delfs–Galbraith, quantum, KLPT/Deuring (supersingular, endomorphism-ring based) |
| P5 | auxiliary computations | dual isogeny, End(E) (ordinary: Kohel; supersingular: quaternions), class-group structure, modular polynomials |
| P6 | higher-dimensional | Richelot, splitting/gluing, theta (2,2)/(ℓ,ℓ), Kani-lemma methods |

## Status table

### P1 — kernel → isogeny

| Algorithm | Ref | Status | Where | Checked by |
|---|---|---|---|---|
| Vélu, odd cyclic kernel | Vélu 1971 | implemented V1 | `kernel/velu.rs` | agrees with Kohel/√élu; homomorphism test |
| Vélu, **any** finite subgroup (even degree, 2-torsion, non-cyclic) | Vélu 1971 | implemented V2 | `kernel/velu.rs::velu_general` | E/E[2] has the j of E; Φ₂; homomorphism; cyclic n = 4…15 vs Kohel |
| Kohel formulas, odd | Kohel 1996 | implemented V1 | `kernel/kohel.rs` | = Vélu |
| Kohel formulas, even / 2-torsion factor | Kohel 1996 | implemented V2 | `kernel/kohel.rs` | = Vélu for degrees 2, 4, 6, 8, 9, 10, 12, 15; kernel E[2] |
| √élu, short Weierstrass | Bernstein–De Feo–Leroux–Smith 2020 | implemented V3 (product trees, remainder-tree multipoint evaluation; three power sums) | `kernel/sqrt_velu.rs` | = Vélu up to ℓ = 100003 |
| √élu, Montgomery (I ± J ∪ K, E_J product tree, codomain from h_S(±1)) | BDFLS 2020 | implemented V3 | `kernel/sqrt_velu_mont.rs` | = Vélu on all 74 CSIDH-512 primes and at ℓ = 1009, 10007 |
| √élu §2: Strassen's modular factorial ℓ! mod n and the deterministic factorisation built on it (gcd(n, ℓ! mod n), binary search over ℓ ≤ √n) | BDFLS 2020 §2; Strassen 1976 | implemented V4 (u64 n, any n incl. even/composite; own remainder tree with precomputed reciprocals) | `strassen.rs` | = naive product on 528 (ℓ, n) pairs; Wilson at 2²⁰; factoring = `factor_u64` on 29 n ≤ 2⁴⁰. Karatsuba trees give O(ℓ^0.79), not Õ(√ℓ) (needs quasi-linear products); crossover ℓ ≈ 10⁴–10⁵ |
| √élu §3, Algorithm 1: h_S(α) = ∏(α − ζˢ) for S = I + J ∪ K (multiplicative version), incl. the cyclotomic case h_S = Φ_n | BDFLS 2020 §3 | implemented V4 (`for_set` index systems for arbitrary S are greedy, not Õ(√#S) for unit groups with many factors) | `kernel/sqrt_mult.rs` | = direct product over Zp, FpM<2>, FpM<8>; Φ_n for 11 n vs integer Φ_n |
| √élu §4.14: irrational kernel generators (P over F_{q^e}, kernel Galois-stable; orbit index systems, F_q-coefficient E_J, norms down) | BDFLS 2020 §4.14 | implemented V4 (Montgomery, F_q input points; slower than Vélu over the extension at the sizes tried, as the paper warns) | `kernel/sqrt_velu_irr.rs` | = Montgomery Vélu over F_{q^e} for e = 2, 3, 4; j = Weierstrass x-only Vélu |
| √élu Appendix A.5: fused Q(±1), Q(α), Q(1/α); remainder trees with precomputed reciprocals | BDFLS 2020 App. A.5 | implemented V4 (fast middle products and merged reductions not done) | `kernel/sqrt_velu_mont.rs` (`fused_codomain_eval`, `RemTreeSqrtVelu`) | = Vélu on CSIDH-512 and at ℓ = 100003 (256 bits) |
| Vélu from x(P) only (no y; kernels over extension fields) | folklore; division-polynomial values | implemented V2, projective + batch inversion V3 | `kernel/xonly.rs` | = point-based Vélu; Galois-stable kernel over F_p² = Kohel over F_p |
| Montgomery x-only Vélu, affine and projective (A24 : C24) | Costello–Hisil 2017 ✔(abstract); Renes 2018; Meyer–Reith 2018 | implemented V2 (affine), V3 (projective, shared (X±Z) for codomain and pushes) | `kernel/montgomery.rs` | j-match with Weierstrass Vélu; CSIDH-512 batched = stepwise = tree |
| ℓⁿ kernel as a chain; naive / balanced / cost-model strategies; composite smooth kernels | De Feo–Jao–Plût 2014; Costello–Longa–Naehrig 2016 | implemented V2 | `kernel/chain.rs` | SIDH key exchange over F_{p²} agrees for all four strategies |
| Twisted Edwards ℓ-isogenies (translate-product definition, explicit map, one-variable X map) | Moody–Shumow 2016 | implemented V3 (odd ℓ; the explicit and X-only maps for a = 1) | `kernel/models.rs` | codomain j = Montgomery Vélu; images on the codomain; homomorphism; three evaluation forms agree |
| Huff-curve ℓ-isogenies | Moody–Shumow 2016 | implemented V3 (odd ℓ) | `kernel/models.rs` | codomain j = Weierstrass Vélu; images on the codomain; homomorphism |
| Radical 3-, 5- and 7-isogenies | Castryck–Decru–Vercauteren 2020; N = 7 derived in this work | implemented V3 (N = 3 flex model, N = 5 Tate normal form, N = 7 via the X₁(7) parameter; unique N-th root when gcd(N, q − 1) = 1) | `kernel/radical.rs` | on CSIDH-512, k steps = CSIDH action 𝔩₃ᵏ / 𝔩₅ᵏ / 𝔩₇ᵏ (k = 1..6); over 40 bits steps are Φ_N-adjacent and non-backtracking |
| └ N = 7 derivation | this work (no published formula used) | X₁(7) order-7 locus b² − bc − c³ = 0, parametrised b = t³−t², c = t²−t; radicand ρ(t) = t(t−1)²; t′ = N(t,α)/D(t) with α = ρ^{1/7}, D(t) = 1+4t−13t²+9t³−t⁴ and an explicit degree-6 numerator in α | `kernel/radical.rs::step7` | radicand fixed by the discriminant ramification + an F₇ split-locus solve; update map by rational reconstruction over four primes; cross-checked three ways — the degree-7 modular correspondence R(t,t′) (CRT over four primes), the Vélu 7-isogeny orbit, and (independently) Φ₇ and the CSIDH-512 𝔩₇ action. Extension cyclic ℤ/7 (every R(t₀,X) splits as 7 or 1⁷) |
| Vélu / Kohel in characteristic 2 | Vélu 1971; Kohel 1996 | implemented V3 (ordinary binary curves y² + xy = x³ + a₂x² + a₄x + a₆) | `binary.rs` | codomain from t = Σx_Q; Kohel x-map; Φ_ℓ mod 2 neighbours = kernel-polynomial neighbours |
| Vélu over F_{p⁴} | — | implemented V3 for the Deuring correspondence (u64 p) | `ext.rs`, `quat/deuring.rs` | j-invariants projected back to F_{p²} and checked |
| Radical isogenies of other degrees (N = 2, 4, 9, 11, 13, …; Onuki–Moriya's CSIDH variants) | CDV 2020; Onuki–Moriya 2022; Castryck–Decru–Houben–Vercauteren 2022 | not implemented – implementation gap (N = 3, 5, 7 done above; N = 11 needs genus-1 X₁(11), N = 13 genus-2 X₁(13), so the single-parameter method does not apply directly) | — | — |
| Montgomery 2- and 4-isogenies (SIKE formulas, projective (A24+ : C24)), 2^e chains with optimal strategies | Costello–Longa–Naehrig 2016; SIKE specification | implemented V3 (over any field, incl. F_{p²} over 434-bit p) | `kernel/two_power.rs`, `fp2.rs` | 4-chain (two strategies) = 2-chain = Weierstrass Vélu chain (j), e = 20; at 434 bits the bench checks the same three |
| Twisted Hessian ℓ-isogenies (product of translates, explicit form; codomain a′ = aˡ, d′ = (ℓd − 6a Σ s_Q)/∏ s_Q derived here) | Bernstein–Kohel–Lange 2015; Dang–Moody 2019 | implemented V3 (odd ℓ prime to 3) | `kernel/hessian.rs` | ℓ = 5, 7, 11, 13, 17 at p = 1009: kernel → identity, images on the codomain, homomorphism, equal point counts; j formula vs a Weierstrass curve |
| General Weierstrass Vélu (any characteristic, incl. char 2, 3 forms) | Vélu 1971 | implemented V3 (odd cyclic kernel) | `weier.rs::gw_velu_cyclic` | y-map by the normalised-differential identity 2Y + A₁X + A₃ = (2y + a₁x + a₃) X′(x) |
| Jacobi-quartic models | various | not implemented – implementation gap | — | — |
| Vélu in characteristic 3 (general Weierstrass, over GF(3ⁿ)) | Vélu 1971 | implemented V3 | `weier.rs`, `gf3n.rs` | = short-form Vélu in large char (j and x-map); over GF(3ⁿ): kernel → O, images on the codomain, homomorphism, equal point counts (ℓ = 5..13) |

### P2 — (E, Ẽ) → isogeny: the BMSS family (Table 1 of the paper ✔)

| Algorithm | Paper complexity | needs σ | Status | Where |
|---|---|---|---|---|
| Linear algebra (undetermined coefficients, §6.1) | O(ℓ^ω) | no | implemented V2 (dense solve, O(ℓ³)) | `find/bmss.rs` |
| Stark 1972 (continued fraction of ℘~ in ℘, §6.2) | O(ℓ M(ℓ)) | no | implemented V2 | same |
| Atkin 1992 (D(℘)=z^{2−2ℓ}exp F, §6.4) | O(ℓ M(ℓ)) | yes | implemented V2 | same |
| Atkin + modular composition (§6.4) | O(M(ℓ)√ℓ + ℓ^{(ω+1)/2}) | yes | implemented V2 with *naive* composition | same |
| Elkies 1992 (µ_{k,j} triangular system, §6.3) | O(ℓ²) | yes | implemented V2 | same |
| Elkies 1998 (recurrence (16), power sums (17), §4.2) | O(ℓ²) | yes | implemented V2 | same |
| fastElkies (ODE for S(x), §4.3) | O(M(ℓ)) | yes | implemented V2; V3 solves the ODE coefficient by coefficient in O(ℓ²) (was O(ℓ³)) with Karatsuba series | same |
| fastElkies′ (rational reconstruction, §4.3) | O(M(ℓ) log ℓ) | no | implemented V2 (dense Padé) | same |
| Padé on the inverted ℘-series (my own formulation) | — | no | implemented V1 | `find/elkies.rs::bmss_isogeny` |
| σ = Σ x(Q) from Φ's second derivatives (Elkies' E2 relation) | — | — | implemented V2; derived here from E2−ℓE2′ and checked against the true σ for ℓ ≤ 23 | `find/bmss.rs::sigma_from_phi` |

Every method returns the kernel polynomial g (D = g²), which is then fed to Kohel's formula, and the
codomain must equal Ẽ. All eight reproduce Kohel's kernel and numerator exactly for ℓ up to 101,
including Galois-stable kernels whose points are not rational; they also agree at 256 and 511 bits
(V3). The arithmetic is Karatsuba, not FFT, so the M(ℓ) bounds are not reached.

Not implemented – implementation gap: Couveignes' "ℓ-isogenies using the p-torsion" (ANTS-II 1996 ✔),
Lercier's characteristic-2 algorithm, Lercier–Sirvent and Lercier–Lubicz–Vercauteren (arXiv 2003.06367).
These solve (E, Ẽ) → isogeny in small characteristic, where the BMSS family fails (it divides by
integers up to ≈ 2ℓ). The binary-curve support of V3 (`gf2n.rs`, `binary.rs`) is the arithmetic they
would need; in characteristic 2 the crate currently finds isogenies only from kernels or by factoring
division polynomials, not from (E, Ẽ).

### P3 — curve + ℓ → ℓ-isogenies

| Algorithm | Status | Where | Notes |
|---|---|---|---|
| Φ_ℓ(j,Y) roots + Elkies normalised codomain + P2 algorithm | implemented V1 (Padé) / V2 (any of 8) | `find/elkies.rs`, `find/bmss.rs::isogenies_via_phi` | |
| Φ_ℓ mod p from q-expansions, dense linear algebra | implemented V1 | `find/modpoly.rs::compute_linear_algebra` | practical to ℓ ≈ 31–43 |
| **Φ_ℓ mod p by Hecke operators and Newton's identities**: Φ(j(q), Y) = (Y − j(q^ℓ)) G(Y), power sums of the ℓ conjugates are ℓ U_ℓ(j^m), each coefficient read back as a polynomial in j from its polar part | implemented V3 (valid for char > ℓ + 1) | `find/modpoly.rs::compute_hecke` | equal to the linear-algebra Φ_ℓ for ℓ ≤ 43 and to the CRT integer coefficients mod 1259; 3.5× (ℓ = 11) to 33.4× (ℓ = 43) faster in the recorded run; ℓ = 127 in 6.1 s at 61 bits (1.47 s after the third speed pass, `results/p12-speed.jsonl`) |
| Integer coefficients of Φ_ℓ by CRT; Φ_ℓ over any field from them (also mod 2) | implemented V3 | `find/modpoly.rs::integer_coeffs`, `Phi::via_crt` | Φ mod 2 roots = kernel-oracle neighbours on binary curves |
| Division-polynomial factoring (Schoof / Elkies / Couveignes ℓ-torsion style) | implemented V1; char 2 V3 | `find/divpoly.rs`, `binary.rs::kernel_polys` | 0.46 ms (ℓ = 5) to 2.10 s (ℓ = 23) at 61 bits, V1 run |
| Frobenius eigenvalue on an Elkies kernel, ± separated via y | implemented V3 | `find/sea.rs::elkies_eigenvalue` | computed in F_q[x, y]/(h, y² − f) |
| **Schoof–Elkies–Atkin point counting** (t mod 2, Elkies primes with isogeny cycles, Atkin candidate sets from the factor degree of Φ_ℓ(j, Y), recombination by baby-step giant-step) | implemented V3 | `find/sea.rs::sea`, `sea_opts`, `atkin_degree_with` | equals BSGS orders at 40 and 61 bits; [#E]P = 0 for random P at 127 bits, same order for every cycle bound. The Atkin factor degree reuses Y^q mod Φ_ℓ(j, Y) from the root test and gets Y^(q^k) by Q-matrix products (Frobenius = composition with Y^q); for squarefree Φ_ℓ(j, Y) without roots all factors have one degree r ∣ ℓ+1 (Frobenius acts freely on P¹(F_ℓ)), so equality at the divisors of ℓ+1 replaces a gcd per k; checked against the gcd definition on 60 Atkin cases |
| Isogeny cycles for Elkies primes: Frobenius eigenvalue mod ℓᵏ on the rational cyclic ℓᵏ-subgroup (pull-back of the next kernel along the non-backtracking chain) | implemented V3 | `find/sea.rs::eigenvalue_cycle` | t mod 3⁴, 5³, 7², 11² = BSGS trace; used by SEA up to kernel degree 40 |
| Atkin primes over F_{p^r} towers (recover t mod l, not just a set) | implemented V3 | `find/sea.rs::atkin_eigenvalue_tower`, `atkin_candidates_tower`, `atkin_trace_tower`, `sea_atkin_tower`; `fpr.rs` (F_{p^r}, runtime r) | over F_{p^d} (d = order of lambda/mu, scanned up to 6) Frobenius^d has an F_l eigenvalue nu = lambda^d; the ordinary Elkies eigenvalue over F_{p^d} returns nu; the refined candidate set contains the true t mod l (checked against the BSGS trace at 39 bits) and is often a singleton — then `sea_atkin_tower` CRTs it like an Elkies congruence. The roots in F_{p^d} come free: F_{p^d} is built as F_p[z]/(f) for an irreducible factor f of Φ_ℓ(j, Y) over F_p (`FpR::from_modulus`), so z is a root |
| Isogeny-graph sampling and integer CRT | implemented V3 | `find/sutherland.rs::phi_mod_p`, `phi_crt` | Original X₁(ℓ)/uniform sampler, rational torsion transport and interpolation; unchanged. Exact agreement with Hecke for ℓ = 3, 5, 7, 11 and integer coefficients for ℓ = 3, 5. Historical performance records apply to this sampler, not the CM path. |
| CM-guided volcano Φ_ℓ (Bröker–Lauter–Sutherland Algorithm 2.1) and explicit CRT | implemented, bounded reference | `find/sutherland/cm.rs::CmPlan`, `find/sutherland/crt.rs::ExplicitCrt`; CLI `modpoly --method cm` | CM orders and norm-equation primes; Hilbert roots with trace/torsion certification; paired surface/floor class actions and projection; interpolation; streaming CRT modulo any positive target. Odd prime ℓ = 3..31 and documented resource limits. See [usage](USAGE.md#cm-implementation-limits-and-api) and [evidence](../../research/cm_volcano_crt_20261016/REPORT.md). Optimized large-degree selection, polycyclic presentations and Enge evaluation remain unimplemented; no performance claim. |
| Hilbert class polynomial / CM method | implemented | `cm.rs::hilbert_class_poly` | Native fixed-point complex evaluation with rounding checks; the CM modular-polynomial path additionally certifies each finite-field seed before use. |

### P4 — two curves → an isogeny

| Algorithm | Ref | Status | Where | Notes |
|---|---|---|---|---|
| Bidirectional BFS in the ℓ-graph | Galbraith 1999 | implemented V1; any neighbour oracle V3 | `path/galbraith.rs`, `path/graph.rs` | oracles: Φ_ℓ over F_p, F_{p²}, GF(2ⁿ), or kernel factoring |
| Random-walk collision (+ volcano normalisation) | Galbraith–Hess–Smart 2002 | implemented V1; any oracle V3 | `path/ghs.rs` | needs the walk primes to connect the two curves |
| Weighted-prime walk, one Φ_ℓ per step | Galbraith–Stolbunov 2013 ✔(abstract) | **partial** V2: the weighted-degree walk only | `path/ghs.rs::galbraith_stolbunov` | their other refinements are not implemented |
| Paths on ordinary binary curves | Galbraith–Hess–Smart 2002 (char 2 was their target) | implemented V3 | `binary.rs` + `path/graph.rs` | verified edge by edge |
| Volcano navigation: height, levels, ascent/descent, crater walk | Kohel 1996 | implemented V1 (partial: same ℓ-volcano only) | `path/volcano.rs` | |
| End(E) of an ordinary curve (conductor via volcano levels) | Kohel 1996 | implemented V2 | `path/endo.rs` | ℓ = 2 dividing the conductor of Z[π] unsupported |
| Class-group action, ordinary, direction by Frobenius eigenvalue; MITM | Couveignes 2006; Rostovtsev–Stolbunov 2006 | implemented V1 | `path/couveignes.rs` | planted exponents |
| CSIDH action (direction by quadratic character): stepwise, CLMPR batched, projective tree strategy | CLMPR 2018; Meyer–Reith 2018 | implemented V2, V3 (generic field, CSIDH-512) | `path/csidh.rs` | commutes, invertible, Φ_ℓ per step; batched = stepwise = tree |
| Group-action inversion by MITM (CSIDH) | — | implemented V2 | `path/csidh.rs::mitm` | classical cost ≈ (2m+1)^{n/2} |
| **CSURF** (CSIDH on the surface: Montgomery⁻ curves y² = x³ + Ax² − x, p ≡ 7 mod 8; radical 2-steps for the ideal above 2, odd steps CSIDH-style with x-only Vélu or √élu) | Castryck–Decru 2019; BDFLS 2020 §6 | implemented V4 (radical 2-step formula derived here; CSURF-512 prime reconstructed, see the module note) | `path/csurf.rs` | every 2-step = Vélu from the kernel point, Φ₂, no backtracking; odd steps = Vélu and √élu, Φ_ℓ (ℓ ≤ 13); actions commute and invert at 62, 111 bits and CSURF-512 |
| Ideal-class order from isogeny cycles | Couveignes 2006; Stolbunov 2010 | implemented V2 | `path/csidh.rs::ideal_order` | order divides h(−4p) |
| **Class-group structure and relation lattice** (binary quadratic forms of disc −4p, h by BSGS around the analytic estimate, discrete logs, LLL basis, Babai reduction of exponent vectors) | Beullens–Kleinjung–Vercauteren 2019; De Feo–Kieffer–Smith 2018 | implemented V3 for 20–50-bit p | `path/relation.rs` | basis vectors act trivially on curves; reduced vectors give the same curve |
| Supersingular over F_p², walk to F_p + F_p-graph search | Delfs–Galbraith 2016 | implemented V1 (j-path only; p ≡ 7 mod 8) | `path/delfs_galbraith.rs` | |
| **Quaternion algebra B_{p,∞}, O₀, ideals (HNF), left/right orders, norms, LLL, Fincke–Pohst** | — | implemented V3 for p ≡ 3 mod 4 | `quat/mod.rs` | |
| **KLPT** (ideal → equivalent ideal of norm 2ᵉ) | Kohel–Lauter–Petit–Tignol 2014 | implemented V3 for left O₀-ideals, ℓ = 2 | `quat/klpt.rs` | N(J) = 2ᵉ, J ⊂ O₀, J = Iξ exactly, primitive; p up to 2¹²⁸ |
| **Deuring correspondence**: ideal → curve (smooth-norm equivalent ideal, kernel over F_{p⁴} by 2D Pohlig–Hellman, Vélu chain) and kernel → ideal | EHLMP 2018; De Feo–Kohel–Leroux–Petit–Wesolowski 2020 | implemented V3 for u64 p (torsion over F_{p⁴} needs N(J) | odd part of p² − 1) | `quat/deuring.rs` | bijection classes → supersingular j checked at p = 1259 and 3499 |
| IdealToIsogeny for 2ᵉ norms at cryptographic size (SQIsign-style, or via higher dimension) | DKLPW 2020; SQIsign 2D | not implemented – implementation gap | — | the KLPT output is not turned into an isogeny at large p |
| **Brandt matrices, class sets of O₀, Mestre's supersingular graph on the j-line** | Mestre 1986; Kohel 1996; Pizer | implemented V3 | `quat/brandt.rs` | Eichler mass formula asserted; tr(Bᵏ) = tr(Aᵏ); B(2) = Φ₂ multiplicity matrix entrywise under the bijection |
| Quantum subexponential (Kuperberg), Childs–Jao–Soukharev; Biasse–Jao–Sankar | — | not implemented – resource limit (quantum model) | — | |

### P5/P6 — auxiliary and higher-dimensional

| Algorithm | Status | Where | Notes |
|---|---|---|---|
| Dual isogeny (from Φ_ℓ; explicit rescaling u = ℓ) | implemented V2 | `find/dual.rs` | dual∘φ = [ℓ] and φ∘dual = [ℓ] on random points |
| Endomorphism ring, supersingular (as an O₀-ideal class) | implemented V3 via the Deuring correspondence (u64 p) | `quat/deuring.rs::ideal_of_kernel` | kernel ↔ ideal round trips |
| **Richelot (2,2)-isogeny** of genus-2 Jacobians: codomain and image of points | implemented V3 | `genus2.rs::richelot` | L-polynomial preserved; images lie on the codomain |
| **Splitting** (Δ = 0: J(C) → E₁ × E₂) and **gluing** (E₁ × E₂ → J(C), Howe–Leprevost–Poonen) | implemented V3 | `genus2.rs::split`, `glue` | L(C) = L(E₁) L(E₂); the published −1 normalisation gives a quadratic twist for these models and is dropped (tested at p ≡ 1 and 3 mod 4) |
| Superspecial Richelot graph over F_{p²}, vertices by Igusa–Clebsch invariants | implemented V3 | `genus2.rs::superspecial_graph` | vertex counts = Ibukiyama–Katsura–Oort for 11 primes 11..83 (tests) and 131, 199 (bench); products h(h+1)/2 |
| **Theta-model (2,2)-isogenies** in level-2 coordinates: codomain from 8-torsion above the kernel (no square roots or inversions), evaluation, gluing E₁ × E₂ → J(C) with the vanishing dual coordinate recovered from x + T′, Rosenhain invariants, split detection with the two j-invariants (formulas re-derived here from the duplication formula) | implemented V3 | `theta.rs` | gluing = Howe–Leprévost–Poonen (Igusa–Clebsch, 12 instances); each theta step is one of the 15 Richelot neighbours (12 steps) |
| **Kani-lemma (2ᵃ, 2ᵃ)-chains**: C × E with kernel {(γP, φP)} for deg φ + deg γ = 2ᵃ, split decision (the core test of the Castryck–Decru / MMPPW / Robert SIDH attacks) | implemented V3: smooth diamonds at u64 p (a ≤ 16), and with γ = u + v·i ∈ End(E₀) (deg 2ᵃ − 3ᵇ = u² + v²) up to 360-bit p (a = 200); optimal strategy with theta doublings | `theta.rs::chain`, `testdata.rs::{kani_instance, kani_endomorphism_instance}` | splits exactly at E₀ × X, X = E₀/(ker φ + ker γ) by Vélu (a = 8, 10, 12, 16); with γ ∈ End(E₀) a factor has j = 1728 (126 to 360 bits); twisted isotropic kernels do not split |
| **Dimension-g theta (2,...,2)-chains** and **dimension-4 Kani embeddings**: the auxiliary isogeny is α ∈ M₂(ℤ[i]) of degree 2ⁿ − 3ᵇ via a four-squares decomposition (so the auxiliary degree can be any positive integer, the point of going to higher dimension), kernel in E₀² × E²; split decision by the vanishing even theta constants | implemented V3 (general g; exercised at g = 2 = the dedicated code, and g = 4) | `theta_g.rs`, `testdata.rs::kani4_instance` | g = 2 reproduces the dimension-2 split (a = 8..16); g = 4 embedding of a 3ᵇ-isogeny splits, twisted kernel does not (n = 12, 16, 20) |
| Reading the secret isogeny *off* the embedding (Robert's evaluation of F on torsion → full key recovery), the Castryck–Decru digit-guessing recovery for 2ᵃ < 3ᵇ, dimension-8 embeddings with a concrete instance, (ℓ,ℓ)-isogenies for odd ℓ, SQIsign2D | not implemented – implementation gap (the split *decision* is implemented up to dimension 4; recovery needs factor-structure extraction from a generic split codomain, and the 2ᵃ < 3ᵇ recovery needs the attack's intermediate-curve bookkeeping) | — | |

### P7 — two curves and a degree d → an isogeny of degree d (V4)

Benčina–Kutas–Merz–Petit–Stopar–Weitkämper, fixed-degree isogenies between supersingular curves
(Isogeny Club 2023 slides; d ≈ p^{1/2+ε}). Module `fixed_degree/`. The full pipeline runs for
u64 p < 2³¹ (Deuring over F_{p⁴}); steps 2–5 are pure quaternion arithmetic and run at any size.

| Step / algorithm | Ref | Status | Where | Checked by |
|---|---|---|---|---|
| 1. End(E) for any supersingular E: bidirectional BFS in the 2-isogeny j-graph from 1728 (Õ(√p)), then the path translated to a left O₀-ideal (sub-ideal whose Deuring image is the next j; reduced after each step) | Eisenträger–Hallgren–Lauter–Morrison–Petit 2018 (classical form) | implemented V4 (quantum Grover variant: resource limit) | `fixed_degree/endring.rs` | Deuring image of the result = j(E); the ℓ + 1 sub-ideals' images = the Φ_ℓ roots (ℓ = 2, 3); right order maximal |
| 2. Connecting ideal I = N O₁O₂ | Kirschmer–Voight 2010 | implemented V4 | `fixed_degree/connecting.rs` | O_L(I) = O₁, O_R(I) = O₂; = class of conj(I₁)I₂; Deuring image of I₁I = E₂ |
| 3. Normalised norm form on an LLL-reduced basis, Q = 2 Nrd / N(I) | slides | implemented V4 | `fixed_degree/connecting.rs::norm_form` | Q(x) N(I) = 2 Nrd(Σ xᵢσᵢ) on random x |
| 4a. Cornacchia with guessing (x₃, x₄): complete the square to X² + |D|Y² = N, factor N, all √−D mod N by CRT, abandon when ω(N) > 11 log log N | slides; Cornacchia 1908; Cornacchia–Smith | implemented V4 (factoring by Pollard rho only: no NFS, so at 127 bits the N do not factor) | `fixed_degree/{cornacchia,factor}.rs` | brute force on boxes; planted solutions; in the pipeline |
| 4b. Coron's bivariate integer Coppersmith (box lattice; certified vectors ‖h(xX)‖₁ < n; Res_y) | Coron 2007 | implemented V4 (+ a simplex lattice variant) | `fixed_degree/coppersmith.rs`, `lattice.rs` (exact integer LLL, any dimension), `mpoly.rs` (multivariate integer polynomials, Bareiss resultants, exact integer roots) | planted roots; reaches Coron's bound at k = 2; steps 2–5 at 61 and 127 bits |
| 4c. Trivariate integer Coppersmith | Bauer–Joux 2007 | **partial** V4: Coron-style trivariate lattice with a run-time algebraic-independence check (resultants mod primes); the rigorous Bauer–Joux Gröbner step is not implemented; reach |xᵢ| ≈ 4 at 30-bit and ≈ 12 at 40-bit Gram entries | `fixed_degree/coppersmith.rs` | planted roots; steps 2–5 at 61 bits (small ε) |
| 5. J = I conj(α)/N(I) of norm d | slides (as in KLPT) | implemented V4 | `fixed_degree/pipeline.rs::solve_degree` | N(J) = d, O_L(J) = O₁, O_R(J) = β⁻¹O₂β exactly |
| 6a. Ideal → isogeny, smooth d: kernel E₀[J + dO₀] pushed through E₀ → E₁′, Vélu chain | slides | implemented V4 | `fixed_degree/represent.rs` | Vélu chain from the caller's E₁ lands on j(E₂) |
| 6b. Ideal → isogeny, any odd d: torsion images χ(S₁), χ(S₂) on E₁[2^(a+2)] via φ_K = [N(K′)⁻¹]φ_K′∘θ (Robert-style representation); **Kani certificate** with an auxiliary isogeny of degree 2ᵃ − d from a second random ideal | Robert 2022; SQIsign2D-style auxiliaries | implemented V4 for d < 2ᵃ (certified) and any d (pairing check only) | `fixed_degree/represent.rs` | (2ᵃ,2ᵃ) split + twisted control does not split; e(χS₁, χS₂) = e(S₁,S₂)^d; Deuring j |
| └ certificate with an endomorphism γ ∈ End(E₁) of degree 2ᵃ − d, or Robert's dimension-4 integer matrix [[u,−w],[w,u]] | Robert 2022 | not implemented – the first does not exist when 2ᵃ − d < √p; with the second the first (2,2,2,2) step is a product of two gluings, which `theta_g::chain_g` cannot start from (checked: integer α fails, Gaussian α works) | — | |
| 7. Hybrid: guess an ℓ^{e₁} prefix on the quaternion side (sub-ideals = End(E) from End(E₁)), solve the ℓ^{e−e₁} rest, compose | slides | implemented V4 | `fixed_degree/hybrid.rs` | Deuring image of the composite = E₂; Kani-certified isogeny (3⁷) |
| Baselines: exhaustive prime-degree search (Φ_d roots or kernels); MITM for d = ℓᵉ; van Oorschot–Wiener claw search | slides; Galbraith et al. 2016 | implemented V4 | `fixed_degree/baselines.rs` | planted paths found; cyclicity checked; counts vs theory |
| Grover O(√d), Tani d^{1/3} | Grover 1996; Tani 2009 | not implemented – resource limit (query-count model only) | `baselines.rs::grover_tani_cost_model` | |

### Field arithmetic underneath (V3)

| Component | Where | Measured (micro bench, this VM) |
|---|---|---|
| Montgomery F_p, N × 64-bit limbs (CIOS, no-carry variant when the top limb has a free bit) | `fpm.rs` | P-256 mul 24.4 ns; CSIDH-512 mul 85.4 ns (portable) |
| MULX/ADCX/ADOX assembly CIOS (generated by `adx_gen.rs`, printed by `cargo run --bin gen_mont_adx`; first written with a Python script, which this native generator reproduces byte for byte) | `fpm.rs::adx` | CSIDH-512 mul 78 vs 86 ns; N = 4 not faster than the portable code, so opt-in (`ISOGENY_ADX4`) |
| Inversion by Pornin's optimised binary GCD | `fpm.rs::bingcd_inv` | P-256 13.1 → 3.2 µs; CSIDH-512 75.5 → 8.0 µs (before/after from the commit message; after = `results/micro-p1-step4.jsonl`) |
| sqrt by one exponentiation when q ≡ 3 mod 4; sliding-window exponentiation | `field.rs` | 512-bit sqrt 370 → 57 µs |
| Karatsuba (threshold measured per field), lazy wide reduction for convolutions, monic division without inversion | `poly.rs`, `fpm.rs` | 256-term poly product over F_{P-256}: 1.05 ms (threshold 4) vs 0.79 ms (16) vs 1.91 ms (schoolbook) |
| GF(2ⁿ), n ≤ 63, PCLMULQDQ | `gf2n.rs` | mul 4.9–5.8 ns; sqrt as a linear map 6.2–8.4 ns; z² + z = c 1.7–5.3 ns (n = 23, 41, 61) |
| GF(3ⁿ), n ≤ 40: packed base-3 digits outside, bitsliced inside (digit 1 / digit 2 bit masks, table conversion 8 digits at a time); shift-and-add multiplication, Frobenius cubing, Itoh–Tsujii inversion; low-weight irreducible (Rabin), searched by exact-weight tails | `gf3n.rs` | mul 98 / 435 → 18 / 42 ns, inv 1.36 / 24.1 → 0.10 / 0.62 µs (n = 5 / 20); n = 40 mul 100 ns, inv 2.2 µs (did not construct before); `results/p10-speed.jsonl` |
| F_{pʳ}, runtime r ≤ 16: stack accumulators with lazy u128 reduction, sparse reduction tail, Frobenius matrix, norm-based inversion | `fpr.rs` | 40-bit p, r = 2 / 3 / 6: mul 54 / 77 / 185 → 37 / 39 / 72 ns, inv 677 / 936 / 1834 → 204 / 250 / 794 ns |
| Powering modulo a polynomial of degree 2..64: table of xᵏ mod m, each reduced coefficient one lazily accumulated dot product, a degree-≤ 1 base multiplied in by a shift; lazy reduction interval sized to p (u64 accumulation when the sum fits) | `poly.rs::TableModulus`, `field.rs::lazy_conv` | x^p mod a degree-24 polynomial (one run, same session): 42 → 23 µs (18-bit p), 333 → 107 µs (61-bit p); in context: division-polynomial factoring median 1.33×, SEA median 1.13–1.39× |
| Φ_ℓ by Hecke/Newton reading only the coefficients used: baby powers S^a and giant powers S^(g√ℓ) of S = q j, each needed [S^m]_N one lazily accumulated dot product; S = (E4 P⁸)³ with P = 1/∏(1 − qⁿ) from Euler's pentagonal recurrence | `find/modpoly.rs::compute_hecke`, `s_series` | 1.7–4.3× for ℓ = 11..127 (`results/p12-speed.jsonl`); Φ_ℓ for ℓ ≤ 89 over the 127-bit field 31.2 → 9.0 s; identical output (= linear-algebra Φ_ℓ, = CRT integer Φ_ℓ) |
| Modular composition h ↦ h(ξ) mod g, Brent–Kung baby-step giant-step with the number of baby steps chosen per caller (√(kn) in EDF, √n then the full Q-matrix in DDF, Q-matrix in the SEA Atkin degree); Frobenius powers x^(q^k) by composition with x^q instead of exponentiation | `poly.rs::Composer`, `ddf`, `edf`; `find/sea.rs::atkin_degree_with` | 127-bit SEA median 4.9× (`results/p11-speed.jsonl`); division-polynomial factoring median 2.1×, EDF for ℓ = 23 at 61 bits 253 → 55–64 ms |
| Signed big integers, Miller–Rabin in Montgomery form, Pollard–Brent | `int.rs`, `bigint.rs`, `field.rs` | — |
| F_{p²} = F_p[i]/(i² + 1) over any prime field (p ≡ 3 mod 4), Karatsuba | `fp2.rs` | 434-bit multiplication 259 ns |

## Mapping to the names in the request

* **Couveignes**: hard homogeneous space / class-group action (V1 ordinary; V2 CSIDH-style + MITM + cycle
  orders; V3 CSIDH-512, tree strategy, relation lattice, radical steps).
  Isogeny cycles (Couveignes–Morain) for Elkies primes are implemented (V3), and Atkin primes are
  now resolved to t mod l over F_{p^r} towers (V3); not covered: the 1996 p-torsion method (small
  characteristic).
* **Kohel**: kernel-polynomial formulas (V1 odd, V2 even, V3 char 2), volcano navigation (V1), ordinary
  End(E) (V2), and the supersingular half of the thesis: quaternion orders, Brandt matrices and the
  Deuring correspondence (V3); KLPT (Kohel–Lauter–Petit–Tignol, V3).
* **Galbraith**: BFS (V1), GHS (V1; char 2 in V3), Galbraith–Stolbunov (V2, weighted walk only),
  Delfs–Galbraith (V1).
* **Kernel-only methods**: Vélu general, Kohel even, x-only Vélu, Montgomery x-only Vélu, ℓⁿ chains
  (V2); √élu (Weierstrass and Montgomery), Edwards, Huff, Hessian, char-2 Vélu, radical isogenies,
  Montgomery 2-/4-isogeny chains, theta (2,2)-isogenies and gluing (V3).

## Duplication note

The main crate already contains related code, which this crate neither imports nor replaces:
`src/pqc/csidh.rs` (toy CSIDH, p = 419), `src/pqc/sqisign.rs` (toy), `src/pqc/fast/isogeny.rs`,
`src/isogeny/*` (small-prime Vélu, CM, volcano), `src/cryptanalysis/isogeny_walk/kernel.rs`
(Elkies kernel polynomial from Φ_ℓ with its own derivation of the E2 relation; **not**
cross-tested against `find/bmss.rs::sigma_from_phi`), `src/cryptanalysis/binary_velu.rs`
(characteristic-2 Vélu; not cross-tested against `binary.rs`).
