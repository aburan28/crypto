# Isogeny-computation algorithms: what exists, what this crate implements, what is missing

Scope: algorithms that *compute* an isogeny (or the data needed to), in the sense of the request
"all the approaches to generating isogenies between curves (Couveignes, Kohel, Galbraith, ...)".
This is a working survey for planning, not a literature review.

**Provenance of the entries.** Marked ✔ in the "checked" column: the statement (algorithm, input,
output, complexity) was read in this session from the source — the Bostan–Morain–Salvy–Schost paper
(arXiv cs/0609020, full text) for everything in §P2, and the abstracts/bibliographic records
retrieved by web search for Couveignes (ANTS-II 1996, LNCS 1122, pp. 59–65), Couveignes–Morain
(ANTS-I 1994), Galbraith–Stolbunov (arXiv 1105.6331) and Costello–Hisil (ASIACRYPT 2017). All other
references are cited from memory; their details have **not** been re-verified.

**Status vocabulary** (following the repository's reporting rules): *implemented* = code + a test that
checks the output independently; *partial* = implemented with a stated restriction; *not implemented
– implementation gap* = nothing prevents it but effort; *not implemented – resource limit* = needs
hardware or a model of computation that is not available here (e.g. a quantum computer). No entry
is withheld for any other reason.

## The problems, and which algorithms solve which

| | Problem | Typical algorithms |
|---|---|---|
| P1 | kernel → isogeny (codomain and map) | Vélu, Kohel, √élu, Montgomery/Edwards x-only Vélu, chains with strategies, radical isogenies |
| P2 | two ℓ-isogenous curves (and maybe σ) → isogeny | the BMSS family: Stark, Elkies 1992/1998, Atkin, fastElkies, … |
| P3 | curve + ℓ → its ℓ-isogenies (kernels / neighbours) | Φ_ℓ roots + Elkies, division-polynomial factoring, isogeny cycles (Atkin primes) |
| P4 | two curves → *some* isogeny between them | Galbraith, GHS, Galbraith–Stolbunov, Kohel (volcanoes), Couveignes (hard homogeneous space), Delfs–Galbraith, quantum, KLPT/SQIsign (supersingular, endomorphism-ring based) |
| P5 | auxiliary computations | dual isogeny, End(E) (ordinary: Kohel; supersingular: quaternions), class-group structure, modular polynomials |
| P6 | higher-dimensional | Richelot, theta (2,2)/(ℓ,ℓ), Kani-lemma methods |

## Status table

`V1` = first delivery (PR #1526 first two commits), `V2` = this extension. Module paths are under `src/`.

### P1 — kernel → isogeny

| Algorithm | Ref | Status | Where | Checked by |
|---|---|---|---|---|
| Vélu, odd cyclic kernel | Vélu 1971 | implemented V1 | `kernel/velu.rs` | agrees with Kohel/√élu; homomorphism test |
| Vélu, **any** finite subgroup (even degree, 2-torsion, non-cyclic) | Vélu 1971 | implemented V2 | `kernel/velu.rs::velu_general` | E/E[2] has the j of E; Φ₂; homomorphism; cyclic n = 4…15 vs Kohel |
| Kohel formulas, odd | Kohel 1996 | implemented V1 | `kernel/kohel.rs` | = Vélu |
| Kohel formulas, **even / 2-torsion factor** | Kohel 1996 | implemented V2 | `kernel/kohel.rs` | = Vélu for degrees 2, 4, 6, 8, 9, 10, 12, 15; kernel E[2] |
| √élu | Bernstein–De Feo–Leroux–Smith 2020 | **partial**: structure (I/J/K index sets, series-ring products) on short Weierstrass; the resultant is a naive O(ℓ) Horner, so no asymptotic saving yet | `kernel/sqrt_velu.rs` | = Vélu, ℓ ≤ 101 |
| Vélu from **x(P) only** (no y; kernels over extension fields) | folklore; division-polynomial values | implemented V2 | `kernel/xonly.rs` | = point-based Vélu; Galois-stable kernel with non-rational points over F_p² = Kohel over F_p |
| Montgomery x-only Vélu (odd ℓ) | Costello–Hisil 2017 ✔(abstract); Renes 2018; Meyer–Reith 2018 | implemented V2 | `kernel/montgomery.rs` | j-match with Weierstrass Vélu along a 12-step chain; the x-map keeps points rational (this is what fixed the pairing of (A±2) with ∏(X±Z) and the sign of A′) |
| ℓⁿ kernel as a chain; naive / balanced / cost-model-optimal strategies; composite smooth kernels | De Feo–Jao–Plût 2014; Costello–Longa–Naehrig 2016 | implemented V2 | `kernel/chain.rs` | SIDH key exchange over F_{p²} (shared j agrees) for all four strategies |
| Montgomery 2-/4-isogeny special formulas; twisted Edwards (Moody–Shumow 2016); Hessian, Jacobi-quartic, Huff models; Montgomery-native √élu | various | not implemented – implementation gap | — | — |
| Radical isogenies | Castryck–Decru–Vercauteren 2020; Onuki–Moriya | not implemented – implementation gap | — | — |
| Vélu over F_{p^k}, k > 2 (arbitrary extension fields) | — | not implemented – implementation gap (the `Field` trait is generic; only F_p and F_{p²} are provided) | — | — |
| Vélu in characteristic 2 / 3 | Vélu 1971 | not implemented here – implementation gap. Characteristic-2 Vélu already exists in the main crate: `src/cryptanalysis/binary_velu.rs` | — | — |

### P2 — (E, Ẽ) → isogeny: the BMSS family (Table 1 of the paper ✔)

| Algorithm | Paper complexity | needs σ | Status | Where |
|---|---|---|---|---|
| Linear algebra (undetermined coefficients, §6.1) | O(ℓ^ω) | no | implemented V2 (dense solve, O(ℓ³)) | `find/bmss.rs` |
| Stark 1972 (continued fraction of ℘~ in ℘, §6.2) | O(ℓ M(ℓ)) | no | implemented V2 | same |
| Atkin 1992 (D(℘)=z^{2−2ℓ}exp F, §6.4) | O(ℓ M(ℓ)) | yes | implemented V2 | same |
| Atkin + modular composition (§6.4) | O(M(ℓ)√ℓ + ℓ^{(ω+1)/2}) | yes | implemented V2 with *naive* composition | same |
| Elkies 1992 (µ_{k,j} triangular system, §6.3) | O(ℓ²) | yes | implemented V2 | same |
| Elkies 1998 (recurrence (16), power sums (17), §4.2) | O(ℓ²) | yes | implemented V2 | same |
| fastElkies (ODE for S(x), §4.3) | O(M(ℓ)) | yes | implemented V2 with schoolbook series | same |
| fastElkies′ (rational reconstruction, §4.3) | O(M(ℓ) log ℓ) | no | implemented V2 (dense Padé) | same |
| Padé on the inverted ℘-series (my own formulation) | — | no | implemented V1 | `find/elkies.rs::bmss_isogeny` |
| σ = Σ x(Q) from Φ's second derivatives (Elkies' E2 relation) | — | — | implemented V2; **derived here** from E2−ℓE2′ and checked against the true σ for ℓ ≤ 23 | `find/bmss.rs::sigma_from_phi` |

Every method returns the kernel polynomial g (D = g²), which is then fed to Kohel's formula, and the
codomain must equal Ẽ. All eight reproduce Kohel's kernel and numerator exactly for ℓ up to 101,
including Galois-stable kernels whose points are not rational. The "paper complexity" column is
what the paper proves for the fast arithmetic; *this code uses schoolbook arithmetic and does not
achieve it*.
Not implemented: Lercier / Couveignes ("ℓ-isogenies using the p-torsion", ANTS-II 1996 ✔) /
Lercier–Sirvent and the characteristic-2 algorithms (Lercier–Lubicz–Vercauteren, arXiv 2003.06367):
these are for small characteristic and need F_{2ⁿ}-type arithmetic and p-adic lifting – implementation
gap, and relevant if binary Koblitz curves are the target.

### P3 — curve + ℓ → ℓ-isogenies

| Algorithm | Status | Where | Notes |
|---|---|---|---|
| Φ_ℓ(j,Y) roots + Elkies normalised codomain + P2 algorithm | implemented V1 (Padé) / V2 (any of 8) | `find/elkies.rs`, `find/bmss.rs::isogenies_via_phi` | needs Φ_ℓ mod p |
| Φ_ℓ mod p from q-expansions (linear algebra) | implemented V1; **practical only to ℓ ≈ 23–31** | `find/modpoly.rs` | checked against the known Φ₂ |
| Division-polynomial factoring (Schoof / Elkies / Couveignes ℓ-torsion style) | implemented V1 | `find/divpoly.rs` | cost grows quickly with ℓ: 0.46 ms at ℓ = 5 and 2.10 s at ℓ = 23 (61-bit p, baseline run) |
| Frobenius eigenvalue on the kernel (± λ) | implemented V1 (support routine) | `path/couveignes.rs::eigen_class` | x-only, so λ and −λ are not separated |
| Isogeny cycles for Atkin primes (Couveignes–Morain 1994 ✔, Couveignes 1996 ✔) | not implemented – implementation gap (needs F_{p^r}, r | ℓ+1) | — | used for t mod ℓᵏ in SEA |
| Complete SEA point counting | not implemented – implementation gap (the Elkies step ingredients exist; the Atkin/Schoof steps and the ± λ resolution via y do not) | — | `order()` is BSGS, not SEA |
| Better Φ_ℓ algorithms (CRT/Sutherland, Bröker–Lauter–Sutherland via volcanoes, Enge quasi-linear) | not implemented – implementation gap | — | needed for ℓ beyond ≈ 30 |
| Hilbert class polynomial / CM method | not implemented here; exists in the main crate (`src/cryptanalysis/hilbert_class_poly.rs`) | — | — |

### P4 — two curves → an isogeny

| Algorithm | Ref | Status | Where | Notes |
|---|---|---|---|---|
| Bidirectional BFS in the ℓ-graph | Galbraith 1999 | implemented V1 | `path/galbraith.rs` | j-line over a prime set |
| Random-walk collision (+ volcano normalisation) | Galbraith–Hess–Smart 2002 | implemented V1 | `path/ghs.rs` | restriction: needs the walk primes to connect the two curves |
| Weighted-prime walk, one Φ_ℓ per step | Galbraith–Stolbunov 2013 ✔(abstract) | **partial** V2: the weighted-degree walk only | `path/ghs.rs::galbraith_stolbunov` | their other refinements are not implemented |
| Volcano navigation: height, levels, ascent/descent, crater walk | Kohel 1996 | implemented V1 (partial: same ℓ-volcano only) | `path/volcano.rs` | |
| **End(E) of an ordinary curve** (conductor via volcano levels) | Kohel 1996 | implemented V2 | `path/endo.rs` | needs Φ_ℓ for every ℓ dividing the conductor of Z[π]; ℓ = 2 unsupported |
| Class-group action, ordinary, direction by Frobenius eigenvalue; MITM | Couveignes 2006; Rostovtsev–Stolbunov 2006 | implemented V1 | `path/couveignes.rs` | planted exponents |
| **CSIDH-style action** on Montgomery curves (direction by quadratic character) | CLMPR 2018; Couveignes | implemented V2 | `path/csidh.rs` | commutes, invertible, every step satisfies Φ_ℓ |
| Group-action inversion by MITM (CSIDH) | — | implemented V2 | `path/csidh.rs::mitm` | classical cost ≈ (2m+1)^{n/2} |
| Ideal-class order from isogeny cycles; independent class-number count | Couveignes 2006; Stolbunov 2010 | implemented V2 (single ideal only) | `path/csidh.rs` | order divides h(−4p) |
| Class-group *structure*, relation lattice | De Feo–Kieffer–Smith 2018; Beullens–Kleinjung–Vercauteren 2019 | not implemented – implementation gap | — | |
| Supersingular over F_p², walk to F_p + F_p-graph search | Delfs–Galbraith 2016 | implemented V1 (j-path only; p ≡ 7 mod 8 here) | `path/delfs_galbraith.rs` | |
| Quantum subexponential (Kuperberg), Childs–Jao–Soukharev; Biasse–Jao–Sankar | — | not implemented – resource limit (quantum model) | — | |
| Endomorphism-ring-based supersingular path finding: quaternion orders, KLPT, ideal→isogeny (Deuring), SQIsign | Kohel–Lauter–Petit–Tignol 2014; EHLMP 2018; De Feo–Kohel–Leroux–Petit–Wesolowski 2020 | not implemented – implementation gap (large). The main crate has only a toy `src/pqc/sqisign.rs` with a precomputed graph and 326-bit F_{p²} arithmetic in `src/pqc/fast/isogeny.rs` | — | |
| Mestre / Pizer / Kohel Brandt-matrix (Hecke-module) supersingular graphs | Mestre 1986; Kohel 1996 | not implemented – implementation gap | — | |

### P5/P6 — auxiliary and higher-dimensional

| Algorithm | Status | Notes |
|---|---|---|
| Dual isogeny (from Φ_ℓ; explicit rescaling u = ℓ) | implemented V2 (`find/dual.rs`) | dual∘φ = [ℓ] and φ∘dual = [ℓ] checked on random points |
| Endomorphism ring, supersingular | not implemented – implementation gap | |
| Richelot, theta-model (2,2)/(ℓ,ℓ)-isogenies, Kani-lemma methods (SIDH attacks, SQIsign2D) | not implemented – implementation gap | |

## Mapping to the names in the request

* **Couveignes**: hard homogeneous space / class-group action (P4: V1 ordinary, V2 CSIDH-style + MITM + cycle orders).
  Not covered: the 1996 p-torsion method (small characteristic) and isogeny cycles (Atkin primes).
* **Kohel**: kernel-polynomial formulas (V1 odd, V2 even), volcano navigation (V1) and End(E) (V2).
  Not covered: the supersingular half of the thesis (Hecke module structure of quaternions / Brandt matrices) and KLPT.
* **Galbraith**: BFS (V1), GHS (V1), Galbraith–Stolbunov (V2, weighted walk only), Delfs–Galbraith (V1).
* **Kernel-only methods**: Vélu general, Kohel even, x-only Vélu, Montgomery x-only Vélu, ℓⁿ chains (V2).

## Duplication note

The main crate already contains related code, which this crate neither imports nor replaces:
`src/pqc/csidh.rs` (toy CSIDH, p = 419), `src/pqc/sqisign.rs` (toy), `src/pqc/fast/isogeny.rs`,
`src/isogeny/*` (small-prime Vélu, CM, volcano), `src/cryptanalysis/isogeny_walk/kernel.rs`
(Elkies kernel polynomial from Φ_ℓ with its own derivation of the E2 relation; **not**
cross-tested against `find/bmss.rs::sigma_from_phi`), `src/cryptanalysis/binary_velu.rs`.
