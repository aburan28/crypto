# Binary-field, Koblitz, ECC2K-130, Weil descent / GHS, genus-2/hyperelliptic transfer and equivariant-F4 index-calculus leads in aburan28/crypto

Scope: local files in `/home/user/crypto` only, read on 2026-10-14. All citations are repo-relative paths with line numbers, written as `[path:lines](path)`. Nothing was run. "Measured" means the cited file reports a measurement. "Derived/model" means the file reports a calculation. Where two files disagree, both are cited.

Directories in the brief that turned out to be out of scope or infrastructure only:
- `cd_attack/` is a Castryck–Decru SIDH reconstruction. It provides genus-2 Cantor arithmetic, Richelot (2,2)-isogenies and Kani gluing, but contains no binary-field IC ([cd_attack/SUMMARY.md:1-47](cd_attack/SUMMARY.md)).
- `eigen_seta/` is an isogeny KEM ([eigen_seta/README.md:1-12](eigen_seta/README.md)).
- `groebner_compare/` and `boolean_representation_audit/` are solver-stage harnesses with explicitly no ECC speedup claim ([groebner_compare/README.md:1-8](groebner_compare/README.md); [boolean_representation_audit/README.md:1-9](boolean_representation_audit/README.md)).

---

## Q1. Inventory: every binary-field and transfer-based IC lead (hypothesis, experiments, numbers, status, why stopped, named next steps)

### Takeaway
The repo built a deep Koblitz/binary IC program in four parts:
1. GGMP Frobenius-invariant bases on small fields.
2. A "compact-orbit" 4-sum pair-table ladder up to n=83.
3. Weil-descent Semaev systems with measured solving degrees.
4. Transfer and cover studies (GHS, hyperelliptic, modular curves, odd-characteristic cubic GHS, JV covers).

Every lead is stopped for one of three reasons:
- It is a generic algorithm, bounded by the Shoup or Corrigan-Gibbs–Kogan (CGK) bound.
- Its measured solving degree grows with the base dimension ℓ.
- 131 is prime and 2 is primitive mod 131, so no Frobenius-stable subspace, intermediate field or small-genus GHS cover exists.

The one symmetry lead with measured upside and no measured ceiling is Frobenius-equivariant F4 (G2). It is rank-only at n ≤ 23, and an internal document mislabels it as closed.

### Cited Findings

**Lead K1: GGMP Frobenius-invariant factor bases (small fields)**
- Hypothesis: a factor base F with π(F)=F collapses the unknowns to orbits. This saves ≈n on relations, ≈n² on linear algebra and m! on the decomposition, against √(2n) for rho ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:23-49](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)).
- Measured savings (|F|/#orbits):
  - n=9: 7.86× on relations and 61.73× on LA, against 4.24× for rho.
  - n=13: 12.96× and 167.99×, against 5.10× for rho ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:87-100](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)).
- Structural kill: Frobenius-stable F₂-subspaces are cyclic codes, so their dimensions are subset sums of the 2-cyclotomic coset sizes. At n=131 the only sizes are 0, 1, 130 and 131, so the construction is "empty exactly at the sizes that matter" ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:244-267](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)).
- Couveignes–Lercier invariant bases were closed twice (2026-10-07 and 2026-10-08). Elliptic periods need an auxiliary curve over F₂ with n | #E'(F₂), and Hasse gives ≤5, so the construction fails for every n ≥ 7 ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:892-906](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)).
- Largest end-to-end subspace-base solve: n=31, `icv1-f2m31-tm90707-c95f16f5`, r=1,439,393, h=1492, m=3, divisor [0,1,6] pruned to 35 orbits, 10.4–10.9 s wall, dominated by the pair table ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:488-525](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)).
- CADO-style logs then descent at the same n=31: 129 probes to full column rank in 15.4 s; each descent takes 1–2 relations ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:555-561](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)).
- Unimplemented items from the GGMP talk ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:907-915](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)):
  - exploit block/homogeneous structure in the Gröbner step;
  - precise complexity estimates across characteristics;
  - multiple invariant factor bases ("the driver uses one at a time").

**Lead K2: decomposition oracles (F4, Buchberger, SAT, Crossbred, symmetrised) on subspace-restricted Semaev systems**
- m=2 is quadratic. Chaining S₃ for m ≥ 3 is cubic, with `unknowns = m·ℓ + (m−2)·n`. At n=63, ℓ=6 that is 633 unknowns against a 64-variable budget ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:134-151](RESEARCH_KOBLITZ_INDEX_CALCULUS.md); [RESEARCH_KOBLITZ_INDEX_CALCULUS.md:308-322](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)).
- Engine comparison at n=9, end-to-end DLP: matrix-F4 21.5 ms against Buchberger 85.5 s ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:174-181](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)).
- SAT on the m=3 cubic system takes 3.44 s against 0.58 ms for search ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:196-208](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)).
- H2 (the 64-variable cap is what binds) was falsified: solve cost saturates at 27 unknowns ([research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md:243-258](research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md)).
- H2 was re-opened for F4 only, after F4 refuted 46 unknowns in 145 s ([research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md:260-265](research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md)).
- H4: F4 beats SAT on refutation by 340× at n=15, m=3. At n=21, m=3 F4 refutes in 50 s while SAT gives no verdict in 31 min ([research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md:294-307](research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md)).
- H5 (no algebraic oracle beats exhaustive search anywhere reachable) still stands. Search costs 0.004–0.02 ms per target; the algebraic oracles cost 0.07–136 ms. The "best candidate rung" to falsify H5, n=63 with ℓ=6, was never reached ([research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md:309-317](research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md)).
- H3 (symbolic S₄ links halve the chaining term, e.g. 190 unknowns down to ≈97 at n=31, ℓ=5, m=7) is listed as "still worth doing". No result was found ([research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md:278-292](research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md); [research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md:415](research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md)).
- Crossbred (Route 1) measured xb/F4 = 0.023 in bit operations and 0.696 in wall time at m=3 ([research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md:72-86](research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md)).
  - The margin shrinks with size (0.040 to 0.696 over one rung).
  - No (D, k) produced a filter.
  - At m=3 exhaustive enumeration beats every algebraic oracle by 1800×.
  - "One more rung decides the route."
- The symmetrised oracle (Route 3) was ~350× faster than the chained x-system at m=3 (10,291 ms to 29 ms). It was closed at its gate because at n=131 the only Frobenius-stable V containing 1 are F₂ and the whole field. Per relation it costs 3.1–14.7× the enumeration it would replace ([research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md:114-126](research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md)).
- Joux–Vitse cover-and-decomposition for ECC2K-130 was closed permanently. It needs a composite-degree tower, and embedding into F_{2^{131m}} gives a GHS genus of ≈2^130 ([research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md:193-200](research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md)).

**Lead K3: compact-orbit 4-sum pair-table ladder (the production Koblitz IC)**
- Landed rungs (a=0): n = 41, 53, 61, 71, 73, all in the ECC2K-130 family. n=83 a=1 also landed ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:24-29](RESEARCH_ECC2K130_IC_FEASIBILITY.md); [RESEARCH_ECC2K130_IC_FEASIBILITY.md:162-167](RESEARCH_ECC2K130_IC_FEASIBILITY.md)).

| rung | rank wall (K=600) | mean probes/relation |
|---|---:|---:|
| n=61, r = 2^47.2 | 13.1 s | ~4·10⁴ |
| n=71, r = 2^52.3 | 1,485 s | ~4.6·10⁶ |
| n=73, r = 2^56.3 | 18,398 s | 5.7·10⁷ |

  Source: [RESEARCH_ECC2K130_IC_FEASIBILITY.md:33-41](RESEARCH_ECC2K130_IC_FEASIBILITY.md).
- Against rho:
  - At n=73 the fully charged precompute is 18,400 s against 46 s for rho online: a ~400× total-work deficit traded for a ~185× online win.
  - In operations (one probe counted as one rho step), the online ratio on a mean target is 23× at n=61, 3.5× at n=71, 0.5× at n=73 and 3.7× at n=83.
  - With precompute charged, IC spends 26–1,200× more operations than one rho solve on every rung ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:113-126](RESEARCH_ECC2K130_IC_FEASIBILITY.md)).
- Parallel target extraction is a negative result: there is no wall-time speedup at n=83 because the first hit sits at ~0.1% of the scan. The online stage is probe-rate bound at ~1.1 M probes/s ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:182-195](RESEARCH_ECC2K130_IC_FEASIBILITY.md)).
- Generic-algorithm kill: the 4-sum oracle is table lookups plus S₃ root solves, so `S·T² ≳ r/|Aut|` (CGK) applies. "No engineering of 4-, 5- or 6-sum tables can beat √(r/2n)." The measured rungs sit 2^25–2^27 above the generic line ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:21-31](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
  - At n=73, (S, T) = (2.6·10⁷, 5.7·10⁷) gives S·T² = 2^76.2 against N = 2^49.1 ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:127-141](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
- Leads ladder L1–L8, measured 2026-10-11 ([docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md:171-363](docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md)):
  - L2, canonicalisation: 2–7% of a probe; ceiling ≈1.04×; closed.
  - L3, residual matching: 4×, 7× and 53× worse at n = 13, 17, 23; closed.
  - L5, ordered τ-slots: 2.8× trials-per-relation gain at n=23 in table form; algebraic form 1.79× (n=23) and 1.34× (n=17) against a predicted 2.
  - L6, phase-aware quotient summation: not implemented; "Two to three hours of Rust."
  - L7 (G1 trace replay): not implemented.
  - L8, conductor gap: on the floor curves τ is lost and nothing is gained.
- Multiplier closure A1 (Z[τ] elements of norm ≤ 16): at n=41 with m=3, 6 columns against 60 cut the one-target total 9.2×. That moves the ladder from 51× to 5.6× behind matched rho, "the first ladder row within one decimal order of rho at 39 bits" ([docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md:424-445](docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md)).
- Transfer-aware IC was measured flat: τ- and integer-multiplier arms are 0.81–1.26× of the Frobenius-only arm at K ≥ 3 at n = 23, 37, 41. C41 is partial ([docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md](docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md) §4.1).
- Conservation lemma: a group-homomorphism transfer cannot change decomposition counts ([docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md](docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md) §1).

**Lead K4: Weil-descent solving degree (the Petit–Quisquater / first-fall-degree line)**
- m=3 chained S₃, refutation degree by ℓ, from four draws per cell ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:77-91](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)):
  - ℓ=2: 5
  - ℓ=3: 5–6
  - ℓ=4: 6–7
  - ℓ=5: exactly 7 at (10,5); ≥7 at (11,5), (13,5) and (15,5)
- The first fall degree is 3 everywhere. That value is Kosters–Yeo's trace equation, not a weakness of the system ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:18-25](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md); [research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md:187-216](research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md)).
- m=4, by ℓ ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:298-316](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)):
  - ℓ=1: 5
  - ℓ=2: 6
  - ℓ=3: ≥7 (at n=9 and 10)
  - The degree is flat in n and rises about one degree per unit of ℓ. The fit is ℓ+4; the m=3 fit is ⌈ℓ/2⌉+4.
- The single cell that discriminates between "levelled at 7" and the semi-regular reference D_reg is (15,5) at degree 7.
  - Its matrix is 3.1M × 2.8M, its dense block about 110–120 GB, so it needs a ~128 GB host ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:281-289](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)).
  - Estimated 1.7–5.6 days per draw. Degree 7 at (13,5) needs 47–50 GB and was not run ([research/notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md:672-673](research/notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md); [research/notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md:812-841](research/notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md)).
- Linearization reach of the bilinear S₃ core: the gap is `70.19 − b_max` bits, independent of arity ([research/linearization_reach_20260930/README.md:1-50](research/linearization_reach_20260930/README.md)). At n=131 linearization resolves ℓ + b_max ≈ 28–31 bits per trial, against 84–99 needed ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:256-261](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
- Gröbner with one summand fixed: deciding degree ≈ ℓ/2 (median 5, 6, 6, 7, 7 over ℓ = 3…7) and measured total 2^{6.26ℓ}. SAT is 5× worse than enumeration at ℓ=5 and 62× worse at ℓ=8 ([research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md:525-545](research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md)).

**Lead G2: Frobenius-equivariant F4 (Faugère–Svartz character blocks)**
- Hypothesis: the barrel system with a free target selector is C_n-invariant, so the Macaulay matrix splits into a trivial block over F₂ plus character blocks over F_{2^d}, d = ord_n(2). At n=131, d=130, giving 1 block over F_{2^130} plus the trivial one ([RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md:10-32](RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md)).
- The gauge-fixing objection was refuted. The symmetric matrix is only 1.5–2× the gauge-fixed one, and its σ-orbits are 4–6× smaller per dimension ([RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md:96-114](RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md)).
- Block decomposition was verified exact at n = 7 and 11 ([RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md:128-144](RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md)).
- Rank-time ratio against the full F₂ rank ([RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md:153-158](RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md); [RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md:180-198](RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md)):

| n | d | ratio |
|---:|---:|---|
| 11 | 10 | 23–41× |
| 13 | 12 | 39× |
| 17 | 8 | 105× |
| 19 | 18 | 105× |
| 23 | 11 | 128× |

  - Against the fair gauge-fixed baseline the implied gain is "still ≥5–10× at n=17–19".
- Model: gain ≈ n²·d / (ρ³·64·pen(d)), giving ≈300× dense at n=131. This is "a target to measure, not a result" ([RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md:115-126](RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md)).
- Caveats: this is a rank comparison, not full F4; the control is the plain kernel (M4RI is ≈1.6× faster); n=11 is the 64-variable engine limit ([RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md:163-170](RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md)).
- Conflict on G2's status:
  - The lead ladder lists "equivariant F4 over gauge-fixing (L7/G2)" as "not an opportunity, measured" and predicts "no gain over gauge-fixing". It cites the §2 prediction, which §3 of the scope note had already refuted ([docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md:306-328](docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md); [docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md:474-477](docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md)).
  - The step-back note calls G2 "the one symmetry lead with measured upside and no measured ceiling" ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:44-47](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).

**Lead T1: GHS/Hess Weil descent of ECC2K-130 and other binary curves**
- For prime n the magic number is m ∈ {1, 2} or m ≥ ord_n(2). The smallest genus above the degenerate case is 2^129 for n=131 and 2^161 for n=163 ([research/cover_route_prime_binary_20261008/README.md:20-42](research/cover_route_prime_binary_20261008/README.md)). The same table gives:
  - n=113: 2^27
  - n=233: 2^28
  - n=283: 2^93
  - n=571: 2^113
- Corrected 2026-10-13 under MMT's magic number: the GHS genus over F_{2^131} is 1, 2^130 − 1 or 2^130 ([research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md:15-31](research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md)).
- Pricing windows:
  - Index calculus on a curve over F₂ beats rho only for genus between 130 and ~290–300.
  - Genus 130 with Jac ~ A (the 130-dimensional abelian variety carrying the subgroup, defined under T2) would cost 2^37.17, which is 2^23.64 cheaper than rho ([research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md:300-324](research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md)).
  - Covers over F_{2^131} itself exist in every genus, but cost ≥ 2^131; genus 2 is 2^70.19× rho ([research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md:63-66](research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md)).
  - The genus-1 branch is the zero map on ⟨G⟩, because Σλ^i ≡ 0 mod r ([research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md:74-77](research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md)).
- N=130 counterfactual (what the field would cost at one bit less), corrected:
  - Originally genus 128 at 2^36.41, 23.9 bits below rho.
  - Under MMT, m=8 does not occur. The least admissible case is genus 255 at 2^55.60 against rho 2^60.31, i.e. −4.71 bits.
  - That margin is "not robust to the unit conversion the note itself estimates" (~2^7 curve additions per trial) ([research/ghs_counterfactual_n130_20261013/README.md:1-40](research/ghs_counterfactual_n130_20261013/README.md)).

**Lead T2: non-GHS covers or curves carrying the ECC2K-130 subgroup**
- The open question: is the simple 130-dimensional abelian variety A/F₂ with #A(F₂) = r isogenous to a Jacobian by an isogeny anyone can evaluate? If so, IC would cost 2^37.17 ([research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md:326-360](research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md)).
- Evidence on that question:
  - Point counts do not forbid it.
  - Torelli codimension is 8128. This is a heuristic.
  - Constructive CM stops at genus 3.
- Best explicit construction: the modular curve X_H(3²·7²·263²) at genus 508,799,809 ≈ 2^28.92, which is 20.69 bits above the window ([research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md:33-42](research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md)).
- Cyclic degree-131 covers of genus ≤ 1 curves need genus ≥ 1300 (cost 2^148.92, +88.11 bits) ([research/notes/ecc2k130/RESEARCH_ECC2K130_CURVE_CONSTRUCTION.md:35-47](research/notes/ecc2k130/RESEARCH_ECC2K130_CURVE_CONSTRUCTION.md)).
- Toy ladder ([research/notes/ecc2k130/RESEARCH_ECC2K130_CURVE_CONSTRUCTION.md:117-136](research/notes/ecc2k130/RESEARCH_ECC2K130_CURVE_CONSTRUCTION.md)):
  - At n=3, genus-3 curves carrying A exist.
  - At n=5 there is nothing at genus 4 and nothing hyperelliptic at genus 5.
  - Non-hyperelliptic genus 5 and genus ≥ 6 are not enumerated.

**Lead T3: odd-characteristic cubic GHS (E/F_{q³} → genus-3 curves over F_q), 2026-10-12 to 2026-10-14**
- Census: the minimal genus over the isogeny class is set by #E(K) mod 4 ([research/ghs_cubic_cover_census_20261012/README.md:30-47](research/ghs_cubic_cover_census_20261012/README.md)):
  - odd: 13
  - 2 mod 4: 5 (with one exception, 9)
  - 0 mod 4: 3
- At q=13 the pre-registered prediction failed in one row: 2 of 1,100 curves have minimum 9 ([research/ghs_cover_census_n5_20261013/README.md:40-48](research/ghs_cover_census_n5_20261013/README.md)).
- The hyperelliptic genus-3 cover (m=3, r=5) was built explicitly, with the transfer φ executed by Cantor arithmetic ([research/ghs_cubic_hyperelliptic_genus3_20261013/README.md:20-40](research/ghs_cubic_hyperelliptic_genus3_20261013/README.md)):
  - The homomorphism check passes 4800/4800.
  - The kernel has order 2 or 4, so φ is injective on the ℓ-part.
- Plane-quartic IC on toys: 90/90 trials correct at q ∈ {5, 7, 11}. The rho arm records step counts only, so no IC-vs-rho ratio exists ([research/ghs_cubic_quartic_ic_20261012/README.md:9-28](research/ghs_cubic_quartic_ic_20261012/README.md); [problems/plane-quartic-jacobian-class-arithmetic-missing.md:15-33](problems/plane-quartic-jacobian-class-arithmetic-missing.md)).
- GHS → Smith → Diem pipeline: open obstruction, because Smith's trigonal construction excludes characteristics 2 and 3 ([problems/ghs-to-smith-pipeline-characteristic-mismatch.md:11-62](problems/ghs-to-smith-pipeline-characteristic-mismatch.md)).
- Bhati–Singh weak genus-2 curves: their factor base generates a subgroup inside J_C0(F_{q^gcd(2,n)}), so the relations never reach the large prime ([research/weak_subfield_hyperelliptic_decomposition_20261012/README.md:35-55](research/weak_subfield_hyperelliptic_decomposition_20261012/README.md)).

**Lead T4: Joux–Vitse cover-and-decomposition on E(F_{p⁶}) (structured odd-characteristic class; the only measured below-rho transfer)**
- The plain route has S/rho falling as p^{−1.92}: ~2,100× at p=53, 8.1–10.2× at p=1009 (ℓ = 2^57.9). Extrapolated parity is p* ≈ 2,940 ([research/notes/index-calculus/RESEARCH_COVER_DECOMPOSITION_LEDGER.md:113-114](research/notes/index-calculus/RESEARCH_COVER_DECOMPOSITION_LEDGER.md)).
- The sieving variant reaches 0.021× rho at p=1009 ([research/notes/index-calculus/RESEARCH_COVER_DECOMPOSITION_LEDGER.md:980-985](research/notes/index-calculus/RESEARCH_COVER_DECOMPOSITION_LEDGER.md)).
- The step-back note summarises this as "0.003–0.08× rho at p ≥ 1009" and says it "lives entirely on E(F_{p⁶}), which has intermediate fields" ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:102-109](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
- It does not apply to ECC2K-130.

**Lead T5: surfaces and genus-2+ tooling (GF(2^83) Koblitz)**
- `binary_surface_check` classifies the binary quadric pencil, verifies quartic transfers, and rejects the opposite product. It is tooling only: "no timing comparison, new IC candidate, catalog promotion, or destination-solver advantage" ([docs/ecdlp-surfaces/REPORT.md:1-12](docs/ecdlp-surfaces/REPORT.md); [.claude/skills/ecdlp-surfaces/SKILL.md:110-147](.claude/skills/ecdlp-surfaces/SKILL.md)).

### Inferences
- Every lead that touches n=131 ends either in a counting identity or in the generic bound. The product law (defined under Q2) and CGK hold whatever the engineering, so only the solving degree D and the A-is-a-Jacobian question remain truly open.
- The lead ladder's L7 "closed" entry looks like a stale citation. It quotes the §2 prediction that §3 of the scope note refuted two days earlier. A synthesizer should treat G2 as open.

### Gaps
- No file I read records G2 beyond rank-only measurements at n ≤ 23. There is no full F4 pass, no inverse DFT and no n=31 wall-clock against `macaulay_profile`.
- I did not open the per-run JSON for the C41 transfer-aware cell (marked partial), nor the ic_leads RESULTS_ROUND2.md tables directly. Those figures are quoted via the lead-ladder doc.

---

## Q2. Concrete sizes solved or attempted, Gröbner/F4/XL/SAT costs, degrees of regularity, and comparison with Pollard rho

### Takeaway
- Algebraic Weil-descent solving was measured only at n ≤ 15 (m=3) and n ≤ 10 (m=4), with ℓ ≤ 5. Refutation degrees are 5–7, rising with ℓ.
- End-to-end Koblitz IC reached n=31 with subspace bases and n=83 with the compact-orbit table.
- No algebraic oracle ever beat exhaustive search (H5 stands).
- No rung beat rho on total work. The closest in operation terms is the A1 closure at n=41, 5.6× behind matched rho.

### Cited Findings
- Small-field F4 Macaulay sizes at n=9, 12 unknowns: degree 2, 8 rows in 33 µs; degree 4, 584 rows in 4.03 ms ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:159-165](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)).
- Divisor bases with m=2 ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:282-287](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)):

| n | dim | vars | search | F4 | SAT |
|---:|---:|---:|---:|---:|---:|
| 7 | 3 | 6 | 11 µs | 198 µs | 144 µs |
| 15 | 10 | 20 | 64 µs | 17 ms | 32 ms |

- n=23, m=2 pipeline wall times after optimisation ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:458-461](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)):
  - enumerate: 0.495 s
  - Gröbner: 20.6 s
  - SAT: 137 s
  - pair table: 6.2 s
- Refutation versus finding, same 27 unknowns: F4 needs 27 ms to find (n=9, m=3) but 75 ms to refute (n=15, m=3), while SAT needs 25.3 s to refute ([research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md:246-252](research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md); [research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md:270-280](research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md)).
- F4 at n=31 ([docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md:57-58](docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md); [research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md:260-265](research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md)):
  - dim 16, m=2: 637 ms median.
  - m=3: refutation in 145 s with 46 unknowns.
- Measured degrees (cells (n, ℓ); N unknowns; four draws each; first fall degree 3 throughout) ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:77-91](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md); [research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:304-308](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)):
  - m=3: ℓ=2 gives 5 (N 11–18); ℓ=3 gives 5–6; ℓ=4 gives 6–7; ℓ=5 gives 7 at (10,5), N=25, and ≥7 up to (15,5), N=30.
  - m=4: ℓ=1 gives 5 (N 12–16); ℓ=2 gives 6 (N 18–24); ℓ=3 gives ≥7 (N 30–32).
- Semi-regular D_reg reference: within one of the measured degree in all 14 exact m=3 cells and equal in 10. It predicted the 7 at (10,5) before the run ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:110-117](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)).
- Literature bounds:
  - Kousidis–Wiemers prove D_ff ≤ m² − m + 1 unconditionally (7 at m=3, 13 at m=4).
  - The conditional complexity is 2^86.0 at n=131 against rho 2^60.8 ("about 2^25 short"), with crossover n ≈ 1250 (Petit–Quisquater's was ≈2000) ([research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md:71-104](research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md)).
  - The rigorous Caminata–Ceria–Gorla bound gives solving degree ≤ 263 at n=131, which is vacuous ([research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md:106-116](research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md)).
- Crossbred at m=3: ratio 0.023 in bit operations, 0.696 in wall time; Crossbred includes `crossbred_covers_the_f4_oracle_on_the_cubic_m3_system` (27 unknowns, D=3, k=12) ([research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md:52-86](research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md)).
- The landed n=83 a=1 rung shows a median 8.2× online win over three paired runs. That is online only; precompute is not charged ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:166-167](RESEARCH_ECC2K130_IC_FEASIBILITY.md)).
- Against preprocessing rho at equal advice (calculation): IC online needs 3,000–6,000× more operations at every rung. Preprocessing-rho online would be ≈490, ≈2,500, ≈9,600 and ≈2,600 steps at n = 61, 71, 73 and 83 ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:127-141](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
- Enumeration-oracle exponent: relation collection costs Θ(2^n) for any ℓ, flat within 1.5× across ℓ = 6…10 ([research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md:20-24](research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md)).
- groebner_compare freezes three planted ECC2K17 three-summand Weil-descended inputs (`ecc2k17_five_way.json`) as correctness controls only. They "cannot establish ordinary relation yield… or a full DLP speedup" ([groebner_compare/README.md:107-128](groebner_compare/README.md)).

### Inferences
- Sizes reached algebraically are tiny (N ≤ 32 unknowns at m=4). The n=131 cost tables therefore extrapolate degree models across roughly 100 to 350 additional unknowns (N=209 at the m=3 optimum; 370 for m=4).
- Every "beats rho" figure is either online-only, against plain rather than preprocessing rho, or on a structured odd-characteristic field (JV E(F_{p⁶})).

### Gaps
- No XL/MutantXL or msolve measurement on the descended Koblitz systems appeared in the files read. The novel-directions note says msolve and WDSat are "already used", but I did not locate their result tables.
- No measured solving degree exists for m ≥ 5 or ℓ ≥ 6 in any file read.

---

## Q3. What the ECC2K-130 feasibility study concluded numerically, and what would have to be true to change it

### Takeaway
- With 4-sum relations, IC needs ~4·10²⁵ probes at K=10⁹. That is ~2·10⁷ × rho's 2^60.9 iterations (≈10 GPU-years on one RTX PRO 6000).
- The total-work crossover would need ~1.8·10¹⁶ columns.
- Letting K float makes it worse: about 8·10⁹ core-years, a 10⁸–10⁹× ratio.
- Through Semaev descent, only m ≥ 4 has a free-oracle floor below rho (2^56.4 against 2^60.81). That leaves the oracle 2^5.71 operations per target, a 2^−72 speedup over search.
- To flip the verdict, one of three things would have to hold:
  - a non-generic, non-fix-summand oracle;
  - a solving degree that stops growing with ℓ (and stays near 2–6 at ℓ ≈ 27–34);
  - an evaluable genus-130…300 curve over F₂ carrying the subgroup.

### Cited Findings
- Curve: #E = 4·l with l = 680564733841876926932320129493409985129 (2^129.0) ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:13-19](RESEARCH_ECC2K130_IC_FEASIBILITY.md)).
- Rho reference ([research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md:103-107](research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md); [RESEARCH_ECC2K130_IC_FEASIBILITY.md:92-94](RESEARCH_ECC2K130_IC_FEASIBILITY.md)):
  - 2^60.8090 group operations with the ⟨−1⟩×⟨π⟩ speed-up; 2^64.83 plain.
  - Measured client: 6.9·10⁹ it/s, ≈10 GPU-years.
- 4-sum extrapolation to n=131 with fitted C=1 ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:67-74](RESEARCH_ECC2K130_IC_FEASIBILITY.md)):

| K | probes/relation | core-years |
|---:|---:|---:|
| 600 | 1.1·10²⁹ | 1.1·10¹⁸ |
| 10⁹ | 4·10¹⁶ | 6.8·10¹¹ |
| 10¹² | 4·10¹⁰ | 6.8·10⁸ |

  - Minimum K for any 4-tuple: ≥ 2·10⁷ columns.
  - A practical yield needs K = 10⁹–10¹², which requires an implicit S₃ index, u256 field words and 130-bit sparse LA ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:76-88](RESEARCH_ECC2K130_IC_FEASIBILITY.md)).
- Total work ≈ r/(n²K): 4·10²⁵ probes at K=10⁹, i.e. ~2·10⁷× rho. Crossover K ≈ 1.8·10¹⁶. The conclusion: "Beating rho on total work requires a relation-shape change (higher-arity sums), not engineering" ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:95-132](RESEARCH_ECC2K130_IC_FEASIBILITY.md)).
- End-to-end optimum with LA included: K* = 4·10¹¹, table 3·10²⁷ states, T* ≈ 2.4·10¹⁷ core-seconds ≈ 8·10⁹ core-years, ratio ≈ 10⁸–10⁹ against rho ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:199-227](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
- 6-sum design: probes/relation ≈ 45·r/(n⁴K⁴). At n=131, K ≈ 2000 gives ~6·10⁷ probes/relation with 5.2·10⁸ states. It is blocked by needing a real S₄ solve per probe ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:196-208](RESEARCH_ECC2K130_IC_FEASIBILITY.md)). The step-back note says it "cannot cross" the CGK line ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:420-422](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
- Product law (the identity tying the required oracle speed-up to rho / (m · 2^131)) and oracle budget ([research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md:478-520](research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md)):

| m | free-oracle floor | oracle budget per target |
|---:|---:|---:|
| 2 | 2^89.25 | — |
| 3 | 2^68.58 | — |
| 4 | 2^56.40 | 2^5.71 |
| 5 | 2^48.44 | 2^13.94 |
| 8 | 2^35.61 | — |

  - The required speed-up over search is `rho/(m·2^131) = 2^−(70.19 + log₂ m)`. "Every tuning knob on this attack is inside the identity."
  - Only a Frobenius-stable base (dividing by n) moves it, and none exists at n=131.
- Priced by degree, at n=131 ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:145-167](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)):

| m | degree model | cost | vs rho |
|---:|---|---:|---:|
| 4 | D=2, w=1 | 2^62.88 | +2.07 bits |
| 4 | D=6 constant | 2^76.31 | +15.5 bits |
| 4 | ℓ+4 fit | 2^173.13 | +112 bits |
| 3 | semi-regular | 2^196.95 | +136 bits |

  - The step-back note: "The 106-bit spread between those two readings **is** the security question" ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:262-266](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
- Conditions that would change the conclusion, as the files state them:
  - (a) "an oracle that is not an elimination of the descended system at all", which is unpriced ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:251-254](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md); [research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md:541-545](research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md)).
  - (b) "A degree that stops growing", decided at (15,5) ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:279-296](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)).
  - (c) a Frobenius-stable base, which does not exist.
  - (d) A isogenous to an evaluable Jacobian of genus 130–300 ([research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md:326-360](research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md)).
  - (e) for the table pipeline, a row with log₂(S·T²) − log₂(r/|Aut|) ≤ 0, "the first non-generic signal; every current row is at +25 or more" ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:229-232](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
- Text-quality caveat: the feasibility note contains merge-duplicated paragraphs with slightly different numbers. Examples are C fitted at 0.33–1.51 against C measured at 0.27 at n=61, and 1.4·10¹⁶ against 1.5·10¹⁶ unordered 5-sums ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:36-65](RESEARCH_ECC2K130_IC_FEASIBILITY.md); [RESEARCH_ECC2K130_IC_FEASIBILITY.md:136-154](RESEARCH_ECC2K130_IC_FEASIBILITY.md)).

### Inferences
- At m=4, the D=2 row (+2.07 bits) and the unconverted units (oracle in monomial-column touches, rho in group operations) leave the result close. If touches per group operation exceed about 2^2, the lowest-degree row would price below rho. The degree, not the counting, is the only thing separating m=4 from parity. This is my inference from the stated "subtract log₂ K" instruction ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:172-180](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)).

### Gaps
- K (column touches per GF(2^131) group operation) "is not measured" ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:173-176](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)).
- No preprocessing-rho measurement on the Koblitz rungs was found.

---

## Q4. Which kill arguments rest on heuristics (first-fall-degree assumption and others) and could be revisited

### Takeaway
- Rigorous kills: the product law, the generic CGK/Shoup bound for table oracles, the cyclotomic-coset classification, the GHS trichotomy at prime n, and transport conservation.
- Heuristic or extrapolated kills:
  - all n=131 solving-degree prices (three-rung fits, semi-regular extrapolation, unconverted units);
  - the 4-sum C-constant extrapolation;
  - the "A is not a Jacobian" rarity argument (Torelli codimension);
  - Gaudry's H1 smoothness model in the genus pricing;
  - the N=130 counterfactual's unit-sensitive margin;
  - the free-relation count's uniformity assumption.
- The first-fall-degree assumption itself is not a kill. It is the optimistic side, and the repo's data contradict it.

### Cited Findings
- First-fall-degree (FFD) assumption: Petit–Quisquater's sub-exponential estimate assumes the gap D_reg − D_ff is bounded. The repo measured the gap widening from 2 at ℓ=2 to 4 at ℓ=5 ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:103-109](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)).
- Literature status: a "stalemate, not a negative". Huang–Kosters–Yeo give a reductio. Their conclusion is scoped "for this system", and "the refutation side is heuristic too" ([research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md:28-69](research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md)).
- The repo's recommendation is to stop measuring FFD, because it is "flat at 3; does not bound D" ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:271-273](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md); [RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:423](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
- Degree extrapolations are explicitly labelled ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:355-364](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)):
  - the constant-degree rows assume the degree stops where n ≤ 15 left it;
  - the semi-regular rows assume one-degree agreement persists to N ≈ 209 (m=4: 370);
  - the fit rows are fitted at ℓ ≤ 5 (m=3) and ℓ ≤ 3 (m=4), the latter with one bound. "None of these is measured."
- Measured degrees come from bounded Macaulay refutation (constant 1 in the row space). They are "neither an F4/F5 solving degree nor a certified first fall degree" ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:335-341](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)).
- Several measured fields (n = 4, 8, 9, 10, 12, 15) have intermediate subfields; GF(2^131) has none. The note says no result depends on a subfield ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:342-351](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)).
- Units: oracle cost in monomial-column touches, not converted to group operations ([research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md:172-180](research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md)).
- 4-sum constant: the structural model predicts C=0.25. The fit is 1.2–1.5 at n=71 and 73. "The 5× gap… is unexplained… should be reconciled… before the fit is extrapolated to n=131" ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:49-57](RESEARCH_ECC2K130_IC_FEASIBILITY.md)).
- A-is-a-Jacobian ([research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md:390-417](research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md)):
  - "The codimension count is a heuristic… needs the class number of a degree-260 CM field."
  - The §5 cost model charges one operation per relation trial. A trial is "perhaps 2^7" curve additions, so the 23.6-bit margin survives with 16 bits to spare.
- Gaudry's H1 heuristic: probability ≈1/g! that a reduced divisor splits, used for the floor ([research/notes/index-calculus/RESEARCH_HYPERELLIPTIC_IC_RHO.md:17-30](research/notes/index-calculus/RESEARCH_HYPERELLIPTIC_IC_RHO.md)). At genus 130 the random-polynomial model and the exact zeta model agree to 0.2 bits ([research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md:313-315](research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md)).
- N=130 counterfactual: its −4.71-bit margin "is not robust to the unit conversion the note itself estimates" ([research/ghs_counterfactual_n130_20261013/README.md:32-34](research/ghs_counterfactual_n130_20261013/README.md)).
- The l ≥ 2 rows use a fixed-genus model that is "pessimistic, not tight" ([research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md:375-378](research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md)).
- Gauge-fixing objection to G2: a heuristic prediction, refuted by measurement ([RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md:96-114](RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md)). The lead ladder still uses it as a kill ([docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md:316-323](docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md)).
- Free relations: the count M²/r assumes pair sums land in the subgroup and that x-values are near-uniform. "Relations found this way are correlated… rank must be checked" ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:177-181](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
- Yokoyama et al. is conditional on unproven assumptions and an unproved conjecture, and should not be cited as a lower bound ([research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md:118-140](research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md)).
- The generic kill is "a theorem, not a heuristic", but only for oracles that use group operations and membership tests. It explicitly does not cover an S₄ solve with an algebraic restriction on the unknowns ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:79-89](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).

### Inferences
- The revisitable kill with the most leverage is the m=4 degree model. Three cells, (4,4), (4,5) and (4,6) (plus (15,5) at m=3), would decide between "levelled at 7" and "ℓ+4". The only cheap enabler named is G2 combined with the FGHR 2-torsion symmetrisation.
- The Torelli-codimension rarity argument is the second revisitable heuristic. It is not testable at genus 130 by search, but a ppav count in A's isogeny class is a well-posed computation the notes name.

### Gaps
- No rigorous D_reg bound for these systems exists in either direction ([research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md:62-69](research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md)).
- I found no measurement of how much FGHR symmetrisation lowers the refutation degree on the m=4 Koblitz chain. The 30–416× FGHR figures quoted refer to other settings ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:283-291](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).

---

## Q5. Next steps and open questions the documents list that were never done

### Takeaway
A long tail of named, cheap falsifiers was never run. The repo's own ranking (step-back note §9) puts five first:
1. preprocessing-rho control at equal advice;
2. free-relation census by sorting existing tables;
3. m=4 solving-degree cells at ℓ = 4–6 with FGHR and G2;
4. G2 taken to full F4;
5. D₄-symmetrised S₅ (prime-field).

I found no evidence that items 1–4 were done. The 2026-10-07 PDP protocols E-1 to E-3 and D-1 to D-3 are all PENDING.

### Cited Findings (status as of the latest file I read)
1. **Preprocessing-rho control** on n=73 and n=83 at the IC table's byte budget. Predicted to flip every online ratio below 1. Status: "has not yet been run"; "should be run before any further rung is landed" ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:143-150](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
2. **Free relations by self-collision** in the K=600 pair-root table:
   - Predicted ≥500 verified relations at n=61 and 83; at n=73 it needs K ≈ 1,500; at n=131 it is dead ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:154-195](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
   - Preregistered as E-3; status PENDING ([research/pdp_preregistrations_20261007/README.md:1-50](research/pdp_preregistrations_20261007/README.md); [research/pdp_preregistrations_20261007/PROTOCOL-E3-index-self-collisions.md:5](research/pdp_preregistrations_20261007/PROTOCOL-E3-index-self-collisions.md)).
3. **Solving-degree cells**:
   - (m, ℓ) = (4,4), (4,5), (4,6) with FGHR and equivariant blocks on a large-memory host, "the only experiment whose outcome bears on ECC2K-130 security" ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:409-412](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
   - (15,5) at degree 7 on a ~128 GB host ([research/notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md:838-841](research/notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md)).
4. **G2 to full F4**:
   - Inverse transform of the reduced rows, integration as a `reduce_system` engine with full-F4 agreement checks, then wall-clock against `macaulay_profile` at n=31, dim 16 (where 10–100× is predicted).
   - Block-side optimisations (log/antilog, carry-less multiplication, blocked elimination) also not done ([RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md:160-170](RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md); [RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md:200-214](RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md)).
5. **G1 F4 trace reuse** (Joux–Vitse F4Remake) at the n=31 cells: 637 ms (m=2) and 145 s (m=3), on 16 fresh targets. Status 2026-10-11: "not implemented" ([docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md:42-65](docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md); [docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md:316-328](docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md)). `groebner_compare` provides a learn/apply harness, via Groebner.jl, only for synthetic systems of 1–12 variables ([groebner_compare/README.md:1-40](groebner_compare/README.md)).
6. **G3, sparse Macaulay / BooleanSolve with block Wiedemann** on the n=41 W8 S6 system (222 variables). Report nnz/row and the needed degree; "nothing runs it" ([docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md:104-128](docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md)).
7. **G4, mixed-field F4** (keep chain intermediates over F_{2^n}). A one-day negative control at n=19, m=3; no result found ([docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md:130-148](docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md)).
8. **A2, cyclotomic-sparse factor bases**: x = Σ_{k∈S} γ^k with γ a 263rd root of unity. The membership test has degree 4 in ONB coordinates. Falsifier at n=23 (ladder n = 11, 23, 83) ([docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md:187-219](docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md)). Status is inconsistent across files:
   - The step-back note groups cyclotomic-sparse bases with others "measured or derived to the generic line" ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:362-365](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
   - The transfer-aware note still lists A2 as a lever "where the next rung comes from" ([docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md](docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md) §5).
9. **N1, group-algebra lift** to F₂[Z/263] (at n=23, F₂[Z/47]): a one-afternoon falsifier; no result found ([docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md:225-247](docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md)).
10. **L6, phase-aware quotient summation** (τ-slots on the T₂-symmetrised system): "not implemented… Two to three hours of Rust" ([docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md:285-304](docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md)).
11. **L8 open items** ([docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md:352-355](docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md)):
    - E2: floor curves at n = 41, 61, 83 under compact-orbit;
    - E3: C37 at four levels;
    - E8(2): odd-N Kani prototype;
    - E9: equal-precompute rho across levels.
12. **Algebraic τ-slots for the no-table regime** (`extract128` with per-summand rotations): "the one item on this list that touches the challenge rung" ([docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md:463-468](docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md)).
13. **Higher-arity extraction and supporting infrastructure** ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:136-174](RESEARCH_ECC2K130_IC_FEASIBILITY.md); [RESEARCH_ECC2K130_IC_FEASIBILITY.md:196-208](RESEARCH_ECC2K130_IC_FEASIBILITY.md)):
    - an m=5 compact-orbit extraction oracle ("main algorithmic lever… unmeasured");
    - an implicit S₃ index;
    - u256 field words and 130-bit sparse LA;
    - the 6-sum S₄-tree design;
    - the remaining rungs n=85 a=0 (K=600) and n=89 a=0 (K ≈ 850).

    The step-back note recommends stopping these as generic ([RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md:420-422](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md)).
14. **H3, symbolic S₄ links** in `koblitz_groebner`, and the scaling target's finish line (primary metric under 64 at n=31 and n=63; H5 falsified at n=63, ℓ=6). Not done ([research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md:415](research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md); [research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md:444-453](research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md)).
15. **Crossbred**: "one more rung decides" Route 1. Route 2 (Crossbred on the symmetrised u-frame system: 13 unknowns at degree 4 against 30 at degree 3, n=15) and Route 4 (chain the symmetrised S₃ for m ≥ 4, `4(ℓ−1)+1+n` unknowns) are listed open ([research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md:84-112](research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md); [research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md:147-164](research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md)).
16. **GGMP open items**: block/homogeneous Gröbner structure for m ≥ 3, complexity across characteristics, multiple invariant bases ([RESEARCH_KOBLITZ_INDEX_CALCULUS.md:907-915](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)).
17. **A-isogenous-to-Jacobian**:
    - Count principally polarised members of A's isogeny class (needs the class number of a degree-260 CM field) ([research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md:392-399](research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md)).
    - Enumerate the toy ladder at n=5 for non-hyperelliptic genus 5 and genus ≥ 6, which is not done ([research/notes/ecc2k130/RESEARCH_ECC2K130_CURVE_CONSTRUCTION.md:135-136](research/notes/ecc2k130/RESEARCH_ECC2K130_CURVE_CONSTRUCTION.md)).
18. **GHS → Smith**: read Diem (2003) and Momose–Chao on odd-characteristic descent to find (p, n) giving genus-3 hyperelliptic targets. The record is open and unassigned ([problems/ghs-to-smith-pipeline-characteristic-mismatch.md:56-67](problems/ghs-to-smith-pipeline-characteristic-mismatch.md)).
19. **Cubic-GHS quartic IC**:
    - Price the rho arm, which currently records steps only, so an S ratio exists ([research/ghs_cubic_quartic_ic_20261012/README.md:25-28](research/ghs_cubic_quartic_ic_20261012/README.md)).
    - Make degree-3 Pic canonical forms unique.
    - The dense factor-base-only line IC does not reach the ℓ-part ([problems/plane-quartic-jacobian-class-arithmetic-missing.md:28-33](problems/plane-quartic-jacobian-class-arithmetic-missing.md)).
20. **PDP protocols D-1 to D-3 and E-1 to E-2**: all PENDING ([research/pdp_preregistrations_20261007/README.md:38-47](research/pdp_preregistrations_20261007/README.md)). They are:
    - D-1: Artin–Schreier coordinates;
    - D-2: four-summand linearisation reach;
    - D-3: filtering against the LA cap;
    - E-1: packed root index;
    - E-2: batch verification.
21. **Koblitz subfield experiment queue**: threads 1–12, among them orbit-coordinate decomposition, Frobenius-aware Macaulay, trace-zero × Frobenius bases, hybrid SAT → Gröbner, and partial Weil descent for n = ab. These are tracked with an "immediate order", but I did not locate completion records ([research/KOBLITZ_SUBFIELD_INDEX_CALCULUS_QUEUE.md:9-46](research/KOBLITZ_SUBFIELD_INDEX_CALCULUS_QUEUE.md)).
22. **Boolean representation audit, next bounded questions**: compare a support-sharing representation; test narrow-interaction preservation before chordal methods; get a quieter host ([boolean_representation_audit/README.md:130-137](boolean_representation_audit/README.md)).
23. **4-sum constant**: reconcile the 5× structural-vs-fit gap in C against first-hit distributions over many targets ([RESEARCH_ECC2K130_IC_FEASIBILITY.md:54-57](RESEARCH_ECC2K130_IC_FEASIBILITY.md)).

### Inferences
- The highest-value never-run items for binary curves are items 3, 4 and 2:
  - Item 3 decides the 106-bit spread.
  - Item 4 is the enabler for item 3 and the only lead with an unbounded measured trend (31× to 128× from n=11 to n=23).
  - Item 2 is a near-zero-cost sort that would make the precompute charge honest.
- Item 1 (preprocessing rho) is essential for honest accounting but is predicted to erase the online wins rather than create a lead.
- Items 13 and 15 are labelled generic or engineering by the repo's own later analysis. Revisiting them only matters if paired with an algebraic restriction on the unknowns.

### Gaps
- I searched by keyword only (grep for "free relation", "self-collision", "equivariant", "preprocessing rho") and in commit messages. Absence of a result file does not prove an item was never run. Some may live in the sibling repos (`cryptanalysis`, `crypto-autoresearcher`), which the docs mention but which are outside this scope.
- The status of A2 (cyclotomic-sparse bases) is contradictory between the step-back note and the transfer-aware note. I could not find a measurement file that resolves it.
