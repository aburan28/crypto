# Prime-field and general ECDLP index-calculus leads in aburan28/crypto: inventory, status, kills, next steps

Scope: the local repo at `/home/user/crypto`, read on 2026-10-14. Sources are the top-level `RESEARCH_*.md` files, `docs/ic/`, `research/notes/index-calculus/`, and dated study directories under `research/`. All source links are repo-relative paths. Git dates are merge dates. The clone's history starts on 2026-09-28, so it is probably shallow, and documents give their own dates, which are quoted below. Binary, Koblitz, ECC2K-130 and Weil-descent leads are covered by another researcher and appear here only where they touch prime-field or general leads.

---

## Q1. Which IC leads exist for prime-field and general ECDLP, and what has each measured? (inventory)

### Takeaway
The repo has a dense prime-field program in three tiers.
- **Tier 1, prime field `F_p`.** It is closed. Every enumeration or table oracle sits on the generic `S·T² ≳ n` line. The only algebraic levers are the PKM tower and Coppersmith small roots, and neither reaches the reach needed. The first arity that could beat rho is `m = 4`, which needs an `S₅` oracle at `δ ≈ 1/5`. The best lattice family reaches `1/80`.
- **Tier 2, extension fields `E(F_{p^k})`, generic curves.** This is the only place the repo has a measured sub-rho component. The `k = 4` linear algebra runs at 0.518× rho, and the handover is extrapolated to about `2^{135}`–`2^{156}` depending on symmetry and on which rho is the reference.
- **Tier 3, weak classes and covers.** The Joux–Vitse cover on `E(F_{p⁶})` is the only measured end-to-end below-rho result: 0.012–0.014× rho at `ℓ = 2^{61.4}`.

Isogeny, fiber, transfer and conductor-gap leads are all closed as constant-factor at most for prime fields.

### Cited findings

**A. The generic ceiling and the decision rule**
- Corrigan-Gibbs–Kogan `S·T² ≳ n` and Shoup `T ≳ √n` cover rho, BSGS, kangaroo and every IC whose decomposition enumerates or table-looks-up factor-base elements. IC-2 direct subtraction has `S·T² = p²/(4B) ≥ p/4`. IC-3 with a pair table has `0.56p²/B² ≥ 0.56p`. IC-`m` is `≥ p`. — [Source](RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md) §1
- With the linear algebra included, the cost is `(m!/2^m)·p·T_dec/B^{m−1} + m·B²`:
  - `m = 2` never beats `√p`, even with `T_dec = 1`;
  - `m = 3` at best ties, and only if `T_dec = O(1)`;
  - `m = 4` beats rho iff `T_dec ≪ p^{1/4}`, at `B ≈ p^{1/5}` (δ ≈ 1/5);
  - `m = 5` beats rho iff `T_dec ≪ p^{1/2}`.
  — [Source](RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md) §2.4
- The 2026-10-10 synthesis restates the bar: "the first arity that could beat rho is `m = 4`, needing an `S₅` small-root oracle at `δ ≈ 1/5` against a measured Coppersmith reach of `1/80`". It also says every `S₃`/`S₄` constant-factor result should be "labelled engineering". — [Source](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md) §0 item 5, §8

**B. Native-word engineering of prime rho and IC (PR #1560, 2026-10-08)**
- Montgomery `u64`, batched inversion, a parallel DP rho, and direct-subtraction or pair-table decomposition give:
  - rho at 28 bits 357×/466× faster than the BigUint baseline;
  - IC-2 at 20 bits 704× faster (352 s → 0.50 s);
  - IC-3 per relation at 16 bits ~10⁷× faster (42 s → 5 µs);
  - IC-3 end to end at 32 bits in 20.5 s, against 0.5 s for rho.
  — [Source](RESEARCH_PRIME_FAST_ECDLP.md) §0, §3–§5
- The verdict is "No vs_rho crossover … the gap grows as `p/B` vs `√p`". The asymptotics are unchanged. — [Source](RESEARCH_PRIME_FAST_ECDLP.md) §5
- §7 lists further speed-ups that were not implemented:
  1. two-level batching of `R ← R + G` (~2×);
  2. Lanczos or Wiedemann over `u64` for `B > 2000`;
  3. rho on an idle host with W = 4096;
  4. a 256-bit path;
  5. GLV/automorphism folding for `j = 0`, plus preprocessing rho;
  6. fixing the false negation-map claim in `examples/eccp79_rho.rs`.
  An idle-host sample is "still owed". — [Source](RESEARCH_PRIME_FAST_ECDLP.md) §6–§7

**C. Root finding in place of Gröbner (2026-10-12)**
- The prime `S₃` root in the third variable is exactly `x(R ∓ F_i)`, so the oracle is one batched inversion plus `|F|` lookups.
- Against the crate's Buchberger at `|F| = 6` (16–24 bits) it is 35,592–46,911× faster per decomposition.
- Against `groebner_f4` on P-192…P-521 at `|F| = 6` it is 3,021–5,802× faster.
- End to end against msolve it is 100× at 16 bits, 105× at 20 bits and 50× at 24 bits.
- The doc says itself: "No security claim … the asymptotic group-operation count of index calculus is unchanged."
— [Source](RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md) §0, §2, §4, §6, §7
- The end-to-end runs of `m = 3, 4`, the NIST curves, and Buchberger above 16 bits are listed as "Not measured end to end". — [Source](RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md) §7
- A GPU port for sm_120 was validated only through host emulation:
  - 21,811 targets and 1,331 relations, 0 mismatches;
  - "The CUDA path has NOT been executed on a GPU", and the throughput of 5×10¹⁰–10¹¹ cells/s is an ESTIMATE;
  - follow-ups: compile with nvcc, run `./bench selftest`, swap in the 32-bit CIOS `fp256.cuh`, and add a global hash table for `|F| > 2048`.
  — [Source](RESEARCH_GPU_ROOTFIND_DECOMPOSITION.md) §0, §3–§5
- ecbench end to end with PKM's Kummer tower base (`prime-kummer:t=T`), 20–28 bits, 896 runs all verified:
  - IC costs 29–141× rho with the negation map;
  - the ratio grows 0.26 ± 0.28 bits per bit of `r`, against the predicted 0.17;
  - the tower base costs the same as a plain base, within 1%: "PKM's structure helps an algebraic solver, not a table".
  — [Source](research/pkm_tower_ecbench_20261013/README.md)
- On binary toys (2026-10-14), complete cold IC still loses to strong rho, about 1.54× (S 21.8 against 14.1). — [Source](research/ecbench_rootfind_vs_rho_20261014/README.md)

**D. Probe A: preprocessing rho on prime curves (2026-10-08, measured)**
- Online steps track `2√(n/S)`. Per-target step ratios against plain rho are 11×, 38×, 42× and 60× at 32, 36, 40 and 44 bits. The break-even target counts are 1, 2, 18 and 116.
- A first run sat 2–7× above the line because the online walk restarted at fruitless cycles. The fix is that "the online and precomputation walks must be identical maps".
— [Source](RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md) §3.1

**E. Coppersmith/lattice small roots for Semaev polynomials**
- The derived reach is `δ_max = 1/(m(m+1)·2^{m−2})`: 1/6 for `S₃`, 1/24 for `S₄`, 1/80 for `S₅`, 1/240 for `S₆`. Against what is needed (1/2, 1/4, 1/6, 1/8) the shortfall is 3×, 6×, 13× and 30×. — [Source](RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md) §2.2
- Probe B measured `S₃` at 40 and 56 bits, levels 1–3:
  - Howgrave-Graham holds to δ = 0.10 and never beyond;
  - the heuristic criterion fades between 0.12 and 0.16;
  - level 3 costs 25–75× per reduction for +0.02–0.04 in δ.
  — [Source](RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md) §3.2
- Probe C: extended Jochemsz–May shifts and the symmetric `(e₁, e₂)` lattice reach "the same threshold". The symmetric lattice is 3–10× cheaper per reduction. — [Source](RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md) §3.3
- W10 (2026-10-11), an exact LP/min-cut optimiser over full-rank triangular single-polynomial lattices:
  - asymptotic δ* is 1/6, 1/24 and 1/80 for `S₃`, `S₄`, `S₅`;
  - the family-wide ceiling is `1/(m·d)` = 1/4, 1/12, 1/32;
  - "The Coppersmith route to the m = 4 oracle is closed" for that family;
  - `S₅` is dense on its `[0,8]⁴` box, all 6,561 monomials.
  — [Source](research/lit2026_w10_s5_coppersmith_20261011/README.md) §2, §6
- W10 corrects a mislabel. The "measured 1/80" in the 10-10 leads note is a derived Jochemsz–May limit. Only `S₃` was ever measured. — [Source](research/lit2026_w10_s5_coppersmith_20261011/README.md) §1

**F. The PKM tower oracle (prime-field algebraic factor base, Petit–Kosters–Messeng 2016)**
- The pre-registered prediction (H0) was that F4's solving degree grows linearly in `N = m·t`. The pilot (2026-09-24) found `D` at 4–5 for `m = 2` through `N = 18`, which it read as inconclusive.
- Measured plateaus:
  - round 2: `D = 6` at `N = 20, 22` for `m = 2`; `D = 7` at `N = 15` for `m = 3`;
  - round 3: `m = 4`, `N = 16` gives `D = 7` (55 min, 8.1 GB);
  - round 7 (2026-10-06): `m = 3`, `N = 18` gives `D = 7`, at 96% of a 14 GB cap.
- "`N = 21` needs a bigger machine." Every plateau is read as A1-inconclusive. A signature-based F5/GVW engine needed rows 1–2 degrees higher in 18 of 18 systems and was not adopted.
— [Source](research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md) status block, §§10–16
- If the solving degree stays bounded (H1), "prime-field index calculus becomes subexponential under a measured heuristic". — [Source](research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md) §0
- From the requirements document (2026-10-12):
  - F4 per-target cost grows 1.6–2.4 bits per tower variable on the prime systems, against 0.33–0.75 for enumeration and the pair table;
  - a dedicated root-finder is "the open problem PKM16 themselves left";
  - candidate mechanisms C1–C8 are listed, and C2 (sparse-spectrum zero finding on the Kummer group algebra) is "the only entry whose success would move β".
  — [Source](research/pdp_root_finder_requirements_20261012/README.md) §0, §8
- C4, the representation technique on Kummer exponents, was falsified on 2026-10-13. Exponent sums of decompositions are indistinguishable from random at `k = 1..6`, while the positive control reads 1.000. "C2 … Not run yet. Solving-degree census with symmetries. Not run yet." — [Source](research/ic_exponent_falsifiers_20261013/README.md)

**G. The P-256-class PKM/Dickson-torus program (24 rounds, 2026-10-04 to 10-06)**
- A depth-18 Dickson-torus factor base for P-256 was reproduced natively: 131,239 negation-folded columns. — [Source](research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_fb1_20261004/RESULT.md)
- Round 24, an `S₁₇` cross-colour collision boundary with the negation quotient. It "reaches ideal collision-constant parity only for K=1", at 0.964× rho (a known-log lower bound, not an attack). The geometric bases stay at 392–395× rho. Fixed-row collection is `2^{137.04}`, 16.96 bits above the `2^{103}` gate. One-list storage is about `2^{142}` bytes. — [Source](research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_negation_cross_colour_round24_20261006/RESULT.md)
- Rounds 20 and 22 carried accounting errors (a double count of 138,031 rows, and `√K` scaling) that were corrected in later rounds. — [Source](research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_factor_base_symmetry_round22_20261006/INVALIDATED.md)

**H. Extension fields `E(F_{p^k})`: the Gaudry/Diem "rho parity" ledger (general ECDLP)**
- `k = 3`: the plain method's best is 528× rho at `2^{33.1}`. Its linear algebra grows as `n^{0.68}`, and it bottoms near 200× rho around `2^{50}`. With double large primes the exponent is `n^{0.424}`, and parity is extrapolated past `2^{230}`. Four `C₃` levers bought 5.5×. — [Source](research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md) §11.15
- `k = 4`:
  - `C₄ = 1.227·10¹²` `F_p` multiplications per `S₅` solve, 17× its floor, with the crossover near `2^{150}` (§11.17);
  - linear algebra measured at `r∞ = 0.518 ± 0.031` of rho on eight curves, all verified (§11.19);
  - the accounting note says that against the negation walk this reads 0.733, a 1.36× edge, with the handover near `2^{156}`.
  — [Source](research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md) §11.17, §11.19
- E18 (`j = 0`, `ψ` fold): `r∞ = 0.138 ± 0.008` against the matched walk (plain 1.232). The handover moves from `2^{156.1}` to `2^{146.0}`. "No end-to-end `k = 4` run below rho exists." — [Source](research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md) §9.3–9.4
- `k = 4` Joux–Vitse three-point decomposition: about 6,945× rho, flat. Parity needs `C′ < 21` against a measured 1.6·10⁵. Verdict: "never". — [Source](research/notes/index-calculus/RESEARCH_RHO_PARITY_PROGRAMME.md) header
- `k = 5`, Joux–Vitse four-point:
  - plain `C″ = 3.0·10¹⁰`;
  - Edwards 2-torsion `C″ = 7.3·10⁷`, a 416× cut;
  - parity extrapolated at `2^{302}` plain and `2^{233}` with 2-torsion; parity at `2^{128}` would need `C″ < 5·10⁴`;
  - the `Q₄` 4-torsion lever was retracted.
  — [Source](research/notes/index-calculus/RESEARCH_K5_TORSION_JOUX_VITSE.md) header
- `j = 0`, `k = 3` `⟨ψ⟩` quotient: exactly ÷3.0 at every size, still 220–540× behind rho. The `Z/3` grading costs 12–30× in F4 multiplications. — [Source](research/notes/index-calculus/RESEARCH_GLV_INDEX_CALCULUS.md) header
- Shape replay (record and replay, F4Remake-style) of the `k = 4` Macaulay solve saves only 2.0–4.7% of multiplications, about 7% extrapolated to `d = 8`. "Record-and-replay does not cut `C₄` 10×." Engineering gave a 6× wall-clock speed-up with identical counts. — [Source](research/gaudry_quartic_replay_20261013/README.md); [Source](research/gaudry_quartic_speed_20261013/README.md)
- C8, the `D₄ = (Z/2)³ ⋊ S₄` symmetry on `S₅` over `F_{p⁴}`, run on 2026-10-13:
  - Bézout 4,096 → 512;
  - Macaulay degree 29 → 14;
  - `C₄` from 1.21·10¹² to 2.54·10¹⁰, a 47.5–48.4× cut (the paper projection was about 630×);
  - decomposition rate falls to 0.133 pooled, about 1/8;
  - the modelled crossover moves from `2^{150}` to `2^{135.1}`–`2^{135.8}` at `h = 2`;
  - one planted target in 20 is degenerate (quotient dimension 517), cause unknown.
  — [Source](research/gaudry_quartic_d4_20261013/README.md)

**I. Weak class and cover: Joux–Vitse cover-decomposition on `E(F_{p⁶})` (general ECDLP, weak class)**
- `S/rho` falls as `p^{−1.92 ± 0.12}`, reaching 10.2× and 8.1× at `p = 1009`.
- With the sieved relation phase the crossover is measured at `p ≈ 430` (`ℓ ≈ 2^{50}`), and the route reads 0.012–0.014× rho at `ℓ = 2^{61.4}`.
- The walk that reaches the class costs more than rho below `p ≈ 8,000`.
— [Source](research/notes/index-calculus/RESEARCH_COVER_DECOMPOSITION_LEDGER.md) lines 113, 142

**J. Isogeny, fiber and transfer leads (general ECDLP)**
- Isogeny fibers (2026-10-10/11). Fibers are torsors, and the summation polynomial of `E₂` pulls back to the product of kernel-twisted ones (Thm 4.2). The ten hypotheses ranked H1–H10:
  - H1 (torsion-symmetrised quotient factor base) is bounded at ×|M^π| ≤ 4;
  - H5–H10 are negative.
  — [Source](RESEARCH_ISOGENY_FIBERS_ECDLP.md) §0, §6
- A1 at `m = 31`:
  - the projected fold halves the unknowns (39 → 19);
  - but the hit rate halves too, so the per-relation metric is 2.95× worse at `m = 3`;
  - the trace-zero base is 1.65× cheaper in S;
  - every arm is 22–223× worse than rho;
  - the 1.5× kill criterion is not met.
  — [Source](RESEARCH_ISOGENY_FIBERS_ECDLP.md) §7.A1-R
- Round 4 (2026-10-11):
  - "isogenies change constants, never exponents" is false over composite `F_{2^n}`;
  - Conjecture 10: over prime fields no class member is weaker in exponent;
  - class-wide multi-target rho at `L = 64` costs 0.998–1.026 of single-curve `L`-target rho;
  - the large-prime conductor gap protects only when the kernel is not rational, with crossover near `ℓ ≈ 2^{21.5}`–`2^{24.7}` at `r = 2^{128}`.
  — [Source](RESEARCH_ISOGENY_FIBERS_ECDLP.md) §17
- Transport conservation lemma: homomorphic transfers cannot change decomposition counts beyond `|κ|^{m−1}`, the multiplier set and `|F|`. The `Z[τ]` multiplier closure is a precompute trade only, M1–M7. — [Source](docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md) §1, §4
- Isogenous neighbours have identical Macaulay ranks. Yield moves by ≤ 1.3× through `|F|`. For the prime-field enumeration PDP the model "does not enter at all". — [Source](RESEARCH_PDP_SPEEDUP_AND_ISOGENY.md) §2.3–2.4
- Isogeny-loop GLV on ECCp-79…ECCp-131:
  - no loop within budget, exhaustively;
  - the cheapest relation costs 20–93× the GLV saving;
  - "Loops are a scalar-multiplication tool (GLV), not a DLP tool".
  — [Source](RESEARCH_ISOGENY_LOOP_GLV.md) §4.2, §4.4

**K. Other closed prime-field probes**
- Symmetric `S₃` rewrite (asymptotic reach 1/9 < 1/6) is closed on paper. — [Source](RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md) §4.6
- Multiplicative-subgroup bases (`x^B = 1`, no exponent change) are closed on paper. — same source
- xedni lifting is closed on paper (JKSST 2000). — same source
- Embedding into `F_{p^n}` is closed on paper (`p^{2−2/n}` against `√p`). — same source
- The zero-minor guess-and-determine method (Mahalanobis, arXiv 2607.09814) has slope 0.905. Its cost ratio to rho grows from 260× to 6,715× over 16 to 26 bits. — [Source](research/lit2026_w12_flagged_claims_20261011/README.md)
- The P-256 backdoor map ranks four hypotheses. "UnknownAlgorithm" is "NotPubliclyTestable". The prime-regime IC ladder "remains above rho at every measured size". — [Source](RESEARCH_P256_BACKDOOR_SPECULATION.md) §4
- The a=−3 deployed-shape IC ladder runs at 16, 20 and 24 bits. The IC/rho gap widens from 19× to 403×. "vs_rho: Not achieved". — [Source](docs/ic/BOUNDARY_TARGETS.md) Regime C rows

### Inferences
- The repo's prime-field IC is closed at three levels:
  1. generic or table oracles, closed by theorem;
  2. lattice small roots, closed by an asymptotic LP for one family;
  3. PKM algebraic towers, inconclusive, stalled by memory at `N ≥ 21`.
- Only level 3 still has an unresolved exponent question (bounded or growing `D`) on real prime-field systems.
- The extension-field `k = 4` line is the closest thing to an exponent-level "positive". Its linear algebra is already below rho, so whether a generic IC route ever beats rho there depends on `C₄`. The `D₄` result shows the symmetry lever is real but about 13× smaller than the `D^{3.1}` projection.

### Gaps
- No prime-field PKM tower point exists at `N ≥ 21` (memory-bound). No solving-degree census with symmetries was found.
- No end-to-end `k = 4` index-calculus run (relations, folded columns, linear algebra) exists in ecbench. Only the model and the linear-algebra panel exist.

---

## Q2. Which leads were killed, how rigorous was each kill, and which assumptions could be revisited?

### Takeaway
Most kills are rigorous in form: a generic-group theorem, an exact LP over a stated lattice family, or a proved lemma on homomorphic transport. But each is scoped. Several carry explicit unexamined assumptions or unrevised numbers that matter:
- the lattice family excludes multi-polynomial and sublattice strategies;
- the PKM "inconclusive" rests on a short plateau;
- the `D₄` crossover uses an `r∞` band that was later revised;
- the 10-10 note's "measured 1/80" is derived, not measured;
- the `k = 4` handovers are against different rho references.

### Cited findings
- **Generic-oracle kill, rigorous (theorem).** It rests on Corrigan-Gibbs–Kogan and Shoup. It covers the `subtract`/`mitm` oracles, residual walks, pair tables and preprocessed variants. — [Source](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md) §1 Q1; [Source](RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md) §1
- W3 maps the newer structured-GGM bound (Corrigan-Gibbs–Henzinger–Wu, ePrint 2026/384) onto decomposition oracles as `δ ≤ c·Σ C(|F|+j−1, j)/r`. That mapping is the authors' derivation: "Unique factorization does not hold for EC sums … *Interpretation, not a theorem of the paper*". The EC section of 2026/384 "could not be extracted". — [Source](research/lit2026_w3_ggm_audit_20261011/README.md)
- **Linear-algebra-aware arity kill, heuristic (cost model).** It assumes `Θ(m·B²)` sparse linear algebra: "nothing sub-quadratic is known for these random-sparse systems". — [Source](RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md) §2.4. *Revisitable:* the `m = 2, 3` "never" verdicts assume the quadratic linear algebra and the `(2B)^m/m!` yield model.
- **Coppersmith kill, rigorous for a stated family only.** W10 §7 lists the routes it does not close:
  - non-full-rank or sublattice strategies (Boneh–Durfee, Herrmann–May unravelled linearisation);
  - cuso Gröbner-optimal shifts that cancel out-of-`M` monomials;
  - more than one polynomial (adding the curve equations, sharing unknowns across targets);
  - unequal or non-box bounds;
  - finite-`p` heuristic effects for `S₄`/`S₅`, which were not measurable;
  - the 2026/1027 LP itself (HTTP 403, not reproduced).
  — [Source](research/lit2026_w10_s5_coppersmith_20261011/README.md) §7
- **Algebraic factor-base kill (§2.3), an argument.** Any subset of size `B` has a membership polynomial of degree `B`. The falsifier for a structured base is the solution count of `{S_{m+1} = 0, L(x_i) = 0}` over `F̄_p` (if `≥ B^{m−1}`, the base is dead). The leads note lists this as an un-run falsifier. — [Source](RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md) §2.3, §4 item 2; [Source](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md) §8 item 2. No run of this census was found in `research/`.
- **PKM tower: not killed.** The pre-registered H0 (linear growth) was not confirmed. `D` reads 6, 6, 7, 7 over `N = 9`–18 at `m = 3`, "one level without a rise, which A1 reads as inconclusive". — [Source](research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md) §16 summary. *Revisitable:* the note's own lower bound (§3.3) says the XL degree is at least `d + δ` with `δ ≈ N` generically. The measured F4 plateau at 7 conflicts with that naive expectation, and that disagreement is unexplained.
- **"Measured 1/80" mislabel.** The 10-10 leads note calls 1/80 "measured" (lines 55–56 and 384–385). W10 says it is derived and that only `S₃` was measured. — [Source](research/lit2026_w10_s5_coppersmith_20261011/README.md) §1
- **`D₄` crossover assumption.** The `D₄` study used `r∞ (φ = 0.73 … 1)` of 0.085–0.158 per summand-normalised term, which is §11.16's derived band. — [Source](research/gaudry_quartic_d4_20261013/README.md). §11.19 later measured `φ = 0.897` and `r∞ = 0.518`, and E18 says the negation walk moves `r∞` to 0.733. — [Source](research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md) §11.19; [Source](research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md) §9.4. *Revisitable:* the `2^{135}` figure should be re-derived with the measured `φ` and the negation-walk reference. Both inputs raise the crossover.
- **`D₄` prediction against measurement.** The `D^{3.1}` scaling projected about 630× on `C₄`. The measured cut is 48×. — [Source](research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md) §9.4; [Source](research/gaudry_quartic_d4_20261013/README.md). The `C ≈ D^{3.1}` heuristic therefore overstates symmetry gains by about 13×. That matters for the projected Faugère–Svartz and `ψ`-grading gains.
- **Record-and-replay kill, measured.** It is bounded at 2–7% of multiplications because zero rows are a small share. Wall-time compilation remains open. — [Source](research/gaudry_quartic_replay_20261013/README.md)
- **Isogeny and transfer kills, proved lemmas.** The transport-conservation lemma and Thm 4.2 are stated with proofs at research-note rigour. Conjecture 10 (no class member is weaker in exponent over prime fields) remains a conjecture. — [Source](docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md) §1; [Source](RESEARCH_ISOGENY_FIBERS_ECDLP.md) §17.7. Prop E (large-prime conductor gap as a reduction) is "the one place where JMV-type equivalence is genuinely open". — [Source](RESEARCH_ISOGENY_FIBERS_ECDLP.md) §9
- **Prop 3.3 correction.** "ℓ ≳ 2⁶⁰ for a 128-bit group" was loose. The generic-kernel crossover is near `2^{22}`–`2^{25}`. — [Source](RESEARCH_ISOGENY_FIBERS_ECDLP.md) §17.3
- **Accuracy of the k=4 rho reference.** §11.19's `r∞` is against `rho4` without the negation map, at S ≈ 1.3. — [Source](research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md) §11.19 accounting note
- **Self-correcting accounting in the P-256 series.** Projections from rounds 20 and 22 were invalidated by a double count and by `√K` scaling. — [Source](research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_negation_cross_colour_round24_20261006/RESULT.md)

### Inferences
- The kill most worth revisiting for exponent potential is the Coppersmith family limit, specifically the multi-polynomial case: adding the curve equations or a shared-target lattice. W10 excludes it explicitly. Second is the PKM tower's unexplained degree plateau.
- The combined-system solution-count falsifier is cheap, named twice, and apparently never run.

### Gaps
- I could not confirm whether any run tested the multi-polynomial Coppersmith variant. Nothing was found under `research/`.
- The `D₄` README does not state which rho variant underlies its `r∞` band. My inference that it is the unrevised band comes from the numbers (0.34–0.63 divided by 4).

---

## Q3. What does each document list as next steps, and was it done later?

### Takeaway
Of the 10-10 synthesis's five "Do" items:
- the `D₄` `S₅` symmetry was done (2026-10-13), but over `F_{p⁴}` rather than as the prime-field lattice probe;
- the Faugère–Svartz equivariant F4 is WIP (Koblitz);
- the preprocessing-rho control on Koblitz, the free-relation census, and the `(m, ℓ) = (4, 4..6)` solving-degree cells were not done as of 2026-10-12/14.

The extension-field `k = 4` program has a clear unfinished queue:
- end-to-end ecbench;
- a compiled straight-line `S₅` kernel;
- the `ψ` `Z/3` grading at `k = 4`;
- investigation of the degenerate `D₄` target;
- a Magma cross-check of the `D₄` system.

### Cited findings
**The 10-10 leads note "Do" list** — [Source](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md) §9
1. *Preprocessing-rho control at equal advice, Koblitz `n = 73/83`.* Not done. On 2026-10-12, BOUNDARY_TARGETS still lists "then a Bernstein–Lange rho reference arm at n=61+" as the next step. Meanwhile an `n = 61` online win of 12.9–13.9× (setup excluded) was landed. — [Source](docs/ic/BOUNDARY_TARGETS.md) line 147. W3 records "Pair-root tables / self-collision (Lead B) … Not measured". — [Source](research/lit2026_w3_ggm_audit_20261011/README.md) line 179
2. *Free-relation census by self-collision.* A protocol was pre-registered on 2026-10-07 (`research/pdp_preregistrations_20261007/PROTOCOL-E3-index-self-collisions.md`). As of 10-11, W3 records it as not measured.
3. *Solving-degree cells `(4, 4..6)`.* Binary, out of scope. Not found as done.
4. *Equivariant F4 to full F4.* "WIP (uncommitted worktree `ict-g2-equivariant-20261012`)". — [Source](docs/ic/IC_TECHNIQUE_INTEGRATION_20261012.md) status table
5. *`D₄`-symmetrised `S₅` lattice probe at `m = 4` (prime field).* Not done as a lattice probe. The `D₄` symmetry was instead applied to Gaudry `k = 4` over `F_{p⁴}` (C8, 48×). — [Source](research/gaudry_quartic_d4_20261013/README.md)

**Prime-field breakthrough program §4/§5**
- Pursue (1), small roots beyond Coppersmith, "only with a concrete new lattice/embedding idea". Keep (2)–(3) closed unless new structure appears. Use (4) to record `S·T²` per ledger row. — [Source](RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md) §5. W10 did (1) for one family.
- Recording `S·T²` per row: the 10-10 note also asks for it (§4, §8 item 3). I found no evidence that it was adopted in BOUNDARY_TARGETS.

**Gaudry `D₄` "Open"** — [Source](research/gaudry_quartic_d4_20261013/README.md)
- The degenerate target (dimension 517): which `Y_R`, and how often.
- The `k = 4` pipeline end to end in ecbench.
- `C₄` against Magma on the `D₄` system, which is not exported yet.

**Scope note's 5-step programme (2026-10-13)** — [Source](research/ic_speed_and_rho_scope_20261013/README.md) §4
- Step 3 (`E(F_{p⁴})` in ecbench) has no evidence of being done.
- Step 4 (compiled `S₅` solver): the replay part is done (2–7%). The compiled straight-line kernel is open.
- Step 5 (the `k = 4` ladder in ecbench, gated on a ≥10× cut in `C₄`) is now unblocked by `D₄`'s 48×, per the `D₄` README's "Decision rule: … the next step is the `k = 4` ladder in ecbench". No ladder was found.

**GLV invariant G2 "next lever, named and not yet registered"** — [Source](research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md) §9.4
1. `ψ` on the system: the `Z/3` grading of `S₅` on `j = 0` curves at `k = 4`, which "needs re-deriving". Not found as done.
2. `D₄` on the system: done 2026-10-13. The E18 caution was that those curves pay about 8 bits back through cofactor and yield.

**Exponent falsifiers (2026-10-13)** — [Source](research/ic_exponent_falsifiers_20261013/README.md)
- C2 (sparse-spectrum zero finding): not run.
- Solving-degree census with symmetries: not run.
- Weighted linear forms `Σ c_i e_i` on Kummer exponents: untested; only the plain sum was tested.

**PKM tower** — [Source](research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md) §16
- `m = 3` at `N = 21` needs a bigger machine.
- `m = 4` beyond `N = 16`.
- Still open.

**Root finding** — [Source](RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md) §5–§7; [Source](RESEARCH_GPU_ROOTFIND_DECOMPOSITION.md) §5
- A `u128` Montgomery field for 62 to 128 bits.
- Solinas fields for P-192/224/384/521.
- `m = 3, 4` end to end.
- Running on a real GPU (`./bench selftest` on sm_120).
- None was found done.

**Prime fast** — [Source](RESEARCH_PRIME_FAST_ECDLP.md) §7
- Two-level IC batching, Wiedemann over `u64`, an idle-host sample, and the `eccp79_rho.rs` negation fix.
- The ecbench root-finder work (10-13/14) built some prime-curve oracles, but the PKM ecbench note says "the prime group has no batched addition (each addition inverts)". So the two-level batching is still undone there. — [Source](research/pkm_tower_ecbench_20261013/README.md)

**Isogeny fibers** — [Source](RESEARCH_ISOGENY_FIBERS_ECDLP.md) §9, §17.8
- All seven §16 leads were run in round 4 (2026-10-11).
- Remaining library moves:
  - quotient-coordinate oracles for saturated bases, re-checking the 2.95× u-frame loss with a pair table on `F/M`;
  - a counted √élu with FFT, then re-measuring the exponents;
  - a pooled multi-target method;
  - a genus-2 `F_p` toolkit.
- Unproved: Lemma A, Prop B, Prop C, Thm D, Prop E and Conjecture F remain "worth writing up properly".

**Transfer-aware IC** — [Source](docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md) §5
- Protocol I-4 (IC yield across a class on prime toys, expected exact flatness) is "preregistered, cheap, expected null". No run was found.
- The C41/C41b cells were partial.

**Sparse algebra follow-up (2026-09-28/29)** — [Source](SPARSE_ALGEBRA_FOLLOWUP.md)
- Unchecked items: frozen sparse, block-structured and fill-in suites; process RSS.
- Chordal-elimination interaction-graph study: "A first general-algebra study should establish the interaction graph and observed elimination width". Not done.
- "m=83 and challenge-specific evidence remain unperformed".

**PDP/isogeny note** — [Source](RESEARCH_PDP_SPEEDUP_AND_ISOGENY.md) §3
- "on the prime side there is none [large-factor candidate] below the exponent barrier".

**Unstarted literature workstreams** — [Source](research/lit2026_speedup_program/workstreams.json)
- W1: solving-degree diagnostics, Hilbert test.
- W2: multi-graded regularity.
- W6: crossbred.
- W7: Wiedemann pruning.
- W9: Las-Vegas `F_p` F4 echelon.
- All are "not started" as of the file.

### Inferences
- The highest-leverage undone items for the general (non-binary) program are:
  1. the `k = 4` end-to-end ecbench ladder with the `D₄` solver, which is now gated open;
  2. stacking `ψ` (`Z/3`) with `D₄` or the negation fold at `k = 4` on `j = 0` curves;
  3. the PKM `N ≥ 21` degree points;
  4. the combined-system solution-count census;
  5. the multi-polynomial Coppersmith lattice.

### Gaps
- Work in uncommitted worktrees (G1 trace replay, G3 Wiedemann, G2 equivariant F4) is invisible from this clone.
- I could not check open PRs (W4, W8, W10, W11 and W12 are "PR open") or any later remote branches.

---

## Q4. Where do the leads connect across threads?

### Takeaway
The docs themselves draw several connections:
- IC with a table is equivalent to preprocessing rho;
- prime `p` and prime `n` are one obstruction (no intermediate field);
- the height parallelogram law blocks lifting;
- Frobenius savings on relations and on F4 blocks are one symmetry;
- `D₄`/FGHR torsion symmetry links the binary FGHR lever to the `m = 4` arity;
- isogeny transport is a torsion fold in disguise.

The cross-lead combinations still open are:
- `ψ`-fold × `D₄` × negation at `k = 4`;
- Faugère–Svartz-style character blocks applied to the `j = 0` `Z/3` grading;
- class-level log tables as precomputation assets.

### Cited findings
- **IC table and preprocessing rho are the same time–memory trade.** Rho achieves `S·T² ≈ 4N`. The Koblitz IC table achieves `2^{25–27}·N`. — [Source](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md) §7(a)
- **Relation collection is a k-list problem.** Wagner does not apply because `E(F_q)` has no visible small quotient. That is why residual walks and large primes cannot help. — [Source](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md) §7(b). The C4 exponent-bias test confirms the same obstruction on the Kummer tower. — [Source](research/ic_exponent_falsifiers_20261013/README.md)
- **Smoothness obstruction.** `ĥ(P+Q) + ĥ(P−Q) = 2ĥ(P) + 2ĥ(Q)` blocks xedni, every lift, and every factor base with membership degree below `B`. — [Source](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md) §7(c)
- **Prime `p` and prime `n` are one obstruction.** Every measured exponent gain (JV cover on `F_{p⁶}`, Gaudry `k = 3, 4`, Diem) lives on fields with a proper intermediate subfield. — [Source](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md) §7(d); [Source](docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md) §2 cover row
- **`D₄` symmetry of `S₅` as the FGHR lever at the first beat-rho arity (`m = 4`).** "If `C₄` fell 630× the k = 4 handover would move by ≈37 bits." Measured: 48×, about 15 bits. — [Source](RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md) §6; [Source](research/gaudry_quartic_d4_20261013/README.md)
- **Isogeny transport and torsion symmetrisation.** The quotient-isogeny factor base equals the FHJRV torsion-symmetrised construction, and kernel saturation is the image law on `E/M` transported back (Thm 7). — [Source](RESEARCH_ISOGENY_FIBERS_ECDLP.md) §0 item 3, §17.7. On `F_{p²}` `j = 0` curves the ζ-fold is free (÷3 exactly on 24 curves), but the repo does not apply it "for this plain case". — [Source](RESEARCH_ISOGENY_FIBERS_ECDLP.md) §17.4
- **GLV/endomorphism and factor-base symmetry.** Folding by any finite endomorphism group equals the order of the subgroup of `(Z/rZ)^*` generated by its eigenvalues: 12 on the `j = 0` GLS twist, 4 on `j = 1728`, 3 on the `j = 0` subfield curve. Degree-2 and degree-3 CM endomorphisms keep zero base points. Isogeny-loop endomorphisms have eigenvalue order about `r`, so they cannot fold. — [Source](research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md) status; [Source](RESEARCH_ISOGENY_LOOP_GLV.md) §4.2
- **Folds stack.** E13 shows the FGHR 2-torsion symmetry and the fold compound (S ÷ 27, then ÷ 2.5). E18 shows the `ψ` fold at `k = 4` gives `r∞` 0.138. Neither has been combined with `D₄`, and `D₄` needs rational 2-torsion, which `j = 0` prime-order-style curves lack. E18 notes those curves "pay ≈ 8 bits back through cofactor and yield". — [Source](research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md) §9.4
- **Class-level precomputation.** A factor-base log table on `E` is a log table for every curve reachable by a cheap isogeny. — [Source](docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md) §3. Class-wide multi-target rho was measured at `√L`. — [Source](RESEARCH_ISOGENY_FIBERS_ECDLP.md) §17.5
- **Conductor gap and automorphisms (crater ascent, H3).** The `drastic_prime_conductor_gap` rounds (2026-10-12) measure √3/√2 automorphism rungs for rho (`D = −3`, `−4`). They are boundary controls, not IC. — [Source](research/drastic_prime_conductor_gap_round3_20261012/PROTOCOL.md)
- **Same requirement in the binary and prime programmes.** "C₄ is exactly the requirements document's problem: an exact solver for the zeros of one fixed-shape system". — [Source](research/ic_speed_and_rho_scope_20261013/README.md) §3

### Inferences
- The one connection with exponent-adjacent value is the `k = 4` extension-field route. Its linear algebra already beats rho, and the symmetry levers (negation, `ψ`, `D₄`, Faugère–Svartz blocks) all act on the same `C₄` system. They have been measured one at a time on different curve families, never stacked on one family with an end-to-end run.
- For pure prime fields, every connection the docs draw reinforces the closure. The only threads that do not close it are the multi-polynomial lattice and the PKM degree plateau.

### Gaps
- No document quantifies whether the `ψ` (`Z/3`) and `D₄` symmetries can coexist on one curve over `F_{p⁴}` (needing `j = 0` plus rational 2-torsion plus `p ≡ 1 mod 3`), or what cofactor that forces.
