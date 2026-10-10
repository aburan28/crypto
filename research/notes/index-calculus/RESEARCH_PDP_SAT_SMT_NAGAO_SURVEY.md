# SAT and SMT solvers, WDSat and Nagao decomposition for the point decomposition problem: what is measured, what the literature adds, what is left

**Date.** 2026-10-10.  **Status.** Survey and gap analysis. Nothing new is
measured here; every number below is quoted from a frozen repository artefact
or a cited paper, and each row says which.  **Class.** Accounting (no
operation count moves).  **Conductor task.** T-1.
**Visuals.** [`research/pdp_sat_smt_survey_20261010/`](../../pdp_sat_smt_survey_20261010/)
(`oracle_ladder.svg`, `product_law.svg`, `REPORT.pdf`).

**The question.** Is there anything in SAT or SMT solving for the point
decomposition problem (PDP) on binary curves that this repository has not
exhausted and that could give a better speed-up, and is Nagao's decomposition
a better route than summation polynomials?  The prompt behind it is the
Trimoska–Ionica–Dequen WDSat line, which reports large speed-ups over Gröbner
bases on exactly these systems.

**Bottom line.**

1. **As an exponent lead the SAT route is exhausted, and the repository has
   measured why.**  On the Weil-descended Semaev system with the summands in an
   `F₂`-subspace `V` of dimension `ℓ`, every SAT, SMT and WDSat back end tried
   here spends about one conflict per candidate tuple of the factor base.  It is
   an enumerator: plain WDSat visits `2^{3ℓ}/3!` nodes to within 0.2 %, the
   native CDCL+XOR solver 0.97–1.02 conflicts per triple from `ℓ = 5` to `8`,
   CryptoMiniSat `36.7×` the pair count at `d = 8` on the ECC2K-130 curve.
   The best oracle in the tree enumerates `m − 1` summands and root-finds the
   last (`2^{(m−1)ℓ}`, with memory `2^{⌈m/2⌉ℓ}`); it beats every solver by
   `10²`–`10⁴` at the largest rung both finish (§2).  WDSat's published
   speed-ups are real and are relative to Magma F4, which on these systems is
   worse than enumeration (deciding degree `≈ ℓ/2`, §3.2); they are not
   speed-ups relative to enumeration.
2. **No solver constant matters for ECC2K-130.**  The product law
   `relations × targets × oracle = m·2^n` (`/n` for a Frobenius-stable base)
   puts any oracle that enumerates sub-tuples `2^{+64…72}` above rho; an oracle
   would have to beat exhaustive search over its own candidate set by
   `2^{−(70.19 + log₂ m − log₂ 131)}` (§1).  A SAT solver cannot do that, and
   neither can an SMT solver, whose Boolean and bit-vector theories bit-blast to
   the same CDCL core (measured: cvc5 and Z3 within `0.44–2.6×` of the native
   solver at `n = 13`, §2).
3. **Nagao's decomposition is not a better route at prime extension degree.**
   Nagao's gain is a *composite-degree* phenomenon: with a subfield factor base
   the condition "the norm polynomial has coefficients in `F_q`" is a quadratic
   Weil-restriction system with no summation polynomial to expand.  At prime
   `n = 131` there is no subfield, the condition becomes "roots in `V`", i.e.
   divisibility by the degree-`2^ℓ` subspace polynomial, and the repository's
   Riemann–Roch encodings reduce to `S′₄` with the linear part hidden.  Measured
   on ECC2K-130: the RR solver resolves `32/32` matched slots where the
   Semaev controls resolve `0/32`, and brute-force pair enumeration then costs
   `0.494×` the RR solver's counted field operations (§5).
4. **What is not exhausted is engineering, and it is listed with falsifiers in
   §6.**  The cheapest items reuse harnesses already in the tree: a matched run
   of the newer SATIC ANF solver and an audit of WDSat's incomplete `-x` mode on
   the frozen 60-input regression; a PB-XOR solver on the one factor base where
   SAT has a structural edge over Gröbner (Hamming-weight in a normal basis);
   ANF preprocessing (Bosphorus) in front of every SAT arm; and completing the
   SMT ladder to `n = 23` so the open SMT study can be classed.  None is
   expected to move a ratio to the floor; each is a stage diagnostic.

---

## 1. Boundary, stated before anything else

Setting: `E/F_{2^n}`, `n` prime (ECC2K-130: `K₀: y² + xy = x³ + 1`, `n = 131`,
`#E = 4r`, `r ≈ 2^129`), factor base `F = {P : x(P) ∈ V}` with `V` an
`F₂`-subspace of dimension `ℓ`, or a Frobenius-stable set (Hamming-weight in a
normal basis), `m` summands.

**Reference.** Pollard rho with the `⟨−1⟩ × ⟨π⟩` automorphism discount:
`2^{60.81}` group operations on ECC2K-130
([`RESEARCH_ECC2K130_DECOMPOSITION.md`](../ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md) §5.2).

**Floor for the whole method (the product law, same note §5.1).**

```text
    relations · targets/relation · oracle/target
  = 2^ℓ · 2^n / C(2^ℓ, m) · C(2^ℓ, m−1)  ≈  m · 2^n          (independent of ℓ)
```

divided by `n` when the base is Frobenius-stable.  Turned around (§5.3 there):
an oracle that reaches rho must beat exhaustive search over its own candidate
set by `2^{−(70.19 + log₂ m)}`, or `2^{−(63.2 + log₂ m)}` with the Frobenius
collapse.  No tuning of `m`, `ℓ` or the solver changes that number; it is the
identity.  The free-oracle floor (charge nothing per decision) sits below rho
from `m = 4` on (`2^{56.40}`), so the method is not structurally impossible —
it is the oracle that protects the curve.

**Floor for the per-target oracle (this note's unit).**  The trivial oracle
enumerates `m − 1` summands and decides the last by root-finding the
univariate `S_{m+1}` against `V`; with the `S_m` symmetry that is

```text
    T_enum(m, ℓ) = 2^{(m−1)ℓ} / (m−1)!    candidate sub-tuples per target,
```

each costing `O(ℓ)` field operations through `gcd(q, L_V mod q)`
([`RESEARCH_SEMAEV_DECOMPOSITION.md`](RESEARCH_SEMAEV_DECOMPOSITION.md)).  With
`2^{⌈m/2⌉ℓ}` memory a pair table (meet in the middle) lowers it to
`2^{⌊m/2⌋ℓ}` lookups.  A solver that reports more than `T_enum` conflicts or
decisions per target is enumerating with overhead; one that reports fewer is
doing algebra the enumerator does not.  **Unit:** the ratio of a solver's
per-target work to `T_enum`, in the solver's own counter where the artefact
reports one and in counted field operations where it does not.  Wall time is
quoted beside it as a practicality note only.

**Falsification target for the "exhausted" verdict.**  The verdict is wrong if
any solver, encoding or preprocessing in §6 produces, on uniform (not planted)
targets at `m = 3`, a per-target count below `0.5 · T_enum` at two consecutive
rungs `ℓ ≥ 8`, with every verdict checked against the exhaustive oracle.  It is
confirmed if every arm stays at or above `T_enum` out to the largest rung it
decides.

---

## 2. The oracle ladder: one table

Every row is a frozen artefact; the right-hand columns are the ratio to the
enumeration floor of §1 on the largest rung the arm decides, and the verdict
the artefact records.  `m = 3` unless stated.  "Planted" rows measure
first-solution cost on a target known to decompose; "uniform" rows measure
exhaustion cost, which is the oracle's real job.

| arm | encoding → solver | base | largest rung decided | work per target vs `T_enum` | correctness | verdict in source | source |
|:--|:--|:--|:--|:--|:--|:--|:--|
| pairs-and-solve (reference) | `S₄` quartic in last variable, `gcd(q, L_V mod q)` | subspace | `ℓ = 12`, `n = 36` (11.5 G triples in 15 s) | `1.0×` by construction; `8 µs → 10 ns` per candidate after one-word field | exhaustive, cross-checked | the oracle the ledger prices | [`RESEARCH_SEMAEV_DECOMPOSITION.md`](RESEARCH_SEMAEV_DECOMPOSITION.md), [`RESEARCH_SAT_SEMAEV.md`](RESEARCH_SAT_SEMAEV.md) §"What was done instead" |
| direct MITM / pair table | group law, no polynomial | subspace, Frobenius orbits | `n = 59, ℓ = 9` (80/80 verdicts, 3.86 core-s median) | `2^{−ℓ}×` time at `2^{2ℓ}` memory | 160/160 | "lowest median core time in every cell" | `research/sat_factor_base_review_20260908/…/STAGE13_RESULTS.md`, `STAGE26_32_PHASE_B_RESULTS.md` |
| matrix-F4 with splitting | Weil-descended `S₃` chain / `S₄`, Macaulay rows | subspace | `n = 21, m = 3` refuted in 50 s; `n = 31, m = 3` refutation 145 s | `8.8·10⁴` GAE at `n = 23, m = 2` against `0.011–25` GAE for the pair table | agrees with enumeration target by target | best *algebraic* oracle; deciding degree `≈ ℓ/2` | [`RESEARCH_IC_BOUNDARY_LEDGER.md`](RESEARCH_IC_BOUNDARY_LEDGER.md) §5, [`RESEARCH_KOBLITZ_SCALING_TARGET.md`](../ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md) |
| native CDCL + Gauss–Jordan XOR (`sat.rs`) | symmetrised `S₄` → native parity rows, free-variable branching, lex symmetry breaking | subspace | `n = 24, ℓ = 8` uniform rejection: 2,877,905 conflicts on 2,829,056 triples | `≈ 1.02` conflicts per *triple*, i.e. `≈ 2^ℓ/3 × T_enum`; wall `62×` brute force and widening (`5×, 8×, 18×, 62×` over `ℓ = 5…8`) | all verdicts match enumeration | "SAT is a strictly worse way to enumerate the factor base" | [`RESEARCH_SAT_SEMAEV.md`](RESEARCH_SAT_SEMAEV.md) §"Does it scale? No." |
| native CDCL + XOR, IC pipeline (`sat_decompose`) | chained `S₃` + Kosters–Yeo trace row + degree-2 Macaulay rows | invariant subspace | `n = 15, m = 3` refuted in 25.3 s (F4: 75 ms); `n = 21, m = 3` no verdict in 31 min | `1.3·10⁷` GAE at `n = 23, m = 2` (`150×` F4) | agrees where it answers | fails 2/16 targets at 20 unknowns, 10/16 at 22 | [`RESEARCH_KOBLITZ_SCALING_TARGET.md`](../ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md) §240–305, ledger §17.11 |
| WDSat `61c6ff3`, plain and symmetry-broken | Trimoska ANF (`.anf`), XORGAUSS off | subspace | frozen corpus `n ∈ {15, 17, 19}`, `ℓ ∈ {5, 6}`, 60 inputs | symmetry on/off: `1,280,281 / 6,246,481` summed median conflicts (`4.88×`, the `3!` of `S₃`) | 30 SAT, 29 UNSAT, 1 algebraic-only rejection, all certified | regression target, "no solver algorithm was changed" | `research/index_calculus_baseline_20260914/regression/RESULTS.md` |
| WDSat, Semaev ANF, `n = 59, ℓ = 9` | same | subspace | none: outer watchdog at 120 s on all 40 targets | `264×` direct-MITM wall, no terminal | — | inconclusive | `…/continuation-05-sota-gates/GATE_STATUS.md` Stage 193 |
| WDSat on the real curve (`solver_16`) | `S′₄`, polynomial basis, `V = {deg x < d}` | subspace, ECC2K-130 | `d = 8` | plain: `2^{3d}/3!` decisions to 0.2 % (`= 2^d/3 × T_enum`); `-x` returned UNSAT on a decomposable target in 1 of 16 cells | `-x` **incomplete** at `n = 131` | "plain WDSat is brute force" | [`RESEARCH_ECC2K130_RR_SOLVER_PANEL.md`](../ecc2k130/RESEARCH_ECC2K130_RR_SOLVER_PANEL.md) §11 |
| CryptoMiniSat 5.14 on the real curve (`solver_16`) | `S′₄` → CNF + native XOR | subspace, ECC2K-130 | `d = 8` (`d = 9, 10` censored at 300 s) | uniform exhaustion: `1` conflict at `d ≤ 6` (the linear certificate), `839` at `d = 7` (`0.31×` pairs), `312,231` at `d = 8` (`36.7×` pairs) | every witness verified | "super-quadratic in `|F|` from `d = 7`" | same, §11 |
| CryptoMiniSat on the ONB Karatsuba circuit | gate-level CNF + XOR, Sinz cardinality for the weight bound | **Hamming-weight `w = 3`** in the normal basis (Frobenius-stable at prime `n`) | `n = 11` (median 2.64 s); `n = 23` nothing in 25 min | `≈ 2^{2.3}` per field bit; encoding untuned | planted | the only SAT-only base; "the shape is already clear enough" | `ecc2k130/README.md` §Measured |
| Kissat 4.0.4 / CaDiCaL 3.0.1 | Tseitin CNF of the five-summand `S₃` chain | window slots, `n = 83` | none with the target free (120 s); half-pinned planted instance 157,941 conflicts in 30 s, no model | — | planted control only | inconclusive | `research/f6_n83_factored_s3_sat_20261005/RESULT.md` |
| cvc5 1.4.2 / Z3 (`QF_UF`, `QF_BV`) | identical rows to the native arm, SMT-LIB 2 | subspace | `n = 13, ℓ = 5` (20/20 each); `n = 19` every arm at the 120 s cap | Z3-bool `0.76×`, Z3-bv1 `1.82×` native conflicts; cvc5-bool `0.44×`, Z3 `2.2–2.6×` native wall | 100/100 | unclassed: parity reasoning is not the gap; `n = 23` not run | `research/pdp_smt_oracle_20261008/README.md` |
| Crossbred (`D = 2, 4`) | Macaulay extraction + `2^k` linear solves | subspace | `n = 9, m = 3` | `xb/F4 = 0.696` wall and rising; **zero filters at any `(D, k)`**; pipeline `1800×` behind enumeration | agree | "enumeration in disguise" | [`RESEARCH_ECC2K130_CROSSBRED.md`](../ecc2k130/RESEARCH_ECC2K130_CROSSBRED.md) |
| Riemann–Roch / Nagao (`quadratic-image`) | function coefficients `(a, b)` as unknowns, norm-polynomial coefficient match, support test `H | L_V` | subspace, ECC2K-130 | `d = 6, 7` (32/32 matched slots; 28 verified 3-summand relations on the real curve) | `Θ(|F|²)` candidates; pair enumeration `0.494×` its counted field ops | every witness re-derived in the group | "the null object wins" | [`RESEARCH_ECC2K130_RR_SOLVER_PANEL.md`](../ecc2k130/RESEARCH_ECC2K130_RR_SOLVER_PANEL.md) §7–8 |
| Riemann–Roch under SAT (WDSat / CryptoMiniSat) | RR-norm and RR-incidence ANF | subspace | RR-norm `n = 9` (8/8); RR-incidence does not fit WDSat (SIGSEGV above `≈ 23 M` table entries) | XORGAUSS `79,151 → 62` conflicts (`1,277×`); RR costs `3.5·10⁴×` the `S′₄` conflicts at `d = 5` | witnesses verified | "the RR encoding is `S′₄` with the linear part hidden" | [`RESEARCH_ECC2K130_WDSAT.md`](../ecc2k130/RESEARCH_ECC2K130_WDSAT.md), panel §11 |
| Hamming-ideal Gröbner (La Scala–Marchesin–Tiwari lifts) | `S₃` + C / FC / QFC Hamming ideals, `MultiSolve` | Hamming-weight, normal basis, `m = 2` | `n = 19` | calls per target grow as `|F_w|^{1.9–2.5}` against a target exponent `0.8`; tame depth within 1–3 of a full summand | all targets correct | "the exhaustive oracle with a truncated Gröbner computation attached" | [`RESEARCH_HAMMING_IDEAL_PDP.md`](RESEARCH_HAMMING_IDEAL_PDP.md) |
| one-hot Frobenius-phase SAT | orbit representative + shift variables | Frobenius orbits | `n = 13` (`≈ 50 s` median); `n = 19` fails at 1,200 s | — | — | no gain | `research/frobenius_quotient_sat_20260930/` |

Three readings.

- **Every arm above the reference line is at or above `T_enum`.**  The two
  that look below it — CryptoMiniSat at `d ≤ 7` on the real curve, and
  matrix-F4 at `m = 2` — are the linear NO-certificate of the `S₄` value set
  (available only for `d ≤ 6`, saturating at `d = 7`; panel §9) and the
  Kosters–Yeo trace fall at `m = 2` (first fall degree 2, a theorem), which
  the enumerator also gets for free.
- **The SAT arms sit a factor `2^ℓ/m` *above* `T_enum`, not at it.**  One
  conflict per *full* tuple means the solver assigns all `m` summands before
  the descended equations bind; the enumerator assigns `m − 1` and root-finds.
  Only WDSat's `-x` mode might close that gap, and the one artefact that tests
  it at full size records a wrong UNSAT (§6, item 1).
- **Gröbner is the arm WDSat beats, and it is below enumeration too.**  The
  deciding degree of F4 on the one-summand-fixed system runs `5, 6, 6, 7, 7`
  over `ℓ = 3…7` (slope `0.500`), so elimination alone is `≈ 2^{3.6ℓ}` against
  `2^{2ℓ}` for the enumerator (decomposition note §5.4).  That is the baseline
  the WDSat papers improve on.

---

## 3. Why every solver is enumeration in disguise

This is the structural reason the §6 list is engineering and not a new lead.
It is a derivation from measured fall degrees, not a theorem.

### 3.1 Where the equations bind

Weil descent of `S_{m+1}(x₁, …, x_m, x_R)` with `x_i ∈ V` gives `n` Boolean
equations in `mℓ` unknowns.  In characteristic 2 the even powers `x^{2^k}` are
`F₂`-linear in the coordinates, so the descended system is quadratic for
`m = 2` and cubic when `S₃` is chained (`m ≥ 3`); the symmetrised `S₄` is
quadratic in the elementary symmetric unknowns with a cubic correspondence
system behind it.  With `m − 1` summands fixed the residual is a univariate
polynomial of degree `2^{m−1}` in the last `x`, whose roots in `V` are found by
`gcd` with `L_V` in `O(ℓ)` field operations: that is the enumerator's leaf.
With only `m − 2` fixed, the residual is the `m = 2` problem for the target
`R − P₁ − ⋯`: `n` quadratic equations in `2ℓ` unknowns.

### 3.2 The `m = 2` residual is not cheaper than `2^ℓ` by any known method

- **Enumeration:** `2^ℓ` candidates for `x₂`, each a group subtraction and a
  membership test.  This is what pairs-and-solve and the pair table do.
- **Gröbner/XL on the residual:** the first fall degree is 2 (Kosters–Yeo
  Cor. 4.11, the trace row), but that fall determines one bit; the measured
  deciding degree grows as `≈ ℓ/2` and the only proved upper bounds on the last
  fall degree are linear in `n`
  ([`RESEARCH_FALL_DEGREE_BOUNDS.md`](RESEARCH_FALL_DEGREE_BOUNDS.md) §2.3).
  Kosters–Yeo also show that "`d_reg ≈ d_ff`" for the split system would decide
  an NP-complete problem in polynomial time, so the subexponential claims that
  rest on it (Petit–Quisquater, Semaev 2015, Kousidis–Wiemers) are unproved
  heuristics, and the toy data here go the other way.
- **Generic MQ solvers on the residual:** `N = 2ℓ` unknowns and `n ≈ 1.5 N`
  equations.  The best generic exponents over `F₂` (Dinur's polynomial method
  `2^{0.6943 N}` for `M = N`, Crossbred and FES around `2^{0.8 N}`, Lokshtanov
  et al. `2^{0.8765 N}`) give at best `2^{1.39 ℓ}`, above the `2^ℓ`
  enumerator.  Crossbred measured on these systems found no filters at any
  `(D, k)` and a margin over F4 that shrinks with `n`.
- **CDCL/DPLL on the residual:** the quadratic rows do not propagate until
  nearly all of `x₂`'s bits are assigned, which is the "tame depth within 1–3
  of a full summand" the Hamming-ideal thread measured and the one conflict per
  leaf the SAT threads measured.

So the leaf cost of any tree search on this system is bounded below by the
`m = 2` residual, and the `m = 2` residual is bounded below by `2^ℓ` on all
measured and all proved evidence.  A solver can only lose constants against
`T_enum`; the table says the measured ones lose `2^ℓ/m`.

### 3.3 What a solver *can* do that the enumerator cannot

Two things, both already in the tree:

- **Cardinality constraints.**  A Hamming-weight factor base in a normal basis
  is Frobenius-stable at prime `n` and has a popcount membership test; its
  membership ideal has degree `≥ 3` and `O(n log n)` auxiliaries (Hamming
  ideals), which Gröbner handles badly and a SAT solver handles natively.  This
  is the one base where the solver family matters.  The product law is
  unchanged by the base, so the gain is confined to the oracle constant.
- **UNSAT certificates.**  CaDiCaL with DRAT-trim produced externally checked
  refutations of five `n = 13` branch queries
  (`research/notes/ecc2k130/n13_oaware_unsat_proof_20260925/`).  That is a
  correctness instrument, not a speed.

---

## 4. The literature, 2019–2026


### 4.1 SAT and SMT solvers on the descended systems

Sources were read in full where the sandbox could fetch them (Trimoska's
thesis, the CP 2020 and SAT 2021 papers, GGMP 2020, Ozdemir et al. CAV 2023,
Hader–Ozdemir 2024, Danner–Kreuzer JAIR 2026); the items that could only be
read as abstracts are marked so.

**WDSat (Trimoska–Ionica–Dequen; AFRICACRYPT 2020, CP 2020, SAT 2021; thesis
2021).**  The solver is a DPLL (no clause learning) with three propagation
modules: CNF unit propagation, parity propagation (XORSET) and Gaussian
elimination on the parity system with substitution tracking read straight from
ANF (XORGAUSS, "XG-ext").  It branches only on the `m·ℓ` coordinate bits, and
breaks the `S_m` symmetry inside the search by refusing the branch that would
violate `x₁ ≤ ⋯ ≤ x_m` (gain exactly `m!`, no overhead).  Two facts from the
thesis decide this note's question:

- **Rigorous bound.**  Conflicts per PDP are at most `2^{mℓ}/m!` (thesis
  eq. 7.14), and the relation-search phase for `S′₄` is `Õ(2^{n+ℓ})` (Thm 7.5.1),
  against rho's `2^{n/2}`.  The authors say so themselves.
- **The bound is attained.**  Thesis Table 7.6, UNSAT instances, `m = 3`:
  `2.80·10⁶` conflicts at `ℓ = 8` against `2^{24}/6 = 2.796·10⁶`; `22.39·10⁶`
  at `ℓ = 9` against `22.37·10⁶`; `179.0·10⁶` at `ℓ = 10`; `1.432·10⁹` at
  `ℓ = 11` against `1.431·10⁹`.  WDSat prunes nothing above the last branching
  level ("the system generally does not become linear until the
  second-to-last branching").  The per-candidate cost at `ℓ = 11, n = 89` is
  `≈ 176,000 s / 1.43·10⁹ ≈ 120 µs`, four orders of magnitude above this
  repository's one-word enumerator.

Measured on Koblitz curves (one thread, Xeon E5-2640): Magma F4 takes 207 s at
`ℓ = 6, n = 17` (3.6 GB), 3,855 s at `ℓ = 7, n = 19` (38 GB) and runs out of
200 GB at `ℓ = 8`; WDSat with symmetry breaking takes 0.22–0.61 s, 2.2–6.9 s
and 30–86 s on the same rungs, and reaches `ℓ = 11, n = 59…89` in 1.5–2 days
per instance in 19 MB.  The headline "up to 1700× over Gröbner" is this
comparison at `ℓ = 6, 7`.  Generic CDCL solvers (MiniSat, Glucose, CaDiCaL,
CryptoMiniSat with the authors' branching patch) stop at `ℓ = 8`.  For `S₃`
(`m = 2`) the XG-ext module with a minimum-vertex-cover branching order
solves `n = 41, ℓ = 20` in 4–14 s; the MVC has `ℓ = 20` elements, i.e. the
solver guesses one summand and solves a linear system, which is the trivial
`m = 2` oracle.  Gaussian elimination *hurts* on `S′₄` (the monomial graph is
complete, so the MVC is every variable).  The thesis' own reality check: at
`n = 59` their parallel collision search takes 0.8 h single-threaded, while
`2⁹` successful decompositions at `ℓ = 9` would take more than 86 h.

**Follow-ups.**  Nothing after 2021 re-ran WDSat on PDP or passed its `ℓ = 11`
instances (Semantic Scholar citation graph, checked 2026-10).  The direct
successor is Blomme, Cherif, Ionica, Dequen, *ANF-based satisfiability for
Weil-descent cryptographic attacks* (CoDIT 2025, HAL hal-05176414, abstract
only: lazy watched structures on the ANF itself, no CNF-XOR translation); its
benchmark sizes are unverified.  Andraschko–Danner–Kreuzer (2-XNF, Math.
Comput. Sci. 2024) and Danner–Kreuzer (CDXCL / Xorricane, JAIR 2026) add
clause learning to XOR-native solving with a proof system (XLIN) polynomially
equivalent to resolution with parity, which is exponentially stronger than
resolution; it was benchmarked on Ascon, Bivium, CTC2 and random MQ, not on
PDP, and loses to Kissat on structured instances.  CryptoMiniSat ≥ 5.8
(Soos–Gocht–Meel, CAV 2020) runs bit-packed Gauss–Jordan at every decision
level; Trimoska's CMS numbers predate it, and this repository's CMS 5.14 runs
are the post-5.8 data point (§2: `36.7×` pairs at `d = 8`).  Bosphorus
(Choo–Soos–Chai–Meel, DATE 2019; ANF/CNF simplification with XL, ElimLin and
learnt facts) is untested on Weil-descent systems anywhere.  Cube-and-conquer
in cryptanalysis exists only for hash inversion (Zaikin, JAIR 2022).
ML-guided SAT benchmarks (SAT4CryptoBench, NeurIPS 2025) carry no ECC
instances.  In-memory XOR–CNF hardware (Im et al., Nat. Commun. 2026) and
GPU clause evaluation (GaloisSAT, 2026) report simulated or generic-instance
gains and cite WDSat without running PDP.  No hybrid Gröbner+SAT or XL+SAT
paper on PDP appeared 2021–2026.

**What the SAT literature never does** is compare against enumeration.  Every
paper's baseline is F4; the MQ literature's baseline is FES (Bouillaguet et
al., CHES 2010: Gröbner `≈ 10⁶×` slower than libFES on `F₂`; Joux–Vitse
Crossbred 2017; Dinur 2021 `2^{0.6943 N}`; Vidal–Delaplace–Ionica, IACR CiC
2025 on Crossbred).  This repository made the comparison (§2) and the answer
is the one the thesis' Table 7.6 already implies.

**SMT.**  cvc5's finite-field theory (Ozdemir–Kremer–Tinelli–Barrett, CAV
2023; split Gröbner bases, CAV 2024) is a Gröbner-basis decision procedure
for **prime** fields; the SMT-LIB theory proposal (Hader–Ozdemir, SMT 2024)
defines extension sorts `(_ FiniteField p n)` but no solver implements them
(cvc5 and Yices support `QF_FFA` for prime fields only).  Bit-vector
encodings of `F_{2^n}` arithmetic are a Weil descent with a Tseitin variable
per multiplication gate, 10–20× larger than the descended ANF (thesis Table
7.1), and bit-blast to the same CDCL cores (cvc5 → CaDiCaL, Z3 → its own);
Ozdemir et al. report bit-vector solvers already inferior on fields of
`≥ 2^{40}` elements.  No peer-reviewed paper evaluates any SMT solver on
ECDLP, PDP, `GF(2^n)` multiplication constraints or Semaev polynomials.  This
repository's own SMT study (§2) is therefore the only data: cvc5 and Z3 land
within `0.44–2.6×` of the native solver at `n = 13` and hit the cap with it at
`n = 19`.  A native `GF(2^n)` theory would be Gröbner bases again.  Verdict:
no unexploited SMT avenue.

**Koblitz specifics (Galbraith–Gebregiyorgis–Murphy–Petit, SAC 2020).**
The `m!` symmetry can be put in the *model* with `π`-shifted factor bases
`F_i = π^{i−1}(F)` and recovered as `R = P′₁ + λP′₂ + ⋯`; it is the same `m!`
WDSat takes in the search, so the two do not stack, but GGMP's form works
under any solver (F4, FES).  Frobenius-invariant subspaces need a factor of
`x^n − 1`; at `n = 131` it is `(x − 1)·Φ₁₃₀` with `ord₁₃₁(2) = 130`, so the
only invariant subspaces are `{0}`, `F₂`, the trace-zero hyperplane and the
field (this repository's §6 of the decomposition note, independently).  GGMP's
conclusion: "our improvements do not lead to index calculus faster than
Pollard rho on instances used in practice."

**Records.**  End-to-end binary ECDLP by Semaev index calculus: `n = 19`
(Huang–Petit–Shinohara–Takagi 2013); relation collection only at `n = 24`
(WDSat, `ℓ = 8`, 18 h).  PDP by Gröbner: `n = 53, ℓ = 6` symmetrised `S₄`
(HPST, 35 s); `m = 2` at `n = 45` needed 126 GB with `D_reg = 5`
(Kosters–Yeo).  PDP by SAT: `ℓ = 11, n = 89` (`m = 3`, WDSat); `ℓ = 20,
n = 41` (`m = 2`).  This repository's panels sit inside those: `n = 59, ℓ = 9`
by MITM in 3.9 core-s where WDSat and CMS time out, and `ℓ = 12` by
pairs-and-solve in 15 s.  Every binary ECDLP record remains a rho record;
ECC2K-130 is unsolved and "no known index calculus attack applies" (Bailey et
al. 2009) is still the state of the art.

### 4.2 Nagao's decomposition and the other un-eliminations of `S_{m+1}`

**Nagao (ePrint 2007/112; ANTS 2010).**  For `Jac(C)/F_{q^n}` with factor base
"`x ∈ F_q`", decompose a divisor `D` by finding `f ∈ L(ng·O − D)` whose norm
polynomial `F(x) = f(x, y) f(x, −y)/u(x)` has **all coefficients in `F_q`**:
the coefficients are quadratic in the unknown coefficients of `f`, and the
`F_q`-rationality condition Weil-restricts to `(n − 1)·ng` quadratic equations
in `(n − 1)·ng` unknowns over `F_q`, generically zero-dimensional of degree
`2^{(n−1)ng}`.  For `g = 1` that is `n(n − 1)` quadrics in `n(n − 1)` unknowns
against Semaev's `n` equations of degree `2^{n−1}` in `n` symmetric unknowns —
**the same Bézout degree `2^{n(n−1)}`, un-eliminated**.  No first-fall or
`D_reg` data for the elliptic case is published.  Nagao 2013/549 claims a
subexponential bound under a first-fall-degree assumption of the kind
Kosters–Yeo and Huang–Kosters–Yeo later refuted experimentally.

**Joux–Vitse (J. Cryptology 2013; EUROCRYPT 2012).**  The people who built
both: "This approach is less efficient than Semaev's in the elliptic case,
but is the simplest otherwise [higher genus]"; Nagao-style tests are feasible
for `(n, g) ∈ {(2, 2), (2, 3), (3, 2)}` and cost 22 ms per core on their
genus-3 cover over `F_{p²}`, where their own sieve on `F(x)` is "about 960
times faster".  Their elliptic record (`E(F_{p^5})`, 130-bit, oracle-assisted
static DH) uses Semaev with the `n − 1` trick (`n` equations, `n − 1`
unknowns, degree `2^{n−2}`) and F4 trace replay, not Nagao.
Faugère–Huot–Joux–Renault–Vitse (EUROCRYPT 2014) then take `S₆` symmetrised
under 2-torsion to 10 s per decomposition on the Oakley curve over
`F_{2^{155}} = F_{2^{31}}^5` ("too slow to seriously threaten the DLP").
Galbraith–Gaudry (DCC 2016): Nagao's approach "is also valid for elliptic
curves … nobody managed to turn this into an algorithm that is faster than an
approach based on summation polynomials."

**The prime-degree obstruction.**  Nagao's and Gaudry's formulations both
presuppose a factor base over a proper subfield `F_q`, `q` large.  Over
`F_{2^131}` the only subfield is `F₂`, so the rationality condition has
nothing to act on.  The replacement, "all roots of `F` lie in the subspace
`V`", is divisibility by the linearised polynomial `L_V` of degree `2^ℓ`,
which is not coefficient-linear; writing it out either re-introduces a root
variable per summand (and lands on the FPPR / HPST symmetrised system) or
keeps the coefficients as unknowns and tests `H | L_V` by iterated squaring —
which is this repository's function-first encoding (`nagaoannihilator`,
`solver_02…09`), exhaustive in `q(q − 1)` coefficient pairs and `Θ(|F|²)`
in its best image-space form (§2).  Galbraith (ECC 2015) states the real
question: "Semaev's approach is to minimize number of variables at the
expense of exponential degree.  Other choices of coordinates can lead to lower
degree but more variables.  Problem: determine the optimal tradeoff."  The
published points on that axis are FPPR's symmetrised descent, HPST's mixed
symmetric/non-symmetric variables (`D_reg` 6–7 → 3–4; whole ECDLP at
`n = 19, m = 3` in 9,963 s), Semaev's 2015 `S₃` chains (cubic Boolean system,
first fall 4 proved, `D_reg = 5` already at `n = 45` by Kosters–Yeo), Karabina's
auxiliary-variable chains (`m = 3: n = 26` in 1,646 s; `m = 4, 5: n = 19`; an
explicitly exponential `2^{34}·2^{n/5}·n^{27}` extrapolation), and
Galbraith–Gebregiyorgis' binary-Edwards `T₂`-invariant coordinate
(`D_reg ≤ 4` at `m = 3` out to `n = 97`, `ℓ ≤ 7`; `m = 4` only `ℓ ≤ 4` under
Gröbner, `ℓ ≤ 7` under MiniSat at 4 % success).  This repository has measured
two of those points (chained `S₃`: cubic, slower in every engine, §2;
symmetrised `S₄` under SAT and F4) and the `T₂` symmetrisation on `K₀`.

**Koblitz structure.**  Galbraith–Granger–Merz–Petit (SAC 2020): Frobenius
`π = [λ]` breaks the `m!` symmetry in the model and a Galois-invariant
subspace `V_f` (roots of a linearised polynomial from an irreducible factor
`f | x^n − 1`) divides relation collection by `n` and linear algebra by `n²`.
The construction needs `deg f = ord_n(2)` useful; checked here by direct
computation: `ord₁₂₇(2) = 7`, `ord₁₃₁(2) = 130`, `ord₁₆₃(2) = 162`,
`ord₂₃₃(2) = 29`, `ord₂₈₃(2) = 94`, `ord₅₇₁(2) = 114`.  So at `n = 131` the
invariant subspaces are `{0}`, `F₂`, the trace-zero hyperplane and the field
(the decomposition note's §6 and `ecc2k130/README.md` say the same), while
K-283 and K-571 happen to have invariant subspaces of dimension `≈ n/3` and
`≈ n/5` — a structural curiosity, far outside any reachable PDP size, noted
and not pursued.  Tsakou–Ionica (ePrint 2021/721) extend endomorphism-reduced
bases to composite-degree Jacobians.  Nothing on prime-degree binary PDP
appeared 2022–2026 in ePrint or arXiv; Courtois' `2^{n/3}` splitting claims
(ePrint 2016/003) have no independent check and the author's own text says
they may "violate the generic group model".

**Refuted and unproved.**  Petit–Quisquater (`2^{c n^{2/3} log n}`), Semaev
2015 (FIPS `n = 409, 571` "theoretically broken") and Nagao 2013/549 all rest
on `D_reg ≈ d_ff`; the measured `D_reg` for `S₃` with a random subspace runs
`3, 4, 4, ≥ 5` at `n = 16, 20, 30, 40` (Kosters–Yeo) and `5` at `n = 45` for
Semaev's own system (126 GB).  FPPR's own extrapolation for `n = 131, m = 2`
is `2^{74.5}` Gröbner operations, "still worse than exhaustive search"
(Renault's EUROCRYPT 2012 slides).  Galbraith (ECC 2015): prime `n > 160`
"seems to be completely immune to point decomposition attacks, even just
trying to get a cube-root algorithm (case `m = 4`)".  Yokoyama–Yasuda–
Takahashi–Kogure (JMC 2020) prove the naive (non-descended) method cannot
beat generic algorithms.  ECC2K-130 is unsolved: the project page still shows
the 2009 `2^{60.9}`-iteration estimate and no completion, Certicom lists the
131-bit challenges as open, and the largest binary ECDLP solved by any method
is 113 bits (sect113r2, FPGA rho, 2016).

---

## 5. Nagao's decomposition against summation polynomials

**Verdict: no, and the reasons are structural rather than a matter of
engineering.**

1. **Same object, un-eliminated.**  `S_{m+1}(x₁, …, x_m, x_R) = 0` is the
   condition that some `f ∈ L((m + 1)·O)` vanishes at `P₁, …, P_m, −R`.  Nagao
   keeps the coefficients of `f` as unknowns (quadratic, many variables);
   Semaev eliminates them (one polynomial, degree `2^{m−1}` per variable, few
   variables).  The Bézout degrees agree.  Which point on that axis is best is
   Galbraith's open tradeoff problem, and every published point (§4.2) and
   every point measured here (§2) is enumeration-bound at the sizes reached.
2. **The gain Nagao offers needs a subfield.**  "Coefficients in `F_q`" is a
   linear condition; "roots in `V`" is not.  At `n = 131` the method has no
   handle, and this repository's attempt to give it one — the function-first
   support test `H | L_V` and its image-space successor — is the best
   evidence: it beats the Semaev encodings under the same SAT solver on the
   real curve (`32/32` against `0/32` at `d = 6`) and then loses to pair
   enumeration by `2.02×` in counted field operations, because it visits
   `≈ |F|²` candidates where the double loop visits `|F|²/2`.
3. **Where Nagao-type relations do apply, the repository already has them.**
   The cover-and-decomposition ledger on `F_{p^6}` (genus-3 cover, Nagao
   relations, Joux–Vitse sieve reproduced at `≈ 10³×` per relation) is the
   composite-degree case; its crossover with rho is measured at `ℓ ≈ 2^{50}`
   on the weak class only.  None of it transfers to prime `n`.
4. **The arity and base-shape limits are intrinsic.**  The RR chart fixes
   `m = 3` and needs an `F₂`-subspace base; a Hamming-weight base or more
   summands are outside the encoding, not merely untested.

What this leaves, from the solver panel's own §11: an encoding that is
linear after a preprocessing cheaper than `|F|²` (a statement about the
group law on `V`, not about Nagao), and the non-subspace base.  Neither is a
Nagao direction.

---

## 6. What is not exhausted: ranked, with falsifiers

Each item names the harness it reuses, the measurement, the falsifier, the
expected class, and the cost.  None is expected to be an advance.  Items are
ordered by (value of the answer) / (cost to get it).  Python is not an option
for any of them (`AGENTS.md`, "Implementation language"); the ecc2k130 codegen
path would need its native replacement first.

| # | item | reuses | measurement and falsifier | expected class | cost |
|--:|:--|:--|:--|:--|:--|
| 1 | **XOR-learning CDCL (Xorricane / CDXCL) on the frozen `S′₄` corpus, `ℓ = 8…11`** — the one SAT route with a proof system stronger than DPLL + Gaussian elimination (XLIN ≡ RES(⊕)), and the open question Trimoska's thesis §9.1 names (WDSat + learning) | upstream benchmark generator (`EC-Index-Calculus-Benchmarks`), `semaev_corpus.rs` parameters, ANF → XNF converter to write (Rust) | conflicts per **uniform** UNSAT target against `2^{3ℓ}/3!` (WDSat's attained bound) and against `T_enum`. *Falsifier for the §1 verdict:* `< 0.5 · T_enum` at two consecutive `ℓ ≥ 8`. *Abandon if* `≥ 0.5 · 2^{3ℓ}/3!` at `ℓ = 9`, i.e. no pruning above the last level | engineering at best; the §3 argument predicts "closed" | 1 week |
| 2 | **WDSat `-x` completeness audit, then a matched SATIC run** on the 60-input regression | `research/index_calculus_baseline_20260914/regression/` (runner is Python and must be replaced natively first, `AGENTS.md`), `wdsat_oracle.rs` | reproduce the wrong UNSAT (`solver_16`, 1/16 cells at `n = 131`) on an instance under the 23 M-entry cap, bisect XG-ext vs plain, report upstream; then SATIC (Blomme et al. 2025) on the same 60 inputs. *Falsifier:* SATIC summed median conflicts `< 0.5×` WDSat symmetry-on with no input `> 1.1×` | accounting (the audit); engineering (SATIC) | 3–5 days after the native runner |
| 3 | **PB-XOR solver on the Hamming-weight normal-basis base** (LinPB: RoundingSat + Gauss–Jordan, Yang–Meel CP 2021; or CryptoMiniSat's BNN/cardinality constraints) — the only base where the solver family is structurally different from Gröbner (§3.3), and Frobenius-stable at prime `n` with no storage | ONB Karatsuba circuit from `ecc2k130/codegen` (Python; needs the Rust emitter first), `cnf.py`'s Sinz counters as the control | median solve and conflicts at `n = 11, 23`, `w = 3`, `m = 3`, planted and uniform. *Falsifier:* no instance at `n = 23` decided where CMS decided none in 25 min; per-bit exponent not below the measured `2^{2.3}` | engineering; the product law is base-independent | 1–2 weeks |
| 4 | **Bosphorus ANF preprocessing** in front of CMS and Kissat on the frozen `S′₄` corpus | corpus + `semaev_sat.rs` ANF writer | conflicts at `ℓ = 8` uniform before/after. *Falsifier:* `< 0.5×`. Expected: it rediscovers the trace fall and degree-2 Macaulay rows the native path already injects | accounting | 1 day |
| 5 | **Complete the SMT ladder to `n = 23`** so `research/pdp_smt_oracle_20261008/` can be classed (its P4 needs `n = 23`) | `pdp_smt_oracle_bench` / `_score` (Rust, exists) | `--n 19 --targets 12`, `--n 23 --targets 8`, seed 20261010, isolated host | accounting | hours of compute |
| 6 | **Yield-law cross-check against thesis Table 7.7** (observed decomposition probability `≈ 2×` the `2^{mℓ}/(m!·2^n)` heuristic at `n = 24`, e.g. `0.372` vs `0.167` at `ℓ = 8`), against this repository's 12-cell confirmation of `C(|F|, m)/#E` | `ic` yield tooling | recount on the thesis' `(n, ℓ)` with signed points and `S₄` root multiplicity separated; a `2×` that survives is a constant in the targets-per-relation box of the product law | accounting | hours |
| 7 | **GPU / bit-sliced pairs-and-solve leaf** (the FES-style enumeration the SAT literature never benchmarks against) | `gpu/fes/` kernels (compiled, cross-checked, not wired; no throughput measured), `s3_x_roots` one-word kernel | candidates per second per device at `ℓ = 12…16` vs one core. *Falsifier:* `< 50×` one core | engineering; moves toy-ladder reach (a PDP record at `ℓ ≥ 14` would be new), not the ratio | 1–2 weeks |

Not worth a round, with the reason: any SMT theory (no `GF(2^n)` solver
exists and a native one would be Gröbner); ML-guided SAT (no ECC benchmarks,
poor size generalisation); in-memory or GPU *SAT* hardware (simulated or
generic instances); Frobenius-invariant subspaces at `n = 131` (do not exist);
cube-and-conquer on the descended system (its natural cubes are the `m − 1`
summand prefixes, which is pairs-and-solve under another name); Semaev's 2015
cubic model (first-fall assumption contradicted at `n = 45`, Kosters–Yeo).

---

## 7. Closed: do not re-run

- **A third SAT solver on the subspace Semaev system.**  `solver_16` §7: "not
  more `d`, not more budget, and not a third SAT solver on the same descended
  system."  Kissat, CaDiCaL, CryptoMiniSat, WDSat, the native solver, cvc5 and
  Z3 are all measured; they differ by constants around `T_enum · 2^ℓ/m`.
- **Riemann–Roch encodings under SAT.**  `3.5·10⁴×` the `S′₄` conflicts at
  `d = 5`; "the coefficient block hides the certificate from the solver."
- **Symmetry breaking inside the IC pipeline's SAT oracle.**  Measured not to
  pay once the degree-2 Macaulay rows are present (5,512 → 6,412 conflicts at
  `n = 9, m = 3`); WDSat's own symmetry breaking is worth exactly the `3!` and
  no more.
- **One-hot Frobenius-phase variables in SAT.**  `n = 19` fails at 1,200 s.
- **Crossbred as the pipeline oracle.**  Loses to the double loop `3.4×`–`1800×`.
- **Hamming ideals under truncated Gröbner.**  H0 stands at every size run.

---

## 8. Visuals and graphs checked

- `research/pdp_sat_smt_survey_20261010/oracle_ladder.svg`: measured per-target
  cost of the SAT arms against the enumeration floor over `ℓ` (data from the
  rows of §2 that report conflicts).
- `research/pdp_sat_smt_survey_20261010/product_law.svg`: where a solver
  constant sits inside the `m·2^n/n` identity, and the `2^{−63}` it would need.
- `research/pdp_sat_smt_survey_20261010/REPORT.pdf`: this note with both
  figures, built by pandoc and headless Chromium (`build.sh`).
- Canonical graphs checked and left unchanged, because no operation count or
  verdict moved: `docs/index-calculus-scoreboard.html`,
  `docs/ic/progress-timeline.json`, `docs/ic/leaderboard.json`,
  `docs/ic/BOUNDARY_TARGETS.md`.  The boundary target "sub-`2^{2ℓ}` oracle at
  fixed `ℓ = 8`" stays open and un-attempted by anything in §6.

---

## References

Repository artefacts are linked inline.  External sources, in the order the
note uses them; items marked † were read as abstracts or secondary
restatements only (ePrint and HAL were behind bot challenges from this
environment).

- M. Trimoska, S. Ionica, G. Dequen, *A SAT-based approach for index calculus on binary elliptic curves*, AFRICACRYPT 2020, LNCS 12174; ePrint 2019/313, https://eprint.iacr.org/2019/313.
- M. Trimoska, S. Ionica, G. Dequen, *Parity (XOR) reasoning for the index calculus attack*, CP 2020, LNCS 12333; arXiv 2001.11229, https://arxiv.org/abs/2001.11229; code https://github.com/mtrimoska/WDSat.
- M. Trimoska, S. Ionica, G. Dequen, *Logical cryptanalysis with WDSat*, SAT 2021, LNCS 12831, https://hal.science/hal-03230392.
- M. Trimoska, *Combinatorics in algebraic and logical cryptanalysis*, PhD thesis, UPJV Amiens, 2021, https://mtrimoska.com/Monika_Trimoska_these.pdf (Tables 7.1, 7.5, 7.6, 7.7, 8.4, 8.5; §9.1).
- M. Trimoska, *EC-Index-Calculus-Benchmarks*, https://github.com/mtrimoska/EC-Index-Calculus-Benchmarks.
- † T. Blomme, M. Cherif, S. Ionica, G. Dequen, *ANF-based satisfiability for Weil-descent cryptographic attacks*, CoDIT 2025, https://hal.science/hal-05176414.
- B. Andraschko, J. Danner, M. Kreuzer, *SAT solving using XOR-OR-AND normal forms*, Math. Comput. Sci. 2024, arXiv 2311.00733.  J. Danner, M. Kreuzer, *Conflict-driven SAT solving using XOR-OR-AND normal forms*, JAIR 86 (2026), https://www.jair.org/index.php/jair/article/view/20298.
- M. Soos, S. Gocht, K. S. Meel, *Tinted, detached, and lazy CNF-XOR solving*, CAV 2020; M. Soos, K. S. Meel, *BIRD*, AAAI 2019.
- D. Choo, M. Soos, K. M. A. Chai, K. S. Meel, *Bosphorus: bridging ANF and CNF solvers*, DATE 2019, arXiv 1812.04580.
- J. Yang, K. S. Meel, *Engineering an efficient PB-XOR solver*, CP 2021.
- O. Zaikin, *Inverting cryptographic hash functions via cube-and-conquer*, JAIR 2022, arXiv 2212.02405.
- Zheng et al., *SAT4CryptoBench*, NeurIPS 2025 D&B, https://github.com/void-zxh/SAT4CryptoBench.
- Im et al., *Accelerating hybrid XOR–CNF SAT natively with in-memory computing*, Nat. Commun. 17:2922 (2026), arXiv 2504.06476.
- A. Ozdemir, G. Kremer, C. Tinelli, C. Barrett, *Satisfiability modulo finite fields*, CAV 2023, ePrint 2023/091.  A. Ozdemir et al., *Split Gröbner bases for satisfiability modulo finite fields*, CAV 2024, ePrint 2024/572.  T. Hader, A. Ozdemir, *An SMT-LIB theory of finite fields*, SMT 2024, arXiv 2407.21169.
- C. Bouillaguet et al., *Fast exhaustive search for polynomial systems in F₂*, CHES 2010.  A. Joux, V. Vitse, *A crossbred algorithm for solving Boolean polynomial systems*, ePrint 2017/372.  I. Dinur, *Improved algorithms for solving polynomial systems over GF(2) by multiple parity-counting*, SODA 2021.  C. Bouillaguet, C. Delaplace, M. Trimoska, SOSA 2022, ePrint 2021/1639.
- S. Galbraith, R. Granger, S. Merz, C. Petit, *On index calculus algorithms for subfield curves*, SAC 2020, ePrint 2020/1315.
- † K. Nagao, *Decomposition attack for the Jacobian of a hyperelliptic curve over an extension field*, ANTS 2010, LNCS 6197; ePrint 2007/112.  † K. Nagao, ePrint 2010/511, 2013/548, 2013/549.
- A. Joux, V. Vitse, *Elliptic curve discrete logarithm problem over small degree extension fields*, J. Cryptology 26 (2013); ePrint 2010/157.  A. Joux, V. Vitse, *Cover and decomposition index calculus on elliptic curves made practical*, EUROCRYPT 2012; ePrint 2011/020.
- J.-C. Faugère, L. Huot, A. Joux, G. Renault, V. Vitse, *Symmetrized summation polynomials*, EUROCRYPT 2014.  J.-C. Faugère, P. Gaudry, L. Huot, G. Renault, *Using symmetries in the index calculus for elliptic curves DLP*, J. Cryptology 27 (2014).
- J.-C. Faugère, L. Perret, C. Petit, G. Renault, *Improving the complexity of index calculus algorithms in elliptic curves over binary fields*, EUROCRYPT 2012 († tables via Renault's slides and HPST).  † C. Petit, J.-J. Quisquater, *On polynomial systems arising from a Weil descent*, ASIACRYPT 2012, ePrint 2012/146.
- Y.-J. Huang, C. Petit, N. Shinohara, T. Takagi, *Improvement of FPPR method to solve ECDLP*, IWSEC 2013; and *On generalized first fall degree assumptions*, ePrint 2015/358.
- † M. Shantz, E. Teske, *Solving the ECDLP using Semaev polynomials, Weil descent and Gröbner basis methods*, 2013, ePrint 2013/596.
- S. Galbraith, S. Gebregiyorgis, *Summation polynomial algorithms for elliptic curves in characteristic two*, INDOCRYPT 2014, ePrint 2014/806; S. Gebregiyorgis, PhD thesis, Auckland.
- I. Semaev, *New algorithm for the discrete logarithm problem on elliptic curves*, ePrint 2015/310, arXiv 1504.01175.  K. Karabina, *Point decomposition problem in binary elliptic curves*, ePrint 2015/319, arXiv 1504.02347.
- M. Kosters, S. L. Yeo, *Notes on summation polynomials*, arXiv 1503.08001.  M.-D. Huang, M. Kosters, S. L. Yeo, *Last fall degree, HFE, and Weil descent attacks on ECDLP*, CRYPTO 2015, ePrint 2015/573.  M.-D. Huang, M. Kosters, Y. Yang, S. L. Yeo, *On the last fall degree of zero-dimensional Weil descent systems*, JSC 2018, arXiv 1505.02532.  S. Kousidis, A. Wiemers, *On the first fall degree of summation polynomials*, JMC 2019, arXiv 1906.05594.
- S. Galbraith, P. Gaudry, *Recent progress on the elliptic curve discrete logarithm problem*, DCC 78 (2016), ePrint 2015/1022.  S. Galbraith, ECC 2015 talk, https://lfant.math.u-bordeaux.fr/ecc2015/documents/galbraith.pdf.
- K. Yokoyama, M. Yasuda, Y. Takahashi, J. Kogure, *Complexity bounds on Semaev's naive index calculus method for ECDLP*, JMC 14 (2020).
- N. Courtois, *On splitting a point with summation polynomials in binary elliptic curves*, ePrint 2016/003.  Y. Tsakou, S. Ionica, ePrint 2021/721.
- D. V. Bailey et al., *Breaking ECC2K-130*, ePrint 2009/541; project page https://ecc-challenge.info/ (fetched 2026-10-10, no completion).  D. J. Bernstein et al., *Faster elliptic-curve discrete logarithms on FPGAs* (sect113r2), ePrint 2016/382.
- R. La Scala, M. Marchesin, S. K. Tiwari, *Hamming ideals and Gröbner bases for ISD-like syndrome decoding*, 2026 (as used in `RESEARCH_HAMMING_IDEAL_PDP.md`).
