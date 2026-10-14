# Requirements for a dedicated root-finder for the point decomposition problem, on prime and binary fields

**Date.** 2026-10-12.  **Status.** Requirements, gates and a baseline protocol.
**Nothing new is measured here.**  Every number is quoted from a frozen
repository artefact or a cited paper, and each row says which; every budget is
a calculation shown with its inputs.  **Class.** Accounting (a requirements
document moves no operation count).
**Origin.** The owner's request of 2026-10-12, with two pages of
Petit–Kosters–Messeng, *Algebraic approaches for the ECDLP over prime fields*
(PKC 2016), §3 "Algebraic attacks on prime curves": Algorithm 1 and the
"back door" remark (§1 restates them).  The request: scope the requirements for
a *new* root-finding algorithm for the point decomposition problem (PDP) of
index calculus, specialised enough to be sped up, tailored to prime fields and
to binary fields, with gates stating how fast it must be against F4/F5 and
general Gröbner-basis computation, and a baseline to compare against.
**Visuals.** [`figures/`](figures/) (three SVGs with their plotting source and
data) and [`REPORT.pdf`](REPORT.pdf) (source `REPORT.typ`).
**Conductor.** The control plane was unreachable in this container (§12), so
no task card is attached; `CLAUDE.md` lets the edit through with a warning.

---

## 0. Summary

1. **The thing asked for is the open problem PKM16 themselves left.**  Their
   slide 19 reads: "Sparse, block-triangular structure, and resemblance to the
   (polynomial time solvable) root-finding problem suggest to build dedicated
   algorithms to solve the systems. Open problem!"  Their conclusion: "only
   practical for very small parameters (Pollard's rho definitely better for
   crypto sizes); open problems: asymptotic complexity, dedicated polynomial
   system solving methods."  No published dedicated solver exists for their
   System (4); the follow-ups the repository knows of (Amadori–Pintore–Sala
   2017, Kudo–Yokota–Takahashi–Yasuda 2018, Yokoyama et al. 2020) change the
   presentation or bound it, and the 2024 ePrint survey by Hu records that no
   subexponential method for prime-field ECDLP is known.
2. **The repository has already built and measured PKM's mechanism**, as the
   `prime-tower` factor base (Kummer, Dickson, 2-isogeny towers) with the
   `pkm-tower` oracle and the engines `f4_fp`, `f4_fp_tower`, `sig_fp_tower`
   ([`RESEARCH_PKM_TOWER_ORACLE.md`](../notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md)),
   and it has the binary analogue as eleven `SystemSolver` engines behind
   `descent-algebraic` ([`GROEBNER_F4F5F6_AUDIT_20261010.md`](../../docs/ic/GROEBNER_F4F5F6_AUDIT_20261010.md)).
   **F4 is measured to be the wrong tool on both tracks** (§3): its per-target
   cost grows at 1.6–2.4 bits per tower variable on the prime systems and at
   8.6–9.3 bits per rung of the base dimension on the binary ones, against
   0.33–0.75 for plain enumeration and the pair table; F5-style criteria made
   it worse; replaying F4's trace across same-shape targets is bounded at
   4.3 % of the rows on the binary cell that was measured.
3. **What a dedicated root-finder must do is therefore not "a faster F4".**
   §4–§5 state the functional and structural requirements: solve for the few
   rational zeros of *one* function on the `2^N`-point grid `V^m` without
   passing through a Macaulay matrix at the solving degree, consume the tower
   (or subspace) as structure, amortise everything target-independent, stack
   the symmetries, and report in the counted units the ledger already uses.
4. **The gates (§6) are slopes before they are constants.**  Against F4/F5 and
   the general Gröbner engines (ours at msolve parity, msolve, Singular,
   PolyBoRi, Magma scripts) the candidate must be at least 10× cheaper at the
   largest cell both finish *and* grow slower with disjoint 95 % intervals over
   four sizes (G2).  Against the enumeration floor its ratio must be below 1
   and falling (G1).  The only outcome that moves an exponent is a growth rate
   below the pair table's, `⌊m/2⌋/m` bits per tower variable, at polynomial
   memory (G3).  End-to-end relevance is a separate gate (G6): at `m = 4` the
   per-target budget admits a degree-2 Macaulay-sized computation on
   ECC2K-130 and a degree-3 one on a 256-bit prime curve, against a measured
   solving degree of 7 on both tracks.
5. **The baseline (§7) is eight arms on frozen corpora**, six of which exist
   in the tree today; the two to build are the F4 trace-replay arm (the
   literature's same-shape lever) and the prime enumeration and pair-table rows
   at the PKM cells, so that "beats F4" can never mean "our F4 was slow".

---

## 1. What the two pages say, and what the repository already has for them

### 1.1 PKM16 §3.1, restated

Replace the `F₂`-linear map of binary Weil descent by any algebraic or rational
map `L` over `F_p` such that, for the chosen `m`,

1. `|{x ∈ F_p : L(x) = 0}| ≈ |{x ∈ F̄_p : L(x) = 0}| ≈ p^{1/m}` (the preimage
   set is rational and has the right size), and
2. `L = L_{n'} ∘ … ∘ L_1` is a composition of low-degree maps.

**Algorithm 1** (index calculus for prime curves): fix `m`; find `L` and its
decomposition; take the factor base `F = {(x, y) ∈ E(K) : L(x) = 0}`; compute
`deg L + Δ` relations by picking `a, b ∈ F_p`, forming `(X, Y) = aP + bQ`, and
solving

```text
    S_{m+1}(x_{11}, …, x_{m1}, X) = 0
    x_{i,j+1} = L_j(x_{i,j})           i = 1 … m,  j = 1 … n' − 1        (4)
    0 = L_{n'}(x_{i,n'})               i = 1 … m
```

(denominators of rational `L_j` cleared); each solution, modulo permutations
of the `x_{i1}`, is a relation `Σ P_i = aP + bQ`; then linear algebra.  `m`
is fixed by the cost of the relation search; `Δ` is a priori 10 to absorb
linear dependencies.  The sentence this document is about: *"At the moment,
we have not investigated the existence of any algorithm better than Gröbner
basis algorithms to solve System (4)."*  The heuristics they list: most
solutions are defined over `F_p` (condition 1); the system has a constant
probability of having a solution; all equations have low degree; System (4)
is "very similar" to the binary Weil-descent systems.

**Constructions of `L`** (PKM16 §3, and the PKC 2016 slides 14–16):
- `p − 1 = S·N'` with `S ≈ p^{1/m}` smooth: `L(X) = X^S − 1`, `V` the subgroup
  of order `S`, `L_j(X) = X^{q_j}` for `S = ∏ q_j`.  The slides remark that NIST
  P-224 has `p − 1 = 2^{96}·N'`.
- `p + 1` smooth: the same through the norm-one torus (Dickson polynomials).
- arbitrary `p`: an auxiliary curve `E'/F_p` with `#E'(F_p) = S·N'`, a
  subgroup `G` of order `S`, isogenies `φ = φ_{n'} ∘ … ∘ φ_1` with kernel `G`,
  and `L_j` the x-coordinate part of `φ_j`.  Finding `E'`: random curves
  (about `|F|` trials) or complex multiplication, which is the second page.
- The second page: construct `N` with a known factorisation and a large
  smooth factor, solve `x² − Dy² = 4N` for small `|D|`, check whether the
  implied `p` is prime (probability `≈ 1/log N`), else try another `D`.  The
  remark: this can produce a sort of "back door" when choosing primes for
  standards, judged unlikely because (i) a user aware of it detects it by
  solving the equation for the given `p` and all small `D`, and (ii) equally
  useful auxiliary curves can be built in time dominated by other steps.

### 1.2 What the repository already holds for each piece

| PKM16 piece | repository component | what is measured | source |
|:--|:--|:--|:--|
| `V = L^{−1}(0)` with a composed `L` | `prime-tower` factor base: `kummer` (`y ↦ y²`), `dickson` (`y ↦ y² − 2`), `isogeny` (Vélu 2-isogeny x-maps); mixed radix as a parameter | all three families built, preimage trees split, `|V| = 2^t` | PKM note §2.1 |
| System (4) | `pkm-tower` oracle: `full` (`m ≤ 3`), `chain` (Semaev 2015), `tree` (Karabina); `raw` or `reduced` presentation; `disjoint` cosets for symmetry | all presentations run at toy sizes | PKM note §2.2 |
| "similar to the binary systems" | Lemma 1: the tower equations are a Gröbner basis with square-free standard monomials; `R = F_p[y]/(towers) ≅ F_p^{V^m}`; Kummer `R` is the group algebra `F_p[(Z/2^t)^m]` | exact | PKM note §2.4 |
| "Gröbner basis algorithms" | `f4_fp` (level with msolve 0.6.5 on generic `F_p` systems), `f4_fp_tower` (sparse, tower-aware F4, the measuring engine), `sig_fp_tower` (F5/GVW criteria), `macaulay-fp` (the XL instrument) | §3.2 | PKM note §§10–16; [`research/f4_fp_msolve_parity_20261011`](../f4_fp_msolve_parity_20261011/README.md) |
| back door / auxiliary curves | `pkm_criterion.rs` (resistance audit: Solinas weight, smoothness of `p ± 1`, auxiliary-curve search) | audit only; "not a full PKM attack" | [`RESEARCH_PKM_CRITERION.md`](../notes/cm-isogeny/RESEARCH_PKM_CRITERION.md) |
| binary Weil descent (PKM16 §2) | `descent-algebraic` with `f4-f2`, `matrix-f4`, `matrix-f5`, `inherited-f4` (default), `crossbred-f2`, `xl-f2`, `sat-cdcl`, `fes-f2`, `wide`; `semaev_decomp` (pairs-and-solve); `mitm` | §3.1 | [`FRAMEWORK.md`](../../docs/ic/FRAMEWORK.md), audit |

The gap PKM16 named, and this document scopes, is one component: the engine
behind `pkm-tower` and `descent-algebraic` that is neither a Gröbner basis
computation nor an enumeration.

---

## 2. The problem, stated as root-finding, on both tracks

**Common form.**  Let `K` be the field, `V ⊂ K` the base set with
`|V| = 2^t` (a power of two in every family the repository runs; mixed radix
changes nothing below), `N = m·t`, and `R = K[y]/(membership) ≅ K^{V^m}` the
coordinate ring of the grid `V^m` (Lemma 1 of the PKM note for towers; the
Boolean ring for subspaces).  Write `S̄ ∈ R` for the normal form of
`S_{m+1}(x_1, …, x_m, x_R)` modulo the membership relations.  The **PDP
root-finding problem** is:

> given `x_R`, return every point of `V^m` (modulo the `S_m` action, and
> modulo the 2-torsion symmetry where it is used) at which `S̄` vanishes, or
> a certificate that there is none.

Expected zeros: `|V|^m / (m!·r)` per target, about one when `|V|^m ≈ m!·r`;
almost every call in relation collection is a refutation.

What a Gröbner basis computes is the ideal `(S̄) ⊂ R`, i.e. the whole
structure of `R/(S̄)`, a `(2^N − σ)`-dimensional algebra, in order to read off
`σ ≤` a few zeros.  Univariate root-finding over `F_q` is polynomial because
`gcd(f, X^q − X)` touches only the rational roots; the multivariate PDP has no
such gcd, and the requirement is an algorithm that finds the rational zeros of
`S̄` on `V^m` with work sub-linear in `2^N` and, per §6, with a growth rate
below the pair table's.

**The two tracks, side by side.**  Every row is from the cited note.

| | binary, `K = F_{2^n}`, `n` prime | prime, `K = F_p` |
|:--|:--|:--|
| `V` | an `F₂`-subspace of dimension `ℓ` (`t = ℓ`); or a Frobenius-stable set (Hamming weight in a normal basis) | `L^{−1}(c)` for a tower `L = L_t ∘ … ∘ L_1`: a coset of `μ_{2^t}` (Kummer), `{w + w^{−1}}` (Dickson), `x(R' + ⟨Q⟩)` on an auxiliary curve (isogeny) |
| membership, as the engine sees it | `n − ℓ` linear relations plus the field equations `v² = v` | `t` quadrics per summand: `y_{j+1} = L_{j+1}(y_j)` with `y_j²` leading |
| the "field-equation half" of Weil descent | present | present: `R ≅ F_p^{V^m}` with square-free standard monomials (Lemma 1) |
| the "conjugate-equation half" | present: `S = 0` over `F_{2^n}` is `n` equations over `F₂` | **absent**: `F_p` has no subfield; `S̄ = 0` is one equation carrying `log₂ p` bits (PKM note §3.2, Prop. 2) |
| unknowns | `m·ℓ` (unsplit), `m·ℓ + (m−2)·n` (chain) | `N = m·t` tower variables (+ `m − 2` free `u_k` in the chain) |
| symmetries usable | negation; `S_m`; 2-torsion (`x ↦ 1/x` on Koblitz, FGHR); Frobenius `C_n` (Faugère–Svartz blocks) | negation; `S_m` (disjoint cosets); 2-torsion of `E` if rational; the group structure of `V` itself (Kummer: `V ≅ Z/2^t`) |
| the reference oracle (enumeration) | pairs-and-solve: `2^{(m−1)ℓ}/(m−1)!` candidates, `O(ℓ)` field ops each via `gcd(q, L_V mod q)` | `tower-exhaustive` on the same base (the PKM note's cell reference) |
| the generic table oracle | `mitm`: `2^{⌊m/2⌋ℓ}` lookups at `2^{⌈m/2⌉ℓ}` memory | `mitm`, `subtract`; family law `S = 0.75·r^{1/6}` |
| best algebraic engine measured | inherited F4 / F6-IC (default), matrix-F4 with splitting | `f4_fp_tower` |
| measured solving degree | 7 at `m = 3` (`ℓ = 4, 5`, unsplit), 4 (chain); Macaulay refutation degree 6 at `m = 3` | 5→6 at `m = 2` (`N = 12…22`); 6, 6, 7, 7 at `m = 3` (`N = 9…18`); 6, 7, 7 at `m = 4` (`N = 8…16`) |
| obstruction each track must beat | no Frobenius-stable subspace at prime `n`; deciding degree `≈ ℓ/2` | one generator over a `2^N`-point grid: XL degree `≥ d + δ_Z`, `δ_Z ≈ N` generically |

Sources: [`RESEARCH_PDP_SAT_SMT_NAGAO_SURVEY.md`](../notes/index-calculus/RESEARCH_PDP_SAT_SMT_NAGAO_SURVEY.md)
§1–2, [`research/pdp_split_dreg_20261011`](../pdp_split_dreg_20261011/README.md),
[`RESEARCH_DREG_MEASUREMENT.md`](../notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md),
PKM note §§2–3 and §§11–16.

---

## 3. Why F4/F5 are the wrong tool here: the measurements

### 3.1 Binary track

- **The solving degree is set by the shape, and it is high.**  At `m = 3`,
  unsplit `S₄`, F4's solving degree is 7 at `ℓ = 4` (four cells, `n = 11…19`)
  and at `ℓ = 5` (`n = 13, 17`); the `S₃` chain reaches 4 but pays `n` extra
  unknowns per cut.  The Macaulay refutation degree is 6 where the first-fall
  degree is 3 (`n = 5, 7`), the Kosters–Yeo decoupling.  The measured deciding
  degree of matrix-F4 on the full descent is `≈ ℓ/2`, which prices it *above
  enumeration* (survey §3.2).
- **Cost per rung.**  Word XORs per target, unsplit `m = 3`, `g = 0`:
  `4.94·10⁷` at `(n, ℓ) = (13, 4)` and `1.93·10¹⁰` at `(13, 5)`, a factor
  391 (`2^{8.6}`); `3.95·10⁷` at `(17, 4)` and `2.43·10¹⁰` at `(17, 5)`,
  a factor 615 (`2^{9.3}`).  Enumeration grows `2²` per rung of `ℓ` at
  `m = 3`, the pair table `2¹`.  Two rungs are a ratio, not a fit; the
  baseline of §7 extends it to four.
- **Guessing beats everything algebraic that was tried**, and is still
  enumeration: the cheapest variant in every `m = 3` cell was unsplit F4 with
  `g = 4` guessed bits; no variant came within 22× of matched rho on the toy
  curves (split-dreg P6).
- **The trace-replay lever is bounded.**  Joux–Vitse's F4Remake saves the rows
  that reduce to zero; on the landed binary cell (`m = 2`, 16-dimensional
  base) 4.3 % of the degree-4 rows are redundant, identically on every target
  ([`RESEARCH_PDP_SPEEDUP_AND_ISOGENY.md`](../../RESEARCH_PDP_SPEEDUP_AND_ISOGENY.md)
  §1.2).  In Joux–Vitse's own experiment (decomposition systems for
  `E(F_{p⁵})`: 5 equations of degree 8 in 4 variables, 29 steps) F4Remake is
  2.1–3.4× their own F4 and 3.4–124× Magma's, with failure probability
  `0.11` at 8-bit `p` falling to `5.8·10⁻⁹` at 32 bits (CT-RSA 2011, §4.1).
- **SAT is enumeration with overhead.**  Plain WDSat visits `2^{3ℓ}/3!` nodes
  to within 0.2 %; the native CDCL+XOR solver 0.97–1.02 conflicts per triple
  at `ℓ = 5…8`; CryptoMiniSat 36.7× the pair count at `d = 8` on ECC2K-130;
  WDSat hit a 120 s watchdog on all 40 targets at `n = 59, ℓ = 9` where the
  pair table answers in 3.86 core-s (survey §2).  WDSat's published speed-ups
  are relative to Magma's F4: 4.2 s / 13.5 s (SAT / UNSAT) against 16.8 s /
  18.7 s on the third summation polynomial at `n = 41, ℓ = 20`, and 0.6 s /
  3.8 s against 229.3 s / 229.4 s on the fourth (Trimoska–Ionica–Dequen,
  arXiv:2001.11229, Tables 4–5), i.e. 60–380× over a Gröbner engine that is
  itself above enumeration.
- **A full basis is the wrong product.**  On the exported `n = 17` Boolean
  decomposition systems (8 targets) the inherited F4 decides a target in
  11–24 ms and F6-IC in 7–16 ms; a full reduced basis takes msolve
  42.7–46.5 s and PolyBoRi 40.4–52.3 s, Singular over 120 s, all verdicts
  agreeing ([`research/gb_cas_comparison_20261011`](../gb_cas_comparison_20261011/README.md)).
- **The one symmetry lever with measured upside** is the Frobenius-equivariant
  block decomposition of the Macaulay matrix: 31×, 39×, 105×, 105×, 128× on
  the rank computation at `n = 11, 13, 17, 19, 23` (rank only, not a full F4
  pass; [`RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md`](../../RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md)
  §3.3).  It shrinks the matrix by `≈ 1/n`; it does not change the degree.

### 3.2 Prime track

Every row below is a Kummer tower at `p₁`, two targets per cell, F4 counts
from the PKM note (§11.7, §12.7, §13.7, §16.9; the round-2/3 build, before the
round-6 full-rank exit, which lowers multiply-adds by up to 39 % and changes no
`D`).

| `m` | `N` | `D` | F4 multiply-adds | width (columns) | peak memory | F4 time (practicality) |
|--:|--:|--:|--:|--:|--:|--:|
| 2 | 12 | 5 | `6.50·10⁷` | 1,231 | 8 MB | 0.05 s |
| 2 | 16 | 5 | `6.18·10⁹` | 4,789 | 36 MB | 1.1 s |
| 2 | 18 | 5 | `4.60·10¹⁰` | 8,644 | 91 MB | 6 s |
| 2 | 20 | **6** | `4.84·10¹¹` | 33,446 | 473 MB | 57 s |
| 2 | 22 | 6 | — | 51,909 | 2.3 GB | 10 min |
| 3 | 9 | 6 | `1.59·10⁸` | 1,186 | 11 MB | 0.1 s |
| 3 | 12 | 6 | `1.06·10¹⁰` | 5,019 | 64 MB | 2.0 s |
| 3 | 15 | **7** | `3.19·10¹²` | 31,397 | 1.1 GB | 6.0 min |
| 3 | 18 | 7 | `8.99·10¹³` | 106,288 | 12.7 GB | 176 min |
| 4 | 8 | 6 | `4.84·10⁷` | 1,049 | 9 MB | 0.1 s |
| 4 | 12 | **7** | `8.55·10¹⁰` | 9,452 | 171 MB | 9.4 s |
| 4 | 16 | 7 | `2.7·10¹³` | 76,837 | 8.1 GB | 55.5 min |

Reading, as calculations on the rows:

- **Slopes.**  `log₂` of the multiply-adds grows per tower variable by
  `(log₂(4.84·10¹¹/6.50·10⁷))/8 = 1.6` at `m = 2`, `(log₂(8.99·10¹³/1.59·10⁸))/9 = 2.1`
  at `m = 3`, and `(log₂(2.7·10¹³/4.84·10⁷))/8 = 2.4` at `m = 4`.  Enumeration
  of `V^{m−1}` grows at `(m−1)/m` bits per tower variable (0.50, 0.67, 0.75)
  and the pair table at `⌊m/2⌋/m` (0.50, 0.33, 0.50).  **F4 grows three to
  seven times faster than plain enumeration on these systems**, which is the
  binary observation again.
- **The degree is not bounded.**  `D` rose from 5 to 6 at `m = 2` (`N = 20`,
  also in the shape-matched null and in the isogeny family), from 6 to 7 at
  `m = 3` (`N = 15`) and at `m = 4` (`N = 12`); the pre-registered rule A1
  reads every line as inconclusive between "bounded" and "slow growth"
  (PKM note §16.10).  The one-generator law (Prop. 2) predicts XL degree
  `≥ d + δ_Z` with `δ_Z ≈ N` for a generic `1/S̄` on the grid; F4's falls
  keep it far below that, but not bounded.
- **It is not the implementation.**  `f4_fp` computes the same reduced bases
  as Singular and is level with msolve 0.6.5 (ratios 0.96–1.42
  single-threaded on the systems above 50 ms; 0.83–1.52 at four threads);
  signature criteria (F5/GVW, `sig_fp_tower`) remove the 0.03 % zero
  reductions but need rows 1–2 degrees higher on 18 of 18 systems and
  1.7–28× the memory (PKM note §13).
- **Replay has a ceiling here too.**  On generic `F_p` systems about 58 % of
  `f4_fp`'s pair rows still reduce to zero (parity note, "Where the time goes
  now"), so a trace-replay arm can save at most that fraction of the
  reduction: a bound of `2.4×`, not an exponent.

### 3.3 The structural reason, in one paragraph

Weil descent turns one equation over `F_{2^n}` into `n` equations over `F₂`
*and* supplies the field equations; the tower supplies the second half for
`F_p` but not the first, so a PKM system is `1 + N` equations in `N` unknowns
of which `N` are rewrite rules, and the single informative equation must cut
`2^N ≈ p` grid points to about one.  Any method whose work is linear algebra
on the monomials of degree `≤ D` pays `Σ_{k≤D} C(N, k)` for its width, and
the measured `D` is 7 and rising.  On the binary side the `n` equations are
there and the deciding degree is still `≈ ℓ/2`.  A dedicated root-finder is
therefore required to get its leverage from somewhere other than the Macaulay
degree: that is requirement S3.

---

## 4. Functional requirements (F)

| id | requirement | acceptance |
|:--|:--|:--|
| F1 | **Instances.**  Input is `(K, E, V, m, presentation, symmetry, x_R)` with `K ∈ {F_{2^n}: n prime, F_p}`; binary `V` a standard subspace (`binary-subspace:dimension=ℓ`, `koblitz-standard-subspace`) or a Frobenius-stable base; prime `V` a `prime-tower` of kind `kummer`, `dickson` or `isogeny` (mixed radix allowed); presentations `full`, `chain`, `tree`; symmetry `none` or `disjoint`. | the solver refuses, with a reason, any instance outside its declared domain; it never substitutes a smaller one |
| F2 | **Plug points.**  It is a `SystemSolver` (binary; `prepare` once per shape, `solve` per target) and/or an `F_p` engine in `fp_solver_registry()` behind `pkm-tower`, with a new registry id; it runs under `ic descent` and `ic bench`, and is registered in ecbench as a method before any end-to-end row (`ecbench-extend`). | listed by `ic bench --list`; a changed algorithm is a new id |
| F3 | **Output contract.**  Verdict ∈ {`Decomposed(witnesses)`, `Refuted`, `BudgetExceeded`}; witnesses are abscissa tuples with lifts verified in the group (some signed sum has abscissa `x_R`); `Refuted` means the solver *proved* the zero set empty; `BudgetExceeded` is never reported as a refutation and never dropped from a cost. | every witness group-verified; the PS-0 verdict and root-set digests reproduce |
| F4 | **Exactness.**  Complete and sound on the frozen corpora (§7.1): identical verdict streams and root sets to the exhaustive reference, planted and uniform targets alike; a randomised solver reports its per-target completeness, measured on ≥ 1,000 uniform targets, and a run with any unexplained miss is void. | zero false refutations; zero false witnesses; digests equal |
| F5 | **Determinism and replay.**  Seeded; every counter deterministic (as `f4_fp_tower`'s multiply-adds are); a replay on the same inputs reproduces every count. | `compare-trace`-style check on every committed row |
| F6 | **Accounting.**  Per target: counted work in the track's unit (binary: 64-bit word XORs and counted field operations as PS-0 counts them; prime: `F_p` multiply-adds / field multiplications with the `ns_per_fp_muladd` calibration), the target-independent `prepare` cost charged separately, peak memory, the largest intermediate object; wall time only as a practicality note with an isolation receipt. | one table, one unit (AGENTS.md §2) |
| F7 | **Language.**  Rust, in-tree, on the existing field and curve types; external baselines only through `gb_cas_export` and the existing worker contracts. | AGENTS.md "no Python" |
| F8 | **Identity.**  Every curve by its ICV1 slug, registered in `docs/curves/registry.json`; every corpus by seed and digest. | AGENTS.md §11 |
| F9 | **Pre-registration.**  Hypothesis, cells, seeds, budgets, decision rule and inadmissible list frozen in a protocol before the first run, as the PKM rounds and the PS tasks do. | protocol file committed before data |
| F10 | **Both tracks from one design.**  The design document states which parts are field-generic (grid algebra, symmetry, batching) and which are track-specific (Frobenius, tower rewriting), and the binary and prime variants share the corpus format and the table. | §6.3 table filled for both tracks |

---

## 5. Structural requirements (S): what the algorithm must exploit, must not assume, and must report

| id | requirement | why (measured or derived) |
|:--|:--|:--|
| S1 | **Membership is consumed as structure, never as one polynomial of degree `|V|`.**  The algorithm works in `R ≅ K^{V^m}` with the tower (or subspace) relations as rewrite rules, or directly with the group structure of `V` (Kummer: `V ≅ Z/2^t`, `R` the group algebra `K[(Z/2^t)^m]`, square-free monomials are characters; isogeny: `V = x(R' + ⟨Q⟩)`; binary: `L_V` linearised). | the naive presentation `⟨S, F_V(x_i)⟩` is the Yokoyama regime, bounded below and known to be no better than generic (PKM note §1.2) |
| S2 | **It targets the rational zeros, not the ideal.**  Nothing it computes is a Gröbner basis or an equivalent-cost object (a triangular set of `R/(S̄)`, a full normal-form table).  Where it uses Frobenius or field equations (binary `v² = v`; prime: rationality is already in the tower, and `X^p` reductions are admissible only when their cost is counted) it says so. | §2: `σ ≤` a few zeros against a `2^N`-dimensional algebra |
| S3 | **No Macaulay-degree dependence in the exponent.**  If the method contains linear algebra on a monomial basis, the basis has a fixed degree `D` that fits the admissible degree of §6.4 (2 on ECC2K-130 at `m = 4`; 3 on a 256-bit prime curve), and the design states the mechanism by which the single `F_p` equation (prime) or the `n` descended equations (binary) are turned into zeros without raising `D`: evaluation/interpolation on the group `V^m`, norm or resultant folding through the tower, a sparse-spectrum zero search, a representation-technique split on the exponent group, or another named mechanism, each with its cost model written before the run. | §3.2 slopes; Prop. 2 of the PKM note; the admissible degrees of §6.4 |
| S4 | **Shape invariance is exploited and charged.**  Everything independent of `x_R` is in `prepare` and charged once per `(E, V, m)`; the per-target program is branch-free or has bounded branching so that `K` targets run in lock-step; cost is reported as `(T_prepare, T_target)` at `K ∈ {1, 16, 256}` and in the `S·T²` form of G5. | F4Remake's gain is 2.1–3.4× and bounded by the zero-row fraction (§3.1–3.2); the lock-step form is where the PS program's 8–15× per pair came from (PS-2) |
| S5 | **Symmetries are stacked, each with its claimed and measured factor.**  `S_m` (symmetrised or disjoint cosets), negation, 2-torsion (FGHR; `x ↦ 1/x` on Koblitz `b = 1`), Frobenius `C_n` (binary; the Faugère–Svartz blocks are the measured 31–128× on rank), `D₄` on `S₅`; the group structure of a Kummer `V` is a symmetry of the grid, not of the equation, and must be declared as such. | the inventory in [`RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md`](../../RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md) §6: every entry is polynomial in `n` or a constant; they compound (E13) and are not to be double-counted with the orbit factor on relations |
| S6 | **Memory is polynomial in `N`**, or the table is declared and charged, which moves the row to the table-oracle domain (`S·T²` with the table in `S`, the Bernstein–Lange control at equal advice). | the IC-with-table rows sit `2^{25}–2^{27}` above the generic line (leads note §2) |
| S7 | **Track-specific obstruction named and addressed.**  Binary: at prime `n` there is no Frobenius-stable subspace, so the Frobenius lever is orbits and blocks, not the base; prime: no subfield, so no conjugate equations and no `L_V`-style gcd below `2^t` per summand. | PKM note §3.2; survey §1 |
| S8 | **Must not**: skip group verification; count `BudgetExceeded` as a refutation; change the base per target; choose seeds after results; charge the auxiliary-curve search (isogeny tower) as free; report a stage ratio as a speedup; compare a restricted-`p` (smooth `p ± 1`) row to a generic-`p` row without saying so. | PKM note §5.4; AGENTS.md §6 |

---

## 6. Gates: how fast, against what, and what counts as a result

### 6.1 Boundaries, stated before anything is built

**Floors (derived).**
- Enumeration, both tracks: `T_enum(m, t) = 2^{(m−1)t}/(m−1)!` candidate
  sub-tuples per target, `O(t)` field operations each
  (`gcd(q, L_V mod q)` binary; the tower analogue prime).  Growth
  `(m−1)/m` bits per tower variable; `m − 1` bits per rung of `ℓ`.
- Pair table (meet in the middle): `2^{⌊m/2⌋t}` lookups at `2^{⌈m/2⌉t}`
  memory.  Growth `⌊m/2⌋/m` bits per tower variable; `⌊m/2⌋` per rung.
- Generic line for the whole method: `S·T² ≳ r/|Aut|` (Corrigan-Gibbs–Kogan),
  `|Aut| = 2n` for Koblitz, `2` for prime curves without CM automorphisms.
- Product law, binary (survey §1): relations × targets × oracle `≈ m·2^n`
  (`/n` Frobenius-stable); an oracle that reaches rho on ECC2K-130 must beat
  exhaustive search over its own candidate set by `2^{−(63.2 + log₂ m)}`.
- Family floor, prime (PKM note §3.7): `Θ(r^{2/(m+1)})` for a per-target
  cost `r^{o(1)}`; at or above rho for `m ≤ 3`.

**References (measured, with their artefacts).**

| id | arm | track | unit | frozen numbers to beat | artefact |
|:--|:--|:--|:--|:--|:--|
| R1 | pairs-and-solve `semaev_decomp::decompose` | binary | ns per pair, counted field ops, gates | 337.7, 375.4, 393.3, 421.4, 452.3 ns/pair at `ℓ = 8…12` (A/A spread 1.6–8.8 %); PS-2 scalar kernel 0.067–0.12 of that | [`research/pdp_oracle_speedup_20261011`](../pdp_oracle_speedup_20261011/PLAN.md) PS-0, PS-2 |
| R2 | pair table / `mitm` | binary | core-s, lookups, `2^{2ℓ}` entries | `n = 59, ℓ = 9`: 80/80 verdicts, 3.86 core-s median | `research/sat_factor_base_review_20260908` STAGE13 (via survey §2) |
| R3 | repository Gröbner engines: inherited F4 / F6-IC, matrix-F4/F5, wide | binary | word XORs, reductions, splits | split-dreg grid (192 runs, exact); `n = 17` dim 6: 185–256 ms per 8 targets; `n = 31` dim 9: 1.6–2.4 s per target at an exhausted 2,000-node budget | split-dreg README; audit §4.2, §4.5 |
| R4 | SAT: native CDCL+XOR, WDSat `61c6ff3`, CryptoMiniSat 5.14 | binary | conflicts / decisions per target | frozen 60-input regression (`n ∈ {15,17,19}`, `ℓ ∈ {5,6}`; symmetry on/off 1,280,281 / 6,246,481 conflicts); `n = 59, ℓ = 9` watchdog | `research/index_calculus_baseline_20260914/regression`; survey §2 |
| R5 | external Gröbner: msolve 0.6.5, Singular 4.3.2, PolyBoRi; Magma scripts (not run) | both | wall (process) | Boolean `n = 17`: 42.7–46.5 s (msolve), 40.4–52.3 s (PolyBoRi), > 120 s (Singular); `F_p`: msolve ÷ `f4_fp` 0.96–1.42 | [`research/gb_cas_comparison_20261011`](../gb_cas_comparison_20261011/README.md); parity note |
| R6 | `tower-exhaustive`, `mitm`, `subtract` on the `prime-tower` base | prime | counted field ops, table entries | `tower-exhaustive`: `2^t` square roots per target at `m = 2` (PKM note §10.6); table family law `S = 0.75·r^{1/6}` | PKM note §5.1; ledger §11.3 |
| R7 | `f4_fp_tower` (and `f4_fp`, `sig_fp_tower`) | prime | multiply-adds, width, `D`, peak memory | the §3.2 table | PKM note rounds 2–7 |
| R8 | F4 trace replay (F4Precomp/F4Remake) in `f4_fp_tower` and `pq_f4_f2` | both | same as R3/R7 | **to build** (§7.2); ceilings: 4.3 % of rows (binary `m = 2`), ≤ 58 % of pair rows (`F_p` generic) | PDP-speedup note §1.2; parity note |

### 6.2 The gates

| gate | what is compared | numeric condition | what passing means |
|:--|:--|:--|:--|
| **G0** exactness | the candidate against the exhaustive reference on every corpus cell | zero false refutations, zero false witnesses, digests equal; every `BudgetExceeded` listed | a row may be read at all |
| **G1** vs the enumeration floor | counted work per target ÷ `T_enum(m, t)` (in the arm's unit converted by the measured factor) | `< 1` at the largest rung both finish, and falling: ratio at rung `t+1` ≤ 0.6 × ratio at `t` over ≥ 4 consecutive rungs, on uniform targets (binary `ℓ ≥ 8, m = 3`; prime `t ≥ 4, m = 3, 4`).  Minimum: the survey's falsification target, `< 0.5·T_enum` at two consecutive rungs `ℓ ≥ 8` | the candidate does algebra the enumerator does not |
| **G2** vs F4/F5 and the general Gröbner engines (the gate the request asks for) | the candidate against R3 (binary) / R7 (prime) and R8, and against R5 on the exported systems, same frozen systems, same host, interleaved | (a) counted work ≤ 0.1 × the best Gröbner arm's at the largest cell both finish; (b) fitted growth `β` (bits per tower variable, or per rung of `ℓ`) lower than that arm's with disjoint 95 % bootstrap intervals over ≥ 4 sizes.  The slopes to beat today: prime 1.6 (`m = 2`), 2.1 (`m = 3`), 2.4 (`m = 4`) bits per tower variable; binary `2^{8.6}`–`2^{9.3}` per rung of `ℓ` at `m = 3` (to be re-fitted over four rungs by §7).  (c) Not slower than msolve on any exported system where msolve finishes | "faster than F4/F5/GB" in both constant and exponent, and not an artefact of our engines |
| **G3** vs the generic line (the only exponent gate) | the candidate's fitted `β` against the pair table's `⌊m/2⌋/m` per tower variable (`⌊m/2⌋` per rung), at memory polynomial in `N` | `β` below the pair-table slope with a disjoint 95 % interval over ≥ 4 sizes; or, with a declared table, `log₂(S·T²) − log₂(r/|Aut|) ≤ 0` on a rung (N2 of the novel-directions note) | the first non-generic signal in the program; an **advance** |
| **G4** vs SAT (binary only) | the candidate against R4 on the frozen regression and at `n = 59, ℓ = 9` | fewer counted decisions/conflicts-equivalent than WDSat and the native solver on every regression input; decides the `n = 59, ℓ = 9` uniform targets that WDSat left at the watchdog | the published SAT line is not a competitor |
| **G5** amortised accounting | `(S, T_prepare, T_target, K)` per cell, `K ∈ {1, 16, 256}` | report `log₂(S·T²) − log₂(r/|Aut|)` per rung; precompute charged; the Bernstein–Lange preprocessing-rho control at equal advice stated beside any "online" claim | the row is comparable with the rest of the ledger |
| **G6** end-to-end admissibility | the candidate's per-target cost against the budget of §6.4 at the arity that can beat rho (`m ≥ 4`) | binary (ECC2K-130): `T_target ≲ 2^{35}` group-operation equivalents at `m = 4`, `2^{41}` at `m = 5`; prime (256-bit): `T_dec ≲ p^{1/4}/32 ≈ 2^{59}` at `m = 4` with `|F| ≈ p^{1/4}/2`, `≈ 2^{76}` at `m = 5`; equivalently an admissible Macaulay degree of 2 (binary) and 3 (prime).  On toy sizes: an ecbench session with every phase priced, matched rho, `claim` verdict, then `challenge`/`verdict` against the frontier | the result bears on a deployed size; without it a G3 pass is a toy exponent |
| **G7** classification | every row labelled per AGENTS.md §3 | advance / engineering / relabelling / accounting | only a G3 pass is a result |

Wall time is a practicality column in every gate; operation counts are the
metric (AGENTS.md §6).  A gate passed by changing the unit, the cell, the
seeds, the base per target, or by dropping refutations is not passed.

### 6.3 The one table (template every round fills)

| arm | track | cell (`m`, `n` or `log₂ p`, `t`) | unit | work / target (uniform) | ÷ `T_enum` | ÷ pair table | ÷ best Gröbner arm | `β` (95 %) | memory | `T_prepare` | verified | class |
|:--|:--|:--|:--|--:|--:|--:|--:|:--|--:|--:|:--|:--|
| R1 … R8, candidate | | | | | | | | | | | | |

Rows for every reference arm and the unmodified baseline, including the rows
that do not finish (as censored), and the same table drawn in
`docs/index-calculus-scoreboard.html` once a row carries an end-to-end `S`.

### 6.4 Where the budgets come from (calculation)

**Binary, ECC2K-130** (`n = 131`, `r ≈ 2^{129}`, matched rho `2^{60.81}`
group operations with `⟨−1⟩ × ⟨π⟩`).  Base `V` of dimension `ℓ`, `|F| ≈ 2^ℓ`
points, `K = 2^ℓ/(2n)` orbit columns, relations needed `≈ K`, decomposition
probability `≈ 2^{mℓ}/(m!·2^n)`.  Choose `ℓ = (n + log₂ m!)/m` so that the
probability is about one; then every target yields a relation and the
relation phase costs `K·T_target`.  The budget is `T_target ≤ 2^{60.81}/K`
with the sparse linear algebra `≈ 5·K²` also below rho:

| `m` | `ℓ` | `K` | linear algebra `5K²` | `T_target` budget | `T_enum` | required ratio to enumeration | Macaulay width `C(mℓ, D)`, `D = 2 / 3 / 4` |
|--:|--:|--:|--:|--:|--:|--:|:--|
| 3 | 44.5 | `2^{36.5}` | `2^{75.3}` **> rho** | — | — | dead at the linear algebra | — |
| 4 | 33.9 | `2^{25.9}` | `2^{54.1}` | `2^{34.9}` | `2^{99.1}` | `2^{−64.2}` | 9,180 / 410,040 / `1.36·10⁷` |
| 5 | 27.6 | `2^{19.6}` | `2^{41.5}` | `2^{41.2}` | `2^{105.8}` | `2^{−64.6}` | 9,453 / 428,536 / `1.45·10⁷` |
| 6 | 23.4 | `2^{15.4}` | `2^{33.1}` | `2^{45.4}` | `2^{110.2}` | `2^{−64.8}` | 9,870 / 457,310 / `1.58·10⁷` |

A binary group addition at `n = 131` is of the order of `2^{10}`–`2^{12}`
word operations (one inversion plus a few multiplications), so the `m = 4`
budget is about `2^{46}` word XORs per target.  Dense elimination of a
Macaulay matrix of width `W` costs about `W³/64` word XORs: `2^{33.6}` at
`D = 2` (fits), `2^{49.9}` at `D = 3` (over, unless the elimination is
sparse), `2^{65}` at `D = 4` (out).  **Admissible degree at `m = 4`: 2,
perhaps 3; measured: 7.**  These agree with the survey's `2^{−(63.2 + log₂ m)}`
to within the rounding of `ℓ` and the orbit size.

**Prime, 256-bit `p`** (`r ≈ p`, rho `≈ 2^{128}`).  From
[`RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md`](../../RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md)
§2.4, with `B = |F|` points and a per-decomposition cost `T_dec`:
`cost(B) ≈ (m!/2^m)·p·T_dec/B^{m−1} + m·B²`, optimum `B = (p·T_dec)^{1/(m+1)}`.
- `m = 3`: total `≈ 3.5·(p·T_dec)^{1/2}`: never below rho, ties only for
  `T_dec = O(1)`.
- `m = 4`: total `≈ 4·(p·T_dec)^{2/5}`; below rho iff `T_dec < p^{1/4}/32`,
  at `B ≈ p^{1/4}/2`.  At 256 bits: `T_dec ≲ 2^{59}` group operations,
  `|V| = 2B ≈ 2^{64}`, `t = 64`, `N = 256`.
- `m = 5`: the formula's optimum `B = p^{1/4}` saturates the decomposition
  probability (`(2B)^5 > 5!·p`), so the regime is one relation per target:
  `B ≈ (5!·p)^{1/5}/2 ≈ 2^{51.6}`, linear algebra `5B² ≈ 2^{105.5}`, and
  `T_dec ≲ 2^{128}/B ≈ 2^{76}`; `t ≈ 53`, `N ≈ 265`.

A prime-field group addition is `≈ 2^{3}`–`2^{4}` field multiplications
with batched inversion, so the `m = 4` budget is `≈ 2^{63}` multiply-adds per
target.  Macaulay widths `C(256, D)`: 32,640 (`D = 2`), `2.76·10⁶` (3),
`1.75·10⁸` (4), `8.8·10⁹` (5), `3.7·10¹¹` (6), `1.3·10¹³` (7); with an
elimination costing between `W²` and `W^{2.8}` the degree-3 step is
`2^{43}`–`2^{60}` (fits), degree 4 `2^{55}`–`2^{77}` (marginal), degree 5
`2^{66}`–`2^{93}` (out at `m = 4`, marginal at `m = 5`).  **Admissible degree:
3 at `m = 4`, 4–5 at `m = 5`; measured: 7 at `m = 4`, `N ≤ 16`, rising.**

Both conversions (group operation to word XORs or multiply-adds, elimination
exponent) are models, stated so that a row can be checked against them; the
ecbench primitive-level axes (`field_muls`, `field_sqrs`, `field_invs`)
replace them once an end-to-end session exists.

---

## 7. The baseline: what to run before any candidate, and how

### 7.1 Frozen corpora

| track | corpus | cells | targets | reference verdicts | status |
|:--|:--|:--|:--|:--|:--|
| binary | PS-0 corpus | `ℓ = 8…12`, `n = 3ℓ`, Koblitz-family `y² + xy = x³ + x² + 1`, ICV1 slugs `icv1-f2m24-tm5983-28fa8b30` … `icv1-f2m36-t473713-f90ff4d7` | 64 uniform + 16 planted per rung, SplitMix64 seed `20261011 + ℓ` | exhaustive, verdict and root-set digests committed, every root re-added in the group | exists ([`baseline_v1/`](../pdp_oracle_speedup_20261011/baseline_v1/README.md)) |
| binary | split-dreg grid | `m:a:n:ℓ ∈ {3:1:11:4, 3:0:13:4, 3:1:17:4, 3:1:19:4, 3:0:13:5, 3:1:17:5, 4:1:11:3, 4:0:13:3}` | 4 seeded + 256 yield targets per cell | exhaustive (`complete_and_exact`) | exists ([`runs-v2/`](../pdp_split_dreg_20261011/README.md)) |
| binary | SAT regression and the `n = 59, ℓ = 9` fixtures | 60 inputs (`n ∈ {15,17,19}`, `ℓ ∈ {5,6}`); 40 targets at `n = 59` | certified SAT/UNSAT; 80/80 MITM verdicts | exists (`research/index_calculus_baseline_20260914`, STAGE13) |
| prime | PKM rounds | Kummer `p₁`: `m = 2`, `N = 12…22`; `m = 3`, `N = 9…18`; `m = 4`, `N = 8…16`; Dickson and isogeny at the pilot sizes; shape-matched null | 2 random targets per cell (rounds), 8 + 8 per the pre-registration | `pkm_tower_check verify`: exhaustive over `V^m` | exists (`research/pkm_tower_*`); extend to 8 + 8 and to four consecutive `t` per `m` |
| prime | new: the four-rung ladder | per `m ∈ {3, 4}` and kind, four consecutive `t` within a 600 s / 8 GiB budget per call, `p` sized `|V|^m ≈ m!·r`, public seeds | 8 planted + 8 uniform | exhaustive (`verify`) | **to freeze** (B-0) |

### 7.2 Arms and what must be built

| arm | exists | to build | effort (planning estimate, person-days) |
|:--|:--|:--|--:|
| R1 pairs-and-solve, R2 `mitm` (binary) | yes | run on the four-rung ladder with counted units | 0.5 |
| R3 repository Gröbner engines (binary) | yes | re-fit the per-rung slope over four rungs of `ℓ` at fixed `n`, both unsplit and chain | 1 |
| R4 SAT arms (binary) | yes (external binaries recorded by hash) | rerun on the ladder where the binaries are present; `unavailable` recorded otherwise | 0.5 |
| R5 external Gröbner (both) | yes for Boolean and generic `F_p` systems | export the **tower** systems through `gb_cas_export`; run msolve/Singular/PolyBoRi; Magma scripts written, not run | 1 |
| R6 `tower-exhaustive`, `mitm`, `subtract` (prime) | yes | rows at every PKM cell in counted field operations, one unit with R7 via the `ns_per_fp_muladd` calibration | 1 |
| R7 `f4_fp_tower` (prime) | yes | nothing; its rows are frozen | 0 |
| R8 trace replay (both) | no | `F4Precomp`/`F4Remake` mode in `f4_fp_tower` and `pq_f4_f2`: record the useful multiples per step on target 0, replay on the rest, detect a divergent step and fall back to a full step (charged); report the saved fraction against the 4.3 % / 58 % ceilings | 3 |
| B-3 the table of §6.3, run report, PDF | — | `pdp_rootfinder_table` native tool: reads the raw records, fits `β` with bootstrap intervals, writes the table and the SVGs | 1 |

### 7.3 Procedure (AGENTS.md §8, §10; research-visuals)

1. **Baseline first.**  Host manifest, clean build hashes, A/A spread and
   ABAB interleaving (≥ 5 rounds) on an isolated core
   (`tools/isolated_bench.py`); the counting build separate from the timing
   build, both hashed.
2. **Same inputs, same verdicts.**  Every arm reads the frozen corpus; a
   verdict disagreement voids the round until explained.
3. **Fits.**  `β` from `log₂(work)` against `N` (or `ℓ`) over ≥ 4 sizes, 95 %
   bootstrap intervals over targets within a size; a two-size ratio is
   reported as a ratio.
4. **Decision rule frozen** before the candidate runs: the conditions of
   G1–G3, the stop conditions of §9, the inadmissible list.
5. **Report.**  Run report with at least one quantitative graph and a PDF
   beside the records; a results report in this directory when a comparison
   lands; `benchmarks/<date>-<slug>/` copies if a speed is reported; scoreboard
   and `progress-timeline.json` rows only for end-to-end `S`.

---

## 8. Candidate mechanisms that meet S1–S8 (for the record; none measured)

Listed so that a design can be checked against the requirements before it is
built.  Each carries its expected class and the cheapest falsifier.  None is
claimed to work.

| id | mechanism | track | expected class | cheapest falsifier |
|:--|:--|:--|:--|:--|
| C1 | **Compiled per-shape straight-line root test with lock-step batching** (the PS program's direction: incremental coefficient update, bit-slicing, GPU) applied to the *last-summand* root test, so the enumeration constant falls by 10–100× | both | engineering (`T/T_enum` → constant `< 1`, slope unchanged) | PS-2/PS-3 gates; G1 passes on the constant and fails on the slope |
| C2 | **Sparse-spectrum zero finding on the Kummer group algebra**: `S̄` has at most `(Σ_{k≤m−1} C(t, k))^m` characters (square-free monomials of degree `≤ m−1` per summand); find the points of `(Z/2^t)^m` where a polynomially-sparse exponential sum vanishes without the full `N·2^N` transform | prime (Kummer, Dickson) | the only entry whose success would move `β` | compute `δ_Z(1/S̄)` on toy grids (the PKM note's instrument): if the degree of `1/S̄` off the zero set is `N − O(1)` for every planted instance, no low-degree certificate exists and the mechanism needs a different handle; if it is `O(1)` for a family, build the search |
| C3 | **Norm/resultant folding through the tower** (the prime analogue of `gcd(q, L_V mod q)`): fold one summand at a time by the Dickson/Chebyshev recurrence of its tower level, `O(t)` products of degree-`2^{m−1}` polynomials per folded variable | prime | engineering: this *is* `tower-exhaustive`'s leaf; folding two variables at once costs `2^t` per fold and returns the enumeration exponent | counted operations against `T_enum` on the ladder; a slope below `(m−1)/m` would be the surprise |
| C4 | **Representation technique on the exponent group** (Kummer: `x_i = g·ζ^{e_i}`, `e ∈ (Z/2^t)^m`): split `e = e' + e''` as in subset-sum algorithms (Howgrave-Graham–Joux; Becker–Coron–Joux) | prime (Kummer) | long shot: `S̄ = 0` is not additive in `e`, and Wagner-type methods need a homomorphism to a small quotient (leads note §7b) | census on planted instances: does any `F_p`-linear functional of the characters of `S̄` split consistently across `e', e''`?  Predicted no; one afternoon |
| C5 | **Guess-and-reduce with the last-fall-degree invariant**: fix `g` bits (binary) or `g` tower levels (prime), stop each branch at the degree where linear forms first appear, use the forms to prune | both | engineering; measured: `g = 4` cheapest at `m = 3`, `D_lm = D` on the prime rounds (PKM note §14) | the split-dreg protocol with `g ∈ {0, 2, 4, 8}` on the ladder; the per-rung slope decides |
| C6 | **Frobenius-equivariant blocks** as the elimination inside whichever method keeps a Macaulay step (binary) | binary | engineering, 31–128× measured on rank, `≈ 1/n` per dimension | the G2 scope note's §4 milestones |
| C7 | **F4 trace replay** (R8) | both | engineering, ceilings 4.3 % / 58 % | the R8 build itself |
| C8 | **Symmetry-stacked `S₅` with the `D₄ = (Z/2)³ ⋊ S₄` invariants** at `m = 4` (Bézout 4096 → 512 on paper, projected ≈ 630× on `C₄`) | both | constant; the cleanest test of whether torsion symmetry is ever more than a constant | the leads note §8 item 2, never run |

A candidate that cannot be placed in this table with a cost model and a
falsifier does not meet S3.

---

## 9. Decision rules and stop conditions (frozen with this document)

- **Success (an advance).**  G0 and G5 reported; G3 passed on ≥ 4 consecutive
  rungs at `m ≥ 3` on either track, on uniform targets; then G6 projected
  with every phase priced, and an ecbench `verdict` sought.
- **Engineering result.**  G1 and G2(a) passed, G2(b) or G3 failed: record the
  constant, label it engineering, no scoreboard row beyond a stage diagnostic.
- **Stop the candidate.**  Any of: ratio to `T_enum` ≥ 1 at the largest rung
  both finish; fitted `β` not below the best Gröbner arm's over four rungs;
  a G0 failure not traced to an instrument bug; memory above the declared
  bound.
- **Instrument first.**  A verdict the exhaustive checker contradicts, or a
  baseline row that does not reproduce its committed digest, voids the round
  until explained.
- **Inadmissible** (PKM note §5.4, AGENTS.md §6): changing the tower, base or
  cell per target; counting `BudgetExceeded` as "does not decompose"; leaving
  refutations out of a cost; charging the auxiliary-curve search as free;
  seeds chosen after results; a per-call cost quoted as a speed; a
  restricted-`p` row compared to a generic-`p` row without saying so;
  extrapolation presented as measurement.

---

## 10. Status of what was requested

| requested item | status | where |
|:--|:--|:--|
| Requirements for a new PDP root-finder (functional, structural) | **verified complete** as a document; nothing built | §4, §5 |
| Read the two PKM16 pages and reflect them | **verified complete**; the restatement checked against the PKC 2016 slides' text (the paper's PDF was not retrievable from this container, §12) | §1 |
| Tailoring to prime fields and to binary fields | **verified complete** | §2 table, S7, §6.4 |
| Gating requirements against F4/F5 and general Gröbner computation | **verified complete**: G2 with the measured slopes to beat; G1, G3, G6 give the floors the Gröbner comparison alone would miss | §6 |
| A baseline to compare against | **specified, not run**: eight arms, five corpora; two arms to build (R8, prime R6 rows); effort table | §7 |
| The novel algorithm itself | **not attempted** (out of scope: the request was the requirements); the candidate table is a checklist, not a design | §8 |
| Measurements | **none made**; every number quoted with its artefact | — |

---

## 11. Figures

- **Figure 1** [`figures/gating_ladder.svg`](figures/gating_ladder.svg):
  the gate ladder with the two tracks' boundaries and reference arms.
- **Figure 2** [`figures/solving_degree.svg`](figures/solving_degree.svg):
  measured F4 solving degrees on both tracks (split-dreg and the PKM rounds)
  against the admissible degrees of §6.4.
- **Figure 3** [`figures/cost_growth.svg`](figures/cost_growth.svg): measured
  growth of F4's per-target work, indexed to the smallest cell of each line
  (so that word XORs and multiply-adds share one dimensionless axis), against
  the enumeration and pair-table slopes.

Data and plotting source: [`figures/data.json`](figures/data.json),
[`figures/plot.py`](figures/plot.py) (lays out recorded values; computes
nothing).  Graphs checked and left unchanged: `docs/index-calculus-scoreboard.html`,
`docs/ic/progress-timeline.json`, `docs/ic-leaderboard.html`,
`docs/bounds/FRONTIER.md`: this document carries no measured `S`.

---

## 12. Environment record

- Conductor: `conductor doctor` reports the control plane unreachable
  (`http://localhost:8080`), no database at `localhost:55432`, Docker daemon
  not reachable; `conductor up` fails for the same reason.  No task card is
  attached; no other session held these paths in any reachable control plane.
- Literature access: `eprint.iacr.org` PDFs return a Cloudflare challenge
  page through this container's proxy (HTTP 403 on the PDF), and the PKC 2016
  paper's host (`christophe.petit.web.ulb.be`) returned 503 and a TLS error
  through the proxy.  Quoted instead: the PKC 2016 slides (text extracted
  with `pdftotext`), the CT-RSA 2011 Joux–Vitse paper (Grenoble mirror), the
  arXiv version of the WDSat work (HTML), and the repository's literature
  notes (`KN-LIT-c137bd`, `KN-LIT-102cdb`, `KN-LIT-99dfe4`, `KN-LIT-7604`,
  `KN-LIT-440ec2`, `KN-LIT-4e8513` in `crypto-autoresearcher/knowledge`).
- Report rendering: Typst and matplotlib installed as Python wheels in a
  scratch virtualenv outside the repository, per the research-visuals skill.

---

## References

- C. Petit, M. Kosters, A. Messeng, *Algebraic approaches for the elliptic
  curve discrete logarithm problem over prime fields*, PKC 2016, LNCS 9615;
  and the PKC 2016 slides (University of Oxford / UC Irvine / Passau).
- A. Amadori, F. Pintore, M. Sala, *On the discrete logarithm problem for
  prime-field elliptic curves*, Finite Fields Appl. 51 (2018); ePrint 2017/609
  (`KN-LIT-99dfe4`).
- M. Kudo, Y. Yokota, Y. Takahashi, M. Yasuda, *Acceleration of index calculus
  for solving ECDLP over prime fields and its limitation*, CANS 2018.
- K. Yokoyama, M. Yasuda, Y. Takahashi, J. Kogure, *Complexity bounds on
  Semaev's naive index calculus method for ECDLP*, J. Math. Cryptol. 14 (2020).
- J. Hu, *Implementation analysis of index calculus method on elliptic curves
  over prime finite fields*, ePrint 2024/1923.
- A. Joux, V. Vitse, *A variant of the F4 algorithm*, CT-RSA 2011; ePrint
  2010/158 (`KN-LIT-c137bd`).
- M. Trimoska, S. Ionica, G. Dequen, *Parity (XOR) reasoning for the index
  calculus attack*, arXiv:2001.11229 / CP 2020 (`KN-LIT-102cdb`); *A SAT-based
  approach for index calculus on binary elliptic curves*, AFRICACRYPT 2020.
- J.-C. Faugère, P. Gaudry, L. Huot, G. Renault, *Using symmetries in the
  index calculus for elliptic curves discrete logarithm*, J. Cryptology 2014;
  J.-C. Faugère, L. Huot, A. Joux, G. Renault, V. Vitse, *Symmetrized
  summation polynomials*, EUROCRYPT 2014.
- M. Kosters, S. L. Yeo, *Notes on summation polynomials* (`KN-LIT-7604`);
  M.-D. Huang, M. Kosters, S. L. Yeo, *Last fall degree, HFE, and Weil descent
  attacks on ECDLP*, CRYPTO 2015 (`KN-LIT-440ec2`, `KN-LIT-4e8513`).
- J.-C. Faugère, J. Svartz, *Gröbner bases of ideals invariant under a
  commutative group: the non-modular case*, ISSAC 2013.
- H. Corrigan-Gibbs, D. Kogan, *The discrete-logarithm problem with
  preprocessing*, EUROCRYPT 2018.
- Repository notes cited inline: the PKM tower oracle, the SAT/SMT survey,
  the split-dreg study, the solving-degree note, the F4/F5/F6 audit, the CAS
  comparison, the msolve parity note, the PDP-speedup and G2 scope notes, the
  prime-field breakthrough program, the IC leads note, the PS program plan.
