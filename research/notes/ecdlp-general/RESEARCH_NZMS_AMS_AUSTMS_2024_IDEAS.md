# ECDLP ideas from the 2024 NZMS–AMS–AustMS programme: a screen, six closures, one live lever

**Status:** position note and preregistration. **No measurement.** Nothing
here sets `S`, a ratio to a boundary, or a class chip. The scoreboard and the
leaderboard are unchanged because no result lands (AGENTS.md §7a applies to
the outcome PR of §5, not to this one).
**Date:** 2026-10-03.
**Input:** the *Schedule and Abstracts* booklet of the 2024 NZMS, AMS and
AustMS Joint Conference (Auckland, 8–13 December 2024), supplied as text in
the session that wrote this note. Talks are cited by speaker, title and
session. No abstract is reproduced.
**Companions:** the R1–R5 admissibility test in
[`RESEARCH_REPRESENTATION_STRUCTURE.md`](RESEARCH_REPRESENTATION_STRUCTURE.md) §3;
[`RESEARCH_HAMMING_IDEAL_PDP.md`](../index-calculus/RESEARCH_HAMMING_IDEAL_PDP.md);
[`RESEARCH_FACTOR_BASE_SHAPE_SEARCH.md`](../index-calculus/RESEARCH_FACTOR_BASE_SHAPE_SEARCH.md);
[`RESEARCH_FACTOR_BASE_SOLVE_COST.md`](../index-calculus/RESEARCH_FACTOR_BASE_SOLVE_COST.md);
the [base-window screen protocol](../ecc2k130/base_window_screen_20261001/PROTOCOL.md);
[`RESEARCH_ECC2K130_RR_SOLVER_PANEL.md`](../ecc2k130/RESEARCH_ECC2K130_RR_SOLVER_PANEL.md);
[`RESEARCH_EDS_RESIDUE.md`](../cm-isogeny/RESEARCH_EDS_RESIDUE.md);
[`RESEARCH_ECC2K130_HYPERELLIPTIC.md`](../ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md);
[`RESEARCH_ECC2K130_EXTENSION.md`](../ecc2k130/RESEARCH_ECC2K130_EXTENSION.md);
[`RESEARCH_FFD_PROOF_COMPLEXITY.md`](../index-calculus/RESEARCH_FFD_PROOF_COMPLEXITY.md).

## TL;DR

- All 43 sessions were read. The talks whose mathematics touches an ECDLP
  mechanism this repository uses (about 40 talks, §3) were each turned into
  a concrete idea. Each idea was checked against the R1–R5 test and against
  what is already measured here.
- **No talk supplies a new licensing datum (R1).** That is what
  `RESEARCH_REPRESENTATION_STRUCTURE.md` §7 predicts. R1 is the scarce
  resource, and a mathematics programme does not publish one for a fixed
  curve.
- Six closures rest on derivations not previously in this repository (§4):
  1. Pairings that respect CM cannot lower the embedding degree.
  2. The zero lattice of a `Z[τ]`-indexed elliptic net has short vectors of
     norm `≪ r^{1/2}`, but finding them is a two-dimensional baby-step
     giant-step, so it stays generic.
  3. Point counts on decomposition varieties either do not depend on `R` or
     are the decomposition problem.
  4. Character-sum certificates of decomposition counts exist only when
     `|F| > 2^{n(k+1)/(2k)} > √r`. Every certified regime costs more than
     rho.
  5. Noncommutative (quantum-strategy) relaxations refute only a subset of
     what the repository's linear certificate already refutes.
  6. Rotation-invariant matroids that are far from uniform require a Singer
     degree `n = (q^k−1)/(q−1)`. Among the repository's primes that is 31
     and 127, not 41, 53, 83 or 131.
- **One live lever (§5): support-complex factor bases.**
  - In a normal basis, any rotation-invariant simplicial complex `Δ` on
    `Z/n` gives a Frobenius-stable factor base. Its membership ideal is the
    Stanley–Reisner ideal of `Δ`, which as clauses is an anti-monotone CNF.
  - The Hamming base is the uniform matroid `U_{w,n}`. The window and
    τ-closure bases are complexes homotopy-equivalent to a circle. The
    repository has measured those two families separately, never at matched
    size under one encoding.
  - The Stanley–Reisner view puts both families on one axis. Cohen–Macaulay
    status, regularity and face dimension become candidate predictors of
    oracle cost.
  - Preregistered stage diagnostic `H_Δ`, on a ladder of degrees where 2 is
    a primitive root (13, 19, 37, 53). Its stop rule is tied to the
    shape-blind null object.
- The expected value is modest, and it is stated as such. A positive `H_Δ` is
  solver engineering until an oracle beats the pair-table null at matched
  size. Only then does it owe the full one-target pipeline and the m = 83
  gate.

---

## 1. The boundary, stated before anything else

Both boundaries come from the leaderboard's Table A (one target, Koblitz;
[`docs/ic/LEADERBOARD.md`](../../../docs/ic/LEADERBOARD.md), ledger §23, as of
`87055a99f`). Units and reference are as defined there.

| quantity | `icv1-f2m41-tm2308219-7f48b14a` | `icv1-f2m53-tm56619371-dac20a85` |
|:--|--:|--:|
| floor: `√(π/4n)` rho steps per `√r`, `A = 2n` | 0.138 | 0.122 |
| reference: one-target strong rho, measured `S` | 0.694 | 0.503 |
| best index calculus, measured `S` | 5.58 | 7.4 |
| best index calculus × reference | 8.04× | 14.7× |
| best index calculus × floor | 40.3× | 60.8× |

At the other degrees named below the floor is 0.159 (n = 31), 0.0973
(n = 83) and 0.0774 (n = 131, ECC2K-130). An idea in this note is an
**advance** only if it lowers the IC-to-floor ratio of the whole one-target
method. Everything §5 can produce on its own is a stage diagnostic, and
leaves `S`, end-to-end cost and speedup unset.

## 2. How the programme was read

A talk was kept if its mathematics bears on one of these: group or curve
structure, CM and endomorphisms, pairings, isogenies and covers, point
counting, polynomial-system solving and commutative algebra, SAT and proof
complexity, factorisation theory, collision and random-walk statistics, or
combinatorial structure on `F_2`-vectors. A kept talk became an idea stated
as a mechanism: what it would compute, from what input, landing where.

The mechanism was then run through R1–R5:

- R1: a consumable published object.
- R2: a homomorphism into another category.
- R3: a non-generic algorithm in the target.
- R4: a payload not already polynomial-time.
- R5: available at cryptographic parameters. R5a fails when the handle
  exists only in a disjoint toy regime. R5b fails when it exists at real
  parameters but lands somewhere more expensive.

Ideas internal to index calculus inherit R1–R4 from Weil descent to `F_2`,
licensed by the subfield. For them the open item is R5, which here is the
first-fall-degree question.

Most of the programme has no such mapping and is not listed: fluids,
geometric PDE, mathematics education, outreach, inverse problems, low-
dimensional topology, operator algebras and most of dynamics. Abstracts are
paraphrased in this note, never quoted.

## 3. The screen

"Closed" means a derivation or an existing measurement shows the idea cannot
beat the boundary of §1. "Live" means it can be falsified cheaply, and §5
preregisters how.

| # | talk(s) — speaker, title, session | the ECDLP idea it suggests | R-test | verdict | why / where covered |
|--:|:--|:--|:--|:--|:--|
| 1 | K. Stange, *Respecting CM on elliptic curves: sesquilinear pairings, elliptic nets, biextensions* (Computational Number Theory, Thu) | a `Z[τ]`-sesquilinear Tate pairing on `E(F_{2^m})[r]` as a MOV transfer for the `E_0` family | R1 none; R5b | closed | §4.1: any such pairing is trivial on the rational `r`-torsion or lands in `F_{q^k}` with the same `k` (`k > 10^7` for ECC2K-130, `RESEARCH_ECC2K130_EXTENSION.md`) |
| 2 | Stange (as row 1) | relations from the zero set of a `Z[τ]`-indexed elliptic net | R2 fails | closed | §4.2: two-dimensional BSGS in `Z[τ]`; computing net values is ECDLP-equivalent (Lauter–Stange; `RESEARCH_EDS_RESIDUE.md` §2.2) |
| 3 | M. Venkatesh (with N. Saxena), *Counting points on surfaces in polynomial time* (Arithmetic Geometry, Fri); M. Kyng, *Computing zeta functions of algebraic curves using Harvey's trace formula* (Computational NT, Wed) | count points on decomposition varieties to decide or estimate relations | R4 / R5a | closed | §4.3 |
| 4 | L. Long, *The Explicit Hypergeometric-Modularity Method*; A. Salerno, *Hypergeometric motives and invertible K3 surface pencils* (Arithmetic Geometry, Wed); K. Kedlaya, *Towards a database of hypergeometric L-functions* (Computational NT, Thu) | character sums or traces as relation counters | R4 (payload is a trace); R5a | closed | §4.4: certified only above `√r` |
| 5 | B. Helton, *Perfect quantum strategies for XOR games* (Functional Analysis and Operator Algebras, Thu) | a noncommutative relaxation of the descended system as a NO-filter for targets that do not decompose | R3 weaker than the existing certificate | closed | §4.5 |
| 6 | A. Kumar, *Generalized Hamming weights and symbolic powers of Stanley–Reisner ideals of matroids* (Computations in AG/CA, Thu); G. G. Smith, *Hodge theory for modular matroids* (Algebraic Combinatorics, Mon); T. Abe, *Solomon–Terao polynomial and Castelnuovo–Mumford regularity of hyperplane arrangements* (Algebraic Combinatorics, Wed); N. Abdallah (with H. Schenck), *Nets in the projective plane and Alexander duality* (Algebraic Combinatorics, Thu); H. Bao (with X. He), *Acyclic matchings on Bruhat intervals…* (Algebraic Combinatorics, Thu); K. Knudson, *Discrete Morse theory on ΩS²* (Applied Topology, Fri) | Frobenius-stable factor bases from rotation-invariant complexes, with Stanley–Reisner membership ideals; commutative-algebra invariants as predictors of oracle cost; discrete Morse theory to compute them | inherits R1–R4 from Weil descent to `F_2`; R5 open (FFD) | **live** | §5; matroid branch narrowed by §4.6 |
| 7 | S. Frengley, *On the geometry of the Humbert surface of square discriminant* (Arithmetic Geometry, Fri); J. Booher, E. Howe, A. Sutherland, F. Voloch, *Doubly isogenous curves of genus two with a rational action of D6* (Computational NT, Thu) | glue `E` into a genus-2 Jacobian through a degree-`N` cover | R5b | closed | Same field `F_{2^131}`: Gaudry's index calculus on genus-2 Jacobians costs `q^{2−2/g} = 2^{131}`, above rho. Covered by `RESEARCH_ECC2K130_HYPERELLIPTIC.md` and `RESEARCH_REPRESENTATION_STRUCTURE.md` §9.5 |
| 8 | D. Perrin (with J. F. Voloch), *Ordinary Isogeny Graphs with Level Structure*; M. Chen (with C. Petit), *Computing the endomorphism ring of supersingular elliptic curve from a full rank suborder* (Computational NT) | walk the volcano, with level structure, to a weaker curve, or exploit the endomorphism ring | R4 | closed | `End(E_0) ⊇ Z[τ]`, the maximal order of `Q(√−7)`, class number 1, and it is public. The class has been searched (`RESEARCH_ISOGENY_CLASS_SEARCH.md`, `RESEARCH_ISOGENY_DEGREE_SEARCH.md`) |
| 9 | Y. Qiao, *Isomorphism problems for some algebraic structures: algorithms, complexity, and cryptography* (Groups, Actions and Computations, Thu) | MinRank / tensor structure of the descended bilinear map `B_V : V × V → F_{2^n}` | attack: known | closed as an attack; tooling only | The HFE link to Weil-descent ECDLP is Huang–Kosters–Yeo (CRYPTO 2015). The closure bound on `V·V` is measured (Kneser, shape-search note). Isotopy invariants of `B_V` could deduplicate future shape sweeps |
| 10 | P. Diaconis, *Computational Polya Theory revisited* (Groups, Actions and Computations, Fri) | Burnside-process sampling or counting of Frobenius orbit classes | R2/R3 fail | closed | It samples orbits and computes no logarithm. Orbit accounting is already exact at prime `n` |
| 11 | J. Coykendall, *Factorization in Monoids and Domains*; F. Gotti, *Divisibility and ascending chains of principal ideals*; H. S. Choi, *Computing elasticity of certain integral domains* (50 years of Comm. Algebra) | smoothness in the exponent ring `Z[τ]` | R2 fails | closed | `Z[τ]` is a PID, but the exponent is not visible in the encoding. The group law carries no multiplicative structure in `x` |
| 12 | T. Ngo Dac, *On multiple zeta values in positive characteristic*; *On special functions and twisted L-series* (Arithmetic Geometry; Special Functions) | a Drinfeld-module-style linear DLP for the `Z[τ]`-module `E(F_{2^m})` | R2 fails | closed | Frobenius is additive on `x`; point addition is not. The Drinfeld DLP is linear because the whole module action is by additive polynomials (Scanlon 2001) |
| 13 | A. Melnikov, *Computable duality theory* (Computability) | compute the character group of `E(F_q)[r]` | restatement | closed | A nontrivial character into `μ_r` is a discrete-logarithm oracle |
| 14 | B. Caldwell et al., *The Douglas–Rachford algorithm for inconsistent problems*; R. Luke, *Convergence Theory for Expansive Markov Chains*; N. Dizon et al., *Wasserstein DRO with Piecewise SOS-Convexity*; A. Bagirov and S. Taheri, DC optimisation (Optimisation) | continuous projection or SOS relaxations as decomposition oracles or refuters | R3 fails | closed | Real relaxations of `F_2` parity are weak. SOS needs linear degree to refute random 3-XOR (Grigoriev 2001). Polynomial-calculus degree is in `RESEARCH_FFD_PROOF_COMPLEXITY.md` |
| 15 | S. Palau (with A. Blancas), *Coalescent point process of branching trees…*; C. Burden (with R. Griffiths), *Coalescence for Feller diffusions*; G. Froyland et al., quenched hitting-time statistics (Probability; Ergodic Theory) | model the merging of distinguished-point trails in parallel rho | no lever | closed | Known (van Oorschot–Wiener 1999). Cross-target sharing is excluded by the one-target rule |
| 16 | D. Maclagan, *Toric Bertini theorems in arbitrary characteristic*; F. Sottile, *The Critical Point Degree of a Bloch Variety*; A. Deopurkar, *How twisty is that orbit?* (Computations in AG/CA) | Newton-polytope (BKK) root counts and coordinate choice for summation systems | R3 fails | closed | The field equations make the solution set Boolean, so BKK counts `F̄_2`-solutions of a different system. Model choice is covered by `RESEARCH_EXOTIC_COORDINATES.md` |
| 17 | B. Creutz, *Quartic del Pezzo surfaces without quadratic points* (Computational NT, Mon); B. Viray, *Number fields generated by points in linear systems on curves* (Arithmetic Geometry, Tue) | splitting types of fibres in Riemann–Roch decompositions; closed points of low degree | no lever | closed | The complete-splitting fraction `1/k!` is the counting factor (§4.3), not a tunable. The Riemann–Roch encoding is measured (RR panel §2, §8) |
| 18 | F. Voloch, *Irreducibility of curves over finite fields* (Arithmetic Geometry, Fri) | decomposition as list decoding of elliptic AG codes | R3 fails | closed | Decompositions of `R` over `F` are exactly the minimum-weight codewords of the elliptic code `C_L(F, (k+1)O − R)`: dimension `k`, weight `|F| − k`. That is beyond the Guruswami–Sudan radius `|F| − √(k|F|)`; decoding elliptic codes is hard by reduction from subset sum on `E` (Cheng 2008) |
| 19 | H. Van Maldeghem, *Weyl substructures, polar kangaroos and uniclass automorphisms of spherical buildings* (Groups and Geometry, Fri) | — | — | name collision | Unrelated to Pollard's kangaroo |

## 4. The derivations behind the new closures

### 4.1 Pairings that respect CM cannot lower the embedding degree

Let `q = 2^m`, `G = E(F_q)[r]` with `r ∤ q − 1`, and let `e : E[r] × E[r] → μ_r`
be any bilinear pairing that commutes with `Gal(F̄_q/F_q)`. Stange's
sesquilinearity over `R = Z[τ]` adds conditions, so it can only restrict this
class.

For `P, Q ∈ G` the `q`-Frobenius fixes both points. So
`e(P, Q) = e(πP, πQ) = e(P, Q)^q`, which gives `e(P, Q)^{q−1} = 1`. Together
with `e(P, Q)^r = 1` and `gcd(r, q − 1) = 1`, the pairing is trivial on
`G × G`.

A nontrivial pairing involving `G` must therefore pair it with the other
eigenline of Frobenius, the one with eigenvalue `q`. That eigenline is
defined over `F_{q^k}`, `k = ord_r(q)`, and the pairing lands in
`μ_r ⊂ F_{q^k}^*`. That is the MOV/Frey–Rück transfer with the same `k`. For
ECC2K-130, `RESEARCH_ECC2K130_EXTENSION.md` bounds `k > 10^7`. The CM
structure refines which pairings exist, not where they land, so the verdict
is R5b.

### 4.2 The zero lattice of a `Z[τ]`-net is a two-dimensional BSGS

A net indexed by `R^n` has `W(v) = 0` iff `Σ v_i P_i = O`. Take `P` and
`Q = [k]P`, and let `π` be the prime of `Z[τ]` above `r`
(`Z[τ]/π ≅ F_r`, with `τ ↦ λ`). The zero set in `R²` is

```text
    Λ = { (α, β) ∈ Z[τ]² : α + kβ ≡ 0  (mod π) },
```

of `Z`-rank 4 and index `r`. Minkowski's theorem gives a nonzero
`(α, β) ∈ Λ` with `N(α), N(β) ≪ r^{1/2}`, and such a pair yields
`k ≡ −α/β (mod π)`.

To find one, let `α` and `β` range over the `≈ r^{1/2}` elements of norm at
most `r^{1/2}`. The condition `αP = −βQ` is then a collision between two
lists of that size: baby-step giant-step in `Z[τ]`, `Θ(r^{1/2})` group
operations. Multiplying by `τ` doubles the norm, so Frobenius classes do not
partition the norm ball and do no better than rho's `√(2m)` fold.

Only the group law and one known endomorphism are used, so generic bounds
with a known automorphism group of order `2m` apply. R2 fails: there is no
change of category. Evaluating net values instead of searching is
ECDLP-equivalent (Lauter–Stange 2008).

### 4.3 Point counts on decomposition varieties carry nothing

Without a factor-base constraint,
`#{(P_1, …, P_k) ∈ E(F_q)^k : Σ P_i = R} = |E(F_q)|^{k−1}` for every `R`. The
decomposition variety in `x`-coordinates counts this up to sign and
permutation symmetry, so it carries no `R`-dependence beyond degenerate
2-torsion and sign cases.

The constraint `x_i ∈ V` can be imposed algebraically as `L_V(x_i) = 0`, where
`L_V` is the subspace polynomial of degree `2^ℓ`. Then the variety is
zero-dimensional, and its `F_q`-point count is the number of decompositions:
computing it is the decomposition problem.

The point-counting algorithms in row 3 run in time polynomial in `log q` for
fixed geometry over growing `q`. Weil descent instead fixes `q = 2` and grows
the dimension with `n`, which is outside their regime (R5a).

### 4.4 Character-sum certificates exist only above `√r`

Take `k` summands and `V` of dimension `ℓ`. Write the indicator of `V` with
additive characters trivial on `V`:
`[x ∈ V] = 2^{−(n−ℓ)} Σ_{a ∈ V^⊥} ψ(ax)`, with `ψ(y) = (−1)^{Tr y}` and `V^⊥`
the trace-dual. Then

```text
    N_R = 2^{−k(n−ℓ)} Σ_{a ∈ (V^⊥)^k}  Σ_{x ∈ X_R(F_q)} ψ(Σ a_i x_i),
```

where `X_R` is the `(k−1)`-dimensional decomposition variety. The trivial
character contributes `|X_R(F_q)| ≈ q^{k−1}`, so the main term is
`2^{kℓ−n}`.

Grant square-root cancellation `≤ C q^{(k−1)/2}` to every nontrivial sum.
That is generous: `C` grows with the degrees, and degenerate sums are worse.
The total error is then at most `C·2^{(k−1)n/2}`.

A certificate `N_R > 0` needs `kℓ − n > (k−1)n/2`, i.e.
`ℓ > n(k+1)/(2k)`. So `|F| ≈ 2^ℓ > 2^{n(k+1)/(2k)} > 2^{n/2} ≈ √r`. For
`k = 2` the condition is `ℓ > 3n/4`, and it approaches `n/2` only as `k → ∞`.
In every certified regime the factor base alone is larger than rho's whole
cost.

Rigorous yield formulas are a proof device, not a lever (R5a). This covers
the hypergeometric and trace-formula sums of row 4 as well.

### 4.5 Noncommutative relaxations are weaker than the linear certificate

A classical solution is a commuting, scalar strategy. So "no perfect quantum
strategy" implies "no classical solution". That is a sound refutation, but it
refutes only a subset of the infeasible instances.

Apply it to the XOR system obtained by linearising the descended equations
(each monomial a fresh variable). The classical decision on that system is
Gaussian elimination, which is the `F_2`-linear NO-certificate of the RR
panel §9. That certificate stops existing once the `S₄` value set spans
`F_2^131`, at `d = 7`; decompositions only begin at `d = 45`. The quantum
relaxation's horizon is no later.

For linear systems of unbounded arity, deciding perfect commuting-operator
strategies is undecidable (Slofstra). The polynomial-time decision in row 5's
line of work covers 3XOR games. On the unlinearised system the relaxation is
never stronger than the classical question, which is the decomposition
problem itself.

### 4.6 Rotation-invariant matroids at prime `n`

A Frobenius-stable, down-closed support set is a rotation-invariant
simplicial complex (§5.1). Matroid complexes are Cohen–Macaulay, which made
them the attractive sub-case. How much room is there?

- **`F_2`-linear.** Cycle spaces of binary matroids are binary codes. The
  rotation-invariant ones are cyclic codes, i.e. ideals of
  `F_2[x]/(x^n − 1)`. When `ord_n(2) = n − 1` only the trivial ones exist.
  This is AGENTS.md §8b restated.
- **Representable by an orbit `{M^i v_0}`.** Over a field containing `ζ_n`
  this is the column matroid of a cyclic code. A `w`-set is dependent iff a
  generalised-Vandermonde (Schur-polynomial) minor vanishes. *Heuristically*
  that happens for a fraction `≲ 1/|field| ≤ 1/(n+1)` of `w`-sets per
  condition.
  - A rotation-invariant matroid is far from uniform only when the orbit
    fills a large part of a projective space: `n = |PG(k−1, q)| =
    (q^k−1)/(q−1)`, with a Singer cycle acting as the rotation.
  - Checked for prime powers `q ≤ 200`, `k ≥ 2`: 31 is Singer twice,
    `|PG(4,2)| = |PG(2,5)|`, and 127 = `|PG(6,2)|`. 41, 53, 83 and 131 are
    not.
  - The least prime power `≡ 1 (mod 131)` is 263.
- **Non-representable paving matroids.** For sparse paving matroids, the
  circuit-hyperplanes form a packing: at most `C(n, w−1)/w` of the `C(n, w)`
  `w`-sets, a fraction `≤ 1/(n−w+1)`.
  - Rank-3 paving matroids built from translates of a Sidon block of size
    `≈ √n` remove a few percent of triples: `n·C(√n, 3)` of `C(n, 3)`, about
    6% at n = 131.

Consequence: at 53, 83 and 131 a matroid factor base is uniform or close to
it, so it is a Hamming base. At n = 31 the Singer identification
`Z/31 ≅ F_32^*` makes `PG(4,2)` a rotation-invariant modular matroid, and its
base differs substantially from the Hamming ball. It has 114,205 independent
supports against 206,368 of weight ≤ 5.

That matroid exists because `ord_31(2) = 5`, the split cyclotomic block of
§8b. Any n = 31 gain from it must be labelled as split-block-dependent, and
cannot be cited at 53, 83 or 131. The middle bullet is a heuristic for
non-Singer degrees, not a theorem.

## 5. The live lever: support-complex factor bases (preregistration of `H_Δ`)

### 5.1 The object

Fix a normal basis `{α^{2^i}}_{i ∈ Z/n}` of `F_{2^n}` and write `supp(x) ⊂ Z/n`.
Frobenius is a rotation. For a rotation-invariant simplicial complex `Δ` on
`Z/n`, let

```text
    X_Δ = { x ≠ 0 : supp(x) ∈ Δ, x the abscissa of a point of E_0 },
    F_Δ = the points over X_Δ,
```

with the relation and cofactor convention of
`research/hamming_ideal_pdp_20260930/PROTOCOL.md`. `F_Δ` is Frobenius-stable
at every `n`.

Membership is `I_Δ + (x_i² + x_i)`, where `I_Δ` is the Stanley–Reisner ideal,
generated by `x^σ` for the minimal non-faces `σ`. As clauses, membership is
one clause `∨_{i∈σ} ¬x_i` per minimal non-face: an anti-monotone CNF that the
repository's CDCL+XOR solver accepts directly. Conversely, every down-closed
rotation-invariant support set is such a `Δ`. Sets that are not down-closed
(trace forms, arbitrary unions of orbits) are outside this family and belong
to the shape-search note.

| family | `Δ` | Cohen–Macaulay? | max face | exists at 53 / 83 / 131 | status here |
|:--|:--|:--|--:|:--|:--|
| Hamming `H_w` | the `(w−1)`-skeleton of the simplex, i.e. the independence complex of `U_{w,n}` | yes (matroid); `reg k[Δ] = w` | `w` | yes | measured, n = 7–19 (Hamming note) |
| window `W_L` | subsets of cyclic intervals of length `L`; `F_Δ` is the τ-closure of one coordinate subspace | no, for `3 ≤ L`, `(L−1)/n < 1/3` (homotopy-equivalent to `S¹`) | `L` | yes | single, unclosed windows measured (window screen, n = 41) |
| window skeleton `W_L^{≤w}` | faces of `W_L` of size `≤ w` | no, for `w ≥ 3` | `w` | yes | new |
| Sidon-block paving `P_{B,w}` | rank-`w` paving matroid: all sets of size `≤ w−1`, plus the `w`-sets not inside a translate of a Sidon block `B` | yes (matroid) | `w` | yes | new; close to Hamming by §4.6 |
| spread closure `⟨A⟩` | generated by the translates of a Sidon set `A` | no | `|A|` | yes | new |
| Singer `PG(4,2)`, `PG(2,5)` | independent sets under `Z/31 ≅ F_32^*` (resp. the cyclic plane of order 5) | yes (modular matroids) | 5 / 3 | **no — n = 31 only** | new; split-block arm |

The window complex is homotopy-equivalent to a circle by the nerve theorem:
the facets are simplices, and their intersections are intervals while
`3L ≤ n`. The nerve is then the clique complex of the cycle power
`C_n^{L−1}`, which is a circle for `(L−1)/n < 1/3` (Adamaszek 2013). Its
`H̃_1 ≠ 0` lies below the top dimension, so Reisner's criterion fails at the
empty face.

The window skeleton keeps that `H̃_1` for `w ≥ 3`. Hamming skeleta are
Cohen–Macaulay with nonzero top homology, hence `reg = w`. Regularity and
depth of the other families are computed before any oracle run (§5.4) and
are not asserted here.

### 5.2 Why it might matter, and why it probably does not

**For.** The oracle sees the factor base only through its membership ideal.
At matched size the families differ in invariants that govern Gröbner and
proof-system behaviour:

- Castelnuovo–Mumford regularity bounds the solving degree (Caminata–Gorla).
- Cohen–Macaulayness of `Δ` is equivalent to a linear resolution of its
  Alexander dual (Eagon–Reiner), the duality in row 6's nets talk.

The two Frobenius-stable shapes this repository has sit at opposite corners:
windows have few large faces and are not Cohen–Macaulay; Hamming balls have
many small faces and are. They have never been compared at matched size
under one encoding.

**Against.** Yield is fixed by `|X_Δ|`, by counting. The shape-blind null,
pair-table enumeration, already beats every algebraic oracle measured on
ECC2K-130 (RR panel §7–§8). The Hamming oracle's calls per target grew as
`|F_w|^{1.9–2.5}`. The prior is therefore that shape reorders solver costs
without beating enumeration. The stop rule of §5.6 is built around exactly
that prior.

### 5.3 Hypotheses

- **H0 (size only).** At matched `log₂|X_Δ|`, the per-target oracle cost does
  not depend on the family beyond A/A noise.
- **H_CM.** At fixed max face size `w`, the Cohen–Macaulay families
  (`H_w`, `P_{B,w}`) and the non-Cohen–Macaulay family (`W_L^{≤w}`) differ,
  with the same sign, at ≥ 3 of the 4 faithful degrees.
- **H_dim.** At fixed non-Cohen–Macaulay status (`W_L` against `W_L^{≤w}`),
  max face size moves cost, with the same sign at ≥ 3 of 4 faithful degrees.
- **H_reg.** Across all (family, degree) cells, Spearman's correlation
  between `reg k[Δ]` and the matched-size cost effect `E_f(n)` of §5.6 is
  ≥ 0.6.

Every predictor (Cohen–Macaulay status by Reisner's criterion, depth,
regularity by Hochster's formula, max face size, minimal-non-face count) is
computed and committed **before** any oracle run.

### 5.4 Frozen design (hash-locked in a source-lock PR before any target exists)

- **Curves.** The `E_0` family in a normal basis.
  - Faithful ladder, where 2 is a primitive root mod `n` so the nontrivial
    cyclotomic block is irreducible as at 131:
    `icv1-f2m13-t181-515ee569`, `icv1-f2m19-t797-b6cf2467`,
    `icv1-f2m37-tm534059-32aad96b`, `icv1-f2m53-tm56619371-dac20a85`.
  - Split-block control: `icv1-f2m31-tm90707-c95f16f5` (`ord_31(2) = 5`),
    for the Singer arms and as a disclosed non-faithful comparison. Nothing
    measured there supports a statement at 53, 83 or 131.
- **Arity** `k = 3` (`S₄`). `k = 2` is excluded because the pair table
  answers it without any oracle.
- **Sizes.**
  - For each family and degree, take the four parameter values whose
    `log₂|X_Δ|` lie nearest `n/3`: two below, two above, ties broken toward
    the smaller parameter. At that scale a constant fraction of targets
    decomposes when `k = 3`.
  - Each family's cost curve is interpolated in `log₂|X_Δ|` at exactly `n/3`.
    No family's size is chosen after any cost is seen.
- **Targets.** 256 per degree, the same for every arm.
  - `Q = [d]G` with `d = SHA-256("support-complex-20261003-v1|n|j") mod r`,
    read big-endian, with `n` and `j` in plain decimal.
  - Reject zero and repeated signed-Frobenius orbits.
  - Scalars are stored as replay fixtures that no oracle reads.
- **Oracle O1.** The repository's native CDCL+XOR solver
  (`src/cryptanalysis/sat.rs`), on one encoding used for every family:
  - all `n` normal-basis bits of each of the three abscissae as variables;
  - the descended `S₄`, built on the existing native encoder
    (`src/cryptanalysis/semaev_sat.rs`) and extended as needed in the
    source-lock PR;
  - membership as the anti-monotone CNF of `Δ`.
  - For `H_w`, a cardinality-network encoding runs as well, as an encoding
    control. The Stanley–Reisner clause count `C(n, w+1)` is why.
- **Oracle O2 (null).** Native exhaustive enumeration with a pair table over
  `F_Δ`. It is shape-blind: its cost is a function of `|F_Δ|` alone.
- **A/A.** The `H_w` arm runs twice, with two fixed solver seeds.
- **Caps.** `10^7` conflicts per target and 60 CPU-minutes per cell.
  Censored targets stay in the record as right-censored. A cell with more
  than 10% censored targets is `CENSORED`, never a win.
- **Verification.**
  - Every reported decomposition is checked with an independent group law in
    the order-`r` subgroup.
  - Every UNSAT answer at `n ≤ 37` is checked by O2's exhaustive enumeration;
    at n = 53, by a deterministic 32-target subsample.
- **Native only.** No Python in any execution path (AGENTS.md). That includes
  the predictor computation: homology of links by discrete-Morse reduction
  (row 6), and regularity by Hochster's formula.

### 5.5 Cost accounting

- **Primary metric: deterministic counts.** Solver propagations and conflicts
  for O1; field multiplications and table probes for O2. CPU time comes from
  `tools/isolated_bench.py` and is a practicality note, never the headline
  (§6, §10).
- **One-time setup**, reported separately from per-target cost and never
  dropped: complex construction, predictor computation, encoding generation.
- **No conversion** between O1's and O2's counts is claimed without a
  measured conversion factor. Their ratio is reported in CPU, and in counts
  only once such a factor is pinned.
- **Classification.** This is a stage diagnostic: `S`, end-to-end cost and
  speedup stay unset.

### 5.6 Outcomes and stop rule (frozen)

Let `E_f(n) = log₂(cost_f / cost_{H_w})` at `log₂|X_Δ| = n/3`, using the
censoring-aware median per-target cost. The noise band is twice the A/A 95%
half-width.

| outcome | condition |
|:--|:--|
| `SIZE_ONLY` | every non-Hamming `E_f(n)` lies within the noise band at ≥ 3 of 4 faithful degrees. H0 stands, and the Stanley–Reisner family is closed as an oracle lever |
| `CM_EFFECT` / `DIM_EFFECT` | H_CM / H_dim holds as stated in §5.3 |
| `REG_PREDICTS` / `REG_FALSIFIED` | Spearman ≥ 0.6 / < 0.3 for H_reg |
| `LEVER_CLOSED` (overrides the rows above) | for every family, `cost_O1 / cost_O2 > 1` at `n/3` at every faithful degree, **and** the fitted slope of `log₂(O1/O2)` in `n` over the four faithful sizes is `≥ 0`, with the fit shown next to its points (§5). A family effect that only reorders solver costs above enumeration is engineering of a losing oracle |
| `PROCEED` | some family has `O1/O2 < 1` at both n = 37 and n = 53, with a negative slope. Only then does a separate protocol price the full one-target pipeline at 41 and 53 (Table A) and run the m = 83 gate (AGENTS.md §8a). Until then nothing reaches the scoreboard |
| `INCONCLUSIVE` | anything else. Censored cells are reported, never imputed |

### 5.7 What would not count

- Choosing a family's sizes, or the comparison point, after costs are seen.
- Dropping encoding or setup cost.
- Reporting O1 against a different oracle but not against O2.
- Citing an n = 31 Singer result at 53, 83 or 131.
- An exponent fit over fewer than four sizes.
- Any multi-target or batch figure as the primary comparison.
- Calling `CM_EFFECT`, `DIM_EFFECT` or `REG_PREDICTS` an advance. They are at
  most engineering of the decomposition stage.

## 6. Watch list: what would reopen a closure

- §4.1: a computable pairing on `E(F_q)[r]` that respects CM but is not
  Galois-equivariant. That would violate §4.1's hypothesis, and no such
  object is known.
- §4.4: a structured family of nontrivial character sums on decomposition
  varieties that cancels far beyond square root. That would move the
  certified regime below `√r`.
- §4.6: a rotation-invariant matroid on `Z/131` that is far from uniform and
  has an efficiently decidable independence oracle, at a non-Singer degree.
- Row 7: a cover landing in an abelian variety over `F_2` rather than
  `F_{2^131}` with small-genus index calculus. For coefficients in `F_2` the
  GHS construction degenerates: its magic number is 1 and the cover is `E_0`
  itself. Here `Res_{F_{2^131}/F_2} E_0` is isogenous to `E_0 × A`, with `A`
  the 130-dimensional trace-zero variety. No small-genus curve covering `A`
  is known.

## 7. Provenance

- The programme text was supplied in-session and is not stored here. The
  screen depends only on speaker, title, session and the mathematics.
- The arithmetic facts were computed with shell integer arithmetic for this
  note: the orders of 2 mod 31, 41, 53, 83, 127 and 131 (5, 20, 52, 82, 7,
  130), `2^65 ≡ −1 (mod 131)`, the Singer-number check, 263, and the
  `PG(4,2)` counts. The source-lock PR of §5.4 recomputes them natively.
- Boundary numbers cite `docs/ic/LEADERBOARD.md` Table A at `87055a99f`.

## 8. References

- M. Adamaszek, *Clique complexes and graph powers*, Israel J. Math. 196 (2013).
- C. Caminata, E. Gorla, *Solving multivariate polynomial systems and an
  invariant from commutative algebra*, WAIFI 2020.
- Q. Cheng, *Hard problems of algebraic geometry codes*, IEEE Trans. Inf.
  Theory 54 (2008).
- J. A. Eagon, V. Reiner, *Resolutions of Stanley–Reisner rings and Alexander
  duality*, J. Pure Appl. Algebra 130 (1998).
- D. Grigoriev, *Linear lower bound on degrees of Positivstellensatz calculus
  proofs for the parity*, Theor. Comput. Sci. 259 (2001).
- M. Hochster, *Cohen–Macaulay rings, combinatorics, and simplicial
  complexes* (1977); G. Reisner, *Cohen–Macaulay quotients of polynomial
  rings*, Adv. Math. 21 (1976).
- M.-D. Huang, M. Kosters, S. L. Yeo, *Last fall degree, HFE, and Weil
  descent attacks on ECDLP*, CRYPTO 2015.
- K. Lauter, K. E. Stange, *The elliptic curve discrete logarithm problem and
  equivalent hard problems for elliptic divisibility sequences*, SAC 2008.
- J. S. Provan, L. J. Billera, *Decompositions of simplicial complexes related
  to diameters of convex polyhedra* (matroid complexes are vertex-
  decomposable), Math. Oper. Res. 5 (1980).
- T. Scanlon, *Public key cryptosystems based on Drinfeld modules are
  insecure*, J. Cryptology 14 (2001).
- W. Slofstra, *Tsirelson's problem and an embedding theorem for groups
  arising from non-local games*, J. Amer. Math. Soc. 33 (2020).
- P. C. van Oorschot, M. J. Wiener, *Parallel collision search with
  cryptanalytic applications*, J. Cryptology 12 (1999).
