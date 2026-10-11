# Index calculus on Koblitz and prime-field curves: a step back

**Date.** 2026-10-10. **Status.** Analysis and calculations over the repo's
measured results. Nothing below is a new measurement. Every number that is a
calculation says so; every number that is a measurement cites its file.
**No security claim about a deployed curve is made or withdrawn here.**

Sources read for this note: `RESEARCH_KOBLITZ_INDEX_CALCULUS.md`,
`RESEARCH_ECC2K130_IC_FEASIBILITY.md`, `RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md`,
`RESEARCH_PRIME_FAST_ECDLP.md`, `research/linearization_reach_20260930/README.md`,
`docs/ic/BOUNDARY_TARGETS.md`, `docs/ic/LEADERBOARD.md`, `docs/bounds/FRONTIER.md`,
and the `research/notes/index-calculus/` notes on GLV-invariant bases, symmetrised
Semaev, factor-base shape, quasi-subfields, exotic coordinates, Hamming ideals,
degree of regularity, fall degrees, descent crossover, cover decomposition, PKM
towers, residual walks, hyperelliptic IC, and the isogeny conductor gap.

---

## 0. Bottom line

1. **Everything the Koblitz ladder runs today is a generic algorithm.** The
   compact-orbit 4-sum pipeline decides "does R decompose" by a table lookup.
   Table lookups and S₃ root solves are group operations plus membership tests,
   so the Corrigan-Gibbs–Kogan bound `S·T² ≳ r/|Aut|` applies to it exactly as
   it applies to rho. With `|Aut| = 2n` for Koblitz, no engineering of 4-, 5- or
   6-sum tables can beat `√(r/2n)` on total work, and the online-after-precompute
   class must be compared against preprocessing rho at the same advice size, not
   against plain rho. The repo's prime-field program states this for `F_p`
   (`RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md` §1); it has not been applied
   to the Koblitz ladder, and when it is (§2 below) the measured rungs sit
   `2^25`–`2^27` *above* the generic line.
2. **The only non-generic step available on any elliptic curve is an algebraic
   decomposition over a factor base whose membership condition has degree far
   below its size.** Binary fields have such bases (any `F₂`-subspace, via Weil
   restriction). Prime fields have none (`RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md`
   §2.3). Prime extension degree `n` removes the Frobenius-invariant ones and the
   quasi-subfield ones (`RESEARCH_QUASI_SUBFIELD.md`), but **not** the
   non-invariant subspace ones. So for ECC2K-130 the question "secure against IC
   or not" reduces to one measurable quantity: the solving degree of
   Weil-descended `S_{m+1}` systems with `m·ℓ ≈ n` (§5).
3. **Symmetry cannot change the exponent; it changes constants, and the repo has
   now measured most of them.** Frobenius gives `n` on relations and `n²` on
   linear algebra (GGMP), the FGHR 2-torsion symmetrisation gives 30–400× on the
   solve, and automorphism folds give `|Aut|` on columns. Each is polynomial in
   `n` or constant. The one symmetry lead with measured upside and no measured
   ceiling is Frobenius-equivariant F4 (Faugère–Svartz character blocks), 31–128×
   on rank at n = 11–23 (`RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md`).
4. **Two cheap, decisive experiments are missing** and both are calculations
   here (§2, §3): the preprocessing-rho control, which I predict erases the
   online-class wins at every rung, and free relations by self-collision in the
   pair-root table, which I predict replaces the 18,400 s rank precompute at
   n = 73 and the whole rank stage at n ≤ 71 and n = 83 with one sort of a table
   that already exists.
5. **For prime fields the map is closed to the level of a falsifier.** The
   first arity that could beat rho is `m = 4`, needing an `S₅` small-root oracle
   at `δ ≈ 1/5` against a measured Coppersmith reach of `1/80`. Nothing short of
   a new algebraic idea moves that; constant-factor engineering there is done.

---

## 1. The three questions that decide "secure or not"

For a curve `E(F_q)` with prime-order subgroup `r` and automorphism group `Aut`
(negation, Frobenius for subfield curves, CM automorphisms), index calculus is a
three-phase algorithm with cost

```
T_total = T_relations + T_linear_algebra + T_descent
T_relations      = (#relations needed) × (attempts per relation) × (cost per attempt)
#relations       ≈ |F| / |Aut|                 (orbit columns)
attempts/relation ≈ 1 / Pr[random point is an m-sum of base points]
cost/attempt     = cost of the decomposition oracle
T_linear_algebra ≈ w · (|F|/|Aut|)²           (sparse, w nonzeros per row)
T_descent        ≈ one decomposition per target
```

Three questions settle the whole picture.

**Q1. Is the decomposition oracle generic?** If it answers "is R an m-sum of
base points" by enumeration, pair tables, direct subtraction, or any membership
test on group elements, then the whole algorithm is a generic-group algorithm
with advice (the table). Shoup gives `T ≥ √(r/|Aut|)` for the uniform case and
Corrigan-Gibbs–Kogan give `S·T² ≳ r/|Aut|` with preprocessing. **Generic IC can
never beat rho on total work and is dominated by Bernstein–Lange preprocessing
rho on online work at equal advice.** This is a theorem, not a heuristic. It
covers: the Koblitz compact-orbit 4-sum, the planned 5- and 6-sum extractions
(unless the S₄ solve carries an algebraic restriction on the unknowns), the
prime-field `subtract` and `mitm` oracles, residual walks, and all preprocessed
variants of these.

**Q2. If the oracle is algebraic, what does it cost?** The only algebraic
oracle known is: restrict the unknown `x_i` to an algebraic set `V` of low
"membership degree" and solve `S_{m+1}(x_1..x_m, x_R) = 0` by Gröbner, XL, SAT or
linearisation. The cost is governed by the solving degree `D` of that system,
not by `|F|`. This is where the exponent could in principle change. Over
`F_{2^n}` any `F₂`-subspace `V` (dimension `ℓ`) works and the Weil restriction
gives `n` Boolean equations in `m·ℓ` unknowns. Over `F_p` no subset of `F_p`
has membership degree below its size (`RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md`
§2.3), and the lattice substitute reaches `δ = 1/6` against the `1/2` needed
(Probes B and C, same file).

**Q3. Does any structure transfer the problem somewhere Q2 is cheaper?**
Weil descent / GHS covers need an intermediate field; isogenies preserve the
field of definition and (for Koblitz, `h(O_K) = 1`) only change `τ`; lifts to
characteristic 0 fail because heights grow quadratically under addition
(xedni). The measured cover win in the repo (JV12 reproduction, 0.003–0.08× rho
at `p ≥ 1009`, `RESEARCH_COVER_DECOMPOSITION_LEDGER.md`) lives entirely on
`E(F_{p⁶})`, which has intermediate fields. Prime `p` and prime `n` are the
*same* obstruction seen from two sides: no intermediate field.

The repo's answers: Koblitz ladder = Q1 generic (not yet stated as such);
Koblitz algebraic solvers = Q2 measured exponential past `n/3`; prime fields =
Q2 closed by lattice reach, Q3 closed by isogeny/cover audits.

---

## 2. Lead A: the Koblitz online-class wins versus preprocessing rho (calculation)

The ledger's `vs_rho` rows compare IC online time to **plain** rho on the same
target and record ratios of 8–2,864 (`docs/ic/BOUNDARY_TARGETS.md` §"Operation
accounting"). The IC online stage uses a materialised table of `K²n` S₃-root
states as advice. The matched control is Bernstein–Lange preprocessing rho with
the same advice: `S` distinguished points from `G`-only walks, walk length
`W ≈ √(N/S)`, online cost `≈ 2√(N/S)` steps, precompute `≈ √(N·S)` steps, with
`N = r/(2n)` after folding negation and Frobenius.

| rung | `N = r/2n` | advice `S` (= IC states) | BL-rho online steps | BL-rho precompute steps | IC online probes (mean) | IC precompute |
|---|---|---|---|---|---|---|
| n=61 | 2^40.3 | 2.2·10⁷ | ≈ 490 | 5.4·10⁹ | ≈ 4·10⁴ | 13 s |
| n=71 | 2^45.2 | 2.6·10⁷ | ≈ 2,500 | 3.2·10¹⁰ | ≈ 4.6·10⁶ | 1,485 s |
| n=73 | 2^49.1 | 2.6·10⁷ | ≈ 9,600 | 1.3·10¹¹ | ≈ 5.7·10⁷ | 18,398 s |
| n=83 (a=1) | 2^45.5 | 2.9·10⁷ | ≈ 2,600 | 3.8·10¹⁰ | 8.8·10⁶ | not charged |

IC numbers are the measured ones from `RESEARCH_ECC2K130_IC_FEASIBILITY.md` §2
and `docs/ic/BOUNDARY_TARGETS.md`; the rho column is the textbook calculation.
At equal advice, preprocessing rho needs **3,000–6,000× fewer online
operations** than the IC extraction at every rung, and its precompute is
seconds to minutes against hours. Stated as the generic bound: at n=73 the IC
pair `(S, T) = (2.6·10⁷, 5.7·10⁷)` has `S·T² = 2^76.2` against `N = 2^49.1`,
so the IC online class sits `2^27` above the line that preprocessing rho sits
on; at n=83 it is `2^25` above.

**Falsifier.** Run `preprocessing_rho.rs` with Frobenius folding on the frozen
n=73 and n=83 fixtures with `S` set to the IC table's byte budget. Prediction:
online wall in milliseconds, so every `vs_rho` online ratio flips to well below
1. The ledger already lists this control as "has not yet been run"
(`docs/ic/BOUNDARY_TARGETS.md` §"Operation accounting",
`docs/ic/PLAN_IC_ACCOUNTING_FIXES_20261007.md` item 5). It should be run before
any further rung is landed, because if the prediction holds the
`single_target_online` class has no content.

---

## 3. Lead B: free relations by self-collision in the pair-root table (calculation)

The rank stage currently finds each relation by targeting `R = [a]G` and
scanning pair-states for a partner, at `r/(4n²K²)` probes per relation
(`RESEARCH_ECC2K130_IC_FEASIBILITY.md` §2). But a relation among factor-base
points does not need a target. The table already holds `x(P_i ± π^k P_j)` for
every pair-state. **Two pair-states with the same x-coordinate (up to Frobenius
rotation) are a 4-term relation** `P_i ± π^k P_j = ±π^s (P_u ± π^t P_v)`, with
the `λ`-power weights the orbit columns already use. Finding all of them is one
sort of the table, zero probes.

Expected count: the table holds `M ≈ n²K²` distinct pair-sum x-classes drawn
from `≈ r/2` x-values, so coincidences `≈ M²/r = n⁴K⁴/r`. With `K` orbit columns
the stage needs `≈ K` independent relations, i.e. `K ≥ (r/n⁴)^{1/3}`.

| rung | `r` | free relations in the K=600 table | relations needed | verdict |
|---|---|---|---|---|
| n=61 | 2^47.2 | ≈ 11,000 | 600 | rank stage free |
| n=71 | 2^52.3 | ≈ 600 | 600 | borderline; K≈650 makes it free |
| n=73 | 2^56.3 | ≈ 40 | 600 | needs K ≈ 1,500 (table 4× larger) |
| n=83 (a=1) | 2^52.9 | ≈ 730 | 600 | rank stage free |
| n=131 | 2^129.7 | 3·10⁻²⁰ | — | needs K ≈ 1.5·10¹⁰, table 4·10²⁴ states: dead |

Caveats before this is believed: the count assumes pair sums land in the
subgroup (cofactor classes cancel) and that x-values of pair sums are close to
uniform; the exact count is a census the existing `relative_pair_stats` backend
can produce in minutes at n ≤ 53. Relations found this way are correlated
(both sides share the table), so rank must be checked, not assumed.

What it does and does not change. It does not change the exponent: the
condition `K³ ≥ r/n⁴` is the same generic time–memory line as before, and at
n=131 the required table is as impossible as the `K = 10⁹–10¹²` scan. It does
collapse the **fully charged precompute** on the rungs the repo actually runs:
at n=83 the rank stage becomes a sort of 3·10⁷ entries instead of 600 scans of
up to 10¹⁰ probes. It also fixes a reporting problem: the ledger's "online
after precompute" rows would then carry a precompute that is honestly small,
and the comparison in §2 becomes the only one left.

**Falsifier.** Sort the n=61 and n=83 K=600 root tables by Frobenius-canonical
x and count collisions; verify each as a group identity; measure rank of the
resulting matrix mod `r`. Prediction: ≥ 500 verified relations at both rungs
and full column rank at n=61.

---

## 4. End-to-end time: the formula the ledger should carry

For any generic-oracle IC on Koblitz, with `c_p` seconds per probe and `c_LA`
seconds per sparse 130-bit multiply-add:

```
T(K) = c_p · r/(n²K)        relation collection (4-sum scan; or n⁴K⁴/r free relations once K³ ≥ r/n⁴)
     + c_LA · w · K²        sparse Wiedemann/Lanczos over Z/rZ, w ≈ 5 nonzeros per row
     + T_descent            one decomposition per target: c_p · r/(4n²K²)
K*   = (c_p r / (2 c_LA w n²))^{1/3}
T*   ≈ 1.9 · (c_p r)^{2/3} (c_LA w)^{1/3} / n^{4/3}
```

The `r^{2/3}` is the whole story: rho is `r^{1/2}`. Plugging the measured
`c_p ≈ 10⁻⁶ s` (1.1 M probes/s single core, `RESEARCH_ECC2K130_IC_FEASIBILITY.md`
§4.6) and an optimistic `c_LA = 10⁻⁷ s` at n = 131:

| quantity | value (calculation) |
|---|---|
| `K*` | 4·10¹¹ columns |
| table `n²K²` | 3·10²⁷ states |
| `T*` | 2.4·10¹⁷ core-seconds ≈ 8·10⁹ core-years |
| rho, measured client | 2^60.9 steps at 6.9·10⁹/s ≈ 10 GPU-years |
| ratio | ≈ 10⁸–10⁹ |

This agrees with the repo's `2·10⁷×` at `K = 10⁹` and shows that letting `K`
float makes it worse, not better, because the linear algebra enters. The same
formula with `w·K²` replaced by the measured `n^{0.68}` Wiedemann scaling
(`RESEARCH_EXTENSION_FIELD_BOUNDARIES.md`) does not change the conclusion.

**What the ledger should record per row** so that end-to-end time is readable
without re-deriving it: `(S, T_online, T_precompute, |Aut|)` and the derived
`log₂(S·T²) − log₂(r/|Aut|)`. A row with that difference at or below 0 is the
first non-generic signal; every current row is at +25 or more.

For prime fields the corresponding formula is already in
`RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md` §2.4: with `T_dec` the oracle
cost, `cost ≈ (m!/2^m) p T_dec / B^{m−1} + m B²`, minimised at
`B = (p T_dec)^{1/(m+1)}`, and `m = 2, 3` never beat `√p` even with a free
oracle. The first arity that can is `m = 4` with `T_dec ≪ p^{1/4}`.

---

## 5. The one quantity that decides binary-curve security

Strip away everything generic and one question remains for every binary curve,
Koblitz or not, prime `n` or not: **how does the solving degree `D` of the
Weil-descended `S_{m+1}` system grow when the base dimension `ℓ` and arity `m`
are chosen so that `m·ℓ ≈ n`** (the regime where a decomposition exists with
probability about 1)?

- The optimistic published answer (Petit–Quisquater 2012, first-fall-degree
  assumption) gives `2^{O(n^{2/3} log n)}` with an estimated crossover against
  rho near `n ≈ 2000`. Even if it held, ECC2K-130 at `n = 131` would be
  unaffected.
- The repo's measurements say it does not hold in the attack regime: first
  fall degree stays at 3 while the solving degree at `m = 3` is 6 at `ℓ ≤ 3`
  and fits `⌈ℓ/2⌉ + 4`; at `m = 4` it fits `ℓ + 4`; the deciding cell
  `(n, ℓ) = (15, 5)` at degree 7–8 was out of memory at 18–50 GB
  (`RESEARCH_DREG_MEASUREMENT.md`). Linearisation resolves `ℓ + b_max ≈ 28–31`
  bits per trial at `n = 131` against the 84–99 bits needed
  (`research/linearization_reach_20260930/README.md`, gap ≥ 56 bits at every
  arity).
- Priced with the measured degrees, `m = 4` at ECC2K-130 is `2^{+15.5}` over
  rho if `D` stayed at 6 and `2^{+112}` if it follows `ℓ + 4`
  (`research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md`). The
  106-bit spread between those two readings **is** the security question.

Recommendation: the deciding measurement is `D` at `m = 4`, `ℓ = 4, 5, 6` on a
machine with enough memory for degree 8, using the equivariant F4 blocks (§6)
to make it affordable, with the FGHR 2-torsion symmetrisation applied (it is
measured to halve the degree in the Gaudry setting and cut 416× at `k = 5`).
Three cells of `D` would settle "levelled at 7" against "follows `ℓ + 4`", and
nothing else in the binary program has that leverage. Stop measuring first
fall degree; it has been 3 everywhere and does not bound `D`.

---

## 6. Symmetry: an exact inventory

The user's question was whether there is a symmetry we are not exploiting.
Here is every symmetry an elliptic curve over a finite field has, what it buys
in each phase, and whether the repo has measured it.

| symmetry | group | relations | linear algebra | rho | solve (algebraic oracle) | measured in repo |
|---|---|---|---|---|---|---|
| negation | Z/2 | ÷2 | ÷4 | ÷√2 | halves degree in `y` (Edwards/FGHR) | yes, everywhere |
| Frobenius (subfield curves) | Z/n | ÷n | ÷n² | ÷√n | Faugère–Svartz character blocks: 31–128× on rank at n = 11–23 | GGMP factors yes; blocks rank-only |
| CM automorphisms (j = 0, 1728) | Z/6, Z/4 | ÷|Aut| | ÷|Aut|² | ÷√|Aut| | Z/3 grading costs 12–70× (does not pay) | yes, `RESEARCH_GLV_INDEX_CALCULUS.md` |
| 2-torsion translation | Z/2 per rational point | none (base fold costs yield) | none | none | FGHR: degree `2^{m-1}` → `2^{m-1}/2`, 30–416× measured | yes, best measured lever on the solve |
| 4-torsion / Klein | (Z/2)² | — | — | — | not `F_q`-valued on subspace bases; Q₄ saturation 0.99× | yes, negative |
| S_m permutation of summands | S_m | ordered tuples ÷m! | — | — | elementary-symmetric rewrite; S₃ 17→≤10 monomials; S₄–S₈ unmeasured | partial |
| D₄ = (Z/2)³ ⋊ S₄ on S₅ (k = 4) | 64 | — | — | — | Bézout 4096 → 512, projected ≈630× on C₄ | **not run** |
| isogenies (horizontal) | class group | none: `h(O_K) = 1` for Koblitz | — | — | ranks identical across the isogeny class | yes, negative |

Two structural facts bound this table. First, every entry is polynomial in `n`
or a constant, so symmetry alone cannot bridge an exponential gap; the generic
lower bound already includes `|Aut|`. Second, the symmetry gains on the *solve*
compound with the base folds (E13 in `RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`,
combined/product ≥ 0.82), so they should be stacked rather than compared.

Connections the notes do not draw:

- **The `n` and `n²` relation/LA savings and the Faugère–Svartz F4 blocks are
  the same symmetry applied in two phases.** The relation matrix over `Z/rZ`
  and the Macaulay matrix over `F₂` are both `F₂[C_n]`-modules. The former is
  decomposed by choosing orbit representatives (done); the latter by the DFT
  over `F_{2^d}`, `d = ord_n(2)`. At `n = 131`, `d = 130`, so the Macaulay matrix
  splits into the trivial character and one block over `F_{2^130}`: the saving
  is the `n`-fold reduction in rows, the arithmetic moves to a 130-bit field.
  That is a real constant-factor win (projected ≈300× dense, 10–100× realistic)
  and the one symmetry item worth the 3–5 days the scoping note estimates.
- **The D₄ symmetry of `S₅` is the FGHR lever at the first arity where
  prime-field IC could tie or beat rho (§4).** It has been computed on paper
  and never run. If `C₄` fell 630× the k = 4 handover would move by ≈37 bits.
  That changes nothing for `F_p` itself (no subfield) but is the cleanest
  remaining test of whether torsion symmetry can ever be more than a constant.
- **The 2-torsion translation on Koblitz is `x ↦ 1/x`** (for `b = 1`, since
  `x(P + T) = √b/x(P)`), which is why the Möbius frame `u = 1/(x+1)` in
  `RESEARCH_EXOTIC_COORDINATES.md` halves the S₄ degrees: it is FGHR
  rediscovered, as that note concludes. There is no further rational
  symmetry on the `x`-line for `E[2](F_{2^n}) = Z/2`.

---

## 7. Connections across the disparate threads

**(a) IC-with-table and preprocessing rho are the same time–memory trade.**
Both store `S` group elements derived from `G` only and answer each target
with `T` online operations. Rho achieves `S·T² ≈ 4N`; the IC table achieves
`S·T² ≈ 2^{25–27} N` (§2). The IC table stores pair-sums indexed for
membership; rho stores distinguished points indexed for collision. The second
is strictly more efficient for the same bytes because a walk covers `W` points
per stored element while a pair-state covers one.

**(b) Relation collection is a k-list sum problem, and Wagner's algorithm does
not apply.** Finding `P₁ + P₂ + P₃ + P₄ = O` over four lists is exactly the
generalised birthday problem. Wagner's `2^{n/3}` method needs a group
homomorphism to a small quotient visible in the representation (low bits of
`Z/2^n`). `E(F_q)` has no such quotient except `E/⟨G⟩` (the cofactor), which is
why §3's sort-and-match at `n⁴K⁴/r` is the best generic rate and why the
problem stays at the plain birthday line. This also explains why "residual
walks" and "large primes" could not help (`RESEARCH_RESIDUAL_WALKS.md`): there
is no quotient to walk in.

**(c) Why smoothness exists for `F_p^*` and not for `E(F_p)`.** NFS lifts
`F_p^*` to `Z` or to a number ring where the size function (absolute value,
norm) is *sub-multiplicative*, so a product of small things is small and
sieving works. The only lift of `E(F_p)` is to `E(Q)` or `E(K)`, where the size
function is the canonical height and `ĥ(P + Q) + ĥ(P − Q) = 2ĥ(P) + 2ĥ(Q)`: the
sum of two small points is not small. That single identity is the obstruction
behind xedni, behind "no factor base with sub-`B` membership degree", and
behind the failure of every lifting idea in the DEFERRED and novel-directions
lists. The Koblitz `Z[τ]`-module structure (`E(F_{2^n}) ≅ Z[τ]/(τ^n − 1)`,
`τ² + τ + 2 = 0`) is a ring with a multiplicative norm, but the point-to-ring
map is the DLP itself, so the norm is not computable from a point. The `Z[τ]`
multiplier-closure idea (A1 in `RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md`) can
only enlarge orbits by `|Z[τ]/(τ^n−1)|`-many small multipliers, a constant.

**(d) Prime `p` and prime `n` are one obstruction.** Every measured exponent
gain in the repo (JV cover on `E(F_{p⁶})`, Gaudry `k = 3, 4` linear-algebra
ceilings, Diem) lives on a field with a proper intermediate subfield, where
"`x ∈ F_q`" is a low-degree condition. `F_p` has no proper subfield;
`F_{2^131}` has `F₂` only, and `F₂`-rationality is 1 bit. The Hamming-weight
(`RESEARCH_HAMMING_IDEAL_PDP.md`), quasi-subfield, and cyclotomic-sparse bases
were all attempts to manufacture a substitute; all measured or derived to the
generic line. This is the structural reason the two deployed families stand
or fall together.

**(e) The 263 coincidence, the conductor `f = 263 · …`, and the trivial Weil
pairing on `⟨G⟩`** (`docs/ic/EXPERIMENTS_ISOGENY_CONDUCTOR_GAP_20261007.md`
E4–E5) are all consequences of `h(O_K) = 1` and `2n + 1 | f` for about half
of the eligible `n`; none is a lever. The volcano has 262 reachable curves,
all Galois conjugates with identical `τ`-structure and identical Macaulay
ranks, so the isogeny program for Koblitz is closed.

---

## 8. Prime fields: what remains

`RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md` is already the right document.
Only three things to add:

1. **The decision rule is sharp enough to stop engineering.** `m = 2, 3` cannot
   beat `√p` with a free oracle (§2.4 there). The first arity that can is
   `m = 4`, needing an `S₅` small-root method at `δ ≈ 1/5` against a measured
   `1/80` reach. Every constant-factor result on `S₃` or `S₄` oracles (PR #1560
   and after) is irrelevant to that question and should be labelled
   engineering.
2. **The two falsifiers are the D₄-symmetrised `S₅` lattice and the
   combined-system solution count.** Run Probe B/C with the `(Z/2)³ ⋊ S₄`
   invariants of `S₅` at `m = 4` on a 40-bit curve; the milestone row is
   `δ ≥ 1/4` at `m = 3`, the target `δ ≈ 1/5` at `m = 4`. Separately, for any
   proposed structured base, count solutions of `{S_{m+1} = 0, L(x_i) = 0}`
   over `F̄_p` at toy `p`; `≥ B^{m−1}` kills it.
3. **Record `S·T²` per ledger row** so preprocessing rho (Probe A, measured on
   the line) is the control for every "online" claim on prime curves too.

---

## 9. What to do, in order, and what to stop

Do:

1. **Preprocessing-rho control at equal advice** on the n=73 and n=83 Koblitz
   fixtures (§2). Half a day. Decides whether the online class exists.
2. **Free-relation census** by sorting the existing K=600 root tables at n=61
   and n=83 (§3). One day. If the prediction holds, the rank precompute
   disappears at every rung the ladder can run, and the ledger's precompute
   charge becomes honest.
3. **Solving-degree cells `(m, ℓ) = (4, 4), (4, 5), (4, 6)`** with FGHR
   symmetrisation and Faugère–Svartz blocks, on a large-memory host (§5).
   This is the only experiment whose outcome bears on ECC2K-130 security
   rather than on accounting.
4. **Equivariant F4 to full F4** (§6). The only symmetry lead with measured
   upside and no measured ceiling; also the enabler for item 3.
5. **D₄-symmetrised `S₅` lattice probe at `m = 4`** (§8). Settles whether
   torsion symmetry can ever be more than a constant on prime fields.

Stop (or label as engineering, not research):

- Further 4-, 5-, 6-sum extraction designs on the compact-orbit base. They are
  generic (Q1) and bounded by `S·T² ≳ r/2n`; the 6-sum yield model in
  `RESEARCH_ECC2K130_IC_FEASIBILITY.md` §4.7 cannot cross that line.
- First-fall-degree measurements. Flat at 3; does not bound `D`.
- Factor-base shape search at `m = 2`; closed by the trace-class artefact.
- Isogeny, volcano, conductor and Kani transport for Koblitz; closed by
  `h(O_K) = 1` and identical ranks across the class.
- `S₃`/`S₄` constant-factor work on prime fields; cannot reach `√p` at any
  constant.

---

## 10. What this note does not establish

It does not prove ECC2K-130 or any prime-field curve secure. It shows that
every method the repo has built is either generic (and therefore provably no
better than rho) or algebraic with a measured solving degree that, at every
measured cell, prices above rho by 15 to 136 bits at `n = 131`. The gap
between "levelled at 7" and "follows `ℓ + 4`" is unmeasured and is the one
place a surprise could still live for binary curves. For prime fields the
corresponding unmeasured place is a small-root method for `S₅` at `δ ≈ 1/5`,
which no published technique approaches.
