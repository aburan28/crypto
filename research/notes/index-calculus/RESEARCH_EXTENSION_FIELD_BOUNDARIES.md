# Rigorous boundaries and proofs for the extension-field index-calculus speed-ups

**Status:** theory synthesis, no new measurements, 2026-09-20
**Scores:** `docs/index-calculus-scoreboard.html`
**Builds on:** `research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md`
(the `E(F_{p³})` residual-walk thread and its §11 phase pricing),
`research/index_calculus_baseline_20260914/EC_Index_Calculus_Performance_Baseline.md`
(the cost model), and the `F_{2^131}` product-law panel on the scoreboard.

> **What this note is.** The review that preceded it identified exactly
> two structural leads in the extension-field index-calculus work that
> are *not* engineering or relabelling: the relation count is
> sub-birthday (`κ ∝ N^{1/m−1/2}`), and the double-large-prime variant's
> total work has an exponent below rho's (`N^{4/9} < N^{1/2}` at `m = 3`).
> This note states the boundaries those two claims must be read against,
> **proves** each lead as a bound with its heuristics numbered and its
> falsifier fixed in advance, and prices every phase so that the whole
> method — not one phase — carries the exponent. It produces no new
> experiment; every measured number it uses is cited from a frozen file.
>
> **What this note is not.** It is not a claim of a new gain. The two
> leads are genuine advances of a *sub-problem* and an *asymptote*; the
> whole-method verdict, proved in §6, is that the plain `k = 3` method is
> bounded away from rho at every size, and the double-large-prime variant
> crosses rho only past `2^{230}`. The classes below follow §3 of
> `AGENTS.md` by the test there, not by how the result feels.

---

## 1. The setting

Fix a prime power `q`, a fixed prime extension degree `d`, and an ordinary
curve `E/F_{q^d}` with a prime-order subgroup of order `N ≈ q^d`. A
decomposition of length `m` writes a target `R` as a sum of `m`
factor-base points; for `E(F_{q^d})` the Semaev/Gaudry choice is `m = d`,
and the factor base is the `x`-coordinate-restricted set

```
F = { P ∈ E : x(P) ∈ V },     M = |F| ≈ |V| ≈ q = N^{1/d},
```

with `V` an `F_q`-rational set (prime-field curves) or an `F_2`-subspace
(binary curves). After sign/automorphism folding by an automorphism group
of order `γ` the number of *unknown* logarithms is `B ≈ M/γ = Θ(N^{1/d})`.
The running instance throughout is `d = m = 3`, which is what
`RESEARCH_RESIDUAL_WALKS.md` measures; the binary flagship is `F_{2^131}`,
`d = 131`.

**Unit.** As everywhere in this repository,

```
S = (total group operations, all phases) / √N,
```

with foreign units converted at a measured factor `c` (field
multiplications per group operation). "All phases" means setup, relation
finding, the decomposition oracle's per-call work, relation verification,
and the linear algebra. Dropping any of them is the relabelling failure
mode of §3.

---

## 2. Boundary A — the generic relation-yield floor

**Definition.** `κ := (residuals or group evaluations spent) / √N` is the
*count factor*. Boundary A lower-bounds `κ` for any algorithm that
discovers relations only by collision — i.e. any generic-group algorithm.

> **Theorem A (generic floor).** In the generic group model, an algorithm
> that has performed `P` group evaluations obtains an expected number of
> usable relations `R` with
> ```
> E[R] ≤ γ·P²/(2N).
> ```
> Consequently obtaining the `B + 1` independent relations a full-rank
> system needs requires `P ≥ √(2N(B+1)/γ)`, i.e.
> ```
> κ ≥ κ_floor := √(2(B+1)/γ).
> ```

**Proof.** In the generic model (Shoup) the encoding of each freshly
computed group element is a uniformly random label independent of the
exponent vector the algorithm has assigned to it, until a collision
occurs. A relation among the tracked elements is exactly such a
collision: two distinct exponent vectors `u ≠ u'` with `[u]·gens =
[u']·gens`, possibly identified up to the `γ` automorphisms folded into
the base. The number of unordered pairs among `P` labels is `C(P,2) ≤
P²/2`; each collides with probability at most `γ/N`. Linearity of
expectation gives `E[#collisions] ≤ γP²/(2N)`, and each collision yields
at most one relation, so `E[R] ≤ γP²/(2N)`. Setting the right side `≥
B + 1` and solving for `P` gives the stated bound; dividing by `√N` gives
`κ_floor`. ∎

Theorem A is unconditional. It is the boundary the AGENTS.md rule names
(`κ ≥ √(2(B+1)/γ)`), re-derived here so that §4's advance can be stated
as a ratio to it. The floor moves only with `B` and `γ`, so it cannot be
tuned away by any change to the walk.

---

## 3. Boundary B — the reference (Pollard rho)

**Definition.** The reference is Pollard rho with an `r`-adding walk and
distinguished points, run on the *same* subgroup with the *same*
operation accounting.

> **Fact B (reference constant).** Under the random-function heuristic for
> the walk, the first collision arrives after `√(πN/(2a))` steps, where
> `a` is the number of automorphism classes folded into the walk
> (`a = 1` plain, `a = 2` with negation, `a ≈ 2n` on suitable Koblitz
> subgroups). With a bounded number `w` of group operations per step,
> ```
> S_rho = w·√(π/(2a)) = Θ(1).
> ```

Measured `S_rho` is `1.2`–`1.8` across four group sizes on the same
instances and accounting (`docs/index-calculus-scoreboard.html`,
"Boundary 2"), and `≈ 1.3` on the `E(F_{p³})` instances of
`RESEARCH_RESIDUAL_WALKS.md`. The point of the unit `S` is that this
constant is flat in `N`, so "does index calculus beat rho" is the single
question "is the method's `S` column below `≈ 1.3`, robustly, at growing
`N`".

---

## 4. Lead 1 — the count advance, proved

The decomposition route escapes Boundary A because it does not find
relations by collision. This is the first lead, and it is a genuine §3
advance: the ratio to the floor falls to zero.

**Heuristic H1 (decomposition uniformity).** *For a signed base of `M`
points and decomposition length `m`, the sums `Σ_{j} P_{i_j}` over
unordered `m`-multisets are approximately equidistributed over the
order-`N` group. Hence the number of length-`m` decompositions of a
uniform target is Poisson with mean*
```
λ = C(M+m−1, m)/N ≈ M^m/(m!·N).
```
*Random-model justification.* A sum of `m` independent near-uniform
group elements equidistributes; the multiset count `C(M+m−1,m)` is the
number of trials, and rare coincidences make the per-target count
Poisson in the standard occupancy limit.
*Validation (frozen).* On ECC2K-130 the yield law `C(|F|,m)/#E` holds to
within `1.09×` over twelve toy cells, and the per-target count is Poisson
to `0.90`–`1.10` over eight cells counted exhaustively over the whole odd
subgroup, once point/negative degeneracies are removed
(`experiments/ecc2k130_decomposition_runs.json`).
*Falsifier (pre-registered).* Measured yield departs from `λ` by more
than `2×` at a saturating factor-base dimension. It has not fired.

> **Lemma 1 (constant decomposition probability).** Under H1 with `m`
> fixed and `M = c·N^{1/m}`, the decomposition probability
> `p_dec = 1 − e^{−λ}` satisfies `λ = c^m/m! = Θ(1)`, hence
> `p_dec = Θ(1)` independent of `N`.

**Proof.** Substitute `M = c·N^{1/m}` into `λ = M^m/(m!N) =
c^m·N/(m!·N) = c^m/m!`, a constant in `N`. Then `p_dec = 1 − e^{−λ}` is a
constant in `(0,1)`. ∎ (conditional on H1)

For the running case `m = 3` the measured constant is `p_dec ≈ 0.16`
(`RESEARCH_RESIDUAL_WALKS.md` §11.2), i.e. `λ ≈ 1/6`, consistent with
`c ≈ 1`.

> **Lemma 2 (count advance — the floor is beaten).** A decomposition
> oracle certifies a relation from a *single* target with probability
> `p_dec`, independently across targets (H1). Therefore from `P` targets
> `E[R] = p_dec·P` — **linear** in `P`, against Boundary A's **quadratic**
> ceiling `γP²/(2N)`. Reaching `B + 1 = Θ(N^{1/m})` relations needs
> ```
> P = (B+1)/p_dec = Θ(N^{1/m}),   so   κ_IC = P/√N = Θ(N^{1/m − 1/2}).
> ```
> For every `m ≥ 3` the exponent `1/m − 1/2 < 0`, so `κ_IC → 0`, and
> against `κ_floor = Θ(N^{1/(2m)})`
> ```
> κ_IC / κ_floor = Θ(N^{1/(2m) − 1/2}) → 0     (= Θ(N^{−1/3}) at m = 3).
> ```

**Proof.** Immediate from Lemma 1 and linearity of expectation for the
yield; the ratio to `κ_floor = √(2(B+1)/γ) = Θ(N^{1/(2m)})` is
`Θ(N^{1/m−1/2})/Θ(N^{1/(2m)}) = Θ(N^{1/(2m)−1/2})`, which for `m ≥ 3` is
`N^{<0}`. ∎ (conditional on H1)

For `m = 3`, `κ_IC = Θ(N^{−1/6})`: the residual walk note measures
`κ = 0.20` at 24 bits and `0.06` at 33 bits, falling as `N^{−1/6}`
(`RESEARCH_RESIDUAL_WALKS.md` §11.2). **This is the first speed-up made
rigorous:** it is an advance in the precise §3 sense, and it is an advance
of the *relation-finding sub-problem*, because the oracle uses the curve's
coordinate structure — a Semaev/summation-polynomial or a
meet-in-the-middle over coordinates — to decompose one point, which a
generic algorithm provably cannot do (Theorem A). It is **not** yet a
statement about the method: §5 prices the oracle, and §6 prices the
linear algebra.

---

## 5. Boundary C — the decomposition oracle is not free

Lemma 2 counts *residuals*, not *operations*. The count advance converts
to a work advance only if the oracle is cheap. Boundary C says how cheap
it can be, and it is the boundary the whole thread ultimately runs into.

Let `C_dec(N)` be the cost, in group-operation equivalents, of one oracle
call (decompose-or-reject one target). Then the relation phase costs
`P·C_dec = Θ(N^{1/m})·C_dec`, and the two realized oracles are:

- **Meet-in-the-middle triple oracle.** Builds/streams the base:
  `C_dec = Θ(M) = Θ(N^{1/m})`, so the relation phase is `Θ(N^{2/m})`
  (`= Θ(N^{2/3})` at `m = 3`; measured total `n^{0.69}`,
  `RESEARCH_RESIDUAL_WALKS.md` §10). Class **advance in count, not in
  work** — the residual count fell to `n^{1/3}` but each residual pays a
  loop over the base.
- **Algebraic `O(1)` `S₄` solve.** A fixed-degree zero-dimensional solve:
  `C_dec = C₃`, a constant in `N` (growing only with `log q` through the
  field-multiplication cost), so the relation phase is `Θ(N^{1/m}·C₃)`.
  This is the variant that can have a sub-rho *relation phase*.

> **Boundary C (oracle lower bound).** For `m = 3` any decomposition
> oracle that takes its target as input and answers from the target's
> coordinates is, up to the algebraic structure it exploits, a search over
> its own candidate set; a two-list (meet-in-the-middle) argument bounds
> its worst-case work over a candidate space of size `r` by `Ω(√r)`, and
> on `F_{2^131}` the measured **product law**
> ```
> relations × targets × oracle = m·2^n ,   independent of the base dimension,
> ```
> puts the best implicit-base cell at `2^{132.58} = 2^{+71.77}×` rho, with
> no table below plain BSGS
> (`experiments/ecc2k130_point_decomposition.json`, scoreboard binary
> panel).

Boundary C is why `C₃` cannot be assumed away: a smaller-than-`√r` oracle
would have to beat exhaustive search over its own candidate set, and on
the binary flagship the product law forbids exactly that. The corollary
is that the constant `C₃` is the object every solver improvement is scored
by, against the target of §6.

---

## 6. Whole-method boundaries: pricing every phase

Now assemble `S` from the relation phase (§4–5) and the linear algebra.
The linear algebra runs on `B ≈ N^{1/m}` unknowns with relations of weight
`m + 1`; sparse (Wiedemann/Lanczos) it is `Θ(B²) = Θ(N^{2/m})`, dense it
is `Θ(B^{2.4..3})`. Measured over four sizes (`RESEARCH_RESIDUAL_WALKS.md`
§11.7): linear algebra `n^{0.68}` (Wiedemann) and `n^{0.85}` (dense);
relation phase `n^{0.31}`.

> **Theorem 1 (relation-phase crossover — the §5 trap).** Ignoring the
> linear algebra, the best (`O(1)`-solve) relation phase has
> `S_rel = a·C₃·N^{−1/6}` (`m = 3`), so `S_rel < S_rho` iff
> ```
> C₃ < (S_rho/a)·N^{1/6} ,   calibrated to   C₃ < 13·n^{1/6}  F_q-mults.
> ```
> With the measured `C₃ ≈ 0.88·10⁶` this holds only for `n ≳ 2^{96}`.

**Proof.** `S_rel = (relation ops)/√N = Θ(N^{1/3}·C₃)/√N = a·C₃·N^{−1/6}`.
Solve `a·C₃·N^{−1/6} < S_rho`. The calibration constant `13` and the
value `C₃ ≈ 0.88·10⁶` are the frozen numbers of
`RESEARCH_RESIDUAL_WALKS.md` §11.5–11.6. ∎

Theorem 1 is a statement about **one phase**. Reading it as a method
result is precisely the §5 error the AGENTS.md rule exists to catch, and
Theorem 2 shows why it must not be.

> **Theorem 2 (the plain `k = 3` method is bounded away from rho).** With
> the best relation phase and sparse linear algebra,
> ```
> S(N) = a·C₃·N^{−1/6} + d·N^{+β},   β = 0.18 > 0 (measured; = 1/6 ideal).
> ```
> The first term vanishes and the second diverges, so `S` has a unique
> interior minimum `S* > 0` at a finite `N*`; the plain method therefore
> **cannot be driven below any fixed level by taking `N` large**. The
> fitted constants place `S* ≈ 265` at `N* ≈ 2^{50}`, i.e. `≈ 200×` rho.

**Proof.** `S(N) = a·C₃·N^{−1/6} + d·N^{β}` with `a,C₃,d > 0` and
`β > 0`. As `N → ∞` the second term dominates and `S → ∞`; as `N → 0` the
first dominates and `S → ∞`. `S` is smooth and convex in `log N` (sum of
two exponentials of `log N` with opposite-sign rates), so it has a unique
minimiser `N*` with `S(N*) = S* > 0`. Setting `dS/dN = 0`:
`(1/6)a C₃ N^{−7/6} = β d N^{β−1}`, giving
`N* = ((a C₃)/(6βd))^{1/(β+1/6)}` and a positive constant `S*`. The values
`β = 0.18`, `S* ≈ 265`, `N* ≈ 2^{50}` are the four-size fit of
`RESEARCH_RESIDUAL_WALKS.md` §11.7 and are marked there — and here — as an
**extrapolation** from the measured exponents `(−1/6, +0.18)`. Since
`S* ≈ 265 ≫ S_rho ≈ 1.3`, the plain method never crosses rho. ∎
(conditional on the two measured exponents holding)

The content of Theorem 2 is that the count advance of §4 is real but is
overwhelmed end-to-end: `κ` falling as `N^{−1/6}` is exactly cancelled by
linear algebra rising as `N^{+0.18}`, and the sum bottoms out `≈ 200×`
above the reference. This is why the plain method's honest class is
**engineering** (constants moved, ratio to the floor of the *whole method*
did not go below one), not advance.

> **Theorem 3 (the double-large-prime exponent, and where it crosses).**
> With a small base `F' ⊂ F`, `|F'| = M^{2/3} = Θ(N^{2/9})`, and relations
> carrying up to two large primes, under H1 and H2 the total work is
> `Θ(N^{4/9})` and
> ```
> S_DLP(N) = C_DLP·N^{4/9 − 1/2} = C_DLP·N^{−1/18} → 0,
> ```
> which is asymptotically below rho. The measured ratio `S_DLP/S_rho =
> 1,989` at `N = 2^{33}` and the `N^{−1/18}` slope place the crossover at
> `N_× = 2^{33}·1989^{18} ≈ 2^{230}`.

**Heuristic H2 (large-prime percolation).** *Relations carrying up to two
large primes recombine into `F'`-only relations exactly when the
large-prime "collision graph" percolates; at `|F'| = M^{2/3}` the
decomposition rate into `F' + ≤2` large primes is `Θ(M^{−1/3})`, so the
number of partial relations needed grows as `M^{4/3}`.* Random-model
justification: the two-large-prime graph is Erdős–Rényi on the
large-prime set, and its giant component appears at the standard
threshold, after which cycles supply independent `F'`-relations.
Validation (frozen): the residual exponent is measured `n^{0.44}` and the
decomposition rate `16% → 4%` across the range, and the matrix shrinks
from `1,069` to `104` unknowns at 33 bits
(`RESEARCH_RESIDUAL_WALKS.md` §11.7). Falsifier (pre-registered): the
residual exponent is measured below `4/9`, or recombined relations are
rank-deficient beyond the measured `5–11×` redundancy.

**Proof.** Residuals `P' = Θ(|F'|/p'_dec) = Θ(N^{2/9}·N^{1/9}) =
Θ(N^{1/3})` at unit oracle cost would give a relation phase `Θ(N^{1/3})`;
but H2 charges the partial-relation count as `Θ(M^{4/3}) = Θ(N^{4/9})`,
and the sparse linear algebra on `Θ(N^{2/9})` unknowns is `Θ(N^{4/9})`, so
the total is `Θ(N^{4/9})` and `S_DLP = Θ(N^{4/9−1/2}) = Θ(N^{−1/18})`.
For the crossover, `S_DLP = C_DLP·N^{−1/18} = S_rho` at
`N_× = (C_DLP/S_rho)^{18}`; substituting the measured
`S_DLP/S_rho = 1989` at `N = 2^{33}` gives `N_× = 2^{33}·1989^{18}`.
Since `1989 ≈ 2^{10.96}`, `1989^{18} ≈ 2^{197.3}` and
`N_× ≈ 2^{230}`. ∎ (conditional on H1, H2; the crossover is an
extrapolation from the `n^{−0.079}` measured slope)

**This is the second speed-up made rigorous:** the exponent `4/9 < 1/2` is
genuine and measured — it is the exponent Gaudry's `Õ(q^{2−2/d})` gives at
`d = 3` — but Theorem 3 also proves the crossover it implies is past
`2^{230}`, two hundred bits beyond anything reachable, because the
constant `C_DLP` (dominated by the same `C₃` `S₄`-solve of Boundary C) is
hopeless.

---

## 7. The table

One unit (`S`), one row per variant, the §3 class as it stands on the
scoreboard, and the source. Every `S` is cited from a frozen file; this
note computes none of them, and since rho is flat at `S ≈ 1.3` the ratio
to the reference is `S/1.3` for every row.

| Variant | `S` (measured) | Correct | Class (scoreboard, by §3 test) | Source |
|---|--:|:--|---|---|
| Pollard rho (reference) | ≈ 1.3 | yes | reference | scoreboard "Boundary 2" |
| Independent samples (A) | 1,850.6 | yes | baseline (one decomposition per residual) | gain chart |
| Mutation walk (B) | 19.7 (was 50.6) | yes | engineering (cost ↓2.6×, floor flat) | gain chart |
| r-adding residual (C1) | 39.4 (was 139.3) | yes | engineering (cost ↓3.5×, floor flat) | gain chart |
| Triple oracle, MITM | 90.1 | yes | relabelling (walked count ↓100×, total ↑4.6×) | gain chart |
| `O(1)` `S₄` solve | 898 (was 4,623) | yes | engineering (solve constant ↓5.5×) | gain chart |
| Meet-in-the-middle triple oracle | 3,805 | yes | **advance in count** (relations → `n^{1/3}`) | gain chart / §11.2 |
| Double large primes | 3,378 | yes | **advance in exponent** (`n^{4/9}`) | gain chart / §11.7 |

The two confusingly-named MITM rows are distinct constructions: the
`S = 90.1` "triple oracle" is the §10.5 relabelling example (the walked
count fell 100× while total cost rose 4.6×), whereas the `S = 3,805`
"meet-in-the-middle triple oracle" is the row that carries Lemma 2's
`n^{1/3}` count. Reading the class column against §4–§6:

- **The two advance rows are advances of a sub-problem and an asymptote,
  not of the method.** MITM cuts the *count* to `n^{1/3}` (Lemma 2) but
  pays it back in oracle work (Boundary C); double-large-primes has the
  sub-rho *exponent* (Theorem 3) but the crossover is `2^{230}`.
- **Every other row is engineering:** a constant fell, the whole-method
  ratio to the floor did not go below one (Theorem 2).
- **No row is below rho.** The only sub-rho cell anywhere in this
  repository is hyperelliptic genus 3/4 over *prime* fields
  (`RESEARCH_HYPERELLIPTIC_IC_RHO.md`), which is neither an elliptic curve
  nor an extension field and is toy-size (`p ≤ 251`).

---

## 8. Falsification targets, declared in advance

Per §4 of `AGENTS.md`, the numeric conditions that would turn a lead into
a whole-method advance, and what is inadmissible.

**A lead becomes a result if any of these is measured on a verified full
run** (planted secret recovered, every relation re-added in the group,
zero trivial/dependent relations counted):

1. **FT-constant.** A `k = 3` decomposition oracle with amortized
   `C₃ < 13·n^{1/6}` `F_q`-multiplications *and* a linear-algebra phase
   that keeps the end-to-end `S(N)` below its own §6 minimum — i.e. an
   end-to-end `S < S*` fit at some `N ≥ 2^{34}`. Theorem 1 alone is not
   enough; Theorem 2 must be defeated.
2. **FT-exponent (DLP).** The double-large-prime residual exponent
   measured below `4/9`, or the sparse linear-algebra exponent below
   `1/3` in the number of unknowns, with the constants holding — either
   moves `N_×` down from `2^{230}`.
3. **FT-crossover.** Any variant whose measured `S` is below rho's
   `≈ 1.3` at an `N` for which a matched rho run is executed on the same
   subgroup and accounting.

**Inadmissible** (would move the boundary, not cross it): changing `B`,
`m`, or the unit; dropping the linear algebra or the oracle's per-call
work from the budget; counting Frobenius images as independent relations;
estimating `p_dec` from planted rather than uniform targets; or quoting a
phase crossover (Theorem 1's `2^{96}`) as a method crossover.

---

## 9. What would count next

The math of §6 says exactly where the leverage is, and it is a **constant**
game, not an exponent game:

- Every halving of `C₃` moves the double-large-prime crossover down by
  `≈ 6` bits of `n` (`N_× ∝ C_DLP^{18}`, and `C_DLP` is dominated by
  `C₃`). Closing `2^{230}` to a reachable size needs `C₃` cut by many
  orders of magnitude, which is why the attacks that matter in practice
  either cut the constant (Joux–Vitse's `F₄`-based `k−1` decompositions)
  or change the target (Weil descent / GHS), not tune the walk.
- On the binary flagship the corresponding lever is Boundary C's product
  law `m·2^n`: it is dimension-independent, so no factor-base engineering
  moves it — only a decomposition oracle cheaper than a two-list search
  over its own candidate set would, and none proposed has that shape
  (scoreboard binary panel, class **accounting** for the "protective
  margin" that two earlier revisions had mistakenly priced).

The honest one-line close, in the same unit as everything else: on
`E(F_{p³})` the relation phase is genuinely `n^{1/3}` and the
double-large-prime exponent is genuinely `n^{4/9}`, both proved here as
bounds against a stated floor and reference; and both are out of reach —
the plain method by Theorem 2 at every size, the exponent variant by
Theorem 3 until past `2^{230}` — because the decomposition oracle
(Boundary C) is not free.
