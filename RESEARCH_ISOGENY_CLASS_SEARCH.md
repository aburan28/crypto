# Searching the ECC2K-130 isogeny class for an easier Gröbner basis

Lever **L5** of `RESEARCH_DEGREE_REDUCTION.md` — the only one that changes
the *curve* rather than the presentation.

The question: the ECDLP transports along an isogeny of degree coprime to the
subgroup order, so an attacker may solve the problem on **any** curve
`F_{2^131}`-isogenous to ECC2K-130 and pull the answer back.  Galbraith,
Hess and Smart turned exactly that freedom into an attack on Weil descent —
walk the isogeny graph until the GHS genus drops.  Does the same walk buy a
lower Gröbner **solving degree**?

**Bottom line.** No, and the reason is not a search that came up empty.
All four boundaries below are *derived* rather than measured, and between them
they leave the lever no room:

- the class has `2^65.06` vertices against ρ's `2^60.81` operations, so
  enumerating it is already `2^4.25×` more expensive than just solving the
  DLP — and still more expensive, by `2^0.24×`, against plain ρ with no
  automorphism speedup at all (**A**);
- exactly **263** of those vertices are reachable by a computable isogeny —
  everything else needs an isogeny of degree `146505763881528721` whose
  kernel polynomial has degree `2^56` (**B**);
- the curve coefficient `a₆ = 1/j` enters the descended Semaev system
  **strictly below the leading form** (as a constant at `m = 2`; in Boolean
  degree `≤ 5` of `6` at `m = 3`), so the degree of regularity is *constant
  on the whole class* (**C**);
- and the one mechanism that does lower `D*` — an L1 subfield factor base,
  mean `D*` 2.04 against 3.53 — needs a proper subfield, which `F_{2^131}`
  does not have and an isogeny cannot create (**D**).

What Boundary C leaves is the **affine** refutation degree `D*`, which the
curve does move because the curve *is* the inhomogeneous part.  EXP-R6
measures that surface exhaustively — every curve over `F_{2^n}`, a superset
of every isogeny class at that size — and finds the variation is a property
of the **(curve, target) pair**, not of the curve: the exact criterion is
`D* = 2 ⟺ a₆ ∉ S^⊥` for a subspace `S^⊥` fixed *by the target*, verified
with zero disagreements over every curve at four targets.  Its consequence is
measured directly: the number of curves that stay on the `D* = 2` floor for
*every* target decays to **zero** — at `n = 8`, from 68 of 255 at `T = 8` to
**0 by `T = 32`** — so a curve that is easy for every decomposition instance,
which is what moving along an isogeny class would have to find, does not
exist.  A curve selected as best on one target set is `2^{+1.13}` **worse**
than the unmodified Koblitz curve on a disjoint one.

---

## 1. Why this thread

`RESEARCH_DEGREE_REDUCTION.md` §2 had four levers, L1–L4, all of which hold
the curve fixed and re-present its decomposition ideal.  L5 is the column
that taxonomy was missing, and it is the one an attacker gets for free:
nothing about the index-calculus pipeline requires the curve you attack to be
the curve you were given.

It is also the lever with the best historical precedent.  Galbraith–Hess–
Smart (EUROCRYPT 2002) is the standard example of isogeny-hopping turning an
intractable descent into a tractable one, and `src/cryptanalysis/binary_isogeny.rs`
already implements that walk for the GHS magic number.  Asking the same
question of the solving degree is the obvious next move, and it had not been
asked here.

The pay-off is two-sided in the usual way.  If some isogenous curve solves
at a lower degree, that is a route into ECC2K-130.  If none does, the reason
*why* is a statement the defensive FFD theorem needs and currently assumes:
that the solving degree is a property of the field and the factor base, not
of the curve.

---

## 2. The boundaries, before measuring

### A — the class is larger than ρ (floor, counting)

`E: y² + xy = x³ + 1` over `F_{2^131}` has `End(E) ⊗ Q = Q(√−7)`,
`d_K = −7`, `h_K = 1`.  With `τ² + τ + 2 = 0` the Frobenius of `F_2` and
`π = τ^131`, writing `τ^131 = a + bτ` gives `Z[π] = Z + Z·bτ`, so the
conductor of `Z[π]` in `O_K = Z[τ]` is `|b|`:

```
  t_131 = −22283658519494248867
  #E    = 2722258935367507707729280517973639940516 = 4 · r,
          r = 680564733841876926932320129493409985129  (prime, 130 bits)
  c     = 38531015900842053623 = 263 · 146505763881528721   (both prime)
```

verified against `t² − 4q = −7c²` and against the challenge order the
repository already stores for `ecc2k-130`.  The class is every curve whose
endomorphism ring is an order between `Z[π]` and `O_K`, so

```
  #class = Σ_{f | c} h(O_f) = 38531015900842054149 = 2^65.06,
  h(O_f) = f · Π_{ℓ | f} (1 − (d_K/ℓ)/ℓ)         (h_K = 1, unit index 1)
```

with the four terms `1 + 262 + 146505763881528722 + 38384510136960525164`.
The formula is cross-validated the hard way: at `n ∈ {5, 7, 9}` the number
of `a₆` whose curve has `|t| = |t_n|`, by exact point counting over every
element of the field, equals `Σ_{f|c} h(O_f)` on the nose
(`class_size_formula_matches_exact_point_counting`).

Pollard ρ on a group of order `r` whose automorphism group has order `m` and
acts freely walks on `r/m` classes, so the birthday bound is
`√(π(r/m)/2) = √(πr/2m)`.  For ECC2K-130 the automorphism group is negation
together with the 131 Frobenius powers, `m = 2 · 131 = 262`:

```
  plain ρ,        m = 1     √(πr/2)    = 2^64.83
  negation only,  m = 2     √(πr/4)    = 2^64.33
  full,           m = 262   √(πr/524)  = 2^60.81   ← what the ECC2K-130 effort targets
```

> **Boundary A.** Visiting every vertex of the isogeny class costs `2^4.25×`
> the automorphism-assisted ρ, at one operation per vertex and before any
> Gröbner work.  A search over more than `2^{−4.25} ≈ 5.2%` of the class is
> beaten by ρ on enumeration alone, whatever it finds.

**The margin is reported against both ρ variants on purpose.** Against *plain*
ρ the class exceeds it by only `2^0.24×` — the same order, not four bits of
room.  So Boundary A's load-bearing content is the **sign**, which holds
against every ρ variant, and not the size of the margin, which depends on how
much automorphism speedup one credits ρ with.  An earlier revision of this
note quoted `2^60.31` and `2^4.75×`: it divided `√(πr/4)` by `√(2·131)`, which
applies negation twice, since `√(πr/4)` is already the `m = 2` form.  That
understated ρ by half a bit and overstated this boundary's margin.  Caught in
review on PR #203 and corrected here, in the module and in the test.

### B — 263 vertices are reachable, and that is the whole search (floor, reachability)

Because `h(−7) = 1`, the crater of every `ℓ`-volcano is a **single vertex,
`E` itself**.  Horizontal `ℓ`-isogenies act through `Cl(O_K) = 1`, so they
return to `E`; descending ones exist only when `ℓ | c` (Kohel).  Hence:

| `ℓ` | `v_ℓ(c)` | `(−7/ℓ)` | horizontal | descending from crater |
|---|---|---|---|---|
| every prime `∉ {263, p}` | 0 | — | back to `E` | **0** |
| `263` | 1 | `+1` | 2, back to `E` | **262** |
| `p = 146505763881528721` | 1 | `−1` | 0 | `146505763881528722` |

> **Boundary B.** For every prime `ℓ` below 263 the `ℓ`-isogeny graph of
> ECC2K-130 is a **single vertex**: there is no small-degree walk to search.
> The reachable class is `1 + 262 = 263` curves — `2^8.0` of `2^65.06`, a
> fraction of `2^{−57.0}`.  Reaching anything else needs a `p`-isogeny:
> kernel-polynomial degree `(p−1)/2 = 2^56.0`, kernel points in an extension
> of `F_{2^131}` of degree dividing `p − 1`, and `2^28.5` field operations
> even by √élu.

So the exhaustive search *does* terminate — over 263 curves — which is what
makes an exhaustive claim possible at all.  It is worth being precise that
this is the second of three nested notions of "exhaustive" in play:

| notion | size | status |
|---|---|---|
| the isogeny class | `2^65.06` | infeasible, and beaten by ρ before it starts (A) |
| the **reachable** class | `263` | **exhaustive, and the search terminates** (B) |
| every curve over `F_{2^n}`, small `n` | `2^n − 1` | **exhaustive, and a superset of every class** (C, §3) |

**Cross-checked independently.** The table above comes from the CM
conductor.  Root-finding the modular polynomial mod 2 at `j(E) = 1` over the
*actual* `F_{2^131}` — a completely different computation — agrees: `Φ₂(X,1)`
has the single root `1` (Frobenius and Verschiebung, both landing back on
`E`), and `Φ₃(X,1)` has **no root at all** (3 is inert in `Q(√−7)` and
`3 ∤ c`).

### C — the curve enters below the leading form (exact, algebraic)

For a binary curve `y² + xy = x³ + a₂x² + a₆`,

```
  S₃(x₁,x₂,x₃) = (x₁x₂)² + (x₁x₃)² + (x₂x₃)² + x₁x₂x₃ + a₆
```

is independent of `a₂`, and `a₆ = 1/j` appears as an **additive constant**.
So `a₆` is the entire curve-dependence of the `m = 2` system, and after Weil
descent it moves only the constant term of each of the `n` Boolean
equations.  Measured, not assumed: at `n ∈ {8, 10, 12}`, 64 curves compared
coefficient-by-coefficient, `a₆` touches Boolean degree `{0}` of `≤ 2`, and
the degree-`≥ 1` part is identical for every curve.

**Beyond `m = 2`, computed rather than estimated (EXP-R6b).** For `m ≥ 3` the
curve no longer enters as a pure constant, so the question is where exactly it
lands.  Two things make this worth computing:

*Boolean degree is a sum of Hamming weights, not a total degree.* After descent
each `Xᵢ` is a vector of `F_2`-linear forms in its bits, and `Xᵢ^{2^k}` is
**also** `F_2`-linear (Frobenius is linear, and `b² = b` on Boolean
coefficients).  So `Xᵢ^e` is a product of `wt(e)` linear forms and

```
  bdeg(Π Xᵢ^{eᵢ}) = Σᵢ wt(eᵢ).
```

Multiplying factor degrees therefore *overestimates*: in `S₄` the product
`(A₂B₁)(B₁C₂)` looks like `3 × 3 = 6`, but `X₁X₂ · X₁X₂ = X₁²X₂²` collapses to
Boolean degree 2.  **The first version of this note quoted `a₆` reaching degree
5 at `m = 3` from exactly that product bound; the exact value is 4.**

Building `S_{m+1}` by Semaev's recursion
`S_{i+j−2} = Res_Y(S_i(…,Y), S_j(…,Y))` and profiling every monomial by its
`a₆`-power gives:

| `m` | symbolic vars | top Boolean degree | max bdeg of `a₆`-carrying terms | gap |
|---:|---:|---:|---:|---:|
| 2 | 2 | 2 | 0 | **2** |
| 3 | 3 | 6 | 4 | **2** |
| 4 | 4 | 12 | 10 | **2** |
| 5 | 5 | 20 | 18 | **2** |

The top Boolean degree is `m(m−1)`, it is **always `a₆`-free**, and every
`a₆`-carrying term sits **exactly two degrees below it**.  Two independent
corroborations of the top figure: at `m = 2` it is 2, which is why the
descended `S₃` system is quadratic (`ffd_harness::F2BoolPoly` encodes exactly
degree `≤ 2`); at `m = 3` it is 6, which is the degree
`binary_semaev_s4`'s own module doc records for the eliminated `3ℓ`-variable
presentation.

**And `m(m−1)` is not a coincidence — the upper bound holds for every `m`.**
Semaev's construction gives `deg_{Xᵢ} S_{m+1} = 2^{m−1}`, and over the
exponents `e ≤ 2^{m−1}` the Hamming weight is maximised at
`e = 2^{m−1} − 1 = 0b11…1`, where `wt(e) = m−1`.  Since the Boolean degree of a
monomial is `Σᵢ wt(eᵢ)` over the `m` symbolic variables,

```
  bdeg(S_{m+1}) ≤ m · max{ wt(e) : e ≤ 2^{m−1} } = m(m−1)     for every m.
```

So the *ceiling* is proved in general; what the computation adds is that the
ceiling is **attained**, that the monomials attaining it are `a₆`-free, and
that `a₆` stops two short of it.  Those three facts are established at
`m ∈ {2,3,4,5}` and are what Boundary C rests on; they are not proved for all
`m`, and §7 says so.

The gap being **constant rather than shrinking** is the load-bearing part.  A
narrowing gap would predict the boundary failing at some larger `m`; a flat one
says it does not erode.

*Verification.* Each `S_{m+1}` is checked three ways before its profile is
believed: symmetry in all `m+1` arguments (on adjacent transpositions, which
generate the symmetric group), degree `2^{m−1}` in each argument as Semaev's
construction requires, and — for `S₄` — agreement with the repository's own
`binary_semaev_s4` at 60 random field points, which validates the resultant
machinery against trusted code.  `S₅` was additionally checked semantically:
over four curves on `F_{2^5}` it vanished on **all 1146** genuine 5-point
decompositions and on no tuple lacking one.

> **Boundary C.** The leading-form ideal of the descended Semaev system is
> independent of the curve, for every `m ∈ {2, 3, 4, 5}` — `a₆` stays exactly
> two Boolean degrees below the leading form at each.  Therefore the degree of
> regularity in the Bardet–Faugère–Salvy sense — a Hilbert-series invariant of
> the leading forms — is **constant on the entire isogeny class**, and so is
> the degree at which any top-degree cancellation first becomes *available*.
> Changing the curve can make a fall's remainder vanish, i.e. **lose** a fall;
> it cannot create one at a lower degree.  The lever's only possible effect on
> `d_reg` is adverse.

This is the answer to the question as literally posed — "curves where the
degree of regularity is easier" — and it is negative by a two-line argument
that needs no search.

**What C does not cover.** The quantity this repository actually measures is
the *affine* refutation degree `D*` (`pc_degree_harness::refutation_scan`:
the degree at which `1` enters the Macaulay row space).  `D*` is a property
of the inhomogeneous system, and the curve controls precisely the
inhomogeneous part, so `D*` does vary with `a₆` — at `n = 8, l = 4` over all
255 curves it takes values `{2, 3, 4}`.  That variation is the lever's whole
remaining surface, and §3–§4 are about it.

### D — the one mechanism that works is unreachable (exact)

The only lever the FFD program has measured to lower `D*` materially is L1,
a **subfield** factor base: mean `D*` 2.04 against 3.53 for a random base,
with `Δ_low` ratios diverging 6.7 → 67.  It requires a proper subfield of
the base field.  `131` is prime, so `F_{2^131}` has none but `F_2`, and an
isogeny is defined over the same field it starts in.

> **Boundary D.** L1 is unreachable from ECC2K-130 by any isogeny, for every
> curve in the class.  The prime extension degree is not a property an
> isogeny can change.

### Reference

The unmodified Koblitz curve itself, and `2^N` enumeration of `V × V`
(`degree_reduction::log2_enumeration_cost`), both measured in the same unit
as every variant — see §6.

---

## 3. What is measured, and why it is exhaustive

Boundary C is a statement about *all* curves, so it is tested the strongest
way available: sweep **every** curve over `F_{2^n}` — all `2^n − 1` values
of `a₆`.  That set is a **superset of every isogeny class** at that size,
the Koblitz one included, so finding no uniformly-easy curve there rules the
lever out for every class at once, without enumerating any of them.  This containment is what
makes an exhaustive statement possible at a size where the `2^65` class is
not enumerable:

> whatever the 263 reachable curves at `n = 131` would do, the sweep at
> small `n` has already done it — and more, since it also covers the
> `2^65 − 263` unreachable ones and every other isogeny class.

Four numbers come out, and the last one decides.

1. **`dim S` and the exact criterion.** Write the descended system as
   `f_i = h_i + c_i`, `h_i` the degree-`≥ 1` part (curve-independent, by C)
   and `c_i ∈ F_2` the constant.  Let `S = {λ : Σ λ_i h_i = 0}` be the left
   null space of the leading parts — also curve-independent.  For `λ ∈ S`,
   `Σ λ_i f_i = ⟨λ, c⟩` is a *constant*, so

   ```
     D* = 2   ⟺   ∃ λ ∈ S with ⟨λ, c⟩ = 1   ⟺   c ∉ S^⊥
   ```

   and `c` is exactly the coordinate vector of `a₆`, because
   `S₃(0, 0, x_R) = a₆`.  So the curves that refute at the `D* = 2`
   Nullstellensatz floor are precisely those whose `a₆` escapes a subspace
   **`S^⊥` fixed by the target `x_R`**, not by the curve.
2. **The between-curve variance** of mean `D*`, against what independent
   sampling from the pooled distribution predicts.  `≈ 1` means the curve
   explains nothing; `≈ T` would mean `D*` is a function of the curve.
3. **The holdout margin.** Take the curve that looks best on target set A;
   re-measure it on a disjoint set B.  A real lever keeps its margin; a
   winner's curse does not.
4. **The uniform-floor survivor count.** How many curves stay on the
   `D* = 2` floor across *every* one of `T` targets.  That is what an attacker
   moving along an isogeny class would need — a curve that is easy for all the
   decomposition instances a relation search throws at it, not for a lucky
   dozen.  If goodness is a curve property the count is flat in `T`; if it is
   a `(curve, target)` property, the criterion in item 1 says it decays.

### The prediction ledger

| # | Prediction | Status | Evidence |
|---|---|---|---|
| **R6** | Some curve in the ECC2K-130 isogeny class has a materially lower solving degree than `E` | **`killed`** | Boundaries C and D, both exact: `d_reg` is constant on the class, and the one `D*`-lowering mechanism needs a subfield `F_{2^131}` does not have. No search required. |
| **R6′** | The residual `D*` variation is a **curve** effect an attacker can move to | **`killed`** | EXP-R6. The exact criterion of §3.1 makes it a `(curve, target)` property; zero disagreements with the solver over every curve at four targets; holdout margin negative at the largest size. |
| **R6″** | The exhaustive search over the class is *feasible* | **`killed`** | Boundary A (`2^65.06` vertices vs `2^60.81` ρ, or `2^64.83` plain ρ) and B (263 reachable). The search that does terminate covers `2^{−57}` of the class. |
| **R6‴** | `D* = 2` density over curves is `1 − 2^{−dim S}` per target | **`supported`** (exact) | `dim S = 1` at `n ∈ {8, 10}`, escape count `128/255` and `512/1023`, matching `(2^n − 2^{n−1})/(2^n − 1)`; mismatches `0`. |
| **R6⁗** | Some curve is on the `D* = 2` floor for **every** target — the uniformly-easy curve an isogeny walk would need | **`killed`** (exhaustively) | EXP-R6: survivor count `68 → 34 → 0` at `n = 8` over `T = 8/16/32`, and `230 → 73 → 14 → 7 → 1` at `n = 10`. **EXP-R6c settles it without extrapolation, at both sizes: over *all* 240 targets at `n = 8`, `0` of 255 curves avoid above-floor targets entirely, and over all 992 targets at `n = 10`, `0` of 1023.** Mean above-floor count is `39.1` per curve (min 26, max 56) — such targets are common, and the `T = 64` zeroes were small-sample. |
| **R6b** | `a₆` reaches the **leading form** at some `m ≥ 4`, making `d_reg` curve-dependent where index calculus is asymptotically interesting | **`killed`** | EXP-R6b, iteration 2. Computed symbolically for `m ∈ {2,3,4,5}`: the top Boolean degree is `m(m−1)`, always `a₆`-free, with `a₆` exactly **2** degrees below at every `m`. Constant gap, not a narrowing one. `S₄` validated against the repo's own implementation; `S₅` against 1146 genuine decompositions. |
| **R6c** | The residual per-curve variation in `D*` statistics is a *solving-degree* property of the curve | **`killed`** | EXP-R6c, iteration 3. It is **decomposition yield**: `ρ_s(decomposable targets, above-floor targets) = −0.9801` over all 255 curves × 240 targets at `n = 8`, and `−0.9648` against the above-floor *rate*, so it is not the mechanical "fewer refutable targets means fewer bad ones". A curve that decomposes more targets has fewer left that can refute above the floor — a relation-yield property, not a `d_reg` one. |

### Pre-registered gates

- **G-R6.** Decided **at the largest measured `n`**, following the thread's
  existing convention (G-R2: "killed if the gap is `≤ 0` at the largest
  size"); smaller sizes are reported but do not vote, because a margin that
  survives at `n = 6` and dies at `n = 10` is the small-size artifact this
  gate exists to catch.  *Supported* if some curve keeps a mean-`D*` margin
  `≥ 0.5` degrees over the Koblitz curve on a **disjoint holdout** target set
  at every measured `n`, **and** the between-curve variance exceeds twice the
  no-effect prediction.  *Killed* if, at the largest size, the holdout margin
  is `≤ 0` while the selection margin is positive (a winner's curse), or if
  the uniform-floor survivor count reaches zero.  *Blocked* otherwise — a
  margin that is positive but shrinking.
- **G-R6‴.** The criterion of §3.1 is *supported* only at **zero**
  disagreements with the solver over the whole sweep.  One mismatch refutes
  it; the density version alone does not count, because filtering on
  refutability conditions the density (a decomposable target makes the
  system satisfiable, and every such curve necessarily has `c ∈ S^⊥`).

---

## 4. The instrument

| piece | what it does |
|---|---|
| `cryptanalysis::isogeny_class_search::koblitz_isogeny_class` | exact class structure for `F_{2^n}`: trace, order, conductor, factorisation, `Σ_{f|c} h(O_f)`, per-prime volcano |
| `…::reach_within_degree` | coverage bookkeeping for a degree-bounded search, with the blocked primes and their √élu cost |
| `…::leading_form_witness` | Boundary C at `m = 2`, coefficient by coefficient across curves |
| `…::syzygy_mechanism` / `leading_part_left_nullspace` | `dim S`, and the exact `D* = 2 ⟺ c ∉ S^⊥` criterion checked against the solver on every curve |
| `…::exhaustive_a6_sweep` | every curve over `F_{2^n}`, with `D*` and first-fall histograms |
| `…::curve_effect_test` | variance decomposition plus the disjoint-holdout winner's-curse control |
| `…::uniform_floor_survivors` | the decisive statistic: curves on the `D* = 2` floor for *every* target, swept over the whole curve space |
| `…::yield_explanation` / `spearman` | the three-way target partition behind the residual variation, and its correlation with decomposition yield |
| `…::all_traces` | exact point counting, for the class-size cross-validation |
| `cryptanalysis::semaev_leading_form` | symbolic `S_{m+1}` over `F_2[a₆]` by the resultant recursion, and the Boolean-degree profile that decides Boundary C at each `m` (4 tests) |
| `examples/isogeny_class_search.rs` | EXP-R6 driver → `experiments/isogeny_class_search.json` |

`D*` is measured by the thread's existing `pc_degree_harness::refutation_scan`,
and costed by `degree_reduction::log2_expected_cost` (Jensen-correct over the
histogram, not the mean), so the unit matches every other lever's.

Run: `cargo run --release --example isogeny_class_search`.

---

## 5. Iteration log

### 2026-09-12 — iteration 3 (EXP-R6c — the residual effect is decomposition yield)

**Task.** The other item iteration 1 left open: a variance ratio of `1.4`
rather than `1.0`, and one curve holding the `D* = 2` floor over 64 targets at
`n = 10`.  Iteration 1 called it "a small residual curve effect that is not
zero" and left it unexplained, which is the kind of remainder that quietly
becomes a claim if nobody chases it.

**Two wrong explanations first, both discarded.**

1. *A decomposability artifact* — the survivor statistic skips decomposable
   targets, so a curve with few refutable targets could be credited vacuously.
   **Wrong:** `a₆ = 13` had 35 refutable targets of 64 and all 35 refuted at
   the floor.  Nothing vacuous about it.
2. *Target correlation* — the exact criterion lets the survivor count be
   computed with no Gröbner work at all (the survivor set is the complement of
   a union of `T` subspaces), and that computation gives `0` survivors by
   `T = 16` for both consecutive and spread targets.  **Also wrong, as an
   answer to this question:** it counts a *decomposable* target as a failure to
   escape.  For an attacker a decomposable target is a success — it is a
   relation — so that statistic answers a question nobody asked.  The measured
   statistic, over non-decomposable targets, is the attack-relevant one and it
   stands.

**Experiment.** The criterion partitions a curve's targets three ways: `a₆ ∉
S^⊥` (refutes at the floor); `a₆ ∈ S^⊥` and decomposable (skipped, a success);
`a₆ ∈ S^⊥` and not decomposable (**above** the floor — the only costly case).
The third is squeezed by the second, so the hypothesis is that the residual
variation *is* decomposition yield.  `yield_explanation` measures both counts
over **every** curve and **every** target above the factor base.

**Result** (`n = 8`, `l = 4`, 255 curves × 240 targets, exhaustive):

| quantity | min | max | mean |
|---|---:|---:|---:|
| decomposable targets per curve | 58 | 92 | 78.9 |
| above-floor targets per curve | 26 | 56 | 39.1 |

```
  ρ_s(decomposable count, above-floor count) = −0.9801
  ρ_s(decomposable count, above-floor rate)  = −0.9648
  curves with zero above-floor targets       = 0 of 255
```

- The correlation is **near-deterministic**, and it survives normalising by
  refutable count — so it is not the mechanical "fewer refutable targets means
  fewer bad ones".
- **No curve is uniformly easy** once every target is used: `0` of 255.  The
  `T = 64` zeroes of iteration 1 were small-sample; above-floor targets are in
  fact common (mean 39 per curve).

**Gate verdict.** **R6c killed.** The residual variation is a *relation-yield*
property, not a solving-degree one — a different quantity, tracked elsewhere in
this repository, and one where more yield helps an attacker for reasons that
have nothing to do with `d_reg`.

**Ledger delta.** R6c `killed`; R6⁗ upgraded from "monotone and reaching zero"
to exhaustive at `n = 8`. Limitation 4 rewritten from an unexplained remainder
to an explained one.

**Class of the change.** **Accounting** — a caveat was resolved, not an attack
improved.

**Replicated at `n = 10`** (1023 curves × 992 targets, exhaustive): `0` of 1023
curves keep `D* = 2` on every non-decomposable target.  And the specific curve
that generated iteration 1's caveat, `a₆ = 13`, **breaks by `T = 256`** — 171
refutable targets, 137 at the floor, 34 above it — against a clean 35-of-35 on
the first 64.  So the caveat was a small-sample artifact of the target count,
not of the decomposability filter (hypothesis 1) and not of target correlation
(hypothesis 2).

**Next.** Nothing on this item.  The one open question in the thread is the
induction for Boundary C at all `m` (§8).

---

### 2026-09-12 — iteration 2 (EXP-R6b — Boundary C holds to `m = 5`, and the ceiling is explained)

**Task.** The queue head from iteration 1, and the one thing that could
reopen L5: does `a₆` reach the leading form at `m ≥ 4`, where index calculus
is asymptotically interesting?  Iteration 1 proved Boundary C only for
`m ∈ {2,3}` by hand expansion.

**Experiment.** Build `S_{m+1}` symbolically over `F_2[a₆]` by Semaev's
resultant recursion and profile every monomial by its `a₆`-power, for
`m ∈ {2,3,4,5}`.  New module `cryptanalysis::semaev_leading_form` (5 tests):
multivariate `F_2` polynomials, the characteristic-2 quadratic resultant, and
a memoised Sylvester determinant for the quadratic×quartic (`S₅`) and
quartic×quartic (`S₆`) steps.

**Result.**

| `m` | symbolic vars | monomials | top bdeg | `a₆` max | gap |
|---:|---:|---:|---:|---:|---:|
| 2 | 2 | 5 | 2 | 0 | **2** |
| 3 | 3 | 24 | 6 | 4 | **2** |
| 4 | 4 | 729 | 12 | 10 | **2** |
| 5 | 5 | 190252 | 20 | 18 | **2** |

- **Boundary C holds at every `m` reached**, with a **constant** gap of 2 —
  not a narrowing one, which is what would have predicted failure further out.
- **The `m(m−1)` ceiling is proved for all `m`**, not merely observed:
  `deg_{Xᵢ} S_{m+1} = 2^{m−1}` and `max{wt(e) : e ≤ 2^{m−1}} = m−1`, so
  `bdeg ≤ m(m−1)` always.  The computation supplies the three facts the bound
  does not: attainment, `a₆`-freeness of the attaining monomials, and the
  two-degree shortfall of `a₆`.
- **A correction to iteration 1.** It quoted `a₆` reaching Boolean degree
  **5** at `m = 3`.  The exact value is **4**.  The 5 came from multiplying
  factor degrees (`3 × 3 = 6`, minus one), which overestimates: Boolean degree
  is `Σ wt(eᵢ)`, and `X₁X₂ · X₁X₂ = X₁²X₂²` collapses to degree 2 under
  Frobenius rather than doubling to 4.  The direction of the claim is
  unaffected — `a₆` was and is strictly below the leading form — but the
  number was loose and is now exact.

**Verification.** Each `S_{m+1}` is checked three ways before its profile is
used: symmetry in all `m+1` arguments (on adjacent transpositions, which
generate the symmetric group), degree `2^{m−1}` per argument as Semaev
requires, and for `S₄` agreement with the repository's own
`binary_semaev_s4` at 60 random field points — which validates the resultant
machinery against trusted code rather than against itself.  `S₅` was checked
semantically as well: over four curves on `F_{2^5}` it vanished on **all 1146**
genuine 5-point decompositions and on no tuple lacking one.

**Gate verdict.** **R6b killed** — no counterexample at `m = 4` or `m = 5`.

**Ledger delta.** R6b `killed`. Boundary C strengthened from `m ∈ {2,3}` to
`m ∈ {2,3,4,5}` plus a general ceiling; limitation 2 narrowed accordingly.

**Class of the change.** **Accounting** again — a boundary got firmer and one
of its numbers got corrected; no attack moved.

**Next.** The induction: prove the gap is 2 for all `m`, or find a
counterexample at `m ≥ 6`.  `S₇ = Res(S₄, S₅)` is the next computable step.

---

### 2026-09-12 — iteration 1 (EXP-R6 — the lever is empty, and every reason is derived)

> This note's iteration 1 is `RESEARCH_DEGREE_REDUCTION.md`'s **iteration 7**;
> the two numberings are separate because this is a separate note.

**Task.** Score L5: search the isogeny class of ECC2K-130 for a curve with a
lower solving degree.

**Experiment.** (i) Derive the class exactly and price its enumeration
against ρ. (ii) Bound what a degree budget reaches, and cross-check with
modular polynomials over the real `F_{2^131}`. (iii) Establish where `a₆`
sits relative to the leading form, at `m = 2` and `m = 3`. (iv) Sweep every
curve over `F_{2^n}` for `n ∈ {6, 8, 10}`, `l = n/2`, with disjoint
selection and holdout target sets.

**Result.**

- **A.** class `2^65.06`, ρ `2^60.81` (`m = 262`) → enumeration costs `2^4.25×` ρ,
  and `2^0.24×` plain ρ.
- **B.** reachable `263` (`2^{−57.0}` of the class); every prime `< 263` has
  a single-vertex graph. `Φ₂(X,1)` root set `{1}`, `Φ₃(X,1)` empty — the CM
  prediction and the modular polynomial agree.
- **C.** `a₆` touches Boolean degree `{0}` of `≤ 2` at `n ∈ {8,10,12}`,
  degree-`≥ 1` identical across curves; the `m = 3` expansion verified on
  300 combinations.
- **Mechanism.** `dim S = 1`; escape count `128` of 255 at `n = 8` and `512`
  of 1023 at `n = 10`; **`mismatches = 0`** against the solver at four
  targets; `c = bits(a₆)` confirmed.
- **Holdout.** selection margin `+0.273 / +0.375 / +0.600` at
  `n = 6 / 8 / 10`; holdout margin `+0.667 / +0.403 / −0.788`. Variance
  ratio `1.345 / 1.382 / 1.409` against a no-effect prediction of 1 and a
  deterministic-curve prediction of 12.
- **No uniformly good curve.** Curves on the `D* = 2` floor for *every* one
  of `T` targets, swept over the whole curve space:

  | `n` | curves | `T=8` | `T=16` | `T=32` | `T=48` | `T=64` |
  |---|---:|---:|---:|---:|---:|---:|
  | 8 | 255 | 68 | 34 | **0** | 0 | 0 |
  | 10 | 1023 | 230 | 73 | 14 | 7 | **1** |

  Monotone in `T` and reaching zero.  At `n = 8` **no curve over the field**
  is on the floor for 32 targets; at `n = 10` one (`a₆ = 13`) survives 64.
  The decay is slower than independent targets would give — which is why the
  variance ratio sits at 1.4 rather than 1.0, i.e. consecutive targets are
  correlated and there *is* a small residual curve effect — but it is
  monotone and it terminates at zero, which is what the criterion of §3.1
  predicts and what the lever needs to be false.

**Gate verdict.** **G-R6 killed.** The selection margin is positive at every
size, as picking the minimum of a noisy statistic always is; the holdout
margin is negative at the largest size, and no size reaches the `≥ 0.5`
threshold at a variance ratio `≥ 2`.  **G-R6‴ supported** at zero mismatches.
**R6⁗ killed** — the uniformly-easy curve does not exist.

**Ledger delta.** R6 `killed` (by C and D, derived). R6′ `killed`. R6″
`killed`. R6‴ `supported`. R6⁗ `killed`.

**Class of the change (per `AGENTS.md` §3).** **Accounting**, not an
advance: the numbers that moved are the boundaries, and the algorithm the
attacker would run is no better than the one they had.  The selection-set
row of §6 is a textbook **relabelling** — it reports `2^{+0.00}`, i.e. the
floor, and is `2^{+1.13}` worse than baseline once scored on fresh targets.

---

## 6. The one table

One unit: `log₂` operations for the decomposition solve,
`cost(D) = rows·cols^{ω−1}` with the expectation taken over the `D*`
histogram (`degree_reduction::log2_expected_cost`), `ω = 2.807`.
`n = 10`, `l = 5`, 12 selection + 12 disjoint holdout targets.

The **correct** column is not decorative.  Every cell counted here is a
refutation, and each one is cross-checked against an independent brute-force
decomposition search over `V × V` (`pc_degree_harness::is_decomposable`):
a target is admitted only when the brute force finds no decomposition, and
the Macaulay certificate then proves the same fact a second way.  A target
that decomposes is *excluded*, not counted as a cheap solve — that is why the
per-curve target counts in the JSON are below 12.

| variant | log₂ ops | ratio to floor | correct | class |
|---|---:|---:|:--:|---|
| **floor** — every cell at the `D* = 2` Nullstellensatz floor | 13.82 | `2^+0.00` | n/a | floor |
| **reference** — Koblitz curve `a₆ = 1`, scored on the holdout | 23.07 | `2^+9.26` | yes | reference |
| best curve selected on target set A, scored on A | 13.82 | `2^+0.00` | yes | **selection (biased)** |
| **the same curve, scored on the disjoint holdout B** | **24.20** | **`2^+10.38`** | yes | accounting |
| **reference** — `2^N` enumeration of `V × V` | 13.32 | `2^−0.49` | yes | reference |

Read the third and fourth rows together: the selected curve sits *exactly*
on the floor where it was chosen and `2^{+1.13}` **above the unmodified
baseline** where it was not.  That gap is the whole of lever L5.

Two honesty notes on this table, both of which follow the thread's own
earlier lessons:

- The `2^N` enumeration row is *below* the floor at these sizes
  (`2^{−0.49}`), i.e. exhaustive search beats the Gröbner solve outright on a
  10-variable system.  That is the degeneracy iteration 1 of the
  degree-reduction thread walked into and it is why the Gröbner column is
  only meaningful here as a ratio to its own floor, never as an attack.
- The ρ reference of Boundary A is in a *different* unit (group operations
  on the real curve) and is deliberately not folded into this table.  It
  bounds the *search*, not the solve; mixing the two is what §5 of
  `AGENTS.md` warns about.

---

## 7. Honest limitations

1. **`D*` is measured at `m = 2`, `n ≤ 10`.** The `m = 3` system is degree 6
   in `3ℓ` variables, where the dense Macaulay tower has essentially no
   multiplier budget — the same structural reach limit that left L2
   `blocked` in iteration 5 of the degree-reduction thread.  `m ≥ 3` is
   covered here by Boundary C (exact, and verified) rather than by a `D*`
   measurement, and that is a genuinely weaker kind of evidence for the
   affine quantity.
2. **Boundary C is now computed for `m ∈ {2,3,4,5}`, still not proved for all
   `m`.**  The gap is exactly 2 at each, and `m = 4, 5` are the sizes that
   matter (index calculus is only asymptotically interesting from `m = 3`), so
   this is much stronger than the `m ≤ 3` the first version rested on.  But it
   remains a computation at four values of `m`, not an induction: the pattern
   `top = m(m−1)`, `a₆ ≤ m(m−1) − 2` is unproved, and `S₇` and beyond were not
   reached (`S₆` already has 190252 monomials).
3. **The 262 curves on the 263-floor are never written down.** Boundary C
   makes their `a₆` values irrelevant to the question (their `d_reg` is
   `E`'s), and the containment argument of §3 covers their `D*` distribution
   — but the explicit enumeration is not done.  Doing it needs
   `H_{−7·263²}` mod 2 (degree 262, coefficients of several thousand bits) or
   a degree-34584 division-polynomial factorisation over `F_{2^131}`; see §8.
4. **The residual curve effect is now explained, and it is not about the
   solving degree** (EXP-R6c, iteration 3).  Iteration 1 flagged a variance
   ratio of `1.4` and one curve (`a₆ = 13`) holding the floor over 64 targets
   at `n = 10`, and left it as an unexplained remainder.  It is decomposition
   yield.  The exact criterion partitions a curve's targets three ways — `a₆`
   escapes `S^⊥` (refutes at the floor); `a₆ ∈ S^⊥` and the target decomposes
   (satisfiable, skipped, and for an attacker a *success*: a relation); `a₆ ∈
   S^⊥` and it does not decompose (refutes **above** the floor, the only costly
   case).  The third is squeezed by the second, and the squeeze is nearly
   deterministic: `ρ_s = −0.9801` between the decomposable count and the
   above-floor count over all 255 curves × 240 targets at `n = 8`, and
   `−0.9648` against the rate.  So the "curve effect" is variation in relation
   yield, a different quantity from `d_reg`, and one where more yield helps an
   attacker for reasons unrelated to the solving degree.  Over the *full*
   target set no curve is uniformly easy at all: `0` of 255 at `n = 8`, and
   **`0` of 1023 over all 992 targets at `n = 10`**.  The curve that prompted
   this caveat, `a₆ = 13`, breaks by `T = 256`: 171 refutable targets, 137 at
   the floor and **34 above it**.  Its clean run over the first 64 was
   small-sample, and the caveat it generated is retired.
5. **Single field representation.** The sweeps use the first irreducible
   polynomial of each degree, and ECC2K-130's own field is a permuted
   type-II ONB, not a polynomial basis.  The class structure is
   basis-independent (it is a statement about `End`), but `dim S` and the
   `D*` histograms are not, and a basis sweep was not run.

---

## 8. What would change the verdict

- ~~**A counterexample to Boundary C at `m ≥ 4`.**~~ **Done — EXP-R6b,
  iteration 2, and it did not find one.**  `S₅` and `S₆` were built and
  profiled: `a₆` stays exactly two Boolean degrees below the leading form at
  `m = 4` and `m = 5`, the same gap as at `m = 2, 3`.  What is left of this
  item is the induction: a proof that the gap is 2 for *all* `m`, or a
  counterexample at `m ≥ 6`.  `S₇` needs a resultant of two `S₄`-sized
  quartics one level up, so the next step is `S₇ = Res(S₄, S₅)` — tractable,
  since `S₆` took seconds.
- **A target-independent good curve.** The criterion of §3.1 says the good
  set is target-dependent.  A curve whose `a₆` escapes `S^⊥(x_R)` for a
  constant fraction of *all* targets, uniformly in `n`, would contradict it
  and is directly searchable with `syzygy_mechanism` inverted.
- **The explicit 263-floor.** If it is ever wanted for another purpose
  (Boundary C says it is not wanted for this one): compute
  `H_{−484183}(X) mod 2` and root-find over `F_{2^131}`.  `h = 262`, so the
  polynomial has degree 262 with coefficients of several thousand bits; the CRT
  method is standard, and `cryptanalysis::hilbert_class_poly` is the place it
  would go.  The alternative — factoring the degree-34584 263-division
  polynomial over `F_{2^131}` into its 264 degree-131 kernel factors — is the
  more expensive route.

---

## References

- **D. Kohel**, *Endomorphism rings of elliptic curves over finite fields*,
  PhD thesis, Berkeley 1996 — the volcano structure of §2B.
- **S. Galbraith, F. Hess, N. Smart**, *Extending the GHS Weil descent
  attack*, EUROCRYPT 2002 — the isogeny-walk move being scored.
- **I. Semaev**, *Summation polynomials and the discrete logarithm problem
  on elliptic curves*, eprint 2004/031.
- **S. Galbraith, S. Gebregiyorgis**, *Summation polynomial algorithms for
  elliptic curves in characteristic two*, INDOCRYPT 2014 — the binary `S₃`
  and `S₄` used here.
- **M. Bardet, J.-C. Faugère, B. Salvy**, *On the complexity of the `F₅`
  Gröbner basis algorithm*, 2015 — the degree of regularity Boundary C
  invokes.
- **D. Bernstein, L. De Feo, A. Leroux, B. Smith**, *Faster computation of
  isogenies of large prime degree*, ANTS 2020 — the √élu bound of §2B.
- **D. Bailey et al.**, *Breaking ECC2K-130*, eprint 2009/541 — the source of
  the curve and of the negation-plus-Frobenius ρ technique whose `√(2·131)`
  speedup factor §2A applies.  The `2^60.81` there is computed from
  `√(πr/2m)` at `m = 262`, not quoted from the paper.
- `RESEARCH_DEGREE_REDUCTION.md` — levers L1–L4, the cost model and the
  baselines this thread reuses.
