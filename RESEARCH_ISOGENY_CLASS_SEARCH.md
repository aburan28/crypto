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
  DLP (**A**);
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

Pollard ρ on the 130-bit subgroup costs `√(πr/4) = 2^64.33` group
operations with the negation map, or **`2^60.81`** after the further
`√131` Frobenius speedup — combined `√(πr/524)`, the figure the
ECC2K-130 effort actually uses.  So:

> **Boundary A.** Visiting every vertex of the isogeny class costs `2^4.25×`
> a full ρ run, at one operation per vertex and before any Gröbner work.  A
> search over more than `2^{−4.25} ≈ 5.3%` of the class is beaten by ρ on
> enumeration alone, whatever it finds.

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

At `m = 3` the same holds with room to spare.  Expanding
`S₄ = Res_X(S₃(X₁,X₂,X), S₃(X₃,x_R,X))` in `a₆`, with
`A₁ = (X₁+X₂)², B₁ = X₁X₂, A₂ = (X₃+x_R)², B₂ = X₃x_R`:

```
  S₄ = [A₁(X₃x_R)² + A₂(X₁X₂)²]²                 a₆-free, Boolean degree 3
     + (A₁B₂ + A₂B₁)(B₁(X₃x_R)² + B₂(X₁X₂)²)     a₆-free, Boolean degree 6
     + a₆  · (A₁B₂ + A₂B₁)(B₁ + B₂)              Boolean degree ≤ 5
     + a₆² · (A₁ + A₂)²                          Boolean degree 1
```

Boolean degrees are counted with squaring free, since Frobenius is
`F_2`-linear and `x² = x` on Boolean coefficients — this is why the `m = 3`
system is degree 6 in `3ℓ` variables rather than 12.  With `x_R` constant that
makes `A₁, A₂, B₂` degree 1 and `B₁` degree 2, from which the four
annotations above follow.  The top degree is 6 and `a₆` reaches at most 5.
The identity itself is verified over 300 `(x₁,x₂,x₃,x_R,a₆)` combinations in
`s4_a6_expansion_is_quadratic_with_subleading_coefficients`.

> **Boundary C.** The leading-form ideal of the descended Semaev system is
> independent of the curve, for `m ∈ {2, 3}`.  Therefore the degree of
> regularity in the Bardet–Faugère–Salvy sense — a Hilbert-series invariant
> of the leading forms — is **constant on the entire isogeny class**, and so
> is the degree at which any top-degree cancellation first becomes
> *available*.  Changing the curve can make a fall's remainder vanish, i.e.
> **lose** a fall; it cannot create one at a lower degree.  The lever's only
> possible effect on `d_reg` is adverse.

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
| **R6″** | The exhaustive search over the class is *feasible* | **`killed`** | Boundary A (`2^65.06` vertices vs `2^60.81` ρ) and B (263 reachable). The search that does terminate covers `2^{−57}` of the class. |
| **R6‴** | `D* = 2` density over curves is `1 − 2^{−dim S}` per target | **`supported`** (exact) | `dim S = 1` at `n ∈ {8, 10}`, escape count `128/255` and `512/1023`, matching `(2^n − 2^{n−1})/(2^n − 1)`; mismatches `0`. |
| **R6⁗** | Some curve is on the `D* = 2` floor for **every** target — the uniformly-easy curve an isogeny walk would need | **`killed`** | EXP-R6. Survivor count `68 → 34 → 0` at `n = 8` over `T = 8/16/32`, and `230 → 73 → 14 → 7 → 1` at `n = 10` over `T = 8/16/32/48/64`. Monotone and reaching zero. |

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
| `…::all_traces` | exact point counting, for the class-size cross-validation |
| `examples/isogeny_class_search.rs` | EXP-R6 driver → `experiments/isogeny_class_search.json` |

`D*` is measured by the thread's existing `pc_degree_harness::refutation_scan`,
and costed by `degree_reduction::log2_expected_cost` (Jensen-correct over the
histogram, not the mean), so the unit matches every other lever's.

Run: `cargo run --release --example isogeny_class_search`.

---

## 5. Iteration log

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

- **A.** class `2^65.06`, ρ `2^60.81` → enumeration costs `2^4.25×` ρ.
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
   `blocked` in iteration 5 of the degree-reduction thread.  `m = 3` is
   covered here by Boundary C (exact, and verified as an identity) rather
   than by a `D*` measurement, and that is a genuinely weaker kind of
   evidence for the affine quantity.
2. **Boundary C is proved for `m ∈ {2, 3}` by direct expansion**, not for
   all `m`.  The resultant recursion makes it plausible that `a₆` stays
   sub-leading at every `m`, but that is a conjecture here, not a theorem.
3. **The 262 curves on the 263-floor are never written down.** Boundary C
   makes their `a₆` values irrelevant to the question (their `d_reg` is
   `E`'s), and the containment argument of §3 covers their `D*` distribution
   — but the explicit enumeration is not done.  Doing it needs
   `H_{−7·263²}` mod 2 (degree 262, coefficients of several thousand bits) or
   a degree-34584 division-polynomial factorisation over `F_{2^131}`; see §8.
4. **There is a small residual curve effect, and it is not zero.** The
   variance ratio sits at `1.4` rather than `1.0`, and the survivor count
   decays more slowly than independent targets would give — at `n = 10` one
   curve (`a₆ = 13`) is still on the floor at `T = 64`.  Consecutive targets
   are correlated and this design does not fully separate that from a genuine
   curve effect.  What is established is that the effect is far too small to
   matter — the holdout margin is *negative* at the largest size, and the
   survivor count is monotone and reaches zero at `n = 8` — not that it is
   exactly zero.  Calling it zero would be overclaiming.
5. **Single field representation.** The sweeps use the first irreducible
   polynomial of each degree, and ECC2K-130's own field is a permuted
   type-II ONB, not a polynomial basis.  The class structure is
   basis-independent (it is a statement about `End`), but `dim S` and the
   `D*` histograms are not, and a basis sweep was not run.

---

## 8. What would change the verdict

- **A counterexample to Boundary C at `m ≥ 4`.** If `a₆` reaches the
  leading form of `S₅`, `d_reg` becomes curve-dependent and the lever
  reopens — at the one `m` where index calculus is asymptotically
  interesting.  This is the single highest-value follow-up, and it is a
  symbolic computation, not a search.
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
  the curve and of the negation-plus-Frobenius ρ technique whose combined
  `√(πr/524)` cost §2A applies.  The `2^60.81` there is computed from
  `√(πr/4)/√131`, not quoted from the paper.
- `RESEARCH_DEGREE_REDUCTION.md` — levers L1–L4, the cost model and the
  baselines this thread reuses.
