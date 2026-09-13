# Reducing the solving degree — the offensive side of the FFD program

**Status:** open research thread, opened 2026-09-11
**Builds on:** `RESEARCH_FFD_PROOF_COMPLEXITY.md` §5 ("the attacker's
corollary" — stated there, never built), `RESEARCH_FFD_WORKFLOW.md`
(ledger + gate discipline, which this thread copies verbatim),
`src/cryptanalysis/{pc_degree_harness,descent_algebraic,descent_lowgamma}.rs`.
**Module:** `src/cryptanalysis/degree_reduction.rs`
**One-line thesis:** the Gröbner solving degree `D*` is not a property of
the ECDLP instance — it is a property of the **presentation** we hand the
solver. The ideal is fixed; the generating set, the variables, and the
determination ratio are ours to choose. This thread enumerates the
presentation levers, scores each one against a baseline it must actually
beat, and asks whether any of them moves `D*` enough to matter.

---

## 1. Why this thread, and why now

The FFD program is *defensive*: it is assembling a theorem that the
first-fall-degree assumption is **false** for generic curves, so that
index calculus stays exponential. Its central measured law is

```
   D*  decreases monotonically with  Δ_low(Σ),
   at ≈ −7.6 degrees of D* per unit early Hilbert defect
   (50 cells, 2n' ∈ {4,…,14}, pooled Spearman ρ_s = −0.79).
```

`RESEARCH_FFD_PROOF_COMPLEXITY.md` §5 notes in passing that this hands the
attacker a search target — *maximise `Δ_low`* — and that every known
speedup (subfield/Koblitz bases, symmetrisation, GHS-amenable bases) is an
instance of "inject low-degree algebraic relations." That paragraph was
never turned into experiments. This thread is that work.

The pay-off is two-sided and both sides are useful:

- **If a lever works**, we have a concrete route to a lower solving
  degree, and — because `Δ_low` is a few low-degree Macaulay ranks — a
  *polynomial-time screen* that says which curve/basis/formulation
  combinations are worth attacking, evaluable at cryptographic `n` where
  the attack itself is unobservable.
- **If every lever fails**, each failure is a measured statement about
  *why* the presentation cannot be improved, which is exactly the content
  the defensive theorem needs and currently assumes.

So no iteration of this thread can end in limbo. That is the point.

---

## 2. The lever taxonomy

`D*` depends on the pair (ideal, presentation). Five levers, ordered by
how much they are allowed to change — L1–L4 re-present one fixed curve, L5
changes the curve itself:

| # | Lever | What it changes | What it costs | Status |
|---|---|---|---|---|
| **L1** | **Factor-base structure** — subfield, Koblitz, sparse normal basis | the *ideal itself* (multiplicative closure injects relations) | only works on special curves/fields | **measured** by the FFD program: Subfield mean `D*` 2.04 vs Random 3.53 at `2n'=n`, and `Δ_low(Subfield)/Δ_low(Random)` diverges 6.7 → 67 |
| **L2** | **Symmetrisation** — solve over the elementary symmetric variables instead of the `Xᵢ` (Faugère–Gaudry–Huot–Renault) | the *variables*: `3ℓ` at degree 6 → `9ℓ−3` at degree 3 | 3× the variables for ½ the degree | **blocked on reach, iteration 5 — but bounded.** The exact Macaulay-width crossover is `ℓ = 6`; below it symmetrisation is the *wider* presentation, and `ℓ ≤ 3` is all a dense Macaulay scan reaches |
| **L3** | **Hybrid slicing** — guess `k` variables, solve `2^k` slices, raising `ρ = #eqs/#vars` | the *determination ratio* | `2^k` multiplicative | **killed, iteration 1** (degenerate optimum) — but one-sided guessing is the cheapest route to the `D*=2` floor, and iteration 3 found it *dominates* the mutant route once guessing is allowed |
| **L4** | **Degree falls (mutants)** — add the *nonzero* low-degree remainders of top-degree cancellations to the generator set, so `x_k · g` rows become available a degree early | the *generating set*; ideal and variables unchanged | the extraction's own climb to degree 3 | **supported on degrees, regime-dependent on cost — iteration 2.** `D*` drops 4.00 → 2.00 on the generic family; net of extraction cost it pays only where the base degree is high |
| **L5** | **The curve** — move along the isogeny class, solving on an isogenous `E'` and transporting the answer back (Galbraith–Hess–Smart) | the *curve*, hence the ideal's inhomogeneous part | computing the isogeny, and finding the target | **killed, iteration 9** — `a₆` enters the descended system strictly *below* the leading form, so `d_reg` is constant on the class; and only 263 of ECC2K-130's `2^65.06` vertices are reachable. `RESEARCH_ISOGENY_CLASS_SEARCH.md` |

L1 is the known part of the map and is not this thread's subject. L2, L3
and L4 are presentation changes that apply to *any* curve, which is what
makes them worth measuring. L5 is the odd one out — it is the only lever an
attacker gets for free, since nothing requires the curve you attack to be
the curve you were given — and the only one whose verdict is settled by a
derivation rather than a sweep.

> **Can L1 be acquired by moving to an isogenous curve?**  No — settled
> in `RESEARCH_ISOGENY_DEGREE_SEARCH.md`, negatively and for three
> independent reasons.  The curve parameter enters the Semaev system
> only as its **constant term**, so the leading forms, and with them the
> degree of regularity, are identical across an entire isogeny class
> (checked exhaustively over every ordinary binary curve at `n ≤ 9` and
> every class member at `n ≤ 17`: zero improvements in 2 994 curve
> comparisons).  For ECC2K-130 specifically the class holds `2^65.06`
> curves — `19×` more work than the rho it would improve — and because
> 131 is prime it contains exactly **one** curve with subfield
> structure: ECC2K-130 itself.

### 2.0 Syzygies are not degree falls (the iteration-1 error, corrected)

Iteration 1 described L4 as "add the relation EXP-J identified." That was
wrong, and the distinction it missed is the one that decides whether L4
exists at all.

EXP-J's relation is `Σ_i ℓ_i · f_i ≡ 0` — a **pure syzygy**. It vanishes
identically, so it yields no polynomial; adding it to the generating set
adds nothing. Pure syzygies cost the solver zero reductions, a constant
factor (exactly what F5's criterion removes), not a degree.

A **degree fall** is the other outcome of the same cancellation: a
combination of degree-`D` rows whose degree-`D` part cancels leaving a
**nonzero** remainder of degree `≤ D−1`. Its value is not that the solver
cannot find it — at degree `D` the remainder is already in the row space —
but that once `g` is a *generator*, the rows `x_k · g` are available at
degree `D`, and those are degree-`D+1` products of the original generators.
The Macaulay tower genuinely accelerates. This is the mutant mechanism of
MutantXL and the degree-fall strategy used against HFE.

The two populations turn out to be wildly different sizes on these systems:
`degree_reduction::extract_degree_falls` counts **0–4** pure syzygies (EXP-J's
population, consistent with its `δ(3)=1`) against **20–112** genuine falls,
per instance. The falls were never measured before.

### 2.1 The objective function, and its calibration limit

The organising hypothesis is that all four levers act through one channel:

> a lever lowers `D*` **iff** it raises the early Hilbert defect `Δ_low`,
> at the rate the FFD program already measured.

If that holds across levers, `Δ_low` becomes the screen — we can rank
reformulations without ever running a Gröbner basis. Iteration 1 found
the first obstruction to this and it is worth stating up front:

> **`Δ_low` is only calibrated within a fixed system shape.** It is
> normalised by `cols(D_low)`, the monomial count, so it is not comparable
> across systems with different variable counts. In EXP-R1 the `Δ_low`
> column moves *non-monotonically* in `k` while `D*` falls monotonically
> — e.g. at `N=10`, balanced, `Δ_low` runs 0.035 → 0.199 → 0.020 while
> mean `D*` runs 2.75 → 2.45 → 2.00. The screen is not wrong; it is being
> read outside the regime it was fitted in.

**Corrected by iteration 4 (EXP-R4′), and the correction cuts deeper than
the original claim.** The diagnosis above — "wrong units, needs a
shape-corrected denominator" — is wrong. Five candidate normalisations were
measured against `D*` on three designed groups, and the pattern is stable
over four seeds:

- Across shapes, *pooled*, the size-normalised variants look excellent
  (`ρ_s` ≈ −0.89 to −0.91). That would say `Δ_low` is already shape-robust
  and iteration 1 merely read it noisily.
- But **hold `vars` fixed and the same variants collapse to `ρ_s` ≈ −0.16
  to −0.34.** Every variant that clears the bar pooled fails it once size
  is controlled.

So the strong pooled figure is a **size proxy**: `D*` and a size-normalised
defect both trend with `N`, and pooling across sizes reads that shared trend
as correlation. The problem was never the denominator — the quantity being
rescaled does not carry enough structure once the size trend is removed.
No renormalisation can fix that, which is why R4′ is killed rather than
solved.

*Worth checking upstream:* the FFD program's own `ρ_s = −0.79` was likewise
pooled across `2n' ∈ {4,…,14}`. It may share this property. That is a
caveat to test, not a refutation — what EXP-G measured (cell means over ten
operating points) is related to but not identical with what EXP-R4′
measures, and this thread has not re-run their sweep.

**Consequence for L2:** do not lean on `Δ_low` to score symmetrisation.
Compare symmetrised against raw at **matched `(vars, eqs, ρ)`** and score on
measured `D*` and total work — the way iterations 2 and 3 did.

---

## 3. The baseline every lever must beat

A lever is not interesting because it lowers `D*`. It is interesting
because it lowers *total work*. Two baselines, and the second is the one
that bites:

1. **The direct solve** (`k = 0`, raw variables, raw generators).
2. **Brute-force enumeration of `V × V`** — `2^N` Semaev evaluations,
   `N = 2n'`, charged generously at `O(1)` amortised per point
   (`degree_reduction::log2_enumeration_cost`).

Baseline 2 exists because of a trap iteration 1 walked into and had to
back out of. At the sizes where `D*` is measurable (`N ≤ 20`; iteration 1
wrote `N ≤ 16`, and EXP-R3b later reached 20 — see §7), `2^N` is *small*,
so any cost model with a `2^k` term will happily slide its optimum to the
largest `k` available and report a large "saving" against the direct
solve. What it has actually found is that exhaustive search is
fastest on a 14-variable system — which says nothing about cryptographic
`n`. Every sweep in this thread therefore reports
`HybridSweep::optimum_is_interior`, and a boundary optimum is recorded as
a **degenerate** result, not a win.

Cost model: `cost(D, N) = cols(N, D)^ω` for the hybrid sweeps (EXP-R1),
and the **row-aware** `rows(D)·cols(D)^{ω−1}` with
`rows(D) = eqs·cols(D−2)` from EXP-R3 onward. The second is not a
refinement for its own sake: L4 buys a lower degree by *adding
generators*, and a `cols^ω` model is blind to the row count, so it would
score that trade as free. Either way the expectation is taken over the
*histogram* of `D*` rather than its mean (Jensen — a spread of degrees
costs strictly more than its average; asserted in
`cost_model_uses_the_histogram_not_the_mean`). Results are reported at
`ω = 2.807` (Strassen, close to the dense bit-packed elimination this repo
actually runs) and `ω = 2.0` (the optimistic sparse limit, which is the
*hardest* case for any degree-reduction lever because it discounts the
collapse most).

---

## 4. The prediction ledger

Status ∈ {`open`, `supported`, `killed`, `blocked`}. "Supported" means
"survived a genuine refutation attempt at the current reach", never
"proved". Same convention as `RESEARCH_FFD_WORKFLOW.md` §2.

| # | Prediction | Status | Last evidence |
|---|---|---|---|
| **R1** | **Hybrid slicing (L3) has an interior cost optimum that beats `2^N` enumeration**, at every operating point, with the collapse fraction not rising in `N` | **`killed`** | EXP-R1, iteration 1. At every one of 9 cells (`N ∈ {10,12,14}` × 3 guess patterns) and both `ω`, the optimum sits at the largest `k` scanned — the model is choosing exhaustive search. The margins it reports (−0.29 → +0.19 bits at `ω`=2.807) are artifacts of where the scan was truncated. And the collapse fraction *rises* with `N` for 2 of 3 patterns. |
| **R1′** | Among guess patterns, **one-sided** guessing (all `k` bits from the `X₁` half) reaches the `D* = 2` floor at the fewest guessed bits, at `c = k₂/N = 1/2` independent of `N` | **`supported`** | EXP-R1: `c(one-side) = 0.500` at `N = 10, 12, 14` — exactly `k₂ = n'` every time — vs balanced 0.800/0.667/0.786 and spread —/0.750/0.714. Flat where the others drift, and identical across seeds 7/11/23 (9 cells, no exceptions). |
| **R2** | **Symmetrisation (L2) lowers `D*`** at matched shape | **`blocked`** (structurally, not incidentally) | EXP-R2, iteration 5. The two presentations are never simultaneously measurable: at `ℓ = 2` the eliminated system has degree-6 generators in 6 variables, so its multilinear Macaulay tower has **zero multiplier budget** and cannot refute at all; at `ℓ ≥ 3` the symmetrised system's `9ℓ−3` variables exceed the elimination budget before it refutes. Raising the budget does not help — see R2′. |
| **R2′** | Symmetrisation's Macaulay width crosses below elimination's only at **`ℓ = 6`** | **`supported`** (exact, no solver) | EXP-R2: `cols(9ℓ−3, 3)` vs `cols(3ℓ, 6)` — 576 vs 64 at `ℓ=2`, 12384 vs 9949 at `ℓ=5`, 22152 vs 31180 at `ℓ=6`, and the gap widens to 7× by `ℓ=10`. Below `ℓ=6` symmetrisation trades 3× the variables for ½ the degree and **loses**. |
| **R3** | **Adding degree falls (L4) as explicit generators lowers `D*`** at matched targets | **`supported`** | EXP-R3, iteration 2. All three families, 3 operating points each, 8 matched targets per cell: `D*` strictly lower wherever there was headroom. Generic (Random) family **4.00 → 2.00** on 8/8 targets at `N = 12, 14`. `worsened = 0` everywhere, as the ideal-membership invariant requires. |
| **R3′** | The `D*` drop **survives its own cost** — extraction must climb to degree 3, so the net saving must still be positive | **regime-dependent** | EXP-R3: net `+1.44` bits mean on Random (positive at every `N`, seeds 7/11/23 give +1.44/+1.45/+1.49); `−0.60` on Coordinate (sign varies); `−6.36` on Subfield (**killed** — the system already solved at `D* ≈ 2.1`, so the climb to 3 is pure overhead). L4 pays where the system is hard and costs where it is easy. |
| **R3b** | L4's **net saving against the raw Gröbner solve keeps growing** at `N = 16, 18, 20` | **`supported`** (Random, seed-robust); **`unstable`** (Coordinate); **`killed`** (Subfield) | EXP-R3b, iteration 7. Random mean `1.44 → 4.44` bits (seed 7), `1.45 → 4.44` (11), `1.49 → 4.16` (23); the reach cells are `+5.11/+3.87/+4.34`, `+5.10/+3.87/+4.34`, `+4.28/+3.87/+4.34` — **`N = 18` and `N = 20` identical to 2 d.p. on all three seeds**. Coordinate's verdict *flips* KILLED/SUPPORTED/KILLED on `N = 18` values `−0.18/+0.89/−0.12`, so it sits on the gate boundary and the seed picks the answer; recorded as unstable rather than resolved. Subfield `−5.3…−9.1` everywhere. |
| **R3c** | The mutant cascade drives the augmented system to the **`D* = 2` floor at every `N`** | **`killed`** | EXP-R3b. At `N ≤ 16` saturation runs several productive rounds (Random: `sat.eqs = 137` at `N = 16` from 32 falls) and reaches `D* = 2`. At `N ≥ 18` it **stops after one round** — `sat.eqs = 54 = 18 + 36` at `N = 18`, `60 = 20 + 40` at `N = 20` — and the augmented system lands at `D* = 4`. L4's degree saving **halves**: `2.25` at `N = 16` → `1.00` at `N = 18, 20`. Identical on all three seeds. The cascade, not the first batch of falls, is what reached the floor, and the cascade dies. |
| **R3d** | L4 **narrows the gap to the `2^N` enumeration boundary** (§3, baseline 2) | **`killed`** | EXP-R3b. **0 of 54 cells beat `2^N`** (3 seeds × 3 families × 6 sizes). On Random — the one family where L4 pays — the augmented margin runs `−10.20, −9.96, −10.12, −9.51, −12.69, −12.09` over `N = 10…20`: 9.5–12.7 bits *behind* brute force, losing **3.18 bits across the `N = 16 → 18` break** alone and `−0.378` bits per size step overall. The pre-registered G-R3b passes on this same data because it scores the saving against the *raw Gröbner solve*: a race between two routes that both lose. **Scope:** the widening is Random-specific. Subfield's augmented margin *narrows* monotonically (`−10.49 → −7.53`) and Coordinate's is non-monotone, so R3d is killed by "no cell beats the boundary", not by a universal widening — see R3e for what the Subfield trend actually is. |
| **R3e** | *(observation, not a pre-registered prediction)* On the **Subfield** family the **raw, unaugmented** solve is converging on the `2^N` boundary: margin `−4.99 → −1.19` bits over `N = 10…20` | **`observed`; explained, and not an advance** | EXP-R3b. At `N = 20` the plain Gröbner solve on a Subfield factor base is only `1.19` bits behind enumeration, and the trend would cross near `N ≈ 22–24`. This is **arithmetic, not cryptanalysis**: Subfield's `D*` is pinned at `2.12–2.25`, i.e. at the Nullstellensatz floor, so its solve is `poly(N)` while enumeration is `2^N` — a crossover is guaranteed and measures the family's *degeneracy*, which is the FFD program's own L1 result (§2), not a lever of this thread. Recorded because a bare "Gröbner overtakes brute force at `N ≈ 23`" would read as a breakthrough and is not one. Note also that **L4 makes this family dramatically worse** (augmented margin `−10.49` vs raw `−4.99` at `N = 10`): on the only family approaching the boundary, this thread's lever actively hurts. |
| **R4** | **Every lever acts through `Δ_low`**: pooled across lever-generated systems, `ρ_s(Δ_low, D*) ≤ −0.6` | **`killed` as stated** | A pooled *instance-level* `ρ_s` over mixed sizes cannot establish this (R4″), so the gate as written measures nothing. This is a kill of the **statistic**, not of the defect: EXP-R4b shows the family-level version of the same correlation is strong and size-robust (R4c). Scoring a lever still needs matched shape and total work — but "score the lever *family*, not the instance" is now an open question, not a closed one. |
| **R4′** | A **shape-corrected** defect exists that is comparable across systems with different variable counts | **`killed`** | EXP-R4′, iteration 4. Five normalisations × three designed groups × four seeds. Nothing clears `ρ_s ≤ −0.6` on all three; the un-normalised variants fail cross-shape as expected (the control that shows the test has power), and the normalised ones pass *pooled* only. **Sharpened iteration 8 (R4e):** within a family at fixed size the defect is not merely weak, it is *positively* correlated with `D*`. Picking the easiest instance in a group by lowest `Δ_low` would systematically pick the hardest. No renormalisation can fix a sign. |
| **R4″** | The strong pooled **instance-level** defect↔`D*` correlation — one point per target — is a **size proxy**, not structure | **`supported`, scope corrected iteration 6** | EXP-R4′: pooled `ρ_s` −0.89…−0.91 collapses to **−0.16…−0.34** once `vars` is held fixed, for every variant that passed pooled. Stable over seeds 7/11/23/41 at 48 targets/cell. EXP-R4b reproduces it on the same snapshot (instance level, size-controlled: **−0.32…−0.35**; first published as −0.12…−0.35, see the iteration-7 correction). **As originally written this row said "the strong pooled defect↔`D*` correlation", with no unit qualifier, and iteration 4 read that as reaching the FFD program's published figure. It does not — see R4b.** |
| **R4b** | The FFD program's own published `Δ_low ↔ D*` law (EXP-G, `ρ_s = −0.79`) shares R4″'s defect and is **also a size proxy** | **`killed`** | EXP-R4b, iteration 6. On EXP-G's own 50 cells the law survives every size control: mean per-block **−0.6998**, blocked rank **−0.6983**, fixed-effects **−0.7499** against pooled −0.7929. In the ECDLP-relevant critical regime (`2n' = n`) size control makes it *stronger*, not weaker: **−0.9124** controlled vs −0.7781 pooled. This thread raised the flag; this thread withdraws it. |
| **R4c** | `Δ_low` is a **family-level** discriminator (it ranks factor-base constructions) but not an **instance-level** one (it does not rank targets within a construction) | **`supported`** | EXP-R4b Panel B, on iteration 4's *own* 2016-cell snapshot, size-controlled both ways: instance level **−0.35** (mean per-block) / **−0.33** (blocked rank, corrected iteration 7 — first published as −0.12); family level **−0.85** (blocked rank), with the three families ordered correctly in **4 of the 4** size blocks where `D*` varies at all (the other two are floored at `D* = 2`). Same cells, same defect, same `D*` — only the unit of analysis differs. |
| **R4d** | Lever **L4 acts *through* `Δ_low`** — the defect mediates the lever's effect, at either the paired or the presentation-family level | **`killed`, with the sign inverted** | EXP-R4d, iteration 8. Three presentations of the *same system at the same variable count on the same target* (raw / one saturation round / saturated), so nothing is cross-shape and R4′'s objection cannot apply. Both statistics come out **strongly positive** where the gate required `≤ −0.6`: paired blocked rank **+0.8385**, presentation-family blocked rank **+0.8704**. This is a harder kill than the gate anticipated — not "no relation" but "the relation runs the other way". |
| **R4e** | *(discovered, not predicted)* `Δ_low` **changes sign with the grouping**: negative between factor-base families, **positive within a family at fixed size** | **`supported`, mechanism identified** | EXP-R4d. Within-cell `ρ_s(Δ_low, D*)` is positive in **10 of 10** decidable cells (`+0.62…+1.00`), while the between-family law at the same sizes is strongly negative (R4c, `−0.85`). Textbook Simpson's paradox, with a cause: `Δ_low` sums degrees `≤ 3`, so a system can only show a degree-3 defect if its tower *reaches* degree 3 deficiently. A target refuting at `D* = 2` never exercises degree 3, so its cutoff-3 defect is ~0 **by construction** — e.g. Coordinate `N = 14`, `D* = 2` targets average `Δ = 0.00213` against `0.03511` for `D* = 4`. Inside a homogeneous group the statistic is partly a proxy for "did this instance need degree 3", which *is* `D*`, positively. |
| **R4f** | *(discovered, not predicted)* `Δ_low` **cannot score a lever's output at all**, independent of any correlation | **`supported` (mechanical, asserted in a test)** | EXP-R4d: **0 of 192** saturated targets retain any cutoff-3 defect — saturation drives `Δ_low` to *exactly* zero, every time. It must: the cutoff-3 defect **is** the space of degree-3 falls, and saturation adds precisely those as generators. So the lever's action is to *consume* the quantity a screen would measure, and "the defect of the post-lever system" is identically 0 however much or little good the lever did. Pinned by `saturation_consumes_the_early_defect`. |
| **R5** | **Levers compose**: one-sided guessing plus the mutant route beats guessing alone | **`killed`** | EXP-R5, iteration 3. The pre-registered gate (collapse fraction `c < 1/2`) is **degenerate** — the composed route hits `c = 0` in every cell, because mutants reach the floor with no guessing at all. Scored on total work instead (G-R5′): composed loses to raw guessing by a **flat −2.87 bits** at every `N` and seed, and neither route beats `2^N`. |
| **R5′** | Mutants and guessing are **substitutes, not complements** — both drive the system to `D* = 2`, and guessing gets there more cheaply per unit work | **`supported`** | EXP-R5: at `k = 0` the mutants are worth `+1.1…+1.75` bits on Random (iteration 2's result), but the moment guessing is allowed the advantage inverts and stays inverted at every `k > 0`. The gap is flat in `N`, so it is structural, not a small-size artifact. |
| **R6** | **Some curve in the ECC2K-130 isogeny class has a materially lower solving degree** than the Koblitz curve (lever L5) | **`killed`** | Boundaries C and D of `RESEARCH_ISOGENY_CLASS_SEARCH.md`, both exact and needing no search: `a₆ = 1/j` is the whole curve-dependence of the descended Semaev system and it enters *below* the leading form (a constant at `m = 2`; Boolean degree `≤ 5` of `6` at `m = 3`), so `d_reg` is constant on the class; and factor-base structure (L1 included) is a property of the *field and subspace*, which an isogeny leaves untouched — so it is available on every curve in the class equally or on none. *(Reformulated after `RESEARCH_QUASI_SUBFIELD.md` landed: the first version said "needs a proper subfield, which `F_{2^131}` lacks", which is too strong, since quasi-subfield polynomials give subfield-like bases at prime `n`. The field-invariance argument does not depend on that question; the census also happens to find no quasi-subfield cell at `n = 131`.)* |
| **R6′** | The residual variation in the **affine** `D*` across curves is a *curve* effect an attacker can move to | **`killed`** | EXP-R6, iteration 9. The exact criterion `D* = 2 ⟺ a₆ ∉ S^⊥` with `S^⊥` fixed *by the target* holds with **0** disagreements over every curve at four targets; the best curve selected on one target set is `2^{+1.13}` worse than the unmodified baseline on a disjoint one. |
| **R6″** | An exhaustive search of the isogeny class is *feasible* | **`killed`** | The class has `2^65.06` vertices against ρ's `2^60.81` operations (`√(πr/2m)`, `m = 262`), so enumeration is `2^4.25×` a full ρ run — and still `2^0.24×` plain ρ, which is the version of the claim that does not depend on crediting ρ its automorphism speedup; and only `263` vertices are reachable without a degree-`146505763881528721` isogeny (kernel polynomial degree `2^56`). |

### Pre-registered gates

- **G-R1.** *Supported* if the best hybrid beats `2^N` enumeration at every
  measured `N` **and** the collapse fraction `c` is non-increasing in `N`.
  *Killed* if the cost optimum is at the scan boundary at every operating
  point (degenerate — the model is choosing exhaustive search), or if the
  margin is negative and not closing. *Blocked* if the margin is negative
  but closing with `N`.
- **G-R1′.** *Supported* if `c(one-side) ≤ c(other patterns)` at every `N`
  **and** `c(one-side)` is non-increasing. *Killed* if one-sided is ever
  worse than balanced, or if `c(one-side)` rises.
- **G-R2.** *Supported* if mean `D*` on the symmetrised system is strictly
  below the raw system at matched `(n, ℓ)` and matched targets, over ≥ 3
  operating points, with the gap not shrinking. *Killed* if the gap is ≤ 0
  at the largest size. **(Iteration 5: blocked — fewer than 3 cells are
  simultaneously measurable, for the structural reason in R2.)**
- **G-R2′** *(registered iteration 5)*. The Macaulay-width comparison
  `cols(9ℓ−3, 3)` vs `cols(3ℓ, 6)` is exact and needs no solver. *Supported*
  if a crossover exists at finite `ℓ` — which bounds where L2 can possibly
  pay, whether or not any solver reaches it.
- **G-R3.** *Supported* if adding the degree-fall generator strictly lowers
  mean `D*` at matched targets over ≥ 3 operating points. *Killed* if `D*`
  is unchanged — which would mean the solver was already finding the
  relation for free, and L4 is empty. (Scored on degrees as written;
  iteration 2 found this necessary but not sufficient, hence G-R3′.)
- **G-R3′** *(registered iteration 2, after G-R3 proved insufficient)*.
  Net log₂ saving `= cost(base at D*_base) − [cost(extraction) +
  cost(augmented solve)]`. *Supported* for a family if positive at every
  measured `N`; *killed* if negative at every `N`; *blocked* if the sign
  varies. Extraction is charged per saturation round on the system as it
  stood going in — booking only the cheap final solve would count the same
  degree twice.
- **G-R3b** *(registered iteration 7)*. A feasibility probe had already
  reported the **degrees** at `N = 16, 18, 20`, so a gate on those would be
  worthless and they are recorded as an *observation* instead (→ R3c). The
  gate is on the **net saving in bits**, which nobody had computed when it
  was written: *supported* if positive at all of `N = 16, 18, 20` **and**
  the mean over `{16,18,20}` is not more than `0.25` bits below the mean
  over `{10,12,14}` (the materiality threshold from G-R5′); *killed* if the
  saving is `≤ 0` at any of the three; *flat* if positive but the mean has
  fallen by more than `0.25` bits. **A verdict that differs across seeds is
  reported as `unstable`, not resolved by majority** — picking the modal
  answer is choosing it after seeing it.
- **The boundary column is not a new gate.** §3's baseline 2 (`2^N`
  enumeration) was fixed in iteration 1. What EXP-R3b adds is only that it
  is now *reported per cell*, next to the saving. Every lever verdict from
  iteration 7 on must carry it, because G-R3b demonstrated that a gate can
  pass while the ratio to the boundary is going backwards (→ R3d).
- **G-R4d** *(registered iteration 8, before any measurement code was
  written; note the ledger already uses `R4c` for a different row, so the
  experiment queued as "EXP-R4c" is named **EXP-R4d**)*. R4 asked whether
  levers act *through* `Δ_low` and was killed as a pooled instance-level
  statistic; iteration 6 showed the defect works at the **family** level
  for factor-base families. This asks whether it also works for
  **lever-generated presentations**, using three presentations of the
  *same system at the same variable count* — `raw`, one saturation round
  (`one-round`), and fully saturated (`saturated`) — so that nothing in
  the comparison is cross-shape. Two statistics, both size-controlled:
  1. **Paired lever effect.** Per target, `dΔ = Δ_low(sat) − Δ_low(raw)`
     against `dD* = D*(sat) − D*(raw)`, blocked on `(family, N)`. If the
     lever acts through the defect, the targets whose defect rose most are
     the targets whose `D*` fell most — a *negative* correlation.
  2. **Presentation-family level.** Aggregate to
     `(operating point, presentation)` means and run the same three
     size-controlled statistics, exactly as EXP-R4b's panel B did for
     factor-base families.
  *Supported* if **both** reach `≤ −0.6`; *killed* if neither does;
  *partial* if exactly one does.
  **Degeneracy clause, stated in advance because it is the likely
  outcome:** EXP-R3 and EXP-R3b both show the augmented system pinned at
  `D* = 2` over wide ranges, so `dD*` may have *no variance* at some
  cells. Where the response is constant no correlation is defined, and the
  cell must be reported as **degenerate** — not scored, and not silently
  averaged in as 0. A gate that passes only because the degenerate cells
  were dropped is not a pass, so the count of degenerate cells is reported
  next to the verdict.
- **G-R4.** *Supported* at pooled `ρ_s ≤ −0.6` over ≥ 30 lever-generated
  cells. *Killed* at `|ρ_s| < 0.2` or a sign flip. **(Retired, iteration 4:
  a pooled `ρ_s` over mixed sizes is not evidence here — see G-R4″.)**
- **G-R4′** *(registered iteration 4)*. A defect variant is *supported* if
  it reaches `ρ_s ≤ −0.6` on all three groups — `within-shape` (vars held),
  `rho-varies` (slicing), `rho-matched` (sizes vary at `ρ ≈ 1`). *Killed* if
  none does.
- **G-R4″** *(registered iteration 4)*. The **size-controlled** statistic —
  the mean of the per-`vars` Spearman — is what decides whether a pooled
  correlation is structure. A variant's pooled figure counts as real only if
  its size-controlled figure also clears −0.6; otherwise the pooled figure is
  recorded as a size proxy. **Amended iteration 6**: the gate must also name
  its **unit of analysis**, because a size-controlled statistic on targets
  and one on family means answer different questions and can disagree by
  0.7 (R4c). A verdict that does not say which unit it used is not a verdict.
- **G-R4b** *(registered iteration 6)*. Applied to a *published* correlation,
  with three size-controlled statistics rather than one, because any single
  one has a failure mode the others do not: mean per-block (equal weight per
  stratum, blind to strata too small to correlate), blocked rank (all cells
  at once, needs strata rank-comparable), fixed-effects (raw magnitudes, not
  ranks). The law is **vindicated** if all three reach ≤ −0.6; recorded as a
  **size proxy** if pooled clears −0.6 while they do not; **inconclusive** if
  pooled does not itself clear the bar on re-analysis, since then there is
  nothing to control for.
- **G-R5.** *Supported* if `c(composed) < c(one-side)` at ≥ 2 operating
  points, seed-robust. **(Retired as degenerate, iteration 3: the composed
  route reaches the floor at `k = 0`, so the gate passes by construction
  and measures nothing.)**
- **G-R5′** *(registered iteration 3, replacing the degenerate G-R5)*.
  Minimise total work over `k` for each route —
  `raw(k) = 2^k·macaulay(N−k, n, D*_raw)` versus
  `composed(k) = 2^k·[extract(N−k) + macaulay(N−k, aug, D*_mut)]` —
  and compare the best of each against the other and against `2^N`.
  *Supported* if composed beats both everywhere; *killed* if it loses to
  the raw route with a flat gap; *blocked* only if the gap closes by
  ≥ 0.25 bits per size step (a materiality threshold, added because the
  first version of the test called a 0.01-bit wobble "closing").

- **G-R6** *(registered iteration 9)*. L5 is *supported* if some curve keeps
  a mean-`D*` margin `≥ 0.5` degrees over the unmodified curve on a
  **disjoint holdout** target set, at every measured `n`, **and** the
  between-curve variance of mean `D*` exceeds twice the no-effect prediction
  (pooled variance over the effective per-curve sample size). *Killed* if the
  holdout margin is `≤ 0` while the selection margin is positive — a winner's
  curse. *Blocked* if the margin is positive but shrinking in the target
  count. The holdout half is not optional: `D*` takes 3 values here, so
  selecting the minimum over a dozen targets manufactures a margin from
  nothing.

---

## 5. Iteration log

> Newest at top. Format mirrors `RESEARCH_FFD_WORKFLOW.md` §7:
> *Task · Experiment · Result · Gate verdict · Ledger delta · Next.*

### 2026-09-12 — iteration 9 (EXP-R6 — L5, the curve-side lever: empty, and derivably so)

**Task.** Score the one lever the taxonomy was missing: change the *curve*.
The ECDLP transports along an isogeny of degree coprime to the subgroup
order, so an attacker may solve on any curve in ECC2K-130's isogeny class —
the move Galbraith–Hess–Smart turned into an attack on Weil descent.

**Experiment.** `RESEARCH_ISOGENY_CLASS_SEARCH.md`, module
`cryptanalysis::isogeny_class_search`, driver
`examples/isogeny_class_search.rs`, snapshot
`experiments/isogeny_class_search.json`.

**Result.** Three of the four boundaries are derived, and they leave the
lever no room.

- **The curve enters below the leading form.** `S₃` depends on the curve only
  through `a₆ = 1/j`, and only as an *additive constant*, so after Weil
  descent `a₆` moves nothing but the constant term of each equation —
  verified coefficient-by-coefficient at `n ∈ {8,10,12}`. At `m = 3` the
  `a₆`-expansion of `S₄` puts `a₆` in Boolean degree `≤ 5` of `6` (identity
  verified on 300 combinations). So `d_reg` — a Hilbert invariant of the
  leading forms — is **constant on the whole isogeny class**, and a curve
  change can only *lose* a degree fall, never create one lower.
- **The class is bigger than ρ.** `Σ_{f|c} h(O_f) = 2^65.06` vertices
  (`c = 263 · 146505763881528721`, `h(−7) = 1`) against ρ's `2^60.81`
  (`√(πr/2m)` at `m = 262`; `2^64.83` for plain ρ). *Corrected in review: an
  earlier revision divided the `m = 2` form by `√262`, double-counting
  negation and understating ρ by half a bit.*
- **263 vertices are reachable.** `h(−7) = 1` makes every crater a single
  vertex, so horizontal isogenies return to `E` and descending ones need
  `ℓ | c`: every prime below 263 has a *single-vertex* graph, cross-checked
  against `Φ₂(X,1)` and `Φ₃(X,1)` over the real `F_{2^131}`.
- **The affine `D*` does vary** — `{2,3,4}` over all 255 curves at `n = 8` —
  and the mechanism is exact: `D* = 2 ⟺ a₆ ∉ S^⊥` where
  `S = {λ : Σ λ_i h_i = 0}` is the left null space of the *curve-independent*
  leading parts. `S^⊥` is fixed by the **target**, so goodness is a
  `(curve, target)` property. Zero disagreements with the solver.
- **Holdout.** selection margin `+0.273/+0.375/+0.600` at `n = 6/8/10`;
  holdout margin `+0.667/+0.403/−0.788`; variance ratio `1.345/1.382/1.409`
  against 1 (no effect) and 12 (deterministic).

**Gate verdict.** **G-R6 killed** — the holdout margin is negative at the
largest size and no size clears `≥ 0.5` at a variance ratio `≥ 2`.

**Ledger delta.** R6 `killed`, R6′ `killed`, R6″ `killed`.

**Class of the change.** **Accounting.** The numbers that moved are the
boundaries; the attacker's algorithm is unchanged. The selection-set row of
the table is a textbook **relabelling** — exactly on the floor where it was
chosen, `2^{+1.13}` above baseline where it was not.

**Next.** The one thing that would reopen L5 is a counterexample to Boundary
C at `m ≥ 4`: if `a₆` reaches the leading form of `S₅`, `d_reg` becomes
curve-dependent at the first `m` where index calculus is asymptotically
interesting. That is a symbolic computation, not a search.
### 2026-09-12 — iteration 8 (EXP-R4d — the defect changes sign with the grouping)

- **Task picked.** The last cheap item: R4 asked whether levers act
  *through* `Δ_low`, and iteration 4 killed it as a *statistic* rather than
  as a claim. Iteration 6 then found the defect strong at the family level.
  So the question was open, not closed. *(Named EXP-R4d, not the queued
  "EXP-R4c" — the ledger already uses R4c for a different row.)*
- **Design that removes the objection which killed R4.** Three
  presentations of the **same system, same variable count, same target**:
  raw, one saturation round, fully saturated. They differ only in how many
  degree-3 falls have been folded in — a *dose* axis for L4 with nothing
  cross-shape in it, so R4′'s objection cannot apply. Gate **G-R4d** was
  registered in the charter before any measurement code was written,
  including a degeneracy clause anticipating that `dD*` might be constant.
- **The gate is killed, and the sign is inverted.** Paired blocked rank
  **+0.8385**; presentation-family blocked rank **+0.8704**. The gate asked
  for `≤ −0.6`. This is a harder kill than "no relation": the relation runs
  *the other way*.
- **My own design flaw, which the degeneracy clause did not catch.** The
  paired statistic turned out to be *mathematically identical* to the raw
  within-cell correlation. Saturation sends `Δ_low → 0` and `D* → 2`, so
  `dΔ = −Δ_raw` and `dD* = 2 − D*_raw`, and the two sign flips cancel. The
  clause I registered covered the case where `dD*` is constant (it caught
  1 of 12 cells); it did not cover the subtler case where the *difference*
  is an affine image of the raw measurement and therefore carries no new
  information. Statistic 1 was not a second reading, it was the first one
  wearing a disguise.
- **What the data actually says (statistic 3, added after seeing the
  sign).** Stratify each cell's targets by their `D*` and the cause is
  plain:

  | family | `N` | `D* = 2` targets | higher-`D*` targets | within `ρ_s` |
  |---|---:|---:|---:|---:|
  | Coordinate | 14 | `Δ = 0.00213` (10) | `Δ = 0.03511` (6, `D*=4`) | **+0.871** |
  | Coordinate | 16 | `Δ = 0.00215` (8) | `Δ = 0.02672` (8, `D*=4`) | **+0.909** |
  | Random | 14 | `Δ = 0.00000` (1) | `Δ = 0.00213` (15, `D*=4`) | **+1.000** |
  | Subfield | 14 | `Δ = 0.07234` (12) | `Δ = 0.07872` (4, `D*=3`) | **+1.000** |

  **Positive in 10 of 10 decidable cells.** Meanwhile the between-family
  law at the same sizes is `−0.85` (R4c). Same statistic, same data,
  opposite signs at the two groupings — Simpson's paradox, with a
  mechanism rather than a shrug.
- **The mechanism.** `Δ_low` sums degrees `≤ 3`, so a system can only
  exhibit a degree-3 defect if its Macaulay tower **reaches** degree 3 in a
  deficient state. A target that refutes at `D* = 2` never exercises degree
  3, so its cutoff-3 defect is ~0 *by construction*. Inside a homogeneous
  group the statistic is therefore partly a proxy for "did this instance
  need degree 3" — which is `D*` itself, positively. Between families the
  structural differences dominate and the sign flips back. **Both are
  real.** The defect does not carry a single sign, and which one you get
  depends on what you hold fixed.
- **And a second, purely mechanical limit (R4f).** **0 of 192** saturated
  targets retain any cutoff-3 defect: saturation drives `Δ_low` to exactly
  zero, every time. It must — the cutoff-3 defect *is* the space of
  degree-3 falls, and saturation adds precisely those as generators. So
  `Δ_low` cannot score a lever's **output** at all, correlation or no
  correlation: the lever's action is to consume the quantity. Pinned in
  `saturation_consumes_the_early_defect` rather than left as prose.
- **What this does and does not touch.** It does **not** touch the FFD
  program's law: that is a between-family statement, it is what the screen
  is used for, and iteration 6 confirmed it size-controlled at `−0.90` in
  the critical regime. What it adds is the **scope boundary** — the screen
  ranks constructions and must never rank instances inside one, not merely
  because it is weak there (R4′) but because it points the wrong way.
  Choosing the easiest target in a family by lowest `Δ_low` would
  systematically choose the hardest.
- **Gate verdicts.** G-R4d: **killed** (sign inverted). R4 is now closed
  rather than "killed as stated": the defect does not mediate L4 on any
  reading tested, and R4f says it cannot in principle score a lever output.
- **Ledger delta.** R4d registered→killed; R4e registered→supported
  (discovered, not predicted — labelled as such); R4f registered→supported
  (mechanical, test-pinned). R4′ sharpened from "weak" to "wrong-signed".
- **Class (AGENTS.md §3): `accounting`.** Numbers changed, no algorithm
  did, and no ratio to the `2^N` boundary moved — this iteration measures a
  predictor, not an attack. The thread is where iteration 7 left it.

### 2026-09-12 — iteration 7 (EXP-R3b — the gate passes and the boundary says no)

- **Task picked.** The queue head: reach for the L4 trend. EXP-R3 measured
  the Random-family net saving *growing* `+0.94 → +1.75 → +1.62` over
  `N = 10,12,14`, and R3′ left it `regime-dependent`. Three sizes is a short
  lever arm for a trend claim.
- **Reach, and a charter correction.** A feasibility probe found the cells
  complete well past where §3 said they would: `N = 16, 18, 20` at 14 s /
  45 s / 167 s per 8-target cell. §3's "`D*` is measurable to `N ≤ 16`" is
  now "`N ≤ 20`". Because the probe reported degrees, **G-R3b was registered
  on the net saving in bits instead** — the one number nobody had computed.
- **The gate passes.** Random, all three seeds: mean saving `1.44 → 4.44`,
  `1.45 → 4.44`, `1.49 → 4.16` bits, positive at every reach point. `N = 18`
  and `N = 20` come out identical to 2 d.p. on every seed. **G-R3b:
  supported.** Coordinate's verdict flips with the seed
  (KILLED/SUPPORTED/KILLED on `N = 18` values `−0.18/+0.89/−0.12`), so it is
  recorded `unstable` rather than resolved by majority. Subfield stays
  killed.
- **And the boundary says the opposite.** §3's baseline 2 has been fixed
  since iteration 1; EXP-R3b is simply the first experiment to *print it*:

  | `N` | 10 | 12 | 14 | 16 | 18 | 20 |
  |---|---:|---:|---:|---:|---:|---:|
  | degree saving | 1.75 | 2.00 | 2.00 | 2.25 | **1.00** | **1.00** |
  | net saving vs raw solve | +0.94 | +1.75 | +1.62 | +5.11 | +3.87 | +4.34 |
  | **margin vs `2^N`** | −10.20 | −9.96 | −10.12 | −9.51 | **−12.69** | **−12.09** |

  **0 of 54 cells beat `2^N`** (3 seeds × 3 families × 6 sizes). The route
  is 9.5–12.7 bits *behind* brute force, and the gap **loses 3.18 bits
  across the `N = 16 → 18` break** — exactly where the gate starts
  reporting its best numbers — for `−0.378` bits per size step overall.

  Two scope limits on that, because "the gap widens" is the kind of
  sentence that travels further than its evidence. It is **Random-only**:
  Subfield's margin narrows monotonically and Coordinate's is non-monotone.
  What kills R3d is the flat fact that *nothing beats the boundary
  anywhere*, not a universal trend.
- **The one thing trending toward the boundary is not ours (→ R3e).** On
  the Subfield family the **raw, unaugmented** solve closes from `−4.99` to
  `−1.19` bits over `N = 10…20`, and would cross near `N ≈ 22–24`. That is
  arithmetic rather than cryptanalysis: Subfield's `D*` sits at
  `2.12–2.25`, at the Nullstellensatz floor, so its solve is `poly(N)`
  against an exponential baseline and *must* cross eventually. It measures
  the family's degeneracy — the FFD program's L1 result — not a lever of
  this thread. It is written down because "Gröbner overtakes brute force at
  `N ≈ 23`" is exactly the sentence someone would quote out of this table,
  and it would be wrong. And L4 **ruins** that family: augmented margin
  `−10.49` against raw `−4.99` at `N = 10`. On the only family approaching
  the boundary, this thread's lever moves it away.
- **Why both are true.** The saving is measured against the raw Gröbner
  solve, and the raw solve got *harder*: base `D*` goes `4 → 5` between
  `N = 16` and `N = 18`. So the headline rose because the reference moved.
  The lever itself got **worse** over the same step — its degree saving
  halved, `2.25 → 1.00` — because the mutant cascade dies: at `N ≤ 16`
  saturation runs several productive rounds (`sat.eqs = 137` at `N = 16`
  from 32 falls) and reaches the `D* = 2` floor; at `N ≥ 18` it stops after
  one round (`sat.eqs = 54 = 18 + 36`) and lands at `D* = 4`. The cost
  structure flips with it: extraction-dominated at `N ≤ 16`
  (`total 29.51 ≈ extract 29.50`), solve-dominated at `N ≥ 18`
  (`total 34.86 ≈ solve 34.83`).
- **Class (AGENTS.md §3): `relabelling`.** The number this thread had been
  tracking improved 3×, while the ratio to the stated boundary went
  backwards by 3.18 bits at the break. AGENTS.md calls this class "not
  hypothetical", and it is not: the gain is attributable to the baseline
  degrading, not to the method improving. **On the family where L4 pays,
  the ratio to the boundary is not flat — it is getting worse with `N`; on
  no family does anything cross it.** That is the result.
- **Gate verdicts.** G-R3b (registered this iteration): **supported** on
  Random, **unstable** on Coordinate, **killed** on Subfield — and
  *superseded in importance* by the boundary column on the same data.
- **Ledger delta.** R3b registered→supported/unstable/killed by family;
  R3c registered→killed (the cascade does not persist); R3d
  registered→killed (nothing beats the boundary at any size, on any
  family, on any seed); R3e registered→observed-and-explained. §3's reach
  and cost-model description corrected, §7's reach claim with it.
- **The methodological point, which is now six of seven.** Iterations 1, 3
  and 5 each caught a metric that counted one resource while ignoring
  another; iteration 4 mismatched units and thresholds; iteration 6 caught
  a statistic that fixed the nuisance variable it had thought of and not
  the one that mattered. This one is the same failure in its purest form,
  and
  **the gate was mine, registered in advance, and still wrong** — because
  pre-registration fixes *when* you choose the metric, not *whether the
  metric is the right one*. The only defence that worked was the boundary,
  and it worked because it was fixed in iteration 1 and is not allowed to
  move. That is what §3 is for, and from now on every lever verdict carries
  the column.
- **A bug in iteration 6's code, found from outside the thread.** While
  this iteration ran, Cursor Agent pushed a fix to `size_control`
  (`e16b67c`): raw 1-based within-stratum ranks have mean `(n_g+1)/2`, so
  pooling them across **unequal** strata reintroduces a between-block term
  `Σ n_g (μ_g−μ)²` — a function of block size, which is the exact nuisance
  `blocked_rank` exists to remove. It is right. Verified and merged, with
  its test.

  What it changes, checked cell by cell rather than assumed: **exactly one
  number.** Panel A's strata are all 5 cells and the family-level panel's
  are all 3, and with equal strata the injected term is zero — so
  `−0.6983`, the critical-regime `−0.9023`, and the family-level `−0.8528`
  are untouched. The instance-level panel is the one with unequal strata
  (144…432), and its blocked rank moves **`−0.1189 → −0.3301`**. The
  instance/family gap is therefore `0.52`, not the `0.73` first published;
  R4c's direction and conclusion are unaffected.

  The part worth keeping is *how it was missable*. My three instance-level
  statistics read `−0.3494`, `−0.1189`, `−0.3248` — one of them a clear
  outlier against the other two, in a function whose whole purpose is that
  the three should agree. I wrote three statistics precisely so that
  disagreement would be informative, then did not read the disagreement.
  The tests did not catch it either, because every block I wrote in them
  was equal-sized. **A test suite that only exercises the balanced case
  cannot see a bias that is defined as a function of imbalance.**
- **Next.** The queue is thin and should be honest about it. L4 is now
  measured to `N = 20` and loses ground; L3 and L3∘L4 are killed; L2 is
  bounded below `ℓ = 6` and out of reach. **No lever in the taxonomy has a
  path to the boundary at the sizes this instrument can see.** The two
  remaining items are EXP-R4c (family-level lever scoring, cheap, cells
  already exist) and EXP-R2b (sparse F4/F5 at ~51 variables — an
  engineering project that should not start without deciding L2 is worth
  it). Neither is likely to change the boundary verdict.

### 2026-09-12 — iteration 6 (EXP-R4b — the thread withdraws its own flag)

- **Task picked.** The queue head: the upstream caveat iteration 4 attached
  to the FFD program's headline `ρ_s = −0.79`. It was the only open item
  where a cheap experiment could overturn something the repo treats as
  settled — and, because this thread raised it, the only one where *not*
  running the experiment leaves a wrong claim standing in someone else's
  document on this thread's authority.
- **Built** (`src/cryptanalysis/degree_reduction.rs`, +5 tests):
  `size_control`, which reports a pooled correlation next to **three**
  size-controlled counterparts — mean per-block Spearman, blocked rank, and
  fixed-effects. Three rather than one because each has a failure mode the
  others do not (small strata, rank-incomparability, magnitude-blindness);
  agreement between them is the evidence, and disagreement localises which
  part of the pooled figure was the nuisance variable. The tests pin both
  directions: a construction with *exactly zero* within-stratum correlation
  and a strong stratum trend must read pooled `+0.895` and controlled `0` on
  all three; a relation identical inside every stratum must survive all
  three at `< −0.99`. Degenerate strata report `None` and are skipped, never
  averaged in as `0`.
- **Result, panel A — EXP-G's own 50 cells.** The law is **not** a size
  proxy:

  | statistic | all cells | critical (`2n' = n`) | over-determined |
  |---|---:|---:|---:|
  | pooled `ρ_s` | −0.7929 | −0.7781 | −0.7334 |
  | mean per-block `ρ_s` | **−0.6998** | **−0.9124** | −0.3809 |
  | blocked rank `r` | **−0.6983** | **−0.9023** | −0.3979 |
  | fixed-effects `r` | **−0.7499** | **−0.8186** | −0.5752 |

  The pooled column reproduces the published −0.79324/−0.77805. In the
  ECDLP-relevant critical regime, holding size fixed makes the law
  **stronger** than pooling it. The over-determined column is weak exactly
  as **P6** predicts — `ρ ≫ 1` floors `D*` at 2, leaving no variance for any
  predictor to explain; the three weakest single blocks are `n'=2` at every
  `n`, so the weakness tracks the determination ratio, not the size.
- **Result, panel B — iteration 4's own 2016 cells, two units of analysis.**
  This is the part that settles *why*, without appeal to EXP-G's data:

  | statistic | instance level (per target) | family level (per `(vars, family)` mean) |
  |---|---:|---:|
  | pooled `ρ_s` | −0.5164 | −0.8734 |
  | mean per-block `ρ_s` | −0.3494 | −1.0000 |
  | blocked rank `r` | **−0.3301**† | **−0.8528** |
  | fixed-effects `r` | −0.3248 | −0.8337 |

  † *Corrected in iteration 7 from the −0.1189 first published here; see
  that entry. The family-level column and all of panel A are unchanged —
  their strata are equal-sized, where the bug has no effect.*

  Same snapshot, same defect, same `D*`, same size control. Only the unit of
  analysis differs, and it is worth 0.52 in blocked rank. The three families
  are ordered correctly in **4 of the 4** size blocks where `D*` varies at
  all; the other two (`vars = 4, 6`) are floored at `D* = 2`.
- **Diagnosis.** EXP-G's within-block contrast is across three structurally
  different factor-base **families** whose defects span ~27× at fixed size.
  Iteration 4 correlated individual **targets** pooled across families,
  where target-to-target noise swamps the family signal — and then reasoned
  from the *shape* of EXP-G's aggregation ("pooled the same way") to the
  conclusion that it shared the artifact. The aggregations are not the same,
  and the analogy was false.
- **What survives and what does not.** R4′'s kill **stands**: no
  renormalisation makes `Δ_low` compare individual instances across shapes,
  and five variants × four seeds is good evidence for that. What does not
  survive is the inference from it — that the published law is therefore a
  size proxy. Those are different claims about different units, and only the
  first was ever tested. The `Δ_low` screen is a **family-level
  discriminator**, which is the use the FFD proposal put it to: ranking
  curve and factor-base choices, not ranking targets within one.
- **Gate verdicts.** G-R4b (registered this iteration): the law is
  **vindicated**. G-R4″: **amended** — a size-control verdict must now name
  its unit of analysis, since two honest ones can disagree by 0.7.
- **Ledger delta.** R4b registered→**killed** (the prediction was that the
  upstream law *is* a proxy). R4c registered→**supported**. R4″'s scope
  narrowed to the instance level, with the original wording and the reason
  it misled recorded in the row. R4's kill re-justified: it kills the
  *statistic*, not the defect. §2.1 corrected; iteration 4's upstream caveat
  marked withdrawn **in place** rather than deleted — a ledger that removes
  its wrong calls stops being a record.
- **Class (AGENTS.md §3): `accounting`.** Numbers changed, algorithm did
  not, and no ratio to the `2^N` boundary moved — this iteration produced
  no attack progress and must not be read as any. What it produced is a
  correction to a claim, which AGENTS.md says is worth committing and not
  worth calling a result. Stated plainly: **the thread is still at zero
  against its boundary**, exactly where iteration 5 left it.
- **Method note, against this thread's own recurring error.** Iterations 1,
  3, 4 and 5 each caught a metric that counted one resource and ignored
  another. This one is the same class with a different resource: a statistic
  that fixed the nuisance variable it had thought of (size) while leaving
  the unit of analysis unstated. The fix is the same shape as the others —
  name the thing being held constant, and report the comparison that would
  falsify you next to the one that supports you.
- **Next.** EXP-R3b (reach for the L4 trend at `N = 16–18`) is now the queue
  head. The newly-open question from R4 is worth a line in the queue too:
  the defect is dead as an instance-level lever score, but **family-level
  lever scoring was never tested** and is not excluded by anything measured
  so far.

### 2026-09-12 — iteration 5 (EXP-R2 — L2 is bounded, and out of this instrument's reach)

- **Task picked.** R2, the last untested lever and the thread's remaining
  upside. Unblocked by iteration 4, which settled that `Δ_low` must not be
  used to score it.
- **Why `S₄`.** Symmetrisation's saving scales like `m!`. The rest of this
  thread runs on `S₃` (`m = 2`), where `m! = 2` — testing L2 there would
  measure noise. `binary_semaev_s4.rs` already carries a symmetrised `S₄`
  descent at `m = 3`, the first `m` where index calculus beats Pollard ρ at
  all. So L2 was testable, on a harness the thread had not used.
- **Built** (`src/cryptanalysis/degree_reduction_anf.rs`, 6 tests): an
  **arbitrary-degree ANF Macaulay path**, since the thread's harness is
  quadratic-only and both halves of the `S₄` model exceed that. It reuses the
  same bit-packed elimination and the same "is `1` in the row space" test, so
  the degrees are commensurable with iterations 1–4. Plus the **eliminated
  presentation** — substitute `eᵢ = σᵢ(x)` — which gives the *same ideal*
  over `3ℓ` variables at degree 6, so the baseline is exact rather than an
  independently re-derived system. A test checks the two presentations agree
  on every one of the `2^{3ℓ}` assignments.
- **Result — G-R2 BLOCKED, structurally.** The two presentations are never
  simultaneously measurable:

  | `ℓ` | symmetrised | eliminated |
  |---|---|---|
  | 2 | `D* = 4` measured | **degenerate** — degree-6 generators in 6 variables leave *zero* multiplier budget, so the tower cannot refute whatever the satisfiability |
  | 3 | **censored** at `D ≥ 5` (24 variables) | `D* = 8–9` measured |

  Two honesty guards were added for exactly these failure modes
  (`elim_degenerate`, `censored_at`), because a degenerate or censored scan
  reports "no refutation" and must never be read as one.
- **Result — R2′ SUPPORTED, and it is the iteration's real deliverable.**
  Whether symmetrisation can pay is an arithmetic question with an exact
  answer, no solver required. Macaulay width at each presentation's
  generator degree:

  | `ℓ` | `cols(9ℓ−3, 3)` | `cols(3ℓ, 6)` | narrower |
  |---:|---:|---:|---|
  | 2 | 576 | 64 | eliminated |
  | 5 | 12 384 | 9 949 | eliminated |
  | **6** | **22 152** | **31 180** | **symmetrised** |
  | 10 | 109 824 | 768 212 | symmetrised (7×) |

  **Symmetrisation is the narrower presentation only from `ℓ = 6` upward.**
  Below that it trades 3× the variables for ½ the degree and loses. Since a
  dense Macaulay scan reaches `ℓ ≤ 3`, **EXP-R2 as designed could not have
  answered the question at any budget** — the block is not a compute
  shortfall, it is that the interesting regime starts past the instrument.
- **What it would take:** a solver handling 51 variables at degree 3 — real
  F4/F5 with sparse linear algebra, not a dense scan. That is a different
  engineering project, and naming it precisely is more useful than a
  censored number.
- **Gate verdicts.** G-R2: **blocked** (< 3 comparable cells, structurally).
  G-R2′ (registered this iteration): **supported**.
- **Ledger delta.** R2 open→blocked; R2′ registered→supported.

### 2026-09-12 — iteration 4 (EXP-R4′ — the defect screen is a size proxy)

- **Task picked.** R4′, promoted to queue head by iteration 3: L2 is the
  only untested lever and it changes the variable count, so `Δ_low` cannot
  score it until we know whether a shape-corrected defect exists.
- **The hypothesis going in was wrong, which is why it was worth running.**
  I expected over-determination to flip the law's sign (both `D*` and the
  defect fall as `ρ` rises), making `rho-varies` the failing group. The data
  says the opposite: `rho-varies` is where the normalised variants score
  *best*.
- **Experiment** (`degree_reduction::{DefectVariant, collect_defect_cells,
  collect_rho_matched_cells}`, `examples/degree_reduction_defect.rs`,
  snapshot `experiments/degree_reduction_defect.json`). Five candidate
  defect summaries — normalised (the incumbent `Δ_low`), raw, generic-rank
  fraction, per-equation, row fraction — scored against `D*` on three
  designed groups: `within-shape` (vars held, family/target vary),
  `rho-varies` (slicing, both shape and `ρ` move), `rho-matched` (sizes vary
  along `2n'=n`, `ρ ≈ 1` held).
- **Result — R4′ KILLED, and the reason is not units.**

  | variant | pooled (`rho-matched`) | size-controlled |
  |---|---:|---:|
  | normalised (`Δ_low`) | **−0.891** | −0.292 |
  | raw | −0.884 | −0.292 |
  | generic-frac | −0.905 | −0.164 |
  | per-equation | −0.905 | −0.292 |
  | row-frac | −0.905 | −0.164 |

  Every variant that clears the bar *pooled* fails it once `vars` is held
  fixed. `D*` and a size-normalised defect both trend with `N`, and pooling
  across sizes reads that shared trend as correlation. **The pooled figure
  is a size proxy at the instance level** (→ R4″). No denominator repairs
  that, because the problem is not the denominator.

  **Correction, iteration 6.** The qualifier "at the instance level" was not
  in the original text, and its absence was not cosmetic — it is what let
  iteration 4 carry the conclusion upstream to a figure it had not measured.
  Aggregate these same cells to `(vars, family)` means and the correlation
  returns at **−0.85** (blocked rank, size-controlled). The defect ranks
  *constructions*, not *targets*. See R4b/R4c and EXP-R4b.
- **Controls that show the test has power.** The un-normalised variants
  (raw, per-equation) score −0.10…−0.21 on `rho-varies` while the
  normalised ones score −0.85…−0.94 — exactly the expected failure of an
  unscaled quantity when the variable count moves. And within a fixed shape
  all five variants rank instances *identically* (they are monotone
  transforms of `Σδ` there), which a test asserts; that is why no
  renormalisation can improve the within-shape column.
- **Method corrections made during the iteration**, both of which changed
  the numbers:
  1. The first grouping pooled `vars = 12` and `vars = 14` into
     "within-shape" — already a cross-shape comparison wearing the wrong
     label. Now one correlation per fixed `vars`.
  2. The −0.6 threshold is calibrated against EXP-G's `ρ_s = −0.79`, which
     was measured on **cell means**. Correlating raw instances against a
     threshold fitted to aggregates compares unlike things; the pooled
     groups are now aggregated to `(family, vars)` means first.
  Statistics were also raised from 6 to 48 targets per cell after the first
  run proved seed-unstable (`ρ_s` swinging −0.24…−0.90 across seeds). At 48
  the figures are stable over seeds 7/11/23/41.
- **A caveat that points upstream, stated as a caveat.** The FFD program's
  own `ρ_s = −0.79` was also pooled across `2n' ∈ {4,…,14}`, so it may share
  this property. This thread has not re-run EXP-G and what it measured is
  related but not identical, so this is a flag to check, not a refutation.
  If it does hold there, the `Δ_low` screen — the FFD proposal's headline
  deliverable — is weaker than advertised, and that matters to the defensive
  side of the program as much as this one.

  > **Withdrawn, iteration 6 (EXP-R4b).** It was checked, and it does not
  > hold there. EXP-G's law survives every size control (−0.70/−0.70/−0.75
  > against pooled −0.79), and in the critical regime size control makes it
  > *stronger* (−0.91). The flag was raised on an analogy — "pooled the same
  > way, so possibly the same artifact" — and the analogy was false, because
  > EXP-G pools family cell means where iteration 4 pooled individual
  > targets. Flagging it was right; what was wrong was reasoning from the
  > shape of someone else's aggregation without re-running it. The caveat is
  > left standing above rather than deleted, because a ledger that quietly
  > removes its wrong calls stops being a record. See R4b, R4c, G-R4b.
- **Gate verdicts.** G-R4′: **killed**. G-R4″ (registered this iteration):
  **supported**. G-R4: **retired** — a pooled `ρ_s` over mixed sizes cannot
  establish it.
- **Ledger delta.** R4′ open→killed; R4 open→killed; R4″ registered→
  supported. §2.1's iteration-1 diagnosis corrected. *(R4″'s scope narrowed
  to the instance level in iteration 6; R4's kill re-justified there.)*
- **Next.** R4′'s kill *unblocks* R2 rather than blocking it: the
  prescription is simply not to use `Δ_low` for scoring. L2 is measured the
  way iterations 2 and 3 measured their levers — matched `(vars, eqs, ρ)`,
  scored on `D*` and total work. `binary_semaev_s4.rs` already carries a
  symmetrised `S₄` descent with its correspondence system (`m = 3`, where
  the `m!` saving is actually meaningful, unlike the `m = 2` harness the
  rest of this thread runs on), so the build is smaller than it looked.

### 2026-09-12 — iteration 3 (EXP-R5 — the levers are substitutes, not complements)

- **Task picked.** R5, the composition test promoted to the head of the
  queue by iteration 2: one-sided guessing reaches the floor at `c = 1/2`
  (R1′) and mutants lower the working degree on the hard family (R3), so
  does composing them reach the floor more cheaply than either alone?
- **The pre-registered gate turned out degenerate — reported as a gate
  failure, not a result.** G-R5 asked for `c(composed) < 1/2`. The composed
  route gives `c = 0.000` in every cell: the mutant system is already at
  `D* = 2` with **no guessing at all**, which is iteration 2's finding
  restated. The gate passes by construction. This is the *third* instance of
  one error in this thread — a metric that counts one resource while
  ignoring what the other costs (iteration 1: boundary optimum; iteration 2:
  uncharged extraction). Replaced by **G-R5′**, scored on total work.
- **Experiment** (`degree_reduction::run_composed_sweep`,
  `examples/degree_reduction_composed.rs`, snapshot
  `experiments/degree_reduction_composed.json`). For each non-decomposable
  target and each `k`, take slices under one-sided guessing and measure each
  slice **twice** — raw, and saturated with its own degree-3 falls — so the
  two routes are compared on identical systems. Then minimise
  `2^k·macaulay(...)` (raw) and `2^k·[extract + macaulay(...)]` (composed)
  over `k`.
- **Result — R5 KILLED.** The composed route is worse at its own optimum,
  by a gap that does not move:

  | `N` | log₂ 2^N | best raw | best composed | composed − raw |
  |---:|---:|---:|---:|---:|
  | 10 | 13.32 | 14.94 | 17.82 | **−2.89** |
  | 12 | 15.58 | 17.20 | 20.07 | **−2.87** |
  | 14 | 17.81 | 19.42 | 22.29 | **−2.87** |

  Identical on Random and Coordinate (at the optimum both families sit at
  the `D* = 2` floor, so the cost is shape-determined), and stable across
  seeds 7/11/23. Neither route beats `2^N` enumeration at these sizes, and
  **both optima sit at the largest `k` scanned** — the iteration-1
  degeneracy again, now flagged automatically by
  `ComposedSweep::optimum_is_interior`.
- **Result — R5′ SUPPORTED, and it explains the kill.** At `k = 0` the
  mutants are worth `+1.1…+1.75` bits on Random (iteration 2, reproduced).
  The instant guessing is allowed the sign flips and stays flipped at every
  `k > 0`. The mechanism is not subtle: **both levers drive the system to
  the same `D* = 2` floor**, and once a slice is already there, the mutant
  route's degree-3 extraction is pure overhead — exactly the Subfield effect
  of iteration 2, reappearing because guessing *manufactures* easy slices.
  So the two levers are **substitutes, not complements**; there is no
  composition gain to find, and the flatness of the gap in `N` says that is
  structural rather than a small-size artifact.
- **Method note.** The first version of G-R5′ called the `−2.89 → −2.87`
  drift "closing". It is not — that is a flat line. The gate now requires
  the gap to shrink by ≥ 0.25 bits per size step before it may be called
  closing; the threshold was added after seeing the wobble, and is recorded
  here rather than quietly applied.
- **Gate verdicts.** G-R5: **retired as degenerate**. G-R5′: **killed**
  (both families).
- **Ledger delta.** R5 open→killed; R5′ registered→supported.
- **Where the thread stands.** Of the three levers that apply to *any*
  curve, L3 is killed, L4 is supported on degrees but pays only where the
  base degree is high, and their composition is killed. **L2
  (symmetrisation) is the only untested lever left**, and it is now the
  thread's whole remaining upside — with R4′ (a shape-corrected defect)
  needed first if `Δ_low` is to score it, since symmetrisation changes the
  variable count.

### 2026-09-11 — iteration 2 (EXP-R3 — degree falls are real, and they pay where the system is hard)

- **Task picked.** R3 (lever L4), per iteration 1's queue: with hybrid
  slicing dead, the remaining levers have to lower `D*` structurally, and L4
  was the cheapest to test.
- **Correction found before running anything.** Iteration 1 said EXP-J had
  "already identified the relation to add." It had not — EXP-J's relation is
  `Σ ℓ_i f_i ≡ 0`, a **pure syzygy**, which contributes no polynomial at all.
  Adding it is adding zero. The object worth adding is a **degree fall**: a
  top-degree cancellation with a *nonzero* remainder (§2.0). The experiment
  was rebuilt around that distinction, which is also what makes it
  non-trivial: the two populations are 0–4 syzygies against 20–112 falls per
  instance, and only the syzygy population had ever been counted.
- **Experiment** (`degree_reduction::{extract_degree_falls, saturate_with_falls,
  run_mutant_cell}`, `examples/degree_reduction_mutants.rs`, snapshot
  `experiments/degree_reduction_mutants.json`). Per non-decomposable target:
  measure `D*`; extract the degree-3 falls; iterate to saturation (≤ 6
  rounds, MutantXL-style, since a new generator creates new products);
  re-measure `D*` on the *same* target. Three families × `N ∈ {10,12,14}` ×
  8 matched targets.
- **Result — R3 SUPPORTED.** `D*` drops wherever there is headroom, and the
  effect is largest exactly where the defensive program says the problem is
  hardest:

  | family | base `D*` | augmented `D*` | working degree | improved |
  |---|---|---|---|---|
  | **Random** (generic) | 4.00 | **2.00** | 3.00 | 8/8 at `N`=12,14 |
  | Coordinate | 2.50 | 2.00 | 3.00 | 2/8 |
  | Subfield | 2.17 | 2.00 | 3.00 | 1/8 |

  `worsened = 0` in every cell, as ideal membership requires — the invariant
  that makes the measurement trustworthy, and a test asserts it.
- **Second correction, made after the first numbers.** The initial cost
  accounting charged only the augmented *solve* (`D* = 2`) and reported a
  saving of +8 to +11 bits. That books the same degree twice: extraction has
  to build the degree-3 rows and take a kernel, every round. Charging it
  (`G-R3′`) cuts the Random saving to **+1.44 bits mean** — still positive at
  every `N` and seed-robust (+1.44/+1.45/+1.49 over seeds 7/11/23) — and
  turns Coordinate break-even and Subfield sharply negative. The honest
  statement is about the **working degree** `max(3, D*_aug)`: L4 converts the
  generic family's degree-4 solve into a degree-3 one. Exactly one degree.
- **Result — R3′ regime-dependent.** Net saving `+1.44` (Random, supported),
  `−0.60` (Coordinate, sign varies), `−6.36` (Subfield, killed). The pattern
  is coherent: extraction costs a degree-3 climb, so a system that already
  solved near degree 2 pays for a degree it did not need. **L4 pays where the
  system is hard and costs where it is already easy** — the mirror image of
  L1, which only helps on special curves.
- **Gate verdicts.** G-R3: **supported** (all three families). G-R3′:
  **supported** on Random, **blocked** on Coordinate, **killed** on Subfield.
- **Ledger delta.** R3 open→supported; R3′ registered→regime-dependent.
- **Next.** Two things follow. (i) The Random-family saving *grows* with `N`
  (+0.94 → +1.75 → +1.62); whether that is a trend or noise needs `N = 16–18`,
  which the degree-3 extraction can reach cheaply even though `D*` cannot.
  (ii) R5 (composition) is now testable and is the interesting one: one-sided
  guessing reaches the floor at `c = 1/2`, and L4 lowers the degree on the
  hard family — composing them asks whether the mutant route reaches the
  floor at `c < 1/2`, which is the only way anything here beats the `2^N`
  baseline iteration 1 established.

### 2026-09-11 — iteration 1 (EXP-R1 — hybrid slicing is degenerate; one-sided guessing is not)

- **Task picked.** R1, the cheapest lever to instrument: P6 already says
  over-determination collapses `D*`, and `ρ` is a dial the attacker can
  turn by guessing variables. Highest refutation-power-per-CPU-hour of the
  four levers, and it establishes the cost baseline the others must beat.
- **Experiment** (`src/cryptanalysis/degree_reduction.rs`,
  `examples/degree_reduction_hybrid.rs`, snapshot
  `experiments/degree_reduction_hybrid.json`). For non-decomposable
  targets over `V` — so the system is unsatisfiable, every slice is
  unsatisfiable, every slice refutes, and the `2^k` slices *partition*
  what the direct solve covered — fix `k` of the `N = 2n'` Boolean
  variables and measure each slice's `D*`, `d_ff` and `Δ_low`. Three guess
  patterns (one-side / balanced / spread) × `N ∈ {10, 12, 14}` at the
  critical operating point `2n' = n`, 8 targets × 8 slices per cell.
- **Result — R1 KILLED, degenerately.** The degree collapse is real and
  large: mean `D*` falls from 2.75 → 2.00 (`N=10`) and 3.00 → 2.00
  (`N=14`) as `ρ` rises. But at **every one of 9 cells and both `ω`** the
  total-cost optimum sits at the largest `k` scanned. The model is not
  finding a hybrid; it is sliding toward exhaustive search, which is
  genuinely fastest at `N ≤ 14`. Against the `2^N` baseline the best
  "hybrid" comes in at −0.29 / −0.03 / +0.19 bits (`ω`=2.807) — i.e. it
  *is* enumeration, to within the accounting noise. The collapse fraction
  also rises with `N` for balanced (0.800 → 0.786 via 0.667) and spread
  (— → 0.750 → 0.714), the wrong direction for the lever to survive.
- **Result — R1′ SUPPORTED.** Guessing is not pattern-neutral. One-sided
  guessing — all `k` bits taken from the `X₁` half — hits the `D* = 2`
  floor at `k₂ = n'` exactly, hence `c = 1/2`, at **all three** sizes and
  **all three seeds** (7/11/23 — 9 cells, no exceptions), while balanced
  and spread need `c ≈ 0.64–0.80` and drift with `N`. The
  mechanism is structural rather than numerical: `n'` guesses determine
  `X₁` completely, leaving a system in `X₂` alone at `ρ = n/n' = 2`. So
  one-sided is the efficient frontier of *pure* guessing, and it is the
  pattern any composed lever (R5) should be built on.
- **Method corrections made during the iteration** (both changed the
  verdict, both are locked in by tests):
  1. The guess patterns were initially compared on **different targets** —
     `Spread` consumes randomness the other patterns do not, so a single
     RNG stream desynchronised the target draws. Split into independent
     target/slice streams; `k0_row_is_pattern_independent` now asserts the
     `k = 0` row is bit-identical across patterns.
  2. The gate was initially keyed on beating the **direct solve**, which
     any `2^k` model passes trivially at toy size. Replaced with the `2^N`
     enumeration baseline plus the `optimum_is_interior` check, which is
     what caught the degeneracy.
- **Incidental finding (→ R4′).** `Δ_low` moves non-monotonically in `k`
  while `D*` falls monotonically, because `Δ_low` is normalised by a
  monomial count that shrinks with the variable count. The `Δ_low` screen
  is calibrated *within* a system shape and cannot rank reformulations
  that change the shape. Registered as R4′; it gates R4.
- **Gate verdicts.** G-R1: **killed** (degenerate). G-R1′: **supported**.
- **Ledger delta.** R1 open→killed; R1′ registered→supported; R4′
  registered→open.
- **Next.** R1's kill removes the "just guess harder" shortcut, which
  means the remaining levers have to lower `D*` *structurally*. Queue
  order below puts L4 (EXP-R3) first: the relation to add is already
  identified by EXP-J, so it is the cheapest structural lever to test, and
  G-R3 is a clean two-sided gate — if `D*` does not move, L4 is empty and
  we learn the solver was finding that relation for free.

---

## 6. Where the thread stands, and the queue

**All five levers have now been tested.** That was the thread's charter, so
this is a natural reporting point rather than a pause.

| # | Lever | Outcome | Class (AGENTS.md §3) |
|---|---|---|---|
| **L1** | factor-base structure | the FFD program's own result; not this thread's subject | — |
| **L2** | symmetrisation | **blocked structurally** — cannot pay below `ℓ = 6`, and a dense Macaulay scan reaches `ℓ ≤ 3` | none (untested at the size where it could pay) |
| **L3** | hybrid slicing | **killed** — the cost optimum is exhaustive search | `accounting` — the apparent +8.83 bits was a boundary artifact |
| **L4** | degree falls (mutants) | **supported on degrees to `N = 16`, then degrades** — the cascade dies at `N = 18` and the degree saving halves; on Random the margin to `2^N` widens, and on Subfield augmenting is far worse than not | `relabelling` — the net saving tripled while the Random boundary margin lost 3.18 bits at the break |
| **L5** | the curve (isogeny class) | **killed** — `d_reg` is constant on the class for **all** `m` (proved, EXP-R6d), and 263 of `2^65.06` vertices are reachable | `accounting` — the boundaries moved, the algorithm did not |
| — | L3 ∘ L4 | **killed** — the levers are substitutes, not complements | none |

**Ratio to the boundary: moving the wrong way.** Per AGENTS.md §3 the only
class that counts as a result is `advance` — a fall in the ratio to the
stated floor. This thread's boundary is `2^N` enumeration (§3, fixed in
iteration 1), and **no lever tested beats it at any operating point
measured — 0 of 54 cells in the widest sweep (EXP-R3b).** Worse, on the
family where L4 actually pays it is *losing ground*: Random's margin goes
from `−9.5` bits at `N = 16` to `−12.7` at `N = 18` and `−12.1` at
`N = 20`, `−0.378` bits per size step overall.

The one quantity in this thread trending *toward* the boundary is the
**raw** solve on the Subfield family (`−4.99 → −1.19` bits over
`N = 10…20`), and it is not a result: `D*` there is pinned at the
Nullstellensatz floor, so a `poly(N)` solve against a `2^N` baseline has
to cross eventually. It measures that family's degeneracy — the FFD
program's L1 finding — and L4 makes it *worse*, not better (R3e).

L4's `+1.44` bits — and the `+4.44` that replaced it at larger `N` — are
gains *against the unaugmented Gröbner solve*, which is a race between two
routes that both lose. Every positive number this thread has produced is of
that kind.

That is the honest bottom line, and the pattern behind it is the thread's
main empirical finding about its own method: **six of seven iterations
found a metric that moved for a reason other than the one being claimed.**
Three levers were scored by a metric that counted one resource and ignored
another (`accounting`); one by a statistic that controlled for the nuisance
variable it had thought of and not the one that mattered; and one — L4 at
reach — by a gate that passed while the boundary went backwards
(`relabelling`). The gate in that last case was **pre-registered**, which
is why the boundary column is not optional: pre-registration fixes *when*
the metric is chosen, not *whether it is the right metric*.

**L5 is a later addition and the counts above do not include it.** It was
measured in iterations 9 (and iterations 1–3 of its own note,
`RESEARCH_ISOGENY_CLASS_SEARCH.md`) and classified `accounting` on the same
test: three of its four boundaries are derived rather than measured, so
nothing about the attacker's algorithm moved. It fits the pattern rather than
breaking it — one of its own numbers was wrong for exactly the reason this
section describes, a ρ cost that double-counted the negation map and made a
boundary margin look four bits wider than it is.

Nothing found here reduces the solving degree at a price worth paying in the
regime the instrument can see. The two positive results are narrow and
specific: mutants buy exactly one degree on the generic family (`+1.44` bits,
iteration 2), and one-sided guessing reaches the `D*=2` floor at `c = 1/2`
(iteration 1). The two most useful results are negative and structural:
levers aimed at the same floor do not compose (iteration 3), and the `Δ_low`
screen does not rank individual instances across shapes, however it is
normalised (iteration 4) — though it does rank **families**, which is what
the FFD program uses it for (iteration 6).

**The thread's third useful output is a correction to itself.** Iteration 4
carried its instance-level kill upstream and flagged the FFD program's
headline law as a probable size proxy. Iteration 6 tested that flag and
withdrew it: the law survives every size control, and strengthens under it
in the critical regime. The flag stood for two iterations and was cited in a
merged PR. The lesson is not "don't flag things" — flagging it is what got it
tested — but that an inference from the *shape* of someone else's
aggregation is a hypothesis, and this thread's own ledger discipline exists
to stop hypotheses being filed as findings.

**Queue, if the thread continues.** It should be said plainly that the
queue is thin and that none of it is likely to change the boundary verdict:

1. **EXP-R2b — reach for L2.** The only lever that is bounded rather than
   killed. Needs a sparse F4/F5 reaching ~51 variables at degree 3
   (`ℓ = 6`). This is an engineering project, not an experiment, and should
   not be started without deciding that L2 is worth that much — a decision
   that should weigh §6's boundary verdict, not just L2's open status.

*(EXP-R3b was run in iteration 7: the L4 trend continues on its own metric
and reverses against the boundary. EXP-R4d — queued as "EXP-R4c" — was run
in iteration 8: the defect does not mediate L4 on any reading, and changes
sign with the grouping. Both queue items that could be answered cheaply have
now been answered, and neither moved the boundary.)*

**So the queue has one item left, and it is not cheap.** Every lever in the
taxonomy has been measured; the only unmeasured claim is L2 above the
`ℓ = 6` crossover, which needs a solver this repo does not have. A thread
that has tested its whole taxonomy and found nothing that beats its boundary
should say that plainly rather than generate further variations of the
measurements it has already made.

*(L5's own follow-ups are all closed: EXP-R6b computed Boundary C to `m = 5`,
EXP-R6c explained its residual as decomposition yield, and EXP-R6d
(iteration 4) settled the induction that was the last open question there —
`a₆` stays strictly below the leading form at **every** `m`, proved rather
than sampled, so `d_reg` is curve-independent on every isogeny class for all
`m`. L5 is now killed by a theorem instead of by four data points. See
`RESEARCH_ISOGENY_CLASS_SEARCH.md` §2C′.)*

## 7. Honest limitations

- **Reach.** `D*` is measurable to `2n' = 20` with the dense single-pass
  `rank_and_refute` — iterations 1–6 asserted `14–16`, and EXP-R3b
  (iteration 7) measured 16, 18 and 20 at 14 s / 45 s / 167 s per 8-target
  cell. The earlier figure was a guess that nobody had tested, which is
  worth recording because it had been quietly limiting the queue: the
  single most informative result in this thread (R3c/R3d, the `N = 18`
  break) lives entirely beyond where the charter said the instrument
  stopped.
  Every number here is still toy-scale; the trends, not the magnitudes,
  are the result. R1's kill is a statement about the *shape* of the cost
  curve (boundary vs interior optimum), which is robust to scale in a way
  a margin in bits is not — but it is still three sizes.
- **A trend over six sizes is still a short lever arm.** R3d says the
  margin to `2^N` widens at `N ≥ 18`, on two operating points past the
  break and three seeds. It is enough to kill "the gap narrows"; it is not
  enough to fit an exponent to how fast it opens, and no exponent is
  claimed.
- **`D*` is a refutation degree.** It is measured on non-decomposable
  targets, i.e. on the failing relation attempts. That is the right cost
  driver for index calculus (most attempts fail), but it is not identical
  to the solving degree on a satisfiable instance.
- **The cost model is a proxy.** `cols^ω` (EXP-R1) and
  `rows·cols^{ω−1}` (EXP-R3) ignore sparsity, and real F4/F5 never builds the
  full Macaulay matrix. They are used only for *comparisons at matched
  shape*. EXP-R3 in particular charges extraction at the same crude rate as
  the solve; a real implementation would share work between the two, so the
  Random-family `+1.44` bits is a conservative floor rather than an estimate.
- **L4's saving is one degree, not two.** The `D* = 4.00 → 2.00` drop is
  real but must be read as `working degree 4 → 3`: the extraction itself
  operates at degree 3. Quoting the augmented `D*` alone double-counts.
- **Nothing here threatens a deployed curve**, and nothing here is
  expected to. The deliverable is a screen and a map of which
  presentations can and cannot help — which is a parameter-selection
  input, not an attack.

---

## References

- J.-C. Faugère, P. Gaudry, L. Huot, G. Renault, *Using symmetries in the
  index calculus for elliptic curve DLP*, J. Cryptology 2014. (L2)
- L. Bettale, J.-C. Faugère, L. Perret, *Hybrid approach for solving
  multivariate systems over finite fields*, J. Math. Cryptol. 2009. (L3)
- J. Ding, J. Buchmann, M. S. E. Mohamed, W. S. A. E. Mohamed, R.-P.
  Weinmann, *MutantXL*, SCC 2008; M. S. E. Mohamed et al., *MXL2*, PQCrypto
  2008. (L4 — the mutant/degree-fall mechanism)
- M.-D. Huang, M. Kosters, S. L. Yeo, *Last fall degree, HFE, and Weil
  descent attacks on ECDLP*, CRYPTO 2015. (what `D*` is)
- S. Galbraith, S. Gebregiyorgis, *Summation polynomial algorithms for
  elliptic curves in characteristic two*, INDOCRYPT 2014. (L1)
- Bardet–Faugère–Salvy, *complexity of Gröbner bases of semi-regular
  systems*. (the generic baseline and the `cols^ω` cost model)
- I. Semaev, *Summation polynomials and the discrete logarithm problem on
  elliptic curves*, ePrint 2004/031.
