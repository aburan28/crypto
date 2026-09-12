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

`D*` depends on the pair (ideal, presentation). Four levers, ordered by
how much they are allowed to change:

| # | Lever | What it changes | What it costs | Status |
|---|---|---|---|---|
| **L1** | **Factor-base structure** — subfield, Koblitz, sparse normal basis | the *ideal itself* (multiplicative closure injects relations) | only works on special curves/fields | **measured** by the FFD program: Subfield mean `D*` 2.04 vs Random 3.53 at `2n'=n`, and `Δ_low(Subfield)/Δ_low(Random)` diverges 6.7 → 67 |
| **L2** | **Symmetrisation** — solve in the elementary symmetric variables `e_1..e_m` instead of `x_1..x_m` (Faugère–Gaudry–Huot–Renault) | the *variables*; the ideal is the same up to a change of coordinates | free (a one-off rewrite) | infrastructure exists (`symmetrized_semaev.rs`) but **has never been measured against `D*` or `Δ_low`** |
| **L3** | **Hybrid slicing** — guess `k` variables, solve `2^k` slices, raising `ρ = #eqs/#vars` | the *determination ratio* | `2^k` multiplicative | **killed, iteration 1** (degenerate optimum) — but one-sided guessing is the cheapest route to the `D*=2` floor, and iteration 3 found it *dominates* the mutant route once guessing is allowed |
| **L4** | **Degree falls (mutants)** — add the *nonzero* low-degree remainders of top-degree cancellations to the generator set, so `x_k · g` rows become available a degree early | the *generating set*; ideal and variables unchanged | the extraction's own climb to degree 3 | **supported on degrees, regime-dependent on cost — iteration 2.** `D*` drops 4.00 → 2.00 on the generic family; net of extraction cost it pays only where the base degree is high |

L1 is the known part of the map and is not this thread's subject. L2, L3
and L4 are presentation changes that apply to *any* curve, which is what
makes them worth measuring.

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
back out of. At the sizes where `D*` is measurable (`N ≤ 16`), `2^N` is
*small*, so any cost model with a `2^k` term will happily slide its
optimum to the largest `k` available and report a large "saving" against
the direct solve. What it has actually found is that exhaustive search is
fastest on a 14-variable system — which says nothing about cryptographic
`n`. Every sweep in this thread therefore reports
`HybridSweep::optimum_is_interior`, and a boundary optimum is recorded as
a **degenerate** result, not a win.

Cost model throughout: `cost(D, N) = cols(N, D)^ω`, with the expectation
taken over the *histogram* of `D*` rather than its mean (Jensen — a
spread of degrees costs strictly more than its average; asserted in
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
| **R2** | **Symmetrisation (L2) lowers `D*`** on the descended system at matched `(n, n')` | `open` | — (needs a symmetrised descent; `symmetrized_semaev.rs` has the algebra but not the descent) |
| **R3** | **Adding degree falls (L4) as explicit generators lowers `D*`** at matched targets | **`supported`** | EXP-R3, iteration 2. All three families, 3 operating points each, 8 matched targets per cell: `D*` strictly lower wherever there was headroom. Generic (Random) family **4.00 → 2.00** on 8/8 targets at `N = 12, 14`. `worsened = 0` everywhere, as the ideal-membership invariant requires. |
| **R3′** | The `D*` drop **survives its own cost** — extraction must climb to degree 3, so the net saving must still be positive | **regime-dependent** | EXP-R3: net `+1.44` bits mean on Random (positive at every `N`, seeds 7/11/23 give +1.44/+1.45/+1.49); `−0.60` on Coordinate (sign varies); `−6.36` on Subfield (**killed** — the system already solved at `D* ≈ 2.1`, so the climb to 3 is pure overhead). L4 pays where the system is hard and costs where it is easy. |
| **R4** | **Every lever acts through `Δ_low`**: pooled across lever-generated systems, `ρ_s(Δ_low, D*) ≤ −0.6` | **`killed`** | Follows from R4″: a pooled correlation over mixed sizes is a size proxy here, so "pooled `ρ_s`" cannot establish the claim however the defect is scaled. Scoring a lever needs matched shape and total work, not a defect correlation. |
| **R4′** | A **shape-corrected** defect exists that is comparable across systems with different variable counts | **`killed`** | EXP-R4′, iteration 4. Five normalisations × three designed groups × four seeds. Nothing clears `ρ_s ≤ −0.6` on all three; the un-normalised variants fail cross-shape as expected (the control that shows the test has power), and the normalised ones pass *pooled* only. |
| **R4″** | The strong pooled defect↔`D*` correlation is a **size proxy**, not structure | **`supported`** | EXP-R4′: pooled `ρ_s` −0.89…−0.91 collapses to **−0.16…−0.34** once `vars` is held fixed, for every variant that passed pooled. Stable over seeds 7/11/23/41 at 48 targets/cell. |
| **R5** | **Levers compose**: one-sided guessing plus the mutant route beats guessing alone | **`killed`** | EXP-R5, iteration 3. The pre-registered gate (collapse fraction `c < 1/2`) is **degenerate** — the composed route hits `c = 0` in every cell, because mutants reach the floor with no guessing at all. Scored on total work instead (G-R5′): composed loses to raw guessing by a **flat −2.87 bits** at every `N` and seed, and neither route beats `2^N`. |
| **R5′** | Mutants and guessing are **substitutes, not complements** — both drive the system to `D* = 2`, and guessing gets there more cheaply per unit work | **`supported`** | EXP-R5: at `k = 0` the mutants are worth `+1.1…+1.75` bits on Random (iteration 2's result), but the moment guessing is allowed the advantage inverts and stays inverted at every `k > 0`. The gap is flat in `N`, so it is structural, not a small-size artifact. |

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
  below the raw system at matched `(n, n')` and matched targets, over ≥ 3
  operating points, with the gap not shrinking in `N`. *Killed* if the gap
  is ≤ 0 at the largest `N`.
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
  recorded as a size proxy.
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

---

## 5. Iteration log

> Newest at top. Format mirrors `RESEARCH_FFD_WORKFLOW.md` §7:
> *Task · Experiment · Result · Gate verdict · Ledger delta · Next.*

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
  is a size proxy** (→ R4″). No denominator repairs that, because the
  problem is not the denominator.
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
- **Gate verdicts.** G-R4′: **killed**. G-R4″ (registered this iteration):
  **supported**. G-R4: **retired** — a pooled `ρ_s` over mixed sizes cannot
  establish it.
- **Ledger delta.** R4′ open→killed; R4 open→killed; R4″ registered→
  supported. §2.1's iteration-1 diagnosis corrected.
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

## 6. The experiment queue

Re-prioritised from the ledger each iteration; this is the current guess.

1. **EXP-R2 — symmetrisation on `S₄` (R2).** *(Queue head after iteration
   4.)* The thread's whole remaining upside, and no longer blocked: R4′
   settled that `Δ_low` should not be used to score it, so L2 is measured on
   `D*` and total work at matched `(vars, eqs, ρ)`.
   `binary_semaev_s4.rs` already builds the symmetrised `S₄` descent plus
   its correspondence system at `m = 3` — the regime where the `m!` saving
   is meaningful, unlike the `m = 2` (`S₃`) harness the rest of this thread
   uses. Two build items: (a) a Macaulay/refutation path that accepts
   arbitrary-degree ANF, since the correspondence half is cubic and the
   current harness is quadratic-only; (b) the eliminated presentation
   (substitute `eᵢ = σᵢ(x)`) as the matched non-symmetrised baseline — same
   ideal, no new algebra needed.
3. **EXP-R3b — reach for the R3′ trend.** The Random-family net saving grows
   +0.94 → +1.75 → +1.62 over `N = 10,12,14`. Degree-3 extraction is cheap
   enough to run at `N = 16–18` even where `D*` itself is out of reach, using
   the *working degree* rather than measured `D*` on the base side.
4. **EXP-R2b — symmetrisation composed with L4.** Only if R2 shows
   symmetrisation moves `D*` at all — and iteration 3's R5′ is a warning
   that two levers aiming at the same floor tend to be substitutes, so this
   should be scored on total work from the start, never on a degree count.

---

## 7. Honest limitations

- **Reach.** `D*` is measurable to `2n' ≈ 14–16` with the dense
  single-pass `rank_and_refute`. Every number here is toy-scale; the
  trends, not the magnitudes, are the result. R1's kill is a statement
  about the *shape* of the cost curve (boundary vs interior optimum),
  which is robust to scale in a way a margin in bits is not — but it is
  still three sizes.
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
