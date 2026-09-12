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
| **L3** | **Hybrid slicing** — guess `k` variables, solve `2^k` slices, raising `ρ = #eqs/#vars` | the *determination ratio* | `2^k` multiplicative | **killed, iteration 1** — see §5 |
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

Any lever that changes the variable count (L3 does; L2 does) therefore
needs a **shape-corrected** defect before `Δ_low` can score it. That is
now ledger row **R4′**, and it is a prerequisite for R4 rather than a
by-product of it.

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
| **R4** | **Every lever acts through `Δ_low`**: pooled across lever-generated systems, `ρ_s(Δ_low, D*) ≤ −0.6` | `open` | blocked on R4′ |
| **R4′** | A **shape-corrected** defect exists that is comparable across systems with different variable counts | `open` | EXP-R1 shows the raw `Δ_low` is not (§2.1) |
| **R5** | **Levers compose**: a structural lever (L1/L2/L4) plus one-sided guessing reaches the floor at `c < 1/2` | `open` | — (the composition test; the only route left to beating `2^N`, given R1) |

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
  cells. *Killed* at `|ρ_s| < 0.2` or a sign flip.
- **G-R5.** *Supported* if `c(composed) < c(one-side)` at ≥ 2 operating
  points, seed-robust.

---

## 5. Iteration log

> Newest at top. Format mirrors `RESEARCH_FFD_WORKFLOW.md` §7:
> *Task · Experiment · Result · Gate verdict · Ledger delta · Next.*

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

1. **EXP-R5 — composition (R5).** *(promoted after iteration 2.)* One-sided
   guessing reaches the `D*=2` floor at `c = 1/2` (R1′); L4 lowers the
   working degree on the hard family (R3). Compose them: does the mutant
   route reach the floor at `c < 1/2`? This is the only remaining path to
   beating the `2^N` baseline iteration 1 established, and both halves are
   now built.
2. **EXP-R3b — reach for the R3′ trend.** The Random-family net saving grows
   +0.94 → +1.75 → +1.62 over `N = 10,12,14`. Degree-3 extraction is cheap
   enough to run at `N = 16–18` even where `D*` itself is out of reach, using
   the *working degree* rather than measured `D*` on the base side.
3. **EXP-R4′ — shape-corrected defect (R4′).** Needed before `Δ_low` can
   score any lever that changes the variable count. Candidates: normalise
   by the *generic* rank rather than the column count, or compare raw
   syzygy counts `Σδ` at matched `(vars, eqs)`. Pure post-processing of
   data the existing harness already emits.
4. **EXP-R2 — symmetrised descent (R2).** The most interesting lever and
   the most work: `symmetrized_semaev.rs` has the symmetric-function
   algebra but the descent is built on raw coordinates, so the
   `e`-variable descent has to be written. Do it after R4′ so the result
   can be scored on both `D*` and a defect that means something.
5. **EXP-R2b — symmetrisation composed with L4.** If R2 shows symmetrisation
   moves `D*` at all, the two generating-set levers should be measured
   together rather than separately.

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
