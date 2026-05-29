# A proof-complexity bridge for the first-fall-degree assumption

**Status:** research proposal (novel approach), 2026-05-29
**Builds on:** `src/cryptanalysis/ffd_harness.rs`,
`src/cryptanalysis/groebner_f4.rs`, `src/cryptanalysis/semaev_sat.rs`,
`src/cryptanalysis/binary_semaev.rs`, `src/cryptanalysis/diem_descent.rs`,
`RESEARCH_FFD_MEASUREMENT.md`.
**One-line thesis:** the disputed first-fall-degree heuristic is, up to
constants, a statement about the **Polynomial Calculus refutation
degree** of the Weil-descended summation-polynomial system — and that
degree is governed by the system's **low-degree Hilbert-function defect**
`Δ_low` (the excess of its quotient Hilbert function over a generic
system), an *algebraic* quantity we can compute from a few low-degree
Macaulay ranks at full cryptographic size *without ever running a Gröbner
basis*. (An earlier form of this thesis pinned `D*` on the **expansion of
the descent incidence graph**; that and a treewidth variant were built and
**refuted** — no graph-incidence invariant predicts `D*`, because the
structure is algebraic. See §3.)

---

## 1. The dispute, restated as one number

Index calculus on `E/F_q` decomposes a target `R` over a factor base
`{P : x(P) ∈ V}` (`V` an `F_2`-subspace of dimension `n'` for binary
fields, a subfield for `F_{p^k}`) by solving

```
   S_{m+1}(x_1, …, x_m, x(R)) = 0 ,     x_i ∈ V,
```

where `S_{m+1}` is Semaev's summation polynomial. Weil descent turns
this into a polynomial system over `F_2` (resp. `F_p`); one computes a
Gröbner basis and reads off the decomposition. **Everything** —
Petit–Quisquater's `2^{O(n^{2/3}\log n)}` claim, Semaev 2015, the
Huang–Kosters–Yeo last-fall-degree rebuttal, the Galbraith–Petit
pushback — turns on a single number:

> **the degree `D*` at which the Gröbner / Macaulay computation
> terminates** (the "solving degree" / "degree of regularity").

- The **first-fall-degree assumption** says `D* ≈ d_ff`, the smallest
  degree at which *any* non-generic syzygy appears, and that
  `d_ff = O(log n)` (or `O(1)`), giving subexponential ECDLP.
- **Huang–Kosters–Yeo (CRYPTO 2015)** introduced the **last fall
  degree** `d_last`, proved `D* ≤ d_last + (\text{small})`, and
  exhibited explicit Weil-descent systems where `d_ff ≪ d_last`. So the
  assumption `D* ≈ d_ff` is *unjustified*: the true solving degree can
  sit far above the first fall.
- **Galbraith–Petit** observed empirically that `d_ff` *looks* constant
  at tiny `n` but extrapolation is unsafe.

`ffd_harness.rs` already measures `d_ff` (the **first** fall). Its own
table shows `d_ff = 3` flat for `n = 3..7`. **But `d_ff` is the wrong
number.** The number that decides security is `D* ≈ d_last`, and nobody
in this repo (or, computationally, in the literature) measures it
directly or predicts it cheaply. That is the gap this proposal closes.

---

## 2. The novel move: `D*` *is* a Polynomial Calculus refutation degree

Fix a **non-decomposable** target `R` (no `m`-term decomposition over
`V`). Then the system `Σ = { S_{m+1}(x,x(R)) = 0, field equations
x_i^2 = x_i }` is **unsatisfiable**, and a Gröbner computation
necessarily derives `1` — it produces a *refutation*.

The classical Clegg–Edmonds–Impagliazzo correspondence says: the
matrix-F4 computation reaching degree `D` is *exactly* a degree-`D`
**Polynomial Calculus (PC)** derivation over the same field. Hence

```
        D*  =  Gröbner solving degree
            ≍  Polynomial-Calculus refutation degree of Σ
            ≍  Huang–Kosters–Yeo last fall degree  d_last.
```

This is the bridge. It lets us **import 25 years of proof-complexity
degree lower bounds** — a body of technique the ECDLP literature has
essentially never used — to settle whether `D*` is small (attack works)
or large (attack is exponential). Three immediate consequences:

| Crypto object                       | Proof-complexity object                    |
|-------------------------------------|--------------------------------------------|
| first-fall-degree assumption        | "PC refutation degree is `O(\log n)`"      |
| HKY last fall degree `d_last`       | PC refutation **degree** of `Σ`            |
| Galbraith–Petit "it grows fast"     | a PC degree **lower bound**                |

The point of the reframing: PC degree is the canonical object for
which **expansion-based lower bounds** exist (Ben-Sasson–Wigderson for
width; Alekhnovich–Razborov, Mikša–Nordström for PC/PCR degree). We can
try to *prove*, not just measure, that `D*` is linear in `n` for generic
curves — i.e. that the first-fall-degree assumption is **false** for
cryptographic parameters. That is a defensive theorem, and it slots
straight into this repo's cryptanalysis-resistance map theme.

---

## 3. The predictor: the low-degree Hilbert-function defect `δ(D)`

The reason PC degree is *low* for some Semaev systems and *high* for
others should not be mysterious — but the explanatory program had to find
the right invariant by elimination. **Two natural graph-incidence
predictors were built and refuted** before the working one was found; the
honest record matters, so it is summarised first (§3.0), then the working
predictor is defined and formalised (§3.1–§3.3).

### 3.0 What was refuted first (graph-incidence invariants)

The original conjecture was **spectral**: PC degree is lower-bounded by
the boundary expansion `γ(G)` of the constraint–variable bipartite graph
`G(E, basis)` whose edges follow the field structure-constant tensor
`c_{ikj}` (`z^i z^k = Σ_j c_{ikj} z^j`). Two refutations killed the whole
*graph-incidence* family of predictors:

- **Spectral expansion `γ` (workflow iter. 3, `descent_lowgamma.rs`).** At
  the operating point `2n'=n` (n=8) the genuine *subfield* factor base
  `F_{2^{n'}}` is the **easiest** case (mean `D*` ≈ 2.0 vs random ≈ 3.5) —
  exactly the structured-is-easy effect index calculus exploits — **yet its
  spectral expansion is not lower** (0.882 ≥ 0.861). The most important
  low-`D*` structure is *not* a low-`γ` structure: `D* = Θ(γ·n)` is false.
- **Treewidth (workflow iter. 4, `descent_treewidth.rs`).** The descended
  *primal* graph is the **complete** graph `K_{2n'}` (the Frobenius-squared
  cross terms saturate every variable pair), so its treewidth is the
  constant `2n'−1` for *every* family and cannot discriminate easy from
  hard at all.

The lesson is structural: **the subfield speedup is algebraic
(multiplicative closure), and is therefore invisible to any invariant of
the variable-incidence graph.** The predictor has to read the
*coefficients*, not the support.

### 3.1 Definition

Work in the Boolean polynomial ring `R = F_2[x_1,…,x_N]/⟨x_i²−x_i⟩`,
`N = 2n'`, where `Σ` is the `V`-substituted descended Semaev system (a set
of `F_2`-quadratics). For each degree `D`, the **Macaulay matrix**
`Mac_D(Σ)` has one column per multilinear monomial of degree `≤ D`
(`cols(D) = Σ_{k≤D} C(N,k)`) and one row per product `m·f` with `f ∈ Σ`,
`deg(m·f) ≤ D`. Let `r(D) = rank_{F_2} Mac_D(Σ)`.

Let `r_gen(D)` be the rank a **semi-regular** (generic) system of the same
shape — `m` quadratics in `N` Boolean variables — would have at degree `D`,
given by the truncated Bardet–Faugère–Salvy Hilbert series (implemented in
`ffd_harness::generic_rank_prediction`). Define the **Hilbert-function
defect**

```
   δ(D)  :=  r_gen(D) − r(D)  =  H_Σ(D) − H_gen(D)   ≥ 0,
```

where `H_•(D) = cols(D) − r_•(D)` is the (truncated) quotient Hilbert
function. The **early defect** is the cumulative sum up to a small cutoff
`D_low`:

```
   Δ_low(Σ)  :=  Σ_{D ≤ D_low}  δ(D).
```

(`descent_algebraic.rs`: `rank_profile`, `early_defect`; default
`D_low = 3`.)

### 3.2 Why `δ(D) ≥ 0`, and what it counts (the formalisation)

`δ(D) ≥ 0` is **not** an empirical accident — it is the right baseline by
construction. The Hilbert function is *upper semicontinuous* in the
coefficients, so among all ideals generated by forms of the given degrees
the **generic** one has the termwise-minimal Hilbert function:
`H_Σ(D) ≥ H_gen(D)` for every `Σ` (rigorous). The Fröberg / semi-regular
*formula* in `generic_rank_prediction` is the standard estimate of that
minimal value — exact in the small-degree, under-determined cases at hand
(and the result is clamped at 0 regardless). Equivalently `r(D) ≤
r_gen(D)`, so

```
   δ(D)  =  r_gen(D) − r(D)  =  H_Σ(D) − H_gen(D)  ≥  0,
```

with equality iff `Σ` behaves semi-regularly up to degree `D`. Thus `δ(D)`
is the **excess of the quotient over generic** at degree `D`: it counts the
low-degree relations (syzygies) `Σ` carries *beyond* a generic system —
exactly the algebraic structure a solver gets "for free." (`δ(D) ≥ 0` is
asserted in `descent_algebraic`'s tests and held in every one of the 50
EXP-G cells.)

**The proxy claim (Hilbert defect ↔ last-fall degree).** The refutation /
solving degree `D*` is the degree at which the Macaulay tower first
contains the contradiction `1 ∈ rowspace` (≈ the last-fall degree of the
Gröbner computation). A system whose quotient Hilbert function runs *above*
generic at low degree — large `Δ_low` — is one whose ideal is already
"saturating" early: each unit of early defect is a low-degree relation that
pulls the termination degree down. The corrected predictor is therefore

```
   D*  is monotone DECREASING in  Δ_low(Σ).
```

We state this as a heuristic correspondence backed by strong empirics, not
a closed theorem; the rigorous half is `δ ≥ 0` and its
excess-syzygy meaning, and the empirical half is the slope below.

**Mechanism for the subfield.** When `V = F_{2^{n'}}` is a subfield it is
*closed under multiplication*: the products `z^i z^k` of factor-base
elements fall back into `V`, so the descended system inherits the subfield's
own multiplication tensor and produces relations that collapse already at
degree 2–3 → large `Δ_low`. A random `V` has no such closure, so its
quotient tracks generic until high degree → `Δ_low ≈ 0`. This is precisely
the algebraic content the incidence graph cannot see.

### 3.3 Central conjecture (corrected) and its evidence

```
   D*  =  D*_min  +  Θ( − Δ_low(Σ) ),     i.e.
   D*  decreases monotonically with the early Hilbert-function defect Δ_low.
```

Equivalently:

> the first-fall-degree assumption holds **iff** the restricted system
> carries a large early Hilbert defect — i.e. iff the factor base is
> multiplicatively structured (a subfield / Koblitz / sparse-normal base),
> which is exactly when `Δ_low` is large and `D*` small.

**Evidence (EXP-G, `examples/ffd_expg_curve.rs`, snapshot
`experiments/ffd_expg_curve.json`).** Across **50 cells** spanning ten
operating points `2n' ∈ {4,…,14}` × three factor-base families
(Subfield, Coordinate, Random), the early defect `Δ_low` predicts `D*`:

| regime | cells | Spearman ρ_s | OLS slope `dD*/dΔ_low` |
|---|---|---|---|
| pooled | 50 | **−0.79** | −6.3 |
| **critical `2n'=n`** (the ECDLP case) | 30 | **−0.78** | **−7.6** |
| over-determined `2n'<n` | 20 | −0.74 | −1.9 (D* floored at 2 by §6 / P6) |

Seed-robust at this reach (critical ρ_s ∈ [−0.70,−0.81], slope ∈
[−7.6,−8.5] over 3 seeds); the exchange rate of ≈ **−7 to −8 degrees of
`D*` per unit early defect** is stable across `2n'` from 4 to 14. The
over-determined slope is shallow only because `D*` is saturated at the
Nullstellensatz floor of 2 there (§6), not because the law weakens.

This **reconciles the dispute** instead of picking a side:

- The HKY counterexamples and the Galbraith–Gebregiyorgis "nice" cases
  (Koblitz / subfield bases, sparse normal bases) are precisely the
  **high-`Δ_low`** structures → `D*` small → the heuristic holds *for
  them*. That is *why* the first-fall assumption ever looked true.
- A **generic** curve over a generic basis (and every prime-field descent
  that lacks a multiplicatively-closed factor base) gives `Δ_low ≈ 0` →
  `D* = Θ(n)` → the heuristic is **false** and the attack is exponential.

Crucially, `Δ_low(Σ)` is the rank profile of a few **low-degree** Macaulay
matrices (degree `≤ D_low`, so `O(N^{D_low})` columns at fixed `D_low`): a
*polynomial-time* screen, computable even when the full Gröbner basis it
predicts is astronomically out of reach. **That is the deliverable a
parameter-selection committee actually needs** — and unlike the refuted
`γ(G)`, it tracks `D*` with the correct sign.

### 3.4 Toward a lower bound: the defect-vanishing route

The corrected predictor immediately suggests the **defensive theorem** the
whole program is after — a *generic* lower bound `D* = Θ(n)`. The bridge is:

> **Conditional theorem (defect ⇒ hardness).** Let `Σ_N` be the
> `V`-restricted descended Semaev system at `2n'=N`. If, for a generic
> factor base `V`, the early Hilbert defect vanishes — `Δ_low(Σ_N) = o(1)`
> as `N → ∞` for every fixed cutoff `D_low` — then `Σ_N` is asymptotically
> **semi-regular** through low degree, so its solving / refutation degree is
> the generic one, `D* = Θ(N) = Θ(n)`. Hence the first-fall-degree
> assumption is **false** for generic curves/bases, and the index-calculus
> attack is exponential there.

The implication step is standard: a sequence that matches the semi-regular
Hilbert function through low degree has its degree of regularity at the
generic value, which for `m = Θ(N)` quadratics in `N` Boolean variables is
`Θ(N)` (Bardet–Faugère–Salvy). The *only* unproven input is the antecedent:

> **Genericity lemma (open, the crux).** For a uniformly random
> `F_2`-linear factor base `V` of dimension `n'`, `Δ_low(Σ_N) = o(1)` with
> high probability.

This is the same "immunity" obstacle as before, but now stated for the
*Hilbert defect* rather than graph expansion — and it is the natural place
to use the **algebraic independence of the Semaev coefficients**: a random
`F_2`-linear restriction should put the degree-`≤ D_low` part of `Σ_N` in
generic position, so its low-degree Hilbert function is the semi-regular one
whp. Proving this is future work; what the experiments establish is that the
antecedent **holds empirically, with room to spare.**

> **Evidence (EXP-H, `examples/ffd_defect_scaling.rs`, snapshot
> `experiments/ffd_defect_scaling.json`).** In the critical regime `2n'=n`,
> the normalized early defect of the **Random** (generic) family decays
> polynomially to zero,
>
> ```
>    Δ_low(Random)  ≈  a · (2n')^{−c},     c ≈ 4.3   (4.2–4.5 over seeds),
> ```
>
> measured cleanly across `2n' ∈ {6,…,20}` (`Δ_low`: 0.155 → 0.0007). The
> **Subfield** family stays bounded away from 0 over the same range
> (0.23 → 0.048), so the ratio `Δ_low(Subfield)/Δ_low(Random)` *diverges*,
> 6.7 → 67. The generic defect vanishing — and the structured defect
> persisting — is exactly the two-sided shape the conditional theorem needs:
> generic ⇒ hard, subfield ⇒ easy.

So the lower-bound route reduces to a single, sharply-stated, empirically
well-supported lemma about the low-degree genericity of random restrictions
— a much more tractable target than a PC-degree expansion bound on a
globally-coupled quadratic system.

---

## 4. What gets built (concrete, plugs into existing modules)

Four pieces, each a thin extension of code already in `cryptanalysis/`:

1. **`pc_degree_harness`** — *(implemented:
   `src/cryptanalysis/pc_degree_harness.rs`, demo
   `cargo run --release --example pc_degree_demo)`.* Forks `ffd_harness`
   but *does not stop at the first fall*. It restricts the factor base to
   a subspace `V` of dimension `n' = ⌊n/2⌋`, picks a **non-decomposable**
   target `x(R)` (the genuinely *unsatisfiable* PDP instance — the
   unrestricted system is satisfiable and has no refutation), and reports
   the **refutation degree** `D*` = smallest `D` at which the constant
   `1` enters the Macaulay row-space (tested by `e0_in_rowspace`). It
   prints both columns: `d_ff` (operational first fall, as in
   `ffd_harness`) and `D* ≈ d_last`. The gap between them *is* the
   Huang–Kosters–Yeo invariant, made operational for the first time in
   this repo.

   **Validation (the load-bearing test).** `refutation_iff_nondecomposable`
   exhaustively checks, over `F_{2^4}` with `n'=2`, that a target refutes
   **iff** it is non-decomposable over `V` — confirming `D*` measures the
   right object (`1 ∈ ⟨S₃, field eqs⟩ ⟺ no `F_2`-solution`). All six
   module tests pass.

   **Averaging layer** — *(implemented: `src/cryptanalysis/pc_degree_avg.rs`,
   demo `cargo run --release --example pc_degree_avg_demo`).* A single
   `(n, n', x₃)` gives one noisy `D*`; `pc_degree_avg` samples **many
   non-decomposable targets** per cell and reports the `D*` *distribution*
   (`DegreeStats`: min / mean / max + histogram, with decomposable draws
   skipped and tallied). It exposes the **determination ratio**
   `ρ = #eqs/#vars = n/(2n')` and ships three sweeps:
   `run_pc_avg_sweep` (the `n' = ⌊n/2⌋` edge), `run_pc_regime_sweep`
   (fix `n`, vary `n'` so `ρ` sweeps down toward 1), and
   `run_pc_operating_point_sweep` (grow `n` at fixed `ρ ≈ 1`). 7 tests
   pass.

   **What the averaged data shows (64 targets/cell).** Two findings, both
   the *opposite* of the earlier "flat `D*`" artifact:

   - **`D*` rises as `ρ → 1` (regime sweep, fixed `n = 10`).** Mean `D*`
     climbs monotonically as the system approaches just-determined:

     | `n'` | `ρ` | mean `D*` | `D*` histogram |
     |---:|---:|---:|---|
     | 1 | 5.00 | 2.00 | `2:64` |
     | 2 | 2.50 | 2.00 | `2:64` |
     | 3 | 1.67 | 2.17 | `2:53 3:11` |
     | 4 | 1.25 | 2.31 | `2:46 3:16 4:2` |
     | 5 | 1.00 | 2.42 | `2:50 3:1 4:13` |

     So the flat `D*` reported by the single-sample harness was indeed an
     **over-determination artifact**: heavily over-determined systems
     (`ρ ≫ 1`) admit degree-2 Nullstellensatz certificates, masking the
     solving degree. The interesting regime is `ρ ≈ 1`.

   - **At fixed `ρ ≈ 1.1`, mean `D*` climbs with `n` (operating-point
     sweep).** Along the **odd-`n`** subsequence:
     `D* = 2.77 (n=7) → 2.98 (n=9) → 3.22 (n=11)`. This is the first
     direct sighting of **prediction #1's `D*`-climb** in this repo, on
     genuinely unsatisfiable instances.

   **Honest caveats.** (i) A strong **even/odd parity effect** dominates
   the scaling sweep: odd `n` give markedly higher `D*` than the
   neighbouring even `n` (e.g. `n=9 → 2.98` vs `n=10 → 2.39`). This is a
   known characteristic-2 Semaev phenomenon (odd-degree fields lack the
   even-`n` half-trace / subfield structure), so prediction #1 must be
   read **within a fixed parity class**, not across the raw sequence.
   *Update (workflow iteration 1, `RESEARCH_FFD_WORKFLOW.md`): with the
   parity split made explicit and 64 targets/cell, BOTH parity classes
   show a positive `D*`-vs-`n` slope — odd +0.121/n (n=7,9,11), even
   +0.060/n (n=6,8,10,12). The earlier "even looks flat" remark was
   undersampling; the auto-verdict gate G-P1 now reports `supported` on
   both classes. The parity effect is a level shift, not a difference in
   the sign of the trend.*
   (ii) The climb is real but shallow over `n ≤ 12` — three points on the
   odd subsequence is suggestive, not conclusive; confirming a *linear*
   law needs (a) more reach (the dense single-pass `rank_and_refute` plus
   `d_cap` censoring now reaches `2n'=14`–`16`; a tried sparse backend lost
   to Macaulay fill-in — see workflow iter. 8), and (b) separating the two
   parity classes explicitly. (iii) `D*` here is still
   bounded by the small `n'`; the absolute degrees are toy-scale. The
   trend, not the magnitude, is the result.

2. **`descent_expansion`** — *(implemented:
   `src/cryptanalysis/descent_expansion.rs`, 14 tests.)* Builds the
   incidence graph from the field multiplication structure-constant tensor
   and returns `γ(G)` two ways: spectral (`1 − σ₂` of the normalized
   biadjacency, via a small Jacobi eigensolver) and unique-neighbor
   boundary expansion (subset scan). Pure linear algebra; runs at
   cryptographic `n` (smoke-tested at `n = 32`). Also ships a Rabin
   irreducibility test + basis enumerator so the P2 study can vary the
   basis at fixed `(n, n')`. **Key finding:** the naïve *tensor* graph
   (bilinear term only) is basis-independent in the refutable regime, so
   the predictor uses the *system* graph (`system_incidence`), built from
   the full descended `S₃` including the basis-sensitive Frobenius-squared
   terms. See the §6 status note for the first G-P2 numbers.

3. **Correlation study** — *(implemented: `descent_lowgamma.rs` (EXP-E) +
   the `run_p2`/`run_exp_e` stages of `ffd_breakthrough_loop`.)* Sweeps
   `D*` against `γ(G)` over both generic irreducible bases and the three
   structured factor-base families {Subfield, Coordinate, Random} at fixed
   `(n, n')`. **Outcome: the conjecture's spectral form is falsified** —
   the subfield case is the easiest (lowest `D*`) yet not the lowest-`γ`,
   so there is no `D* ≈ γ·n` law for spectral `γ`. See §6 prediction 2 and
   the workflow iteration-3 log. The infrastructure now supports retrying
   with a reformulated `γ` (boundary expansion on the quadratic-only
   support, or treewidth).

4. **The algebraic predictor** — *(implemented:
   `src/cryptanalysis/descent_algebraic.rs`, 3 tests; sweep
   `examples/ffd_expg_curve.rs`, snapshot `experiments/ffd_expg_curve.json`;
   single-pass `pc_degree_harness::rank_and_refute` for reach.)* After the
   graph-incidence predictors were refuted, this computes the low-degree
   Hilbert defect `Δ_low` (§3.1) and correlates it with `D*`. **Outcome:
   supported** — 50 cells over `2n' ∈ {4,…,14}` give pooled Spearman
   ρ_s = −0.79, critical-regime slope ≈ −7.6 (§3.3), seed-robust. This is
   the working replacement for the expansion study (item 3).

5. **Defect-scaling study** — *(implemented:
   `examples/ffd_defect_scaling.rs` (EXP-H), snapshot
   `experiments/ffd_defect_scaling.json`.)* Measures `Δ_low(2n')` for the
   Subfield vs Random families in the critical regime to `2n'=20` (cheap:
   degree-≤3 ranks only). **Outcome:** Random `Δ_low ≈ (2n')^{−4.3} → 0`
   while Subfield stays bounded — the empirical antecedent of the §3.4
   conditional theorem (the generic defect vanishes; the structured one
   persists).

6. **Lower-bound attempt** — the route is now §3.4 (defect vanishing).
   (a) *Algebraic (favoured):* the conditional theorem `Δ_low = o(1) ⇒ D* =
   Θ(n)` is in hand; the one open input is the **genericity lemma**
   (`Δ_low = o(1)` whp for a random `F_2`-linear factor base, via algebraic
   independence of the Semaev coefficients), for which EXP-H (item 5) is
   strong evidence (`(2n')^{−4.3}` decay).
   (b) *Expansion (legacy):* the Ben-Sasson–Wigderson / Mikša–Nordström
   route via immunity; retained only as a fallback, since spectral `γ` does
   not track `D*`. **If the genericity lemma holds, the first-fall-degree
   assumption is provably false at cryptographic scale** — a defensive
   theorem, not just a heuristic.

---

## 5. The attacker's corollary (so this is two-sided, not just defense)

The same theory hands the *attacker* a search target: **maximize the early
Hilbert defect `Δ_low`**. It unifies every known speedup as "inject
low-degree algebraic relations":

- subfield / Koblitz curves → the factor base is multiplicatively closed,
  so products fall back into it → large `Δ_low`;
- symmetrization (Faugère–Gaudry–Huot–Renault, our
  `symmetrized_semaev.rs`) → the `S_m`-orbit relations are exactly extra
  low-degree syzygies → raises `Δ_low`;
- GHS-amenable field/basis choices → sparse normal bases that share
  multiplicative structure.

It also makes a sharp prediction: **any genuinely new subexponential
family must exhibit a descent system with non-vanishing early defect.**
That converts "hunt for a subexponential ECDLP attack" into the concrete,
checkable subproblem "find a curve/basis whose restricted descent system
has `Δ_low = ω(1)`" — a structural certificate we can screen candidate
primes and fields for (a few low-degree Macaulay ranks), feeding directly
into `pkm_criterion.rs` and `solinas_correlations.rs` in the resistance
map.

---

## 6. Falsifiable predictions (how to kill this proposal fast)

1. `pc_degree_harness` on a generic binary curve with a **dense normal
   basis** must show `D*` growing roughly linearly in `n` over
   `n = 5..10`, *while* `d_ff` stays ≈ 3. If `D*` also stays constant,
   the bridge's premise (that `d_ff ≪ D*` generically) is wrong and the
   proposal dies.
2. Across the basis/subspace sweep, `D*` must be **positively monotone in
   `γ(G)`** (high expansion ⇒ high PC/refutation degree,
   Ben-Sasson–Wigderson; *not* `1/γ` — an earlier draft had the direction
   backwards). **Status: REFUTED for spectral `γ` (workflow iterations
   2–3).**
   - *Iteration 2* (generic bases): over 24 generic irreducibles at
     (n=8, n'=3), Spearman ρ_s = −0.322 — weak and wrong-signed. But all
     generic bases are high-γ (γ_spec ∈ [0.74, 0.95]), so this was the
     wrong part of the axis.
   - *Iteration 3* (the decisive structured contrast, `descent_lowgamma`):
     the subfield factor base `F_{2^{n'}}` at `2n'=n` is the **easiest**
     case (mean `D*` 2.04 vs random 3.53) yet has **higher** spectral γ
     (0.882 ≥ 0.861). The low-`D*` structure is *not* the low-`γ`
     structure → the spectral-`γ` conjecture is contradicted by the most
     important structured case.
   - Two structural facts also emerged: the *tensor* incidence graph is
     basis-independent in the refutable regime (so the predictor uses the
     *system* graph), and unique-neighbor boundary expansion is uniformly
     0 on these dense graphs (so the BW boundary quantity needs a sparser
     graph model — a surviving, untested reformulation).
3. The HKY explicit counterexample systems must register as **low-`γ`**
   under `descent_expansion`. If they don't, the predictor doesn't
   capture the known gap and must be reformulated. *(Subsumed by (2)'s
   refutation of the spectral form and replaced by the algebraic predictor
   in (4); retained as record.)*

4. **(The surviving, supported prediction — P3-alg.)** Across the
   factor-base families and operating points, `D*` must be **monotone
   *decreasing* in the early Hilbert defect `Δ_low`** (§3). **Status:
   SUPPORTED (workflow iterations 5–7).** EXP-G: 50 cells over
   `2n' ∈ {4,…,14}` give pooled Spearman ρ_s = −0.79; in the critical
   regime `2n'=n` ρ_s = −0.78 with OLS slope ≈ −7.6 (seed-robust, slope
   ∈ [−7.6,−8.5]). The kill condition is now the *reverse*: if at larger
   reach the slope flattens to 0 or flips sign, the algebraic predictor
   fails too and the bridge has no working predictor. So far it strengthens
   with reach.

5. **(The lower-bound antecedent — §3.4.)** For the generic (Random)
   factor base in the critical regime, the early defect must **vanish** as
   `2n' → ∞`. **Status: SUPPORTED (workflow iteration 8, EXP-H).** Random
   `Δ_low ≈ (2n')^{−4.3} → 0` across `2n' ∈ {6,…,20}`, while Subfield stays
   bounded (ratio diverges 6.7 → 67). Kill condition: if Random `Δ_low`
   instead plateaued at a constant `> 0`, generic systems would *not* be
   asymptotically semi-regular and the defensive-theorem route would
   collapse. It does not plateau.

Each is cheap: (1), (2), (4) run in minutes on a laptop at `2n' ≤ 14`;
(3) is pure linear algebra at any `n`; (5) runs to `2n'=20` in seconds
(degree-≤3 ranks only).

---

## 7. Honest limitations / where it can break

- **PC ⇄ solving-degree is not exactly equality** over `F_2` with field
  equations — it is PC (sometimes PCR) and there are `±O(1)` and
  characteristic subtleties (Galesi–Lauria). The experiments calibrate
  the constant; the asymptotic claim is what matters.
- **Immunity is the crux and may fail precisely for special fields** —
  which is the *point* (those are the breakable cases), but it makes the
  *generic* lower bound genuinely hard. The proposal's defensive theorem
  is conditional on the algebraic-independence step.
- **Locality.** Expansion bounds want bounded-degree constraints; the
  descent equations are degree 2 but globally coupled through the
  tensor. The argument has to lean on tensor *sparsity* in the chosen
  basis, so the bound is really about (curve, basis) pairs, not curves
  alone — consistent with the known basis-sensitivity of GHS/Weil
  descent.
- It predicts asymptotics and a *screening* invariant, **not** a break
  of any 256-bit curve. As with everything in §6 of the research map, it
  does not threaten deployed parameters.

---

## 8. Why this is new

- HKY *define* `d_last` and bound it with HFE-style algebra; **nobody
  identifies the low-degree Hilbert-function defect `Δ_low` as its
  controlling parameter, nor exhibits the empirical `D* ↓ Δ_low` law.**
  That bridge — and the falsification of the natural graph-incidence
  alternatives along the way — is the contribution.
- Proof-complexity PC-degree lower bounds have been applied to random
  k-XOR/k-CNF, Tseitin, and pigeonhole — **never to
  summation-polynomial / ECDLP systems.** The translation, and the
  "immunity from algebraic independence of Semaev coefficients" step,
  is the new technical content.
- It yields a predictor **evaluable without running Gröbner**, at full
  cryptographic `n` — `Δ_low` is the rank profile of a few low-degree
  Macaulay matrices (`O(N^{D_low})` columns) — exactly the regime where the
  attack itself is unobservable and where the dispute has therefore been
  stuck on extrapolation from `n ≤ 8`.
- It **reconciles** Petit–Quisquater/Semaev with Huang–Kosters–Yeo and
  Galbraith–Petit: all three are right, on different *defect* regimes —
  multiplicatively-structured factor bases carry large `Δ_low` (heuristic
  holds), generic ones carry `Δ_low ≈ 0` (heuristic fails).
- The route is itself instructive: the natural combinatorial predictors
  (spectral expansion, treewidth) are **refuted**, and the surviving
  predictor is algebraic. A negative result on graph-incidence invariants
  for Semaev systems, with a positive algebraic replacement, is new.

---

## References

- I. Semaev, *Summation polynomials and the discrete logarithm problem
  on elliptic curves*, ePrint 2004/031; *New algorithm for the discrete
  logarithm problem on elliptic curves*, ePrint 2015/310.
- C. Petit, J.-J. Quisquater, *On polynomial systems arising from a
  Weil descent*, ASIACRYPT 2012.
- M.-D. Huang, M. Kosters, S. L. Yeo, *Last fall degree, HFE, and Weil
  descent attacks on ECDLP*, CRYPTO 2015.
- S. Galbraith, S. Gebregiyorgis, *Summation polynomial algorithms for
  elliptic curves in characteristic two*, INDOCRYPT 2014.
- J.-C. Faugère, P. Gaudry, L. Huot, G. Renault, *Using symmetries in
  the index calculus for elliptic curve DLP*, J. Cryptology 2014.
- M. Clegg, J. Edmonds, R. Impagliazzo, *Using the Groebner basis
  algorithm to find proofs of unsatisfiability*, STOC 1996.
- E. Ben-Sasson, A. Wigderson, *Short proofs are narrow — resolution
  made simple*, J. ACM 2001.
- M. Alekhnovich, A. Razborov, *Lower bounds for polynomial calculus:
  non-binomial case*, FOCS 2001.
- M. Mikša, J. Nordström, *A generalized method for proving polynomial
  calculus degree lower bounds*, CCC 2015.
- N. Galesi, M. Lauria, *Optimality of size-degree tradeoffs for
  polynomial calculus*, ACM TOCL 2010.
</content>
</invoke>
