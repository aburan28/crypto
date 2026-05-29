# A proof-complexity bridge for the first-fall-degree assumption

**Status:** research proposal (novel approach), 2026-05-29
**Builds on:** `src/cryptanalysis/ffd_harness.rs`,
`src/cryptanalysis/groebner_f4.rs`, `src/cryptanalysis/semaev_sat.rs`,
`src/cryptanalysis/binary_semaev.rs`, `src/cryptanalysis/diem_descent.rs`,
`RESEARCH_FFD_MEASUREMENT.md`.
**One-line thesis:** the disputed first-fall-degree heuristic is, up to
constants, a statement about the **Polynomial Calculus refutation
degree** of the Weil-descended summation-polynomial system — and that
degree is governed by the **expansion of the descent incidence graph**,
a combinatorial quantity we can compute at full cryptographic size
*without ever running a Gröbner basis*.

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

## 3. The predictor: expansion of the descent incidence graph

The reason PC degree is *low* for some Semaev systems and *high* for
others should not be mysterious — proof complexity tells us exactly what
controls it. For a bounded-degree polynomial system over `F_2`, the PC
refutation degree is lower-bounded by the **boundary expansion** of the
constraint–variable bipartite (Tanner) graph, *provided* the
constraints are "immune" (no small subset is satisfiable by a
low-degree assignment that a clever derivation could exploit).

Concretely, after the Weil descent in `ffd_harness::weil_descend_s3`,
each of the `n` output equations is a degree-2 `F_2`-polynomial in the
`2n` bit-variables, and *which* bit-variables appear in equation `j` is
dictated by the **structure-constant tensor** `c_{ikj}` of the chosen
`F_2`-basis of `F_{2^n}`:

```
   z^i · z^k  =  Σ_j c_{ikj} z^j        (the field-multiplication tensor).
```

Build the bipartite graph

```
   G(E, basis)  :  equations  ⟷  bit-variables,   edge iff variable
                   occurs in equation (with multiplicity from c_{ikj}).
```

> **⚠ STATUS (workflow iteration 3): the spectral form of this conjecture
> is REFUTED.** EXP-E (`descent_lowgamma.rs`) measured `(γ, D*)` for a
> genuine *subfield* factor base `F_{2^{n'}}` at the operating point
> `2n'=n` (n=8): the subfield is the **easiest** case (mean `D*` 2.04 vs
> random 3.53) — exactly the "structured ⇒ easy" effect index calculus
> exploits — **yet its spectral expansion is not lower** (0.882 ≥ 0.861).
> So the *most important* low-`D*` structure is **not** a low-`γ`
> structure, contradicting `D* = Θ(γ·n)` for the spectral `γ`. The
> conjecture is retained below as originally written for the record; see
> §6/§8 and `RESEARCH_FFD_WORKFLOW.md` iteration 3 for the refutation and
> the surviving (weaker) possibilities.

**Central conjecture (falsifiable — spectral form now falsified, see
above).**

```
   D*  =  Θ( min( m·n' ,  γ(G)·n ) ),
```

where `γ(G)` is the (small-set) boundary/spectral expansion of
`G(E, basis)`. Equivalently:

> the first-fall-degree assumption holds **iff** the descent incidence
> graph has expansion `o(1)` — i.e. iff the basis/field structure-
> constant tensor is sparse and clustered.

This **reconciles the dispute** instead of picking a side:

- The HKY counterexamples and the Galbraith–Gebregiyorgis "nice" cases
  (Koblitz / subfield bases, sparse normal bases) are precisely the
  **low-expansion** structures → `D*` small → the heuristic holds *for
  them*. That is *why* the first-fall assumption ever looked true.
- A **generic** curve over a generic basis (and every prime-field
  descent that lacks a subfield) gives a **high-expansion** tensor →
  `D* = Θ(n)` → the heuristic is **false**, and the attack is
  exponential.

Crucially, `γ(G)` is the spectral gap of an `n × 2n`-ish graph: it is
computable in milliseconds at `n = 256` even though the Gröbner basis it
predicts is astronomically out of reach. **That is the deliverable a
parameter-selection committee actually needs.**

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
   law needs (a) a sparse-F4 backend to reach `2n' ≳ 14`, and (b)
   separating the two parity classes explicitly. (iii) `D*` here is still
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

4. **Lower-bound attempt** — adapt Ben-Sasson–Wigderson /
   Mikša–Nordström expansion ⇒ PC-degree to `Σ`. The one genuinely hard
   step (and the honest risk) is **immunity**: Semaev equations are not
   random XORs, they carry the group law, so we must show that no small
   subset of descent equations is "prematurely refutable" at low degree.
   The proposed route: prove that the summation-polynomial coefficients
   are **algebraically independent** functions of the curve (a
   Schwartz–Zippel / generic-coordinates argument), so that under a
   random `F_2`-linear change of factor-base basis the system is, with
   high probability, an expander-supported quadratic system to which the
   expansion bound applies. **If this step holds for generic curves, the
   first-fall-degree assumption is provably false at cryptographic
   scale** — a defensive theorem, not just a heuristic.

---

## 5. The attacker's corollary (so this is two-sided, not just defense)

The same theory hands the *attacker* a search target: **minimize
`γ(G)`**. It unifies every known speedup as "an expansion reduction":

- subfield / Koblitz curves → block-structured, low-expansion tensor;
- symmetrization (Faugère–Gaudry–Huot–Renault, our
  `symmetrized_semaev.rs`) → collapses the `S_m`-orbit, shrinking the
  *effective* variable set and hence `γ·n`;
- GHS-amenable field/basis choices → sparse normal bases.

It also makes a sharp prediction: **any genuinely new subexponential
family must exhibit an explicit low-expansion descent.** That converts
"hunt for a subexponential ECDLP attack" into the concrete, checkable
subproblem "find a curve/basis whose multiplication-tensor Tanner graph
has `o(1)` small-set expansion" — a structural certificate we can screen
candidate primes and fields for, feeding directly into
`pkm_criterion.rs` and `solinas_correlations.rs` in the resistance map.

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
   capture the known gap and must be reformulated. *(Not yet run; given
   (2)'s refutation of the spectral form, this should target a
   reformulated `γ` — boundary expansion on the quadratic-only support, or
   treewidth — not spectral `γ`.)*

Each is cheap: (1) and (2) run in minutes on a laptop at `n ≤ 10`; (3)
is pure linear algebra at any `n`.

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
  bounds it via combinatorial expansion of the descent graph.** That
  bridge is the contribution.
- Proof-complexity PC-degree lower bounds have been applied to random
  k-XOR/k-CNF, Tseitin, and pigeonhole — **never to
  summation-polynomial / ECDLP systems.** The translation, and the
  "immunity from algebraic independence of Semaev coefficients" step,
  is the new technical content.
- It yields a predictor **evaluable without running Gröbner**, at full
  cryptographic `n` — exactly the regime where the attack itself is
  unobservable and where the dispute has therefore been stuck on
  extrapolation from `n ≤ 8`.
- It **reconciles** Petit–Quisquater/Semaev with Huang–Kosters–Yeo and
  Galbraith–Petit: all three are right, on different expansion regimes.

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
