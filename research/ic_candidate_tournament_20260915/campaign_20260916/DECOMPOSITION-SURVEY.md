# Decomposition survey: what could move the exponent against a matched rho

A scoping note for the IC-vs-rho tournament on binary Koblitz curves
`K_a : y² + xy = x³ + a x² + 1` over `F_{2^n}`, `n` prime (tournament `n = 13…61`,
cryptographic relevance `n ≥ 131`). Read-mostly: **no tournament ran, nothing was
promoted, no tracked file was changed.** The new files are this note and three small
model scripts, each with its output next to it (§8 lists exactly what ran).
Classification: **scoping / stage diagnostic**. It computes no `S` and no new rho ratio.

> **Addendum 2 (2026-10-02), appended at the end of this file:** §3.2's symmetric-group
> and Frobenius items and §3.3's three algebraic forms are now measured or closed by
> structure; pointers only.
>
> **Erratum 1 (2026-09-30), appended at the end of this file.** §0 item 2, §2 and §3.3
> say `m = 3` cannot beat rho even with a free oracle. That holds for the no-large-prime
> model only. With two large primes, the free-oracle exponent at `m = 3` is `4/9`. No
> `m ≥ 4` bar is loosened. Apart from this note, the sections below are unchanged.

The question: is there a way of decomposing `R = [a]G + [b]Q` into factor-base points
that could change the **cost exponent** against rho, rather than the constant in front
of it? Or is none known, and if so, why not?

## 0. Bottom line

1. **Nothing the tournament can currently enter can change the exponent. This is a
   theorem, not an engineering gap.** Both IC arms (the pair table, round 0023's
   `scaled`, and the triple table, `counted`) find relations with the group law, equality
   of encodings and the Frobenius. So they are generic algorithms. Their cost is
   `r^{j/(2j−1)}`, where `j` is the arity of the stored sum table: `r^{2/3}` for pairs,
   `r^{3/5}` for triples, and above `r^{1/2}` for every `j`. Shoup's `Ω(√r)` generic lower
   bound rules out anything below `1/2` for the whole family: any table, fold, walk or
   summand count. The measured losses to the matched rho (RESULTS.md round 0024: pair
   1.9–9.3×, triple 1.2–3.1× past `r ≈ 4·10⁶`) are this law showing up. Constants are
   all that is left there (§3.5).
2. **The only known route to an exponent below `1/2` uses the coordinate
   representation.** That means algebraic decomposition over an `F_2`-subspace factor
   base (Semaev summation polynomials, Weil descent, or Nagao's function-first form of
   the same variety) with **`m ≥ 4` summands**. With the per-trial solve cost written
   as `2^{c·n}`, the optimum exponent is
   `E(c, m) = 2(1+c)/(m+1)` for `c ≤ 1/m`, and `1/m + c` above that. It is below `1/2`
   only if **`c < 1/4` at `m = 4`**, `c < 0.30` at `m = 5`, and `c < 1/2 − 1/m` beyond.
   `m = 3` cannot do better than tie rho, even with a free oracle (§2).
3. **Every measurement on record points against `c` being that small.** None of them
   is a proof.
   - The `m = 3` refutation degree grows by about `ℓ/2` (5, 6, 6, 7, 7 over `ℓ = 3…7`).
   - The only `m = 4` per-target costs on record grow by **0.82–0.98 bits per unit
     `n`**, against a threshold of 0.25 (`decomp-m4-readout.txt`, read from
     `research/chain_split_order_20260924/tables.md`: two sizes, confounds listed).
   - The literature is a stalemate. Every subexponential claim rests on a
     first-fall-degree-type assumption that Huang–Kosters–Yeo turn into a reductio. No
     rigorous bound exists in either direction.
4. **The tournament range cannot host a win for any algebraic arm, even on the most
   optimistic published assumption.** The model (`decomp-exponent-model.txt` Part C)
   puts an `m ≥ 4` algebraic arm this far above the matched rho at `n = 61`:
   - **+34 to +41 bits** under Semaev's "degree ≤ 4" assumption. It first crosses rho
     at `n ≈ 190–260`.
   - **+114 bits** at `m = 4` under the Kousidis–Wiemers first-fall bound read as a
     solving degree. It crosses at `n ≈ 715–1750`.
   - Under the measured degree trend it **never crosses**: slope 1.0 per unit `n`,
     worse than enumeration.
5. **Top recommendation.** Do **not** open a tournament round for a new decomposition
   yet. Pre-register one **`m = 4` exponent audit** instead: a stage diagnostic that
   measures `ĉ`, the per-trial refutation-cost growth rate of the chained `m = 4`
   system at `ℓ ≈ n/4`, over `n = 9…19` on tournament-style curves. It has a
   random-system null and an enumeration null, and the decision threshold
   `c* = 0.25` is fixed now (§5). It reuses existing code and can be scoped to under
   an hour of CPU. Only if it comes back "alive" does an `m = 4` tournament arm become
   worth registering, and then as a **slope instrument**, not a promotion attempt (§6).
   A "closed" result closes `m = 4` for this solver family at these sizes. It does not
   close the route (§5.6).

## 1. Setting, unit, and numbered heuristics

**Matched rho.** The signed-Frobenius walk (Wiener–Zuccherato, SAC 1998;
Gallant–Lambert–Vanstone, CRYPTO 2000), with normal-basis orbit naming, as in round 0024
(`research/ic_triple_counted_20260923/rho-normal-basis.patch`). Its expected cost is
`√(πr/(4n))` group operations. Every exponent below is compared with rho's `1/2`. The
`√n` it gains is also available to every IC arm through the orbit fold, so it cancels
out of exponent comparisons.

**Heuristics** used throughout, so each claim can name the ones it rests on:

- **H1 (yield law).** A uniform target has on average `C(|F|, m)/#E` decompositions over
  an `m`-summand base `F`. Measured to within 1.09× on twelve toy cells
  (`research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md` §2). `yield/ceiling` is
  close to 1 on every base not confined to a subgroup
  (`research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md` §3.5).
- **H2 (relations).** About one independent relation per column is needed. Columns are
  base points up to sign, and up to Frobenius when the base is folded. Measured in every
  tournament round as `K` plus a small surplus.
- **H3 (linear algebra).** Sparse elimination costs about `m·K²`.
- **H4 (solve cost).** The per-trial decomposition oracle costs `T = 2^{c·n}` at the
  chosen `ℓ/n`. **This is the unknown on which the whole question turns.**
- **H5 (Macaulay model).** An algebraic oracle costs about
  `(number of monomials of degree ≤ D)^ω` with `ω = 2`, where `D` is the refutation
  (solving) degree. This is optimistic: it ignores the rows and the constant.
- **H6 (generic accounting).** Every table entry and every probe costs at least one
  group-operation equivalent. The tournament's instruction count only adds to this.

## 2. The exponent audit: two formulas that decide the ranking

`decomp_exponent_model.py` evaluates both. They are derived here, not measured.

**Generic `j`-table decomposition (H1, H2, H6).**
- Store all `j`-sums of an `F`-point base and loop `m − j` summands per trial.
- A trial succeeds with probability `F^m/r` and costs `F^{m−j}` probes. `F` relations
  are needed.
- Total cost: `F^j + r/F^{j−1}`, minimised at `F = r^{1/(2j−1)}`.
- **Result: `r^{j/(2j−1)}`, independent of `m`.** The IC/rho ratio then grows as
  `r^{1/(2(2j−1))}`:

| `j` | table | exponent | IC/rho grows as |
|--:|:--|--:|--:|
| 2 | pair (the tournament incumbent) | 0.667 | `r^{0.167}` |
| 3 | triple (`counted`) | 0.600 | `r^{0.100}` |
| 4 | 4-sum | 0.571 | `r^{0.071}` |
| ∞ | — | → 0.5 from above | → 1 |

The round-0021 two-point measurement read the pair arm at `r^{0.153}` against the model's
`r^{1/6}`. Running `research/ic_triple_table_20260923/model.py` gives the triple optimum's
local slope as 0.51 → 0.57 between `r = 2^26` and `2^42`, rising toward `3/5`.

**Algebraic subspace decomposition (H1–H4).**
- Base: points with `x ∈ V`, `dim V = ℓ`, `x = ℓ/n`.
- Relation collection costs `2^{ℓ}·2^{n−mℓ}·m!·T` while `mℓ ≤ n`.
- Linear algebra costs `2^{2ℓ}`.
- Minimising over `x` gives `E(c, m)`:

| `m` | `c*(m)`: largest `c` still below `1/2` | `E(0, m)` (free oracle) | `γ*`, chained S₃ (`N ≈ (m−1)n` unknowns) | `γ*`, direct `S_{m+1}` (`N ≈ n`) |
|--:|--:|--:|--:|--:|
| 3 | 0 (ties at best) | 0.500 | 0 | 0 |
| 4 | **0.250** | 0.400 | **0.083** | 0.250 |
| 5 | 0.300 | 0.333 | 0.075 | 0.300 |
| 6 | 0.333 | 0.286 | 0.067 | 0.333 |
| 8 | 0.375 | 0.222 | 0.054 | 0.375 |

`γ*` is the largest per-unknown exponent (`T = 2^{γN}`) a solver may have on the
descended system and still win asymptotically. For comparison:
- exhaustive search is `γ = 1`;
- the chained systems the repository actually builds need `γ < 0.09`, i.e. far below
  any generic Boolean solver;
- a constant or logarithmic degree would give that, and a degree linear in `ℓ` would
  not.

In the H5 model, a refutation degree `D = s·ℓ` keeps `c < c*` only for
`s ≲ 0.05` (`m = 4`) to `≈ 0.10` (`m = 6`), against the `s ≈ 0.5` measured at `m = 3`.

**The oracle budget at tournament sizes** (Part D).
- This is the most a per-target solve may cost with every other phase free except the
  linear algebra, maximised over `ℓ`, with the optimistic `n`-fold.
- It is **17–161 group operations at `m = 3`** and **10²–10⁶ at `m = 4, 5`** for
  `n = 23…61`.
- For comparison, the cheapest measured algebraic oracle prices are:
  - 8.8·10⁴ GAE per target (F4, `n = 23`, `m = 2`, 22 unknowns);
  - 2.4·10⁶ GAE (F4, `n = 15`, `m = 3`, 30 unknowns; `BOUNDARY_TARGETS.md` Regime B);
  - a pair-table probe costs 0.5–25.
- The rows of the budget table need systems of 47–253 unknowns.

## 3. The ranked shortlist

Ranked by how plausible an exponent change is and how cheaply that can be tested. Only
item 1 can change the exponent in principle. Items 2–4 are levers or reformulations of
item 1. Item 5 is the closed family, kept as the null object.

### 3.1 Algebraic `m ≥ 4` subspace decomposition (chained S₃ or symmetrised `S_{m+1}`), Frobenius-stable `V` where one exists

**Mechanism.**
- The factor base is `F = {P : x(P) ∈ V}`, where `V` is an `ℓ`-dimensional `F_2`-subspace.
  It is a union of Frobenius orbits when `V` is stable.
- A target decomposes iff the summation polynomial `S_{m+1}(x₁, …, x_m, x(R)) = 0` has a
  root with every `x_i ∈ V`. In practice this is the chain
  `S₃(x₁, x₂, e₁), S₃(e₁, x₃, e₂), …, S₃(e_{m−2}, x_m, x(R))`.
- Weil descent turns that into a Boolean system with `mℓ + (m−2)n` unknowns (chained)
  or `mℓ` unknowns (direct `S_{m+1}`), which F4, F5, Crossbred or SAT decides.
- Yield per trial is `min(1, 2^{mℓ}/(m!·2^n))` (H1). Per-trial cost is the solve, `T`
  (H4, H5).
- **Exponent: `E(c, m)`, below `1/2` iff `c < c*(m)`** (rests on H1–H5). This is the
  only family on this list whose cost law allows an exponent below rho at all.

**Known literature.** Repo-frozen sources are marked [F]; the others are cited by this
repository's notes and not re-read here.
- Semaev, ePrint 2004/031 (summation polynomials).
- Gaudry, J. Symbolic Comput. 2009, and Diem, Compositio Math. 147 (2011): the
  composite-degree subexponential results. They do not apply at prime `n` (§4).
- Faugère–Perret–Petit–Renault, EUROCRYPT 2012, and Petit–Quisquater, ASIACRYPT 2012:
  subexponential for `F_{2^n}` **under a first-fall-degree assumption**. The turning
  point against generic is quoted as `n ≈ 2000` (PQ) and `n ≈ 1250` (Kousidis–Wiemers,
  JMC 2019, whose `D_ff ≤ m² − m + 1` is a theorem, while its complexity reading is
  assumption-gated).
  - Source for these figures: `research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md`
    §2.
  - This model's `L_const` reproduces the same order independently: `n ≈ 715–1750`
    for `m = 4–6`.
- Semaev, ePrint 2015/310 [F: `crypto-autoresearcher/inputs/SEMAEV-2015-310/`]:
  `2^{c√(n ln n)}` under Assumption 1 (step degree ≤ 4 for chained S₃).
- Karabina, ePrint 2015/319 [F: `inputs/KARABINA-PDP-2015/`].
- Huang–Kosters–Yeo, CRYPTO 2015, ePrint 2015/573: a chain of `m = O(n)` S₃'s has
  `D_ff ≤ 5`, so the assumption would give a polynomial-time ECDLP. They call that
  "highly improbable".
  - The repository's literature note gives this paper's title inconsistently. I cite it
    only by venue and ePrint number. **Title uncertain.**
- Galbraith–Gaudry, DCC 2016: "no consensus". Caminata–Ceria–Gorla give a rigorous but
  numerically vacuous bound.
- Galbraith's 2015 summary: data only to `n ≤ 26` for `m ≥ 3`, and to `n ≈ 40–45` at
  `m = 2` [F: `inputs/ELLIPTICNEWS-CHAR2-2015/`, receipt only].

**What this repository measured.**
- m = 3 degrees.
  - The refutation degree reads 5, 6, 6, ≥7 over `ℓ = 2…5`
    (`research/dreg_ell_grid_20260925/RESULTS.md`), with a surplus confound shown by
    `(7, 4)` = 7 7 6 7.
  - The deciding degree with one summand fixed reads 5, 6, 6, 7, 7 over `ℓ = 3…7`, a
    slope of 0.50.
  - The first fall degree is 2–3 and says nothing about the solving degree
    (Kosters–Yeo's trace identity).
- The only m = 4 costs on record are the chain-split-order ladder: 0.82–0.98 bits per
  unit `n` (`decomp-m4-readout.txt`).
- End to end, every algebraic oracle loses to enumeration: by 3.4–1800× in
  `RESEARCH_ECC2K130_CROSSBRED.md` §4, and by 10³–10⁶× against rho at
  `n = 17–31` in `research/koblitz_symmetrised_e2e_20260927/RESULTS.md`.
- At `n = 131`, `RESEARCH_ECC2K130_DECOMPOSITION.md` §5.3 derives that the oracle must
  beat exhaustive search over its own candidates by `2^{70.19 + log₂ m}`.
- The campaign's own F4/F5 bench (`docs/ic/README.md`) gives m = 2 ≈ rho steps, and
  m = 3 at best a poly-factor tie. F5 prunes nothing at degree 3. This is consistent
  with `E(0, 3) = 1/2`.

**Cheapest pre-compute falsification (proof search map).**
- **Baseline reproduction.** Re-run two frozen cells (`K_0/2^9` m=4 and `K_1/2^15` m=4)
  through `groebner_stage_bench` and match `tables.md`'s word-operation totals exactly.
  The counts are deterministic.
- **Observation collision.** The observable is the per-target refutation cost. If random
  systems with the same unknown count, degree and term density
  (`koblitz_bench::random_control_system`, already the dreg ladder's control) show the
  same growth, the Semaev structure carries no weight.
- **Quantifier order.** The claim is "there exists a solver family such that, for all
  large `n`, `T ≤ 2^{cn}` with `c < 1/4`". Finite data can only falsify within the
  tested family and sizes. A closed verdict is scoped to that.
- **Method ceiling.** Even `c = 0` gives only `E = 0.4`. At `n = 61` the free-oracle
  margin is 2^9.4 (Part C `L_free`), so no tournament-size win is possible for any
  oracle costing more than ≈10⁵ group operations per target.
- **Nearby object.** Run the same ladder on a random-`b` binary curve of the same `n`.
  The route is not Koblitz-specific, so equal growth is expected. A difference would
  point at the orbit structure.
- Runtime: minutes for the baseline reproduction. The full audit is §5.

**Null-object control.** Two nulls: the `Enumerate` oracle (`c = (m−1)x = 0.75` at
`x = 1/4`, the product law) and the random-system control above.

**Crossover with the matched rho.** Model only:
- none under the measured-degree law;
- `n ≈ 190–260` under Semaev's Assumption 1;
- `n ≈ 715–1750` under the Kousidis–Wiemers bound read as a solving degree;
- none at tournament `n` under any law (+34 bits at best, `n = 61`).

**`sota_delta` / `dominated_by`.**
- `sota_delta`: none demonstrated at any size.
- `dominated_by`: the matched rho at every measured size, and the enumeration oracle
  end to end at every measured size (`n ≤ 31`).
- Honest status: the only live exponent route, open in the literature, with the measured
  trend against it.

### 3.2 Symmetry-enlarged descents (2-/4-torsion translations, Frobenius-orbit coordinates) as a lever on the degree slope of 3.1

**Mechanism.**
- Same base, equation and yield as 3.1.
- The system is rewritten in invariants of a larger group acting on the summands:
  - the symmetric group (Faugère–Gaudry–Huot–Renault, J. Cryptology 2014);
  - translation by 2- or 4-torsion (`K_0` has `#E(F_2) = 4`; Galbraith–Gebregiyorgis,
    INDOCRYPT 2014 [F: `inputs/GG-2014-806/`]);
  - Frobenius shifts (orbit-coordinate decomposition, thread 1 of
    `research/KOBLITZ_SUBFIELD_INDEX_CALCULUS_QUEUE.md`).
- This lowers the degree and the unknown count by constant factors.
- **It changes the exponent only if it changes the slope `s` of the degree in `ℓ`.**
  A constant cut in `D` leaves `E(c, m)` alone.
- Galbraith–Gaudry (2016, §9.2) call larger symmetry groups in the chain "an open
  problem". The repository notes that 2-torsion symmetrisation removes the linear trace
  equation. That is absence in ten sources read, not evidence of absence
  (`RESEARCH_ECC2K130_IC_LITERATURE.md` addendum).

**Measured here** (`research/koblitz_symmetrised_e2e_20260927/RESULTS.md`):
- `sym/x` reads 0.84, 0.77 and 0.76 at `m = 2` (`n = 17, 23, 23`), then **4.01 at
  `n = 31`**;
- `sym-m3` sits at 66,027× rho at `n = 17`;
- the same torsion-closed base under the combinatorial oracle reads `mitm-u/mitm-x` =
  0.89–1.28.

So no slope evidence exists yet, and the one trend is the wrong way.

**Cheapest falsification.**
- The `m = 3` refutation-degree ladder (the dreg-grid cells `ℓ = 2…5`) run on the
  symmetrised chain beside the unsymmetrised one. If both slopes read ≈ 0.5, the lever
  is a constant.
- This needs the symmetrised system wired into `ladder_measure`. I did not check how
  much code that is.

**Null-object control.** The existing `x` arm: the same unknown budget on a base not
closed under the translation.

**Crossover.** Inherits 3.1's.

**`dominated_by`.** The matched rho, and the unsymmetrised `x` arm at `n = 31`.

### 3.3 Nagao / Riemann–Roch function-first decomposition, extended to `m ≥ 4`

**Mechanism.**
- Instead of the summation polynomial, parametrise a function in a Riemann–Roch space
  vanishing at `−R`. Require its norm polynomial to split with every root in `V`: a
  remainder condition after division by the subspace polynomial.
- Weil descent is then applied to the function's coefficients (Nagao, ePrints 2013/548,
  2013/549, 2015/984 [F: `crypto-autoresearcher/inputs/NAGAO-*`]).
- It is the same variety as 3.1 with different coordinates, so it has the same yield
  (H1).

**Exponent.**
- Nagao's claims are subexponential (2013/549) and polynomial `O(n^{8w+1})` (2015/984).
  Both rest on a first-fall-degree assumption.
- In 2015/984, Proposition 5 is justified through the *fake* first fall degree plus the
  unproved Lemma 4. Source: `crypto-autoresearcher/coordination/review/sembin-20260916-propquant/report.md`,
  and `knowledge/open-problems/KN-OPEN-7f0511.md`, which lists four inequivalent "first
  fall degree" definitions in use.

**Measured here** (`research/notes/ecc2k130/RESEARCH_ECC2K130_RR_SOLVER_PANEL.md`, on
the real ECC2K-130 curve):
- The RR solvers beat the Semaev controls: 32/32 against 0/32 at `d = 6`.
- **Brute-force pair enumeration then exhausts the same instances in 0.494× the counted
  field operations.** The RR solver visits about `|F|²` candidates, the same order as
  the double loop.
- A linear NO-certificate exists only for `d ≤ 6`.
- The encoding's arity is fixed at 3.
- By §2, **no `m = 3` method can have an exponent below `1/2` even with a free oracle**,
  so an exponent claim needs this encoding extended to `m ≥ 4`.

**Cheapest falsification.**
- The autoresearcher `GOAL-SEMBIN-*` lane already measures the fake-vs-true
  first-fall-degree question on 2015/984's EQS4. Do not duplicate it.
- For the tournament, the check is §5's audit with an RR-encoded `m = 4` arm, once one
  exists.

**Null-object control.** Pair enumeration, already run.

**`dominated_by`.** Pair enumeration (0.494×), and the matched rho.

### 3.4 Frobenius-stable and quasi-subfield factor bases (a poly(`n`) lever and an availability table, not an exponent)

**Mechanism.**
- Take `V = ker g(σ)` for a divisor `g` of `tⁿ − 1`. Then `F` is a union of Frobenius
  orbits, and the relation count, columns and linear algebra all fold by `n`.
- In HKPY's framework (J. Math. Cryptol. 2020), quasi-subfield polynomials generalise
  this.

**Why it cannot change the exponent.**
- It is a factor of `n`. The matched rho already gets the same `√n` from the same
  Frobenius.
- In HKPY's own complexity (Remark 3.1), the barrier is the solving constant
  `κ ≈ 4.876`, not `deg λ`:
  - the optimum is `min_α (1 − α + κα²) = 0.9487`, whatever `deg λ` is;
  - Euler–Petit (FFA 2021) prove `β ≥ 3/4` for linearised polynomials, against a needed
    `β < 0.103`;
  - source: `research/notes/index-calculus/RESEARCH_QUASI_SUBFIELD.md` §5, §5c.

**Availability** (`decomp-frobenius-stable.txt`, exact):
- **None** of the tournament degrees 13, 19, 29, 37, 53, 59 or 61 has a non-trivial
  stable subspace. Nor do 131 or 163: `2` is a primitive root.
- 17, 23, 41, 43 and 47 have them only at dimension about `n/2` or `n/3`.
- **Only `n = 31`** (`ord = 5`: dimensions 5, 6, 10, 11, …) has one near the `m = 4`
  optimum `n/5 … n/4`.
- Among the NIST Koblitz degrees:
  - 233 has 29/30 and 58/59;
  - 571 has 114/115, which is exactly `n/5`, the `m = 4` balance point;
  - 283 has only 94/95; 409 only 204/205; 163 has none.

**Cheapest falsification.** None needed for the exponent claim: the fold is a factor by
construction. The yield of stable against random subspaces is already measured at about
1× the ceiling (ledger §3.4, §3.5).

**Null-object control.** A random subspace of the same dimension.

**`sota_delta`.** At most a factor `n`. **`dominated_by`:** the matched rho, which has
the same quotient.

### 3.5 Generic higher-arity table decomposition (the current tournament architecture): closed as to the exponent

**Mechanism.** §2's `j`-table family: the pair arm is `j = 2`, the triple arm `j = 3`.
Folding by Frobenius and sign, walked probes and base-size rules change constants only.

**Exponent.** `r^{j/(2j−1)} > r^{1/2}` for every `j`. The bound below `1/2` is Shoup,
EUROCRYPT 1997: the Frobenius is `[λ]`, which a generic algorithm simulates in
`O(log r)` operations. **Proved, not heuristic** (up to log factors).

**Measured.**
- The pair arm's lead over the matched rho ends at `r ≈ 4·10⁶` (`n23a1`, 1.126×).
  It wins by constant factors (0.61–0.82× in instructions) on six smaller cells.
- The triple arm is 1.22–3.14× rho in instructions at every cell. Natively it is just
  under 1 on three small cells.

**Falsification.** None needed. It stays in the survey as the **null object** for 3.1:
an algebraic arm has to beat this family's cost on the same base before its
exponent means anything.

**Crossover.** Measured, and in the direction the law predicts. `dominated_by`: the
matched rho for `r ≳ 4·10⁶` (pair) and everywhere in instructions (triple).

## 4. Excluded by structure (not ranked)

| route | why it does not apply at prime `n` on a curve defined over `F_2` | source |
|:--|:--|:--|
| GHS Weil descent to a hyperelliptic cover | the cover genus is `≈ 2^{n−1}` at prime `n` | `research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md`; literature note §4 |
| Joux–Vitse cover-and-decomposition | needs a tower `F_{q^d}/F_q/F_p` with composite degree | ePrint 2011/020; literature note §4 |
| Base change to `F_{2^{nk}}` with a subfield-`x` base | for odd `n`, a point with `x ∈ F_{2^k}` has `y ∈ F_{2^k} ∩ F_{2^{2k}}`, so it lies in the subgroup `E(F_{2^k})`. Its sums never reach the order-`r` part. On an isogenous curve that is not defined over the subfield, Diem's descent has extension degree `n`, and Diem's regime needs `n² ≲ k`, i.e. a base of `≥ 2^{n²}` points | derived here; literature note §4 ("collapses both towers") |
| Koblitz structure beyond Frobenius | `End(E) = Z[τ]`, class number 1, `Aut(E) = {±1}`: all of `⟨−1⟩ × ⟨τ⟩` is already used by rho | literature note addendum |
| Summand count or table arity inside the generic family | §3.5 | Shoup 1997 |

Transfer attacks (MOV/FR, anomalous) and Xedni-type lifting are not decomposition
strategies and are out of scope here. Nothing above is claimed about them.

## 5. Top recommendation: a pre-registrable `m = 4` exponent audit

A stage diagnostic, not a tournament round. What would be registered, before any cell
is run:

1. **Object.** The chained `m = 4` Semaev system:
   - `4ℓ + 2n` unknowns, `ℓ = round(n/4)`;
   - on `K_0` and `K_1` at prime `n ∈ {11, 13, 17, 19}`, plus `n = 9` and `15` as the
     frozen anchors;
   - uniform (natural) targets, 16 per cell;
   - the base's subspace drawn as in `dreg_ladder` (seeded per cell). At `n = 31`, add
     one stable-subspace cell as a side row.
2. **Engine.** The best in the tree: the inherited F4 with interleaved order and linear
   elimination (`RESEARCH_CHAIN_SPLIT_ORDER.md`). Plus the enumeration oracle as a
   second arm.
3. **Metric.** Median per-target word operations, deterministic.
   - `ĉ` is the least-squares slope of `log₂ T` against `n`, with a bootstrap-over-targets
     95% band.
   - Refuted-only and satisfiable-only medians are reported separately.
   - Word operations are not converted to group operations; the conversion does not move
     a slope.
4. **Decision** (fixed now, from §2):
   - **alive** if the band's upper end is below `c* = 0.25`;
   - **closed for this engine at these sizes** if its lower end is above `0.25`;
   - otherwise **inconclusive**.
   - Secondary check: the refutation degree `D(ℓ)` on the small cells, with `s* ≈ 0.05`.
5. **Controls.**
   - Baseline reproduction first: the `tables.md` totals for the two frozen cells must
     match exactly.
   - The random-system null at identical unknown count and degree must grow *faster*
     than the Semaev systems, or "structure" is not what is being measured.
   - The enumeration null gives `c = (m−1)x ≈ 0.75`.
6. **Stop rules and censoring.**
   - Hitting the node budget or a row/column cap is reported as censored. It is never
     treated as negative evidence.
   - A cell that censors on more than half its targets drops out of the fit, and that is
     recorded.
7. **What it cannot show.**
   - An asymptotic statement: the quantifier order in §3.1.
   - Anything about another solver family.
   - Anything at `n ≥ 131`. Any later ECC2K-130 transfer claim owes the `m = 83` gate of
     `AGENTS.md` §8a.
   - An "alive" verdict is a reason to scale, not a result.

**Successor or revisit conditions if the verdict is closed.**
- `m = 5`, where `c* = 0.30`, is the next rung, but the chain adds `n` unknowns per
  summand, which makes `γ*` smaller.
- Any engine that measures below 0.25 at `m = 4` reopens this.
- So does a published degree bound for Weil-descended chains, in either direction.
- So does 3.2's symmetry lever, if it lowers the slope `s` rather than `D` by a constant.

**Implementation.**
- The chained-`m` Gröbner decomposition already runs for `m = 4`:
  `koblitz_index_calculus::groebner_decompose`, `examples/groebner_stage_bench.rs`, and
  the chain-holdout ladders.
- What is missing is a registered `m4-exponent` ladder with `ℓ = round(n/4)` across the
  sizes above, plus the refuted/satisfiable split.
- The degree side needs `koblitz_bench::ladder_measure` generalised from its fixed
  `m = 3` chain.
- Cost estimate (advisory): the `n = 15` `m = 4` cells ran at about 4·10⁶ word operations
  per target. At the readout's roughly 2^{0.8} per unit `n`, `n = 19` is about 4·10⁷.
  That is seconds per target, so the full audit is under an hour of CPU.

## 6. What it would take as a tournament arm

- **Worker.**
  - `examples/ic_tournament_worker.rs` already maps `solver` ∈ {`f4`, `f5`,
    `inherited_f4`, `sat_xor`, `sat_cnf`, `enumerate`} to `KoblitzIcOptions` with
    `m = summands`.
  - The single-word `koblitz_tiny_ic` path covers only `pair_table` with three summands,
    so an algebraic arm runs on the general path. A registry entry such as
    `{"solver": "inherited_f4", "summands": 4, …}` plus a subspace base recipe is the
    whole configuration.
  - `tournament.py` still defaults to `BASE_CONFIG['solver'] = 'pair_table'`.
- **Checker (`oracle.py`).**
  - `verify(..., summands=m)` takes `m` as a parameter. It re-adds every relation's base
    points in the group and checks each descent `[a]G + [b]Q = Σ P_i`. It does not care
    how a relation was found, so algebraic provenance needs **no** checker change.
  - Three things do bind:
    - every relation and every descent must have **exactly** `summands` points. An arm
      that also finds shorter decompositions must drop them. Admitting `≤ m` would be a
      pre-registered checker amendment.
      - A concurrent, uncommitted lane already drafts one:
        `round25-oracle-mixed-summands.patch` in this directory, which I did not
        author or review. It lets a report declare any count `k` from 3 up to the
        configured `m`, but every relation in that report must then have exactly `k`
        points. So an algebraic arm still cannot mix lengths within one report.
    - `factor_base_orbits` requires Frobenius closure, so a non-stable subspace base must
      list `factor_base` explicitly: about 2^{15} points at `n = 61`, fine.
    - the round-0024 note records that main's harness admits only three-summand
      pair-table sources, so a 4-summand arm needs round 0023's frozen evaluator or an
      admission amendment.
- **What it would show.**
  - It cannot win at `n ≤ 61`. The margins are the §0 item 4 numbers, and measured
    algebraic arms sit 10³–10⁶× rho at `n ≤ 31`.
  - So register it as a **slope instrument**: its IC/rho growth in `r` compared with
    `E(ĉ, 4) − 1/2`, and only after §5 says alive.
  - It must not be registered as a promotion attempt, and its loss is not negative
    evidence about the route.

## 7. Scope and honesty

- **Measured versus modelled.**
  - Everything in §2 and every crossover `n` comes from a closed-form model with the
    stated heuristics (H1–H6) and `ω = 2`. Constants below 2× and unit conversions are
    dropped, which moves no exponent.
  - The measured numbers quoted are other lanes' frozen results, cited by path.
  - The only numbers this note computed are listed in §8.
- **Premature closure is avoided on purpose.**
  - §3.1 is **not** closed. The literature is a stalemate. The repository's
    degree-growth evidence is at `m = 3` and small `ℓ`, and its `m = 4` evidence is two
    sizes with confounds.
  - The model's `L_meas` law extrapolates the `m = 3` slope to `m ≥ 4`. That is
    optimistic if the `m = 4` degrees are higher, and it is unmeasured.
- **Overclaiming is avoided too.**
  - No curve is weakened. No crossover is measured for any algebraic arm.
  - The only proved statement is §3.5's generic lower bound, and it is a bound on a
    family, not on the problem.
- **Citations.**
  - Sources marked [F] are frozen in `crypto-autoresearcher/inputs/`.
  - Shoup 1997, Wiener–Zuccherato 1998 and GLV 2000 are standard and cited from
    memory.
  - All others are cited as the repository's notes cite them and were not re-read here.
  - The HKY CRYPTO 2015 title is marked uncertain.
  - `research/notes/index-calculus/RESEARCH_FFD_MEASUREMENT.md` attributes the CRYPTO
    2015 paper to "Huang–Kiltz–Petit". This looks like an error; the author list elsewhere
    in the repository is Huang–Kosters–Yeo. It is recorded here, not corrected.

## 8. What was run

All runs were on this container, single-threaded, about 5 s in total. No build and no
cargo; no tracked file modified.

| command | output | runtime |
|:--|:--|:--|
| `python3 decomp_frobenius_stable.py` | `decomp-frobenius-stable.txt`: `ord_n(2)` and stable dimensions for the tournament and cryptographic degrees | < 1 s |
| `python3 decomp_exponent_model.py` | `decomp-exponent-model.txt`: Parts A–D (generic family, `c*(m)`, four solve laws against matched rho, oracle budgets) | 0.2 s |
| `python3 decomp_m4_readout.py` | `decomp-m4-readout.txt`: `m = 4` per-target growth read from `research/chain_split_order_20260924/tables.md` | < 1 s |
| `python3 research/ic_triple_table_20260923/model.py` (existing, unmodified; output to terminal only) | the triple optimum's local slope, 0.51 → 0.57 over `r = 2^26…2^42` | < 5 s |

## Erratum 1 (2026-09-30): double large primes lower the `m = 3` ceiling

§0 item 2, the `m = 3` row of §2's table, §3.3 ("no `m = 3` method can have an exponent
below `1/2` even with a free oracle") and `decomp_exponent_model.py` line 93 all state the
same `m = 3` ceiling. Nothing above is rewritten.

**What was wrong.** The claim holds for §2's model, which has no large primes:
`E(c, m) = 2(1+c)/(m+1)` gives `E(0, 3) = 1/2`. It is not a ceiling for every method. The
double-large-prime variant works on a small base `F′ ⊂ F` and recombines relations carrying
up to two large primes through the collision graph. It costs `q^{2−2/m}` in place of
`q^{2−2/(m+1)}`.

**The repository already has this result, in its own setting.**
[`RESEARCH_EXTENSION_FIELD_BOUNDARIES.md`](../../notes/index-calculus/RESEARCH_EXTENSION_FIELD_BOUNDARIES.md)
Theorem 3 derives total work `Θ(N^{4/9})` at `k = 3`. It rests on H1 and on H2, the
large-prime percolation heuristic.
[`RESEARCH_RHO_PARITY_PROGRAMME.md`](../../notes/index-calculus/RESEARCH_RHO_PARITY_PROGRAMME.md)
measures it end to end, for the relation phase, at `n^{0.44 ± 0.01}` on `E(F_{p³})` with a
subfield base. So this survey contradicted the repository's own theorem and measurement.

**The corrected exponents.** In the unit `r = q^m`, the double-large-prime free-oracle
exponent is `2(m−1)/m²`. With per-attempt oracle cost `2^{cn}`, not re-optimised,
`E_dlp(c, m) = 2(m−1)/m² + c` ([`decomp_dlp_erratum.py`](decomp_dlp_erratum.py), output in
[`decomp-dlp-erratum.txt`](decomp-dlp-erratum.txt)):

| `m` | `E(0, m)`, §2 | `E_dlp(0, m)` | `c*`, §2 | `c*`, double large primes | best bar |
|--:|--:|--:|--:|--:|--:|
| 3 | 0.500 | **0.444** | 0 | **0.056** (`1/18`) | 0.056 |
| 4 | 0.400 | 0.375 | 0.250 | 0.125 | 0.250 (unchanged) |
| 5 | 0.333 | 0.320 | 0.300 | 0.180 | 0.300 (unchanged) |
| 6 | 0.286 | 0.278 | 0.333 | 0.222 | 0.333 (unchanged) |

**What changes.**
- **`m = 3` has an exponent bar, but a tight one:** `c < 1/18`, so the oracle must be
  nearly free. Every measured `m = 3` algebraic oracle is far above it. Examples: the
  refutation-degree slope of about 1 per unit `ℓ` in `ic_symmetry_lever_slope_20260929`, and
  the F4 prices in §2's budget paragraph. No `m = 3` exponent audit against `c < 1/18` has
  been run.
- **No `m ≥ 4` bar is loosened.** At `m ≥ 4` the double-large-prime bar is *stricter* than
  §2's `c*(m)`, because §2's model already rebalances the base against `c`. So every `m = 4`
  verdict (`ic_m4_exponent_audit_20260928`, `ic_m4_head_engine_20260929`), read against
  `c* = 0.25`, stands.
- **§3.3's inference is weakened.** It argued that "an exponent claim needs this encoding
  extended to `m ≥ 4`". An `m = 3` encoding with a nearly free oracle plus double large
  primes would also qualify.

**What this rests on.**
- **Heuristics.** Carrying Theorem 3 to a subspace base on binary Koblitz curves assumes H1
  (the decomposition rate) and H2 (percolation of the large-prime graph) hold there. That
  has not been validated on binary curves.
- **Memory.** The graph holds `Θ(N^{2/9})` large-prime vertices. It is not charged here,
  and neither is it in §2.

## Addendum 2 (2026-10-02): §3.2 and §3.3 measured

Pointers to records written after this survey; the sections above are unchanged.

- **§3.2, symmetric-group action:** measured by the `m = 3` Riemann–Roch norm-form ladder,
  which is fully symmetric in the summands:
  [ic_rr_norm_ladder_20260930](../../ic_rr_norm_ladder_20260930/RESULTS.md), `s̄ = 1.167`,
  constant lever.
- **§3.2, Frobenius-orbit coordinates:** closed by structure,
  [NOTE-20260930-frobenius-orbit-coordinates.md](NOTE-20260930-frobenius-orbit-coordinates.md);
  the orbit-union reading is the coset-typed base of `H-SEMBIN-c59e50`.
- **§3.3, extended to `m ≥ 4`:** the search form's ceiling is §3.5's; the support form's
  system degree rises 2 per unit `ℓ`
  ([support-degree.txt](../../ic_rr_norm_ladder_20260930/support-degree.txt)); the norm form
  with `V`-typed roots is a constant lever (above). The `m = 4` norm form (`4ℓ + n + 1`
  unknowns, cubic) is an engineering arm for §5's audit, not an exponent candidate. The
  first-fall-degree claim stays with the SEMBIN lane, as this section already says.

