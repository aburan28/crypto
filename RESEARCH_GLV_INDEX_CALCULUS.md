# GLV / `C₃` levers on the `E(F_{p³})` index-calculus harness

**Module:** `src/cryptanalysis/glv_gaudry.rs` (engine changes in `f4_fp.rs`: field-op counter, weight-block reduction; `gaudry_cubic.rs`: `Curve3::new`, `SymmetrisedS4::terms`)
**Bench:**  `cargo run --release --example glv_gaudry_bench -- --exp {quotient,canonical,graded,invariant} --sizes 271,541,1051[,2113] --seeds 2 --json experiments/22_glv_<exp>.json`
**Data:**   `experiments/22_glv_{quotient,canonical,graded,invariant}.{json,log}`
**Tables:** `python3 scripts/glv_gaudry_tables.py experiments/22_glv_*.json` (every number below is printed by it from the frozen files)
**Setting:** §11 of `RESEARCH_RESIDUAL_WALKS.md` — Gaudry's subspace base `{P : x(P) ∈ F_p}` on `E(F_{p³})` with the `O(1)` symmetrised-`S₄` solve — on `j = 0` curves, where the order-3 automorphism `ψ(x, y) = (ωx, y)`, `ω ∈ F_p`, preserves the base.

> **Result in one line.**  Of the four GLV experiments, one moves `S`:
> quotienting the factor base by `⟨ψ⟩` divides the columns, the relations,
> the residuals, the solver calls, the non-zeros and the total cost by
> `3.0` at every size from `2^24` to `2^33` — with the count sitting on
> its (fold-aware) floor before and after, so this is the automorphism
> group's `3`, *engineering*, not an advance, and the method stays
> `270–700×` behind rho.  Canonical generation saves no solver call on
> the harness's own residual stream (the collision floor `3T²/n` is below
> `0.02` at every size, and `0` were found) but is what makes a
> sieve-style pair generator usable on a `j = 0` curve at all: `95 %` of
> its decompositions are `ψ`-conjugates of the pair itself and fold to
> `0 = 0`.  The `Z/3` grading exists — `S₄` is `ψ`-invariant term by term
> — but only on the *orbit* system, which has three times the solutions,
> `50–70×` the F4 multiplications, and does not close under a single-degree
> Macaulay truncation by degree `19`; block-aware elimination cuts its
> dense footprint by `3×` and its multiplications by `1.4 %`.  The
> `C₃`-invariant formulation of the per-residual PDP is the ordinary
> Semaev system term for term, and the function-first (Nagao `L(4O)`)
> formulation is `> 30×` larger in every matrix dimension and unsolved at
> a budget where the ordinary system solves in seconds.

## 1. Setting and what `ψ` does to the harness

`E : y² = x³ + b` over `F_{p³}`, `p ≡ 1 (mod 3)`, prime order `n`, so
`n ≡ 1 (mod 3)` and `ψ = [λ]` with `λ² + λ + 1 ≡ 0 (mod n)`
(`generate_j0_instance3`; the six sextic-twist orders are fixed by `p`,
and `p ∈ {73, 271, 541, 1051, 2113}` are primes with a prime twist,
`n ≈ 2^{18.6}, 2^{24.2}, 2^{27.2}, 2^{30.1}, 2^{33.1}`).  Because
`ω ∈ F_p`, `ψ` maps the subspace base `F` to itself: `x ↦ ωx` keeps
`x ∈ F_p` and `x³ + b`, so `F` is a union of `|F|/3` orbits
`{P, ψP, ψ²P}` (no fixed point: `x = 0` is not on a prime-order curve,
`b` being a non-square).  Every base point is `ψ^k` of its orbit
representative, i.e. `[λ^k]` times it, and every base point already
carries the `⟨−1⟩` fold (canonical `y`), so the control is the
`⟨−1⟩`-folded base of §11 and the experiment is the `⟨−1, ψ⟩` quotient.

The one algebraic fact everything else rests on is checked at run time:
the symmetrised `S₄(e₁, e₂, e₃, x₄)` of a `j = 0` curve has all `49`
terms of `ψ`-weight `a + 2b + d ≡ 0 (mod 3)` (`s4_weight_histogram =
[49, 0, 0]` on every instance), i.e. `S₄(ζx₁, ζx₂, ζx₃, ζx₄) = S₄(x)`.
The module's S₃ identity in `ec_index_calculus_j0` (`S₃(ζx) = ζ·S₃(x)`)
is the weight-`1` version of the same statement.

## 2. Boundaries, stated before measuring

- **Relation count (experiment 1).**  A square system needs at least
  `columns` relations, hence `columns / ρ` residuals at decomposition
  rate `ρ`; both stores see the same `ρ`.  The quotient divides
  `columns` by `3`, so the floor divides by `3` with it.  By §3 of
  `AGENTS.md`, an *advance* would be `residuals / floor < 0.9` on the
  quotient; `≈ 1.0` on both is the count moving with its bound.
- **Reference.**  Plain rho on the same group (`run_rho3`, measured per
  instance, two seeds).  A rho that folds by `⟨−1, ψ⟩` gains the same
  `√3` (`RESEARCH_RESIDUAL_WALKS.md` §10.3 measured it on prime-field
  `j = 0`); that column is an extrapolation and is marked so.
- **Orbit duplicates (experiment 2).**  `T` uniform residuals collide in
  a `⟨−1, ψ⟩`-orbit `C(T, 2)·6/n` times in expectation; the canonical
  generator cannot save more solver calls than that on a uniform stream.
- **The orbit system (experiments 3, 4).**  Any `C₃`-invariant
  formulation of the PDP solves the *orbit* `{R, ψR, ψ²R}`, whose three
  relations are `λ`-multiples of one another (the dependence the
  `ec_index_calculus_j0` module documents).  Its quotient dimension is
  `3 × 64 = 192` against `64`, so a block-aware graded solve can at best
  tie the ordinary per-residual solve, never beat it: the boundary is
  `1.0×` the ordinary solve's field multiplications.
- **Falsification targets.**  Experiment 1: an advance needs
  `residuals / floor < 0.9`.  Experiment 2: solver calls saved above
  `3T²/n` on the uniform stream.  Experiment 3: orbit-block F4 below
  `1.0×` the ordinary solve's multiplications.  Experiment 4: an
  invariant or function-first system whose Macaulay solve closes at
  `D < 10` or with fewer than `286` columns.  Inadmissible: changing the
  base, the decomposition size, the operation accounting, or dropping
  a phase; a run whose answer is not the planted logarithm.

## 3. Experiment 1 — the `⟨ψ⟩` factor-base quotient

One residual stream per instance feeds two relation stores.  Each
residual `R = aG + bQ` goes through the harness's symmetrised-`S₄`
solve once; every verified decomposition `R = Σ c_i P_i` is inserted
into the control as `(i, c_i)` and into the quotient as
`(orbit(i), c_i λ^{k(i)})`, same-orbit terms merged, rows that fold to
`0 = 0` or to a unit multiple of a stored row dropped (none occurred on
this stream).  A store is solved by singleton filtering to a square
core and sequential Wiedemann (`square_core`, `wiedemann_u64`) as soon as
its own core determines `d`, and the residual count at that moment is
what a stand-alone run on the same stream would have needed; the
solver cost per residual is identical by construction.  Costs are the
harness's: group operations, `F_p` multiplications at the measured
`63` per addition, `9` per multiplication modulo `n` in the linear
algebra.

| p | n | base | orbits | rate ρ | control S | quotient S | ratio | control res/floor | quotient res/floor | rho S (plain, measured) | quotient S / rho S | folded rho S (÷√3, extrapolation) | quotient S / folded rho |
|---:|:--|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 271 | 2^24.2 | 140 | 46 | 0.201 | 2,209.9 | 746.5 | 2.96 | 1.00 | 1.00 | 1.79 | 418× | 1.03 | 724× |
| 541 | 2^27.2 | 268 | 90 | 0.165 | 1,854.0 | 647.2 | 2.86 | 1.00 | 1.04 | 1.09 | 592× | 0.63 | 1,025× |
| 1051 | 2^30.1 | 556 | 186 | 0.200 | 1,190.3 | 397.6 | 2.99 | 1.00 | 1.01 | 0.56 | 714× | 0.32 | 1,237× |
| 2113 | 2^33.1 | 1034 | 344 | 0.159 | 991.5 | 329.6 | 3.01 | 1.00 | 1.00 | 1.22 | 271× | 0.70 | 469× |

Means over two seeds; every run recovered the planted logarithm.  The
per-seed table (columns, relations, residuals, non-zeros, bytes, the
square core's rows and non-zeros, Wiedemann time, mod-`n`
multiplications) is the first table of `scripts/glv_gaudry_tables.py`;
at `2^33.1`, seed 1:

| store | columns | relations | residuals | NNZ | bytes | core rows × NNZ | Wiedemann ms | mod-n mults | total ops | S |
|:--|---:|---:|---:|---:|---:|:--|---:|---:|---:|---:|
| control `⟨−1⟩` | 1000 | 1000 | 6806 | 3994 | 95,904 | 727 × 2904 | 71.7 | 8,979,901 | 100,072,225 | 1,030.3 |
| quotient `⟨−1, ψ⟩` | 334 | 334 | 2343 | 1333 | 32,016 | 221 × 883 | 6.9 | 830,739 | 34,059,481 | 350.7 |

What moved and what did not:

- Columns, relations, residuals and solver calls all fall by `3.0×`
  (`2.86–3.01` in `S`), non-zeros and store bytes by `3.0×`, the
  Wiedemann core by `3.3×` in rows and non-zeros, its time and
  multiplications by `≈ 10×` (quadratic in the core).  The linear
  algebra is `0.4 %` of the control's `S` and `0.2 %` of the quotient's,
  so its `10×` is invisible in the total; the `3×` in `S` is the solver
  running on a third of the residuals.
- Peak process RSS is `5–41 MB` for both stores together and is not the
  metric; the exact store sizes are.
- `residuals / floor` is `1.00` on the control and `0.98–1.09` on the
  quotient at every size: the count moved *with* its floor.  **Class:
  engineering**, exactly as the prime-field sixfold fold of §10.3 was.
- Against rho the quotient is `271–714×` behind on the measured plain
  rho and `469–1,237×` behind a folded rho (extrapolated `÷√3`), where
  the tuned `O(1)` solve of §11.6 was `528×` behind at `2^33`.  The two
  seeds' rho `S` range from `0.16` to `2.62`, so the rho column is noisy
  at two seeds; the control-to-quotient ratio is not.

The `3` is the size of the automorphism group acting on the base.  It
cannot be tuned upward, it applies to the relation phase only through
the number of unknowns, and rho takes the same `√3`.

## 4. Experiment 2 — canonical relation generation

`GlvInstance3::canonical` reduces a residual to the lexicographically
smallest of `{(ω^k x, ±y)}` and returns the unit `u ∈ {±λ^k}` with
`R = u·C`; the generator solves `C` and scales the rows by `u`.  A
residual whose canonical key has been seen is a saved solver call (up
to `100` of them per stream are solved anyway and their rows checked to
be unit multiples of the stored ones).  Decompositions are folded to
orbit columns, and a row is dropped when it folds to `0 = 0` or to a
unit multiple of a stored row.  Two streams per instance: `T` uniform
residuals `aG + bQ` (`T ≈ p`, about what the quotient pipeline needs)
and the first `3,000` base pairs `P_i + P_j`, `i < j`, as a sieve-style
generator.

| p | n | stream | residuals | distinct orbits | duplicates | expected (uniform) | solver calls saved | verified / mismatches | rows produced | rows zero (ψ-trivial) | rows unit-duplicate | informative rows |
|---:|:--|:--|---:|---:|---:|---:|---:|:--|---:|---:|---:|---:|
| 271 | 2^24.2 | uniform | 270 | 270 | 0 | 0.0110 | 0 | 0 / 0 | 56 | 0 | 0 | 56 |
| 271 | 2^24.2 | pairs | 3000 | 2163 | 837 | 1.3567 | 837 | 100 / 0 | 6812 | 6468 | 245 | 99 |
| 541 | 2^27.2 | uniform | 540 | 540 | 0 | 0.0055 | 0 | 0 / 0 | 85 | 0 | 0 | 85 |
| 541 | 2^27.2 | pairs | 3000 | 2870 | 130 | 0.1704 | 130 | 100 / 0 | 9021 | 8598 | 252 | 171 |
| 1051 | 2^30.1 | uniform | 1050 | 1050 | 0 | 0.0028 | 0 | 0 / 0 | 219 | 0 | 0 | 219 |
| 1051 | 2^30.1 | pairs | 3000 | 2972 | 28 | 0.0232 | 28 | 28 / 0 | 9459 | 8911 | 285 | 263 |
| 2113 | 2^33.1 | uniform | 2112 | 2112 | 0 | 0.0014 | 0 | 0 / 0 | 302 | 0 | 0 | 302 |
| 2113 | 2^33.1 | pairs | 3000 | 2991 | 9 | 0.0029 | 9 | 9 / 0 | 9395 | 8970 | 234 | 191 |

Seed 1 of two; seed 2 is within a few percent on every column (full
table from the script).

- **On the harness's own stream the canonical generator saves nothing.**
  Zero orbit duplicates in `270–2,112` residuals at every size, against
  an expectation of `0.001–0.011`; zero rows dropped.  This is the
  boundary, not a shortfall: uniform residuals do not repeat orbits at
  these counts, and the `3×` of experiment 1 is already the whole
  content of the fold.  **Class: accounting** — the mechanism is
  correct (`0` mismatches in `337` verified duplicates across the pair
  streams) and changes no cost.
- **On a pair generator it is decisive, and for a reason specific to
  `j = 0`.**  Every point satisfies `P + ψP + ψ²P = O`, so the `S₄`
  solve finds, for `R = P_i + P_j`, "decompositions" such as
  `R = −ψP_i − ψ²P_i + P_j` — the pair itself in `ψ`-conjugates — and in
  the quotient these fold to `0 = 0`: `95 %` of all rows
  (`6,468 / 6,812` … `8,970 / 9,395`), with another `3 %` unit multiples
  of rows already stored.  Only `1.5–2.8 %` of what the naive generator
  would insert carries information.  The orbit duplicates among the
  residuals themselves (`837` of `3,000` at `p = 271`, falling to `9` at
  `p = 2,113` because a `3,000`-pair prefix covers fewer whole
  `ψ`-orbits of a larger base; the structural saving over a full
  enumeration is `2/3`) are the "solver calls saved" the experiment
  asked for.  Without canonicalisation a pair-based sieve on a `j = 0`
  curve inserts a store that is `97 %` zero and dependent rows.

## 5. Experiment 3 — the `Z/3`-graded Macaulay matrix

Two instruments on the same residuals (those the harness decomposes:
`11–23` undecomposable residuals per size were skipped).  The Macaulay
instrument builds, at one degree `D`, every shift `m·f_i` with
`deg m ≤ D − deg f_i` over the grevlex-ordered monomials of degree
`≤ D`, reduces it with the harness's sparse-aware elimination and its
multiplication count, and calls the system solved at the first `D`
where every standard monomial times every variable is a pivot or
standard (the multiplication matrices exist).  Block-aware means the
columns are partitioned by weight `Σ wᵢeᵢ (mod 3)` and each block is
reduced on its own, which is legitimate only when every row lies in one
block.  The second instrument is the repository's F4 (`f4_fp`) with the
same block option added at every step, plus its substitution solver, so
solution counts are checked.

**The ordinary system has no grading to exploit.**  With `x_R` fixed,
the three Weil components of `S₄(e; x_R)` mix all three weight classes:
`100 %` of the Macaulay rows are coupled on every residual.  The system
solves at `D = 10` with a `252 × 286` matrix, quotient dimension `64`,
`6.1·10⁵` multiplications (the harness's own solve, with its row
selection and the eigenvalue step, spends `0.88–1.03·10⁶`).

**The grading exists on the orbit system.**  Making `X = z·x_R` an
unknown with `z³ = 1` (and `z^d` reduced modulo `z³ − 1`) gives four
equations in `(e₁, e₂, e₃, z)` of weights `(1, 2, 0, 1)` that are
weight-homogeneous — the computational content of `S₄`'s
`ψ`-invariance — so the Macaulay matrix splits into three blocks of
equal size and F4 blocks `14` of its steps.

| system | coupled rows | Macaulay `D_solve` | quotient dim | peak matrix | dense cells | multiplications |
|:--|---:|:--|:--|:--|---:|---:|
| ordinary `(e₁, e₂, e₃)` | 1.00 | 10 | 64 | 252 × 286 | 72,072 | 6.1·10⁵ |
| orbit, plain | 0.00 | none ≤ 19 | 576 and rising | 11,985 × 8,855 | 106,127,175 | 4.2·10⁹ summed to 19 |
| orbit, block-aware | 0.00 | none ≤ 19 | 576 and rising | 3 blocks ≤ 4,002 × 2,954 | 35,375,823 | 4.2·10⁹ summed to 19 |

Identical on all twelve residuals (four per size).  The orbit system
has `192` affine solutions but Bézout degree `6·6·6·3 = 648`: `456`
solutions at infinity, and the single-degree affine truncation keeps
finding new standard monomials (`348` at `D = 10`, `576` at `D = 19`)
instead of settling at `192`.  Its degree-by-degree F4 does solve it:

| system | F4 `D_solve` | max matrix | field mults (top run + p substitution runs) | solutions | blocked steps | mults vs ordinary | mults block / plain | ms block / plain |
|:--|---:|:--|---:|---:|---:|---:|---:|---:|
| ordinary | 10 | 221 × 277 | 6.2–13.3·10⁶ | 1–3 | 0 | 1.0× | — | — |
| orbit, plain | 10 | 1,668 × 1,632 | 4.3–4.6·10⁸ | 3× ordinary | 0 | 34–70× | — | — |
| orbit, block-aware | 10 | 1,668 × 1,632 | 4.3–4.5·10⁸ | 3× ordinary | 14 | 34–69× | 0.984–0.988 | 0.50–1.10 |

Ranges over the twelve residuals; the "vs ordinary" ratio falls with
`p` only because the substitution solver's `p` sub-runs are counted on
both sides.  `orbit = 3 × ordinary` solutions held on every residual.

- **Block-aware elimination is a footprint lever, not a field-operation
  lever.**  Sparse-aware elimination never multiplies across blocks, so
  the counts agree to `1.4 %`; what the split buys is `3.0×` fewer dense
  cells (Macaulay) and `0.5–0.8×` the wall time in F4 (one residual at
  `p = 1051` ran `1.10×`).  **Class: engineering** on the orbit system.
- **Against the boundary it is far below parity.**  The orbit system
  costs `34–70×` the ordinary solve in F4 multiplications for three
  dependent relations, i.e. `> 100×` per useful relation, and its
  `d_reg` under a single-degree truncation is beyond `19` against `10`.
  The `Z/3` grading is real, but it is a symmetry of the *orbit*, and
  quotienting the per-residual problem by it is experiment 1, not a
  Gröbner lever.  **Class: relabelling** as a solver lever — the graded
  matrix moves work from three residuals into one system that is larger
  than the three solves it replaces.

## 6. Experiment 4 — the `C₃`-invariant formulation against the baselines

**Derivation.**  `invariant_generators` enumerates the weight-`0`
monomials and their minimal generators.  For the diagonal action on the
symmetrised variables and the residual abscissa, weights `(1, 2, 0, 1)`
on `(e₁, e₂, e₃, X)`, the invariant ring has the eight generators
`e₃, e₁e₂, e₂X, e₁²X, e₁X², e₁³, e₂³, X³` (`1, 3, 8, 11` invariant
monomials in degrees `1–4`); for the unsymmetrised diagonal action,
weights `(1, 1, 1, 1)` on `(x₁, x₂, x₃, X)`, it is the cubic Veronese:
`20` cubic generators and nothing below degree `3`.  Because `S₄` has
weight `0`, `S₄(u₁X, u₂X², u₃, X)` with `u₁ = e₁/X`, `u₂ = e₂/X²`,
`u₃ = e₃` is a polynomial in `u` and the orbit invariant `c_R = X³`
alone (`invariant_uses_orbit_invariant_only = true` on every residual):
that is the PDP written in invariants of the diagonal `C₃`, the residual
entering only through its orbit.  The Weil restriction then pins the
`F_p`-rationality of `e₁, e₂`: `u₁ ∈ X⁻¹F_p`, `u₂ ∈ X⁻²F_p`, and
substituting `u₁ = ẽ₁/X`, `u₂ = ẽ₂/X²` gives back the three Weil
components of the ordinary system **term for term**
(`invariant_identical_to_ordinary = true` on every residual, `12` of
`12`; same Macaulay profile, same F4 profile to the multiplication).
The invariant formulation of the *per-residual* PDP has no separate
existence; the only genuinely invariant object is the orbit system of
experiment 3.

**Baselines on the same residuals.**  The strongest ordinary Semaev
baseline is the harness's symmetrised `S₄` solve (`0.9–1.0·10⁶` `F_p`
multiplications per residual, Macaulay `252 × 286` at `D = 10`).  The
repository's function-first baseline (`research/nagao_relations`) is
binary-field Python evidence with `S` unmeasured, so its prime-field
analogue was built in the harness: for `f = y + c₂x² + c₁x + c₀` through
`R`, the three further zeros lie in the base iff the cubic
`N(f)/(x − x_R)`, `N(f) = (c₂x² + c₁x + c₀)² − x³ − b`, is
`c₂²(x³ + t₁x² + t₂x + t₃)` with `t ∈ F_p` — nine `F_p` unknowns (the
coordinates of `c₁, c₂`, and `t`), nine equations of degree `≤ 3`; `c₀`
is eliminated by `f(R) = 0`, which makes the division exact (checked
symbolically), and `e = (−t₁, t₂, −t₃)`.  Every harness triple with
distinct abscissae yields a point of this system
(`function_first_witnessed = true`, `1–2` witnesses per residual), so
it is the right system; solving it is another matter:

| formulation | unknowns | equations | Macaulay | F4 (bound 14, budget 120 s) |
|:--|---:|---:|:--|:--|
| ordinary symmetrised Semaev `S₄` | 3 | 3 | `D = 10`, `252 × 286`, dim `64`, `6.1·10⁵` mults | `D_solve 10`, `221 × 277`, solved, `0.4–3.3 s` |
| `C₃`-invariant `(ẽ₁, ẽ₂, ẽ₃)` | 3 | 3 | identical to ordinary | identical to ordinary |
| orbit `(e, z)`, `z³ = 1` | 4 | 4 | no closure `≤ 19`, `11,985 × 8,855` | `D_solve 10`, `1,668 × 1,632`, `3×` solutions, `7–17 s` |
| function-first (Nagao `L(4O)`) | 9 | 9 | no closure `≤ 6`, `1,980 × 5,005`, standard `1,099` and rising, `3.2·10⁶` mults | reaches `D = 7–8` at `10,641 × 8,789`, `2.5–4.6·10¹⁰` mults, basis `132`, **unsolved at 120 s** |

- The invariant system is the ordinary system: **accounting**.
- The function-first system needs `> 17×` the columns of the ordinary
  Macaulay matrix by degree `6` without closing (`> 30×` by the degree
  `8` its F4 reaches), and its F4 spends `> 3,000×` the ordinary
  solve's multiplications without producing a univariate element; the
  ordinary solve finishes in seconds.
  On `E(F_{p³})` with the subspace base the symmetrised Semaev system
  is the function-first system with the six function coefficients
  eliminated, and eliminating them is what makes it solvable.
  **Class: relabelling** for function-first as a prime-field PDP
  formulation (its binary-field completion counts in
  `research/nagao_relations` are not contradicted; they were never
  costed in this unit).

## 7. Verdict

| lever | what moved | class | test against the boundary |
|:--|:--|:--|:--|
| `⟨ψ⟩` factor-base quotient | `S ÷ 3.0` at every size; columns, relations, residuals, solver calls, NNZ, LA all `÷ 3`; `residuals / floor` `1.00 → 0.98–1.09` | engineering | count moved with its floor; `271–714×` rho (plain), `469–1,237×` (folded, extrapolated) |
| canonical residuals and decompositions | `0` solver calls saved on the uniform stream (`0` found, `≤ 0.011` expected); `97 %` of a pair generator's rows removed as `0 = 0` or unit-dependent | accounting (uniform) / engineering (pair sieve) | at the collision floor; no `S` changes |
| `Z/3`-graded block-aware F4 | exists only on the orbit system: `3.0×` fewer dense cells, `1.4 %` fewer multiplications, `0.5–1.1×` wall; the orbit system costs `34–70×` the ordinary solve | engineering (footprint) / relabelling (as a solver) | boundary is `1.0×` the ordinary solve; measured `34–70×` |
| `C₃`-invariant PDP | identical to the ordinary system on `12/12` residuals | accounting | — |
| function-first `L(4O)` formulation | `> 30×` columns, unsolved at `> 3,000×` the multiplications | relabelling | — |

The one `3` on this page is the order of the automorphism group acting
on the base, and it sits where §10.3 found it on prime fields: in the
number of unknowns.  The scoreboard's verdict — every index-calculus
variant on `E(F_{p³})` is hundreds of times behind rho at the sizes
that fit — is unchanged; the best cell moves from `898` (§11.6 tuned)
to `309–351` at `2^33`, against a rho that folds by the same group.

## 8. What was not done

- The block-aware F4 was tested on the orbit system only; on any system
  with `x_R` fixed there is no block to split (`100 %` coupled rows), so
  there is nothing to run.
- The `20`-generator Veronese (unsymmetrised diagonal) formulation was
  derived (generators and Hilbert counts) but not solved: it embeds the
  orbit system into `20` unknowns with quadratic Veronese relations and
  is bounded below by the orbit system's cost.
- Two seeds per size; the rho reference at two seeds spans `0.16–2.62`
  in `S`.  The control-to-quotient ratio, measured on one stream, does
  not depend on it.
- `F_p` multiplications of the F4 substitution solver include its `p`
  sub-runs; the top-level run's counts are in the JSON per formulation.
