# Shape search at `m ≥ 3`: goals G1–G6, protocol

**Stage diagnostic.** Refutation degrees of descended Semaev systems on toy
binary fields. `S`, end-to-end cost and speedup are **unset**; nothing here
moves a ratio to rho. Lane `shape-search`, 2026-10-15.

Committed before any measured cell ran. The smoke runs used to size the
cells are disclosed in §6. Conductor's control plane was down in this
container (`connection refused` on `localhost:8080`), so no task was reserved.

## 1. The frozen goals, verbatim

Source: `research/notes/index-calculus/RESEARCH_FACTOR_BASE_SHAPE_SEARCH.md`
§7. The goals table was frozen in commit **`110117d8c`** (2026-10-01,
"Factor-base shape sweep at m=2: trace class explains the parity split,
closure is Kneser-bounded"), found through the GitHub commits API for that
path. The local clone is shallow. Its boundary commit `474584339` carries the
file unchanged. The rows, copied without edits:

| id | axis | hypothesis | instrument | kill condition | status |
|:--|:--|:--|:--|:--|:--|
| G1 | shape at `m = 3` | on class-0 targets, `gp` and `fp` have mean refutation degree ≥ `0.5` below random at two of `n ∈ {13, 17, 19}`, `l = ⌈n/3⌉` | chained-`S₃` ladder of `RESEARCH_DREG_MEASUREMENT.md` (random `ℓ`), extended with the shape builders of `factor_base_closure_sweep.rs` | no cell with the gap | **pending** |
| G2 | τ-twisted ordered decomposition | `R = P₁ + τP₂ + τ²P₃`, `x(τʲP) = x(P)^{2ʲ}` is `F₂`-linear, so the chain costs no degree and the tuple is ordered: yield `× m!` at equal `D`; the summands live in `V, σV, σ²V`, so a Frobenius progression keeps them nearly in one space | the same chained ladder with the second summand squared `s` times; `s ∈ {1, 2, l}` | refutation degree rises by `≥ 1` against the plain chain at equal `l` | **pending** (A13 is the `s = d` normal-basis case; no algebraic solve exists for it) |
| G3 | ratio minimal polynomial | at `m ≥ 3` the ratio's minimal polynomial (weight, tail degree) moves `D` beyond its trace effect | G1's instrument over the modulus panel | class-0 `D` flat in the panel to `0.3` | **pending**; at `m = 2` it is flat (T1) |
| G4 | Frobenius-invariant non-subspace bases | `F = {x : Tr(c_j x^{2^{e_j}+1}) = 0, j ≤ k}`, size `2^{n−k}`, Frobenius-stable at prime `n` (no subspace is), membership by `k` quadratic forms; the orbit fold `÷n` applies to a base of any size | `ic search` census for yield; F4 on the `S₃` system plus `2k` quadratic constraints for `D` | per-relation cost above the pair table at matched yield | **proposed**, not built |
| G5 | Möbius pole | `F = {P : (x+β)^{−1} ∈ V}` for a random pole `β` is a non-subspace set in `x`; only `β = 1` carries the 2-torsion symmetry (A14), and the `D*` of other poles is unmeasured | `coordinate_quotients.rs` engine with a free pole | class-0 `D*` flat in `β` to `0.3` | **proposed**; expected null |
| G6 | the gate | any shape surviving G1–G3 at `n ≤ 19` runs at `n = 31` with the real F4 engine, then `n = 53`, before the `m = 83` gate of AGENTS §8a | `koblitz_groebner` / `dreg_sweep` | per-relation cost above enumeration at `n = 31` | **pending**; the `m = 83` gate is not discharged by anything here |

## 2. Instrument

`examples/factor_base_shape_ladder.rs` (new; Rust, no Python in any computed
value). It builds the chained `m = 3` system of
`koblitz_groebner::build_decomposition_system`: `S₃(x₁, x₂, u)` and
`S₃(u, x₃, x_R)`, with `b = 1` (`E₀`). It uses `3ℓ + n` unknowns, the sparse
modulus `find_irreducible_sparse(n)`, and the ladder's semantics.

- **Exact count first.** Each draw's solutions are counted exactly before it
  is measured, and only draws with no solution are measured.
- **Degree scan.** `solving_degree_from(eqs, n_vars, 1, d_max)` is the same
  call `ladder_measure_m` makes. It runs with `KIC_SPARSE_DENSE_FINISH=1`,
  `KIC_SPARSE_F5=1` and `KIC_SPARSE_DENSE_BUDGET_MB=1200`. All three are
  documented as changing cost and never an outcome.
- **One basis per summand.** The system is assembled summand by summand, so
  each summand can carry its own basis (G2) and the modulus can vary (G3).
- **`selftest`** checks three things:
  - the exact count against an exhaustive evaluation of the built
    equations, at `m = 3` and for the G5 system;
  - the G5 count against a direct count over `x = β + 1/u`;
  - the system against `build_decomposition_system`, term for term, for a
    single basis.

Changes from the ladder, all needed by the goals:

- **Shapes.** The subspace comes from a shape builder. The builders are
  ported to `u64` words from `factor_base_closure_sweep.rs`:
  - `random`: the ladder's own `random_subspace_basis`.
  - `random2`: the same builder on an independent stream. This is the
    negative-control arm.
  - `gp`: `⟨1, g, …, g^{l−1}⟩` with a random `g ∉ {0, 1}`.
  - `fp`: `⟨α, α², …, α^{2^{l−1}}⟩` with a random `α ∉ {0, 1}`.
  - `coord`: `⟨1, z, …, z^{l−1}⟩`, which is `gp` with `g = z`. It is
    deterministic, so it has one instance.
  - `subfield`: `F_{2^l}` when `l | n`.
- **Trace class.** Each target `x_R` is tagged with
  `t = Tr(√b/x_R) = Tr(1/x_R)`, the note's §2.2 class. Only class-0
  targets are measured. Class-1 draws are counted and skipped.
- **Sampling.** Targets are drawn uniformly from `F_{2^n}^*` **without
  replacement**, per shape instance. This matters at `n ≤ 13`, where a
  stream with replacement repeats targets.
- **Per-instance profile.** Each instance records `dim V·V`,
  `dim V ∩ σV` and whether `V ⊂ ker Tr`.

**Outcome per draw:**

| outcome | meaning |
|:--|:--|
| `resolved D` | the refutation (or pinning) degree; `refuted` flags which |
| `at_least d_max+1` | the degree-`d_max` matrix was built in full and did not resolve; a mathematical lower bound |
| `caps_hit` | a resource limit; never evidence |

Resource rules (lane PREAMBLE): `nice -n 10`, one thread, address space
capped at 2.5 GB (`prlimit --as=2684354560`), one heavy run at a time. A
cell that would need more is recorded as out of budget and not run. Seed
`20261015` for every measured cell.

## 3. Cells and counts

### G1: shape at `m = 3`

**Frozen cells** `(n, l) = (13, 5), (17, 6), (19, 7)`:

- **`(17, 6)` and `(19, 7)` are out of budget** (derived). They have 35 and
  40 unknowns. The `ℓ = 5` cells at 28–30 unknowns already read ≥ 7
  (`RESEARCH_DREG_MEASUREMENT.md` Results 7–9). Their degree-7 matrices need
  47–50 GB at `(13, 5)` and more here, and even degree 6 has
  `C(35, 6) ≈ 1.6·10⁶` and `C(40, 6) ≈ 3.8·10⁶` monomial columns. This is
  far past the 2.5 GB cap. They are not run.
- **`(13, 5)` at `d_max = 6`** runs one-sided and last, only if time
  remains: `gp` and `fp`, one instance each, one class-0 target each.
  - The random arm is the disclosed feasibility probe of §6 together with
    the four committed ladder draws in
    `research/dreg_fixed_surplus_20260923/runs/cell-13-5-6.u*.jsonl`, all
    `≥ 7`, with their classes recomputed from their `x_r`.
  - A `gp` or `fp` draw resolving at ≤ 6 shows a gap of ≥ 1 on that draw.
    A draw at `≥ 7` leaves the cell undetermined, because both arms are
    censored.
  - Degree 7 does not fit the cap.

**Scaled cells (disclosed deviation).** Prime `n` with `l = ⌈n/3⌉`, as in
the frozen rule, at sizes this host can run. They stand in for the frozen
cells and are scoped as such. The verdict on the frozen cells is reported
separately (§4).

| cell | unknowns | `d_max` | arms | instances × class-0 unsat targets |
|:--|--:|--:|:--|:--|
| `(5, 2)` | 11 | 7 | random, random2, gp, fp, coord | 4 × 2 (coord 1 × 8) |
| `(7, 3)` | 16 | 7 | random, random2, gp, fp, coord | 4 × 2 (coord 1 × 8) |
| `(11, 4)` | 23 | 6 | random, random2, gp, fp, coord | 3 × 2 (coord 1 × 6) |

**Positive control** `(9, 3)`, `d_max = 7`. Here `l | n`, so the subfield
`F_8 ⊂ F_{2^9}` exists. The arms are subfield (1 × 8) and random (4 × 2).
The subfield is the shape the `m = 2` sweep found easiest. If it does not
sit below random at `m = 3` either, the instrument's sensitivity to shape is
in doubt, and the G1 nulls are read as weaker. The prediction is that the
subfield sits below random by at least 0.5.

**Negative control.** `random2` against `random` in every cell. Its gap
estimates the null spread of the statistic.

### G3: ratio minimal polynomial

- **Shape.** The `coord` window `⟨1, z, …, z^{l−1}⟩` under each modulus
  of the panel. The modulus is the minimal polynomial of the ratio `z`.
- **Panel.** Every trinomial, the lightest pentanomial, and the two
  heaviest moduli, in `enumerate_irreducibles` order. That is the closure
  sweep's panel. When there are at most 8 irreducibles, all of them are
  used.

| cell | panel | class-0 unsat targets per modulus |
|:--|:--|:--|
| `(5, 2)`, `d_max = 7` | all 6 irreducibles | 6 |
| `(7, 3)`, `d_max = 7` | 4 trinomials + 1 pentanomial + 2 heaviest = 7 | 6 |
| `(11, 4)`, `d_max = 6` | 2 trinomials + 1 pentanomial + 2 heaviest = 5 | 4 |

The noise reference is the spread of instance means in G1's `random` arm at
the same cell.

### G2: τ-twisted ordered decomposition

- **Cells.** `(5, 2)`, `(7, 3)` and `(9, 3)`, all at `d_max = 7`.
- **Shapes.** `random` and `fp`: 3 instances × 2 class-0 unsat targets per
  arm.
- **Twists.** `s ∈ {1, 2, l}`.
- **Arms.** Every arm of an instance sees the same target stream.

| arm | bases of `(x₁, x₂, x₃)` | what it is |
|:--|:--|:--|
| `plain` | `(V, V, V)` | the reference chain |
| `second_s` | `(V, σ^s V, V)` | the instrument as written: "the second summand squared `s` times" |
| `ordered_s` | `(V, σ^s V, σ^{2s} V)` | the hypothesis's `R = P₁ + τ^s P₂ + τ^{2s} P₃` |

Each arm records the class-0 satisfiable fraction (the yield side) and the
refutation degree of the unsatisfiable draws. An arm with no unsatisfiable
class-0 target within the field is recorded as such, not as a degree.

### G5: Möbius pole (`m = 2`)

- **System.** With `u = 1/(x + β)` and `F = {P : u(P) ∈ V}`, the instrument
  descends `(u₁u₂)²·S₃(β + 1/u₁, β + 1/u₂, x₃)`. It is written so that
  every term has Boolean degree ≤ 2 (see the source).
- **Cells.** The closure sweep's odd points: `(7, 3, d=8)`, `(9, 4, d=8)`,
  `(11, 5, d=7)` and `(13, 6, d=6)`.
- **Draws.** Random `V`, 4 instances, and 4 class-0 unsat targets per arm.
- **Arms.**
  - `x_subspace`: the plain `V` in `x`, the reference.
  - `β = 0`.
  - `β = 1`: the A14 frame.
  - Three random `β ∉ {0, 1}` per instance, pooled as `β_rand`. Their
    individual spread is the noise reference.
- **Engine (deviation).** The frozen instrument names the
  `coordinate_quotients.rs` engine. That engine computes invariant rings of
  point-map groups and measures no refutation degree. The descended system
  above, measured by the same `solving_degree` as every other cell, is the
  `D*` instrument for this base. The substitution is disclosed.

### G4, G6

- **G4.** Not built in this lane. Building it is a new base and a new
  system, and it is the last item of the brief. It is run only if time
  remains, which is not expected.
- **G6.** Runs only for a survivor of G1–G3.

## 4. Decision rules

**Statistic.** The class-0 mean refutation degree per arm and cell. A
censored draw (`at_least d`) enters at its lower bound `d`, and the
censored count is reported. `caps_hit` draws are excluded and reported.

**G1, per shape `X ∈ {gp, fp}` and cell.** Let
`gap = mean(random) − mean(X)`.

| gap status | condition |
|:--|:--|
| established | `gap ≥ 0.5`, `X` has no censored draw, and censored random draws count at their lower bound (which understates `gap`) |
| excluded | both arms fully resolved and `gap < 0.5` |
| excluded by censoring | `X` has censored draws, so `mean(X)` is only a lower bound, and the lower-bound gap is `< 0.5` with `random` fully resolved |
| undetermined | anything else |

Verdicts:

| verdict | condition |
|:--|:--|
| **survives** | a shape has the gap established at two or more cells |
| **killed** (frozen kill) | no cell has the gap established for `gp` or `fp` |
| **not killed, not supported** | anything in between |

The negative control (`random2` vs `random`) is reported beside the
verdict. If its `|gap| ≥ 0.5` at some cell, a G1 gap at that cell is not
distinguishable from noise and is not counted as established.

Scope: the verdict on the scaled cells is a verdict on those cells. On the
frozen cells `(17, 6)` and `(19, 7)` it is **not testable (resource)**.
`(13, 5)` is as §3 says.

**G3, per cell.** The panel spread is `max − min` of the class-0 means over
the moduli.

| verdict | condition |
|:--|:--|
| **killed** (frozen kill: flat to 0.3) | spread ≤ 0.3 at every cell |
| **survives** | spread > 0.3 at some cell, and the spread exceeds the G1 `random` instance-mean spread at that cell |
| **not killed, within noise** | spread > 0.3 but within that noise spread |

**G2, per cell, shape and twisted arm.** Let `rise = mean(arm) − mean(plain)`.

| verdict | condition |
|:--|:--|
| **killed** (frozen kill) | some twisted arm has `rise ≥ 1` with neither arm censored at the compared draws |
| **survives at these sizes** | every arm with measured draws has `rise < 1` |

The yield factor (satisfiable fraction, `ordered_s` vs `plain`) is reported
as *derived*, against the hypothesis's `× m! = × 6`.

**G5, per cell.** The `β`-spread is `max − min` of the class-0 means over
`{β = 0, β = 1, β_rand}`.

| verdict | condition |
|:--|:--|
| **killed** (frozen kill, the expected null) | spread ≤ 0.3 at every cell |
| **not killed** | spread > 0.3 at some cell; read against the spread among the three individual random-`β` arms |

**G6.** No survivor means no G6 run.

## 5. Predictions (mine, before data)

- **G1.** Killed on the scaled cells. The `ℓ`-dependence dominates, and
  shape is worth under one degree.
- **Positive control.** Subfield below random.
- **G3.** Flat.
- **G2.** No rise, with a yield gain for `ordered_s`.
- **G5.** Flat.

## 6. Smoke runs used to size the cells (disclosed; not results)

All smoke runs used seed `1`, not the measured seed, and none is counted.

- `selftest`: ok, 104 checks.
- **`(7, 3)`**, `d_max = 7`, random, gp and fp, 1 × 2: all resolved at 6,
  about 0.5 s a draw. Subfield is absent there, as `3 ∤ 7`.
- **`(11, 4)`**, `d_max = 6`, random, 1 × 1: resolved at 6, 68 s, peak RSS
  0.34 GB.
- **G5 `(13, 6)`**, 1 instance × 2 targets per arm: all resolved at 4,
  under 0.1 s each.
- **G2 `(7, 3)`**, `s = 1`, 1 instance per shape: all measured draws 6.
  - For random `V`, `ordered_s1` found no class-0 unsatisfiable target.
    The smoke loop then drew with replacement, which is what prompted
    sampling without replacement.
- **Feasibility probe `(13, 5)`**, `d_max = 6`, random, 1 × 1: launched
  before this commit, with its result recorded in RESULTS. It decides
  whether the one-sided `(13, 5)` arm is affordable.
