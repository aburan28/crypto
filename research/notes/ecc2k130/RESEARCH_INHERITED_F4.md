# F5 for the Koblitz Gröbner stage: the criterion, and the reduction it should have predicted

**Modules:** `src/cryptanalysis/matrix_f5_f2.rs` (the F5 criterion over the
Boolean ring), `src/cryptanalysis/inherited_f4.rs` (the reduced basis
specialised down the splitting tree), `src/cryptanalysis/koblitz_groebner.rs`
(`SolverEngine::{MatrixF4, MatrixF5, InheritedF4}`).
**Harness:** `cargo run --release --example groebner_stage_bench -- --out DIR`,
the frozen suite of [`RESEARCH_GROEBNER_STAGE.md`](RESEARCH_GROEBNER_STAGE.md);
`cargo run --release --example inherited_f4_probe` for the node-level probe.
**Frozen evidence:** `research/inherited_f4_20260922/`.
**Predecessor:** [`RESEARCH_GROEBNER_STAGE.md`](RESEARCH_GROEBNER_STAGE.md)
and the four Autolab rounds under
`research/sat_factor_base_review_20260908/autolab_groebner_*_20260921/`, which
took the same stage from `2,475 M` to `781 M` word XORs by engineering the
matrix build and the elimination kernel.  This round leaves both alone and asks
what F5 has to say about the *rows*.

**The question.**  The decomposition oracle reduces a Macaulay matrix at every
node of a DPLL tree.  F5's contribution to Gröbner-basis computation is to
predict, before any elimination, which rows will reduce to zero and not build
them.  Two versions of that question were asked here: the textbook one — *which
products `t·f_i` does the F5 criterion remove from one Macaulay matrix?* — and
the one the splitting solver actually poses — *which of the rows a child node
reduces has its parent already reduced?*

**Bottom line.**  The F5 criterion, correctly stated for the Boolean ring,
removes **nothing** from the quadratic degree-3 matrices the oracle solves (the
trivial syzygies first appear at degree 4) and `0–12%` of the rows where linear
equations are present; net of its own cost the stage moves by `1.001×`.  The
second question has a complete answer: **every** row of a child's matrix is
the specialisation of a row its parent has already reduced, and all but the
`≈ 3/n` of them whose pivot contains the assigned variable are *still reduced*
after specialisation.  Building nothing below the root and re-reducing only
those rows cuts the stage's counted work by **`6.25×`** across the frozen ladder
(`7.16×` on the rung that dominates it), deciding all 176 targets identically
and the 36-target holdout identically, at `4.4×` less wall time.  It pays on
the **quadratic** `m = 2` systems, whose reduced rows stay sparse, and loses
(`0.60×`) on the chained cubic `m = 3` systems, whose 40%-rank-deficient
matrices make the reduced rows dense — so the default engine inherits on
quadratic systems only (§3.4).  Class: **engineering**, per §3 of `AGENTS.md`
— the floor below has no term this work can move.

## 0. The boundaries, stated before anything is measured

**The floor.**  Unchanged from the predecessor note: `RESEARCH_ECC2K130_DECOMPOSITION.md`
§5 bounds the whole attack at `m·2^131` group operations independently of `ℓ`,
and §5.3 restates it as a demand that the oracle beat exhaustive search over its
candidate set by `2^{70.19 + log₂ m}`.  It counts candidate tuples; it contains
no term for how one Macaulay matrix is reduced or how many of them are built.
**No result below can move it**, and its ratio column is flat by construction.

**The reference.**  `SolverEngine::MatrixF4 { max_degree: 3 }` — the engine
that was the default before this round, on the tree as merged at `2ac24346`,
run on the frozen ladder of `RESEARCH_GROEBNER_STAGE.md` with the same targets:
`781,399,171` word XORs over 176 targets, `7.04 s` of stage wall on this host.
Both sides of every comparison are the *same binary*; the reference is selected
with `KIC_F4_INHERIT=0`, so nothing but the engine differs.

**The unit.**  64-bit word operations.  For the reference and the F5 variant
that is the elimination's XORs, as the suite has always counted them.  The
inherited engine is charged for its elimination XORs **and** for its
specialisation, at one word operation per word read and per word written when
a basis is specialised — a cost the from-scratch path never pays, because its
matrix build has never been in the unit (§1 of the predecessor note).  The
comparison is therefore conservative against the candidate.  Wall time is a
practicality note (§6 of `AGENTS.md`).

**Falsification target**, inherited from the predecessor note and applied
unchanged: a variant is an improvement only if (a) `compare.py` accepts it —
identical verdict digest, reductions, infeasibility certificates,
propagations, splits, F4 calls and oversize count on every rung — and (b) the
counted word operations fall on every rung, not merely in total.  A rung with
no splitting tree cannot fall and is reported flat.  Inadmissible: changing the
ladder, targets, node budget, degree or size caps; dropping the specialisation
from the count; reporting the stage as an attack cost.

## 1. The F5 criterion over the Boolean ring

### 1.1 What the criterion says, and why `F_2` needs its own

Order the generators `f_1, …, f_m`, `d_i = deg f_i`.  Over a polynomial ring
the row `t·f_i` is redundant at degree `d` whenever `t = LM(g)` for some
`g ∈ ⟨f_1, …, f_{i−1}⟩` of degree `d − d_i`, because
`t·f_i = (t + g)·f_i + g·f_i`: the first summand is rows of `f_i` with smaller
multipliers, the second lies in the span of the rows of `f_1, …, f_{i−1}`.
Induction over the signature order `(i, t)` makes every pruned row a combination
of rows kept.  Over `F_2` with the field equations there is a second trivial
syzygy, the Frobenius one: `f_i² = f_i`, so `(f_i + 1)·f_i = 0`.

The semi-regularity literature folds it in as "prune `t·f_i` when
`t ∈ LM⟨f_1, …, f_i⟩`", with `f_i` itself standing for its Frobenius syzygy.
**In the Boolean ring that is unsound.**  Multiplying by a monomial is not
order-compatible there — `x₁·(x₁x₂ + x₁ + x₂) = x₁`, a *smaller* leading
monomial — so a witness `g = g' + h·f_i` may have its `f_i`-part `h` carrying
monomials *above* `t`, and the induction runs the wrong way.  The test
`naive_f2_criterion_is_unsound_in_the_boolean_ring` keeps the counterexample:
one generator `f = x₀x₁ + x₀ + x₁` in three variables at degree 5, where the
naive rule prunes `x₀·f` and `x₁·f` and the kept rows span a space of rank 4
against the full 6.

What is sound is to use the Frobenius syzygies as the polynomials they are.  With

```text
    V_j(e)  = span{ s·f_k : k ≤ j, deg s ≤ e − d_k }
    W_i(d)  = V_{i−1}(d − d_i)  +  span{ s·(f_i + 1) : deg s ≤ d − 2·d_i }
```

a witness `g = g' + σ ∈ W_i(d)` with `LM(g) = t` gives
`t·f_i = (t + g)·f_i + g'·f_i + σ·f_i = (t + g)·f_i + g'·f_i`, since `σ·f_i = 0`
identically; every monomial of `t + g` is below `t` and so of degree
`≤ d − d_i`, and `(s·f_i)·f_k` expands into multiples of `f_k` of degree
`≤ d`.  **Prune `(t, i)` iff `t ∈ LM(W_i(d))`.**  The leading-monomial sets
come from echelonising the lower-degree matrices `V_j(e)` incrementally in
generator order, one pass per distinct `e = d − d_i`; that work is charged in
the unit.  `f5_rows_span_the_f4_row_space_on_random_systems` pins row-space
equality with matrix-F4 over 40 systems at degrees 2–4.

### 1.2 What it removes here

For a quadratic generator at `d = 3`, `W_i(3) = V_{i−1}(1)` is the span of the
*linear* generators before it, and the Frobenius span is empty (`3 < 2·2`).  On
a system with no linear equations the criterion prunes nothing at all; the
Koszul and Frobenius syzygies of quadratics live at degree 4.  The measurement
agrees to the row:

| rung | rows built (reference) | rows pruned | criterion XORs | word XORs ref → F5 | ratio |
|:--|--:|--:|--:|--:|--:|
| `K_0/2^9`, m=2 | 43,954 | 2,071 (4.7%) | 3,284 | 1,352,350 → 1,348,747 | 1.003× |
| `K_0/2^9`, m=3 | 72,521 | 1,899 (2.6%) | 4,294 | 4,164,666 → 4,149,441 | 1.004× |
| `K_0/2^13`, m=2 | 150,798 | 18,281 (12.1%) | 117,745 | 20,668,245 → 19,776,482 | 1.045× |
| `K_1/2^15`, m=2 | 3,584 | 0 | 0 | 107,470 → 107,470 | 1.000× |
| `K_1/2^17`, m=2 | 537,336 | 0 | 0 | 55,594,932 → 55,594,932 | 1.000× |
| `K_1/2^23`, m=2 | 2,718,600 | 0 | 0 | 699,511,508 → 699,511,508 | 1.000× |

The three `K_0` rungs prune because their node systems acquire linear
equations after substitution; the three `K_1` rungs never do.  This is the
expected result, and the point of recording it is that the F5 criterion is
now *available and priced* for the degree-4 frontier (`koblitz_f4_frontier_probe`),
where the Frobenius and Koszul rows it removes are `m(m+1)/2` per matrix, not a
speculation about what F5 "would" do at degree 3.

## 2. What the splitting solver actually repeats

The oracle's tree has thousands of nodes per target at `n = 23` (`4,934`
reductions over 16 targets, `2,478` splits, `2,438` refutations).  Every node
builds its degree-2 and degree-3 Macaulay matrices from its current system and
reduces them from scratch.  But the current system is `S|_{v = c}` for the
parent's `S`, and substitution is a **ring homomorphism** `φ_c` of the Boolean
ring.  Hence, with multipliers over every variable,

```text
    V(S|_{v=c})  =  φ_c( V(S) )  =  span{ φ_c(ρ) : ρ a reduced row of the parent }
```

up to two corrections handled below.  The child's row space is spanned by the
`rank` reduced rows of its parent, specialised.  And specialisation preserves
most of the echelon structure:

- under `v := 0` every monomial containing `v` is deleted; a row whose leading
  monomial does not contain `v` keeps it;
- under `v := 1` every monomial `m ∋ v` folds into `m ∖ v`, which is *smaller*
  in any degree-compatible order, so a leading monomial without `v` survives
  and nothing folds onto it.

Only the rows whose **pivot contains `v`** lose their place.  Those are
re-reduced against the rest by leading term; everything else is left as it is.
No Macaulay matrix is built below the root.  This is F5's organising idea —
do not redo a reduction you can predict — applied to the axis the solver
actually iterates along, the tree, rather than to the generator index.

**The two corrections.**  (i) A generator whose degree *drops* under `φ_c`
acquires multipliers of higher degree that no parent row maps to; those
products are built and inserted (*completion rows*), with multipliers over the
variables occurring in the specialised system.  (ii) Multipliers containing the
assigned variable, or a variable occurring in no generator, are not images of
parent rows.  They add only rows `x·g` that reach the linear tail when
`g ∈ {0, 1}`, and `g = 1` is already a refutation; so the inherited basis lies
between `V_occurring(S')` and `V_all(S')`, every space in that sandwich has the
same tail, and the solver — which reads nothing but the tail — behaves
identically.  `specialised_basis_spans_the_from_scratch_row_space` pins the
sandwich and the tail over 60 systems, every variable, both values, two levels
deep; `inherited_and_f5_engines_walk_the_same_tree_as_matrix_f4` pins identical
roots and identical `(reductions, refutations, propagations, splits)` on Semaev
systems for `m = 2` and the chained cubic `m = 3`.

**Which root.**  The root builds its Macaulay matrix with multipliers over the
variables occurring in the system (the from-scratch step's own
active-multiplier policy) and reduces it as the from-scratch step reduces that
shape: full RREF below 24 variables, echelon-only at 24 and above.  A fully
reduced root makes every later displaced row reduce in exactly as many XORs as
it has pivot bits, which pays back over a deep tree; an echelon-only root is
cheaper where the tree is shallow.  Both are retained as controls
(`KIC_F4_INHERIT_ROOT=rref|ref`) and tabled below.

## 3. The table

One unit, one table, every variant a row.  Medians of three repetitions; word
operations are deterministic and identical across all three on every side.
Ratios are reference / candidate.  `compare.py` accepted every comparison shown.

| rung | `m` | `ℓ` | targets | decomposed | reference | F5 criterion | inherited F4 | ratio | ratio to floor | wall | class | correct |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|:--|:--|
| `K_0/2^9`  | 2 | 6 | 40 | 39 | 1,352,350 | 1,348,747 | **656,325** | **2.06×** | flat | 1.50× | engineering | ✓ identical |
| `K_0/2^9`  | 3 | 6 | 16 | 15 | 4,164,666 | 4,149,441 | **3,844,010** | **1.08×** | flat | 1.35× | engineering | ✓ identical |
| `K_0/2^13` | 2 | 12 | 40 | 40 | 20,668,245 | 19,776,482 | **12,526,392** | **1.65×** | flat | 1.34× | engineering | ✓ identical |
| `K_1/2^15` | 2 | 4 | 32 | 0 | 107,470 | 107,470 | 107,470 | 1.00× | flat | 0.99× | — (no tree) | ✓ identical |
| `K_1/2^17` | 2 | 8 | 32 | 7 | 55,594,932 | 55,594,932 | **10,050,078** | **5.53×** | flat | 3.78× | engineering | ✓ identical |
| `K_1/2^23` | 2 | 11 | 16 | 9 | 699,511,508 | 699,511,508 | **97,756,456** | **7.16×** | flat | 5.13× | engineering | ✓ identical |
| **total**  |   |   | 176 | 110 | 781,399,171 | 780,488,580 (1.001×) | **124,940,731** | **6.25×** | flat | **4.36×** | engineering | ✓ |

The inherited-F4 column is the engine as shipped, `KIC_F4_INHERIT=1` forcing
it on the one cubic rung (`K_0/2^9`, `m = 3`) so that every rung prices the
inherited path; the shipped default routes cubic systems to the reference
engine (§3.4), so on that rung the default costs exactly the reference's
`4,164,666`.  `K_1/2^15` refutes every target at its root (`32` reductions,
`0` splits): there is no tree to inherit along, so the engine does exactly the
reference's work and the row is flat rather than improved.  The falsification
target's clause (b) is met on every rung that has a tree.

### 3.1 The root controls

| variant | total word ops | ratio | min rung ratio | wall |
|:--|--:|--:|--:|--:|
| inherited F4, root as the from-scratch shape policy (**default**) | 124,940,731 | 6.25× | 1.00× (`K_1/2^15`, no tree) | 4.36× |
| inherited F4, echelon-only root everywhere | 142,095,482 | 5.50× | 1.08× | 4.18× |
| inherited F4, fully reduced root everywhere | 131,002,042 | 5.96× | 1.00× | 4.40× |

The fully reduced root loses `46%` on `K_0/2^13` and is flat on the shallow
`n_vars ≥ 24` rungs, where the reference itself only echelonises; the
echelon-only root loses `16%` on `K_1/2^23`, where displaced rows cascade
against unreduced pivots over a deep tree.  The default takes each regime's
better half.

### 3.2 Holdout

Two instances the tuning never saw, from the predecessor note's holdout ladder:

| instance | targets | decomposed | reference | inherited F4 | ratio | wall | correct |
|:--|--:|--:|--:|--:|--:|--:|:--|
| `K_0/2^19`, m=2 (`ℓ = 18`, 36 unknowns, the widest matrices) | 12 | 12 | 116,912,104 | 73,382,095 | **1.59×** | 1.32× | ✓ identical |
| `K_1/2^17`, m=2, divisor 1 | 24 | 10 | 33,578,834 | 6,460,770 | **5.20×** | 3.60× | ✓ identical |

`K_0/2^19` has 18 F4 calls per target — a shallow tree with a very wide root —
and shows the floor of this method: when almost all the work is the root, the
root is all there is to pay.

### 3.3 Where the inherited engine's operations go

The node-level probe (`research/inherited_f4_20260922/probes/`) captures every
system the production solver reduces, reduces each from scratch, specialises
it on its lowest occurring variable both ways, and compares against a
from-scratch reduction of each child.  Every child's decisive rows agree with
the from-scratch step's through the solver's view (a refutation is a
refutation; otherwise the same forced variables):

| cell | `m` | nodes | mean rank | children from scratch | children inherited | of which re-reduction | of which specialisation | ratio | mismatches |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `K_1/2^23` | 2 | 1,022 | 338 | 130,016 / node | 11,060 / node | 2,248 | 8,812 | **11.8×** | 0 of 2,044 |
| `K_1/2^17` | 2 | 407 | 198 | 33,141 / node | 3,544 / node | 593 | 2,951 | **9.35×** | 0 of 814 |
| `K_0/2^13` | 2 | 111 | 219 | 71,557 / node | 39,884 / node | 27,652 | 12,232 | **1.79×** | 0 of 222 |
| `K_0/2^9`  | 2 | 124 | 73 | 4,778 / node | 1,022 / node | 265 | 756 | **4.68×** | 0 of 239 |
| `K_0/2^9`  | 3 | 120 | 237 | 36,553 / node | 27,104 / node | 13,246 | 13,858 | **1.35×** | 0 of 236 |
| `K_0/2^15` | 3 | 304 | 532 | 195,254 / node | 79,701 / node | 34,297 | 45,404 | **2.45×** | 0 of 608 |

One level below a freshly reduced root, fill-in is not the story on the
quadratic cells: the re-reduction of displaced rows is `2,248` of `11,060`
operations per node at `n = 23`, because only `≈ 70` of `338` rows are
displaced per child and each reduces in a handful of pivot hits against a
fully reduced basis.  The ladder run, which specialises the same basis up to
22 levels deep, distributes the cost differently, and the profile's
`specialise_word_ops` counter splits it exactly:

| `K_1/2^23`, whole ladder run | word operations | share |
|:--|--:|--:|
| root reductions (16 targets, degrees 2 and 3) | 6,719,658 | 7% |
| re-reduction of displaced rows | 60,351,941 | **62%** |
| specialisation (word reads and writes) | 30,585,461 | 31% |
| tail RREF | 128,323 | 0.1% |
| **inherited F4 total** | **97,756,456** | |

Re-reduction costs `12.2 k` per reduction on the ladder against `2.2 k` per
child in the depth-one probe: a fold `m ∋ v ↦ m ∖ v` can land on another row's
pivot column, so the basis drifts away from reduced form as the tree deepens
and each displaced row cascades through more pivots before it settles.  That
drift is an engineering lever bounded above by the `7%` the roots cost.

### 3.4 Where it loses: the cubic systems

The depth-one probe says the cubic `m = 3` cells inherit at `1.35–2.45×`.  The
whole-solve measurement says otherwise, and the whole solve is the measure.
Priced through `ic boundary --oracles` on the frozen 13-cell oracle ladder
(8 targets per cell, seed `123212651130`, both engines on the same binary,
identical verdicts and solver counters on every target):

| cell | unknowns | equations | reference word XORs | inherited F4 (forced) | ratio | class |
|:--|--:|--:|--:|--:|--:|:--|
| `K_0/2^9`, m=3 | 27 | 18 | 2,134,142 | 1,997,844 | 1.07× | engineering |
| `K_0/2^15`, m=3 | 30 | 30 | 1,653,034,947 | 2,746,829,713 | **0.60×** | **regression** |

The `K_0/2^15` cell is 8,920 reductions deep, `3,617` of them refutations, and
its degree-3 matrices are `≈ 1,020` rows of rank `≈ 600`: some `40%`
rank-deficient, against `< 5%` on the quadratic cells.  The reduced rows of a
rank-deficient matrix are dense combinations, so every displaced row and every
completion row (a cubic generator drops to quadratic under one substitution in
three) re-reduces through a hundred-odd dense pivots, `300 k` operations per
reduction against the from-scratch echelon's `187 k` on the sparse matrix.  A
fully reduced root does not help (`0.67×`); neither did two rejected policies —
dropping a basis whose specialisation cost more than a fixed multiple of the
root's per-unit reduction cost (`0.45–0.48×`: by the time the estimate fires,
most of the subtree has paid) and switching the subtree to the from-scratch
engine at the same trigger (`0.45×`, same reason).  What decides it is not a
per-node quantity but the **degree of the system**, known before the first
reduction: [`SolverEngine::effective_for`] inherits on quadratic systems and
reduces cubic ones as `MatrixF4` always did.  The `1.07×` on `K_0/2^9`, `m = 3`
is forgone by that rule; `KIC_F4_INHERIT=1` reclaims it where a caller knows
better.

## 4. The phases

Stage wall on the `K_1/2^23` rung, the one that dominates:

| phase | reference | inherited F4 |
|:--|--:|--:|
| build (reference: Macaulay build; inherited: root build **and every specialisation**) | 2,518 ms | 1,008 ms |
| reduce (reference: elimination; inherited: tail RREF of the low block) | 3,160 ms | 20 ms |
| readback | 66 ms | 0 ms |
| **stage wall** | **5,890 ms** | **1,148 ms** |

The elimination phase, which the unit has always measured and every earlier
round optimised, is `0.3%` of what it was: there is almost nothing left to
eliminate at a child.  What remains is the specialisation, booked under build
because that is what it replaces.  The `27%` readback the predecessor note
found is gone with the matrices it read from.

## 5. What this does not establish

Per `AGENTS.md` §8, and because the number is large enough to tempt:

- **This is a stage diagnostic.**  The Gröbner stage is one phase of relation
  collection, which is one phase of the attack.  No crossover, exponent or rho
  comparison follows; none is claimed.  The oracle's cost at `n = 131` is set
  by the counting argument of `RESEARCH_ECC2K130_DECOMPOSITION.md` §5, and
  `6.3×` off the stage is `2^{2.6}` off `2^{131}`.
- **The algebra is unchanged.**  Same degree, same row space at every node
  (up to the tail-preserving sandwich of §2), same splitting rule, same node
  budget, same verdicts.  Nothing here changes the first fall degree, the
  solving degree, or how many nodes a target needs.
- **The F5 criterion is a negative result at this degree**, stated in numbers
  (§1.2), and a positive capability: `MatrixF5` is the engine to run at the
  degree-4 frontier, where it prunes `m(m+1)/2` rows per matrix by
  construction.  It has not been measured there in this round.
- **The size-cap semantics differ in one corner.**  The reference re-checks the
  row and column caps at every node; the inherited engine checks them where it
  builds, i.e. at the root and at completion.  A child whose from-scratch
  matrix would exceed the caps while its parent's did not (a generator dropping
  a degree at the cap's edge) is reduced by the inherited engine and skipped by
  the reference.  It did not occur on the ladder or the holdout (`oversize = 0`
  everywhere); if it did, `compare.py` would refuse the comparison.
- **The accounting is asymmetric against the candidate**, deliberately: the
  reference's matrix build is not in the unit, the candidate's specialisation
  is — `30.6 M` of its `97.8 M` operations at `n = 23`.  Dropping that charge
  would read `10.4×` on that rung; it is not dropped, and the `7.16×` stands.
- **It is a quadratic-system result.**  On the chained cubic systems the
  shipped engine reduces exactly as the reference does (§3.4), so those cells
  are flat by construction; the one measured attempt to inherit on them
  regressed `0.60×` and is retained as a rejected variant.
- **Wall time is a practicality note.**  `4.4×` on this host, single thread,
  no pinning.  The `ic run` end-to-end effect is not measured here.

## 6. Reproducing

```bash
cargo build --release --example groebner_stage_bench --example inherited_f4_probe

# reference (the pre-round default engine), candidate (the default), F5 criterion
KIC_F4_INHERIT=0 ./target/release/examples/groebner_stage_bench --label reference --out /tmp/gs/baseline/rep1
./target/release/examples/groebner_stage_bench --label inherited --out /tmp/gs/inherited_f4/rep1
KIC_F4_INHERIT=0 KIC_F4_CRITERION=f5 ./target/release/examples/groebner_stage_bench --label f5 --out /tmp/gs/f5_criterion/rep1
python3 research/groebner_stage_20260915/compare.py /tmp/gs/baseline /tmp/gs/inherited_f4

# holdout, root controls, node-level probe (args: a n targets degree m)
KIC_F4_INHERIT=0 ./target/release/examples/groebner_stage_bench --ladder holdout --out /tmp/gs/holdout_baseline/rep1
./target/release/examples/groebner_stage_bench --ladder holdout --out /tmp/gs/holdout_inherited_f4/rep1
KIC_F4_INHERIT=1 KIC_F4_INHERIT_ROOT=ref  ./target/release/examples/groebner_stage_bench --out /tmp/gs/ref_root/rep1
KIC_F4_INHERIT=1 KIC_F4_INHERIT_ROOT=rref ./target/release/examples/groebner_stage_bench --out /tmp/gs/rref_root/rep1
KIC_F4_INHERIT=1 ./target/release/examples/inherited_f4_probe 1 23 2 3 2
KIC_F4_INHERIT=1 ./target/release/examples/inherited_f4_probe 0 15 8 3 3

# the oracle ladder, both engines on one binary (§3.4 and the scoreboard rows)
KIC_F4_INHERIT=0 ./target/release/ic boundary --regime koblitz --koblitz-degrees 11 --repeats 1 \
    --seed 123212651130 --no-fold-max-degree 31 --s4-max-degree 31 --oracles --oracle-targets 8 \
    --json --out /tmp/pricing-reference.json
./target/release/ic boundary --regime koblitz --koblitz-degrees 11 --repeats 1 \
    --seed 123212651130 --no-fold-max-degree 31 --s4-max-degree 31 --oracles --oracle-targets 8 \
    --json --out /tmp/pricing-inherited.json

cargo test --release --lib -- cryptanalysis::inherited_f4 cryptanalysis::matrix_f5_f2 cryptanalysis::koblitz_groebner
```

`KIC_F4_NODE_DUMP=path` appends every node system the solver reduces as JSON
lines, which is how the probe gets real systems rather than synthetic ones.

## References

- J.-C. Faugère, *A new efficient algorithm for computing Gröbner bases without
  reduction to zero (F5)*, ISSAC 2002.
- M. Bardet, J.-C. Faugère, B. Salvy, *On the complexity of Gröbner basis
  computation of semi-regular overdetermined algebraic equations*, ICPSS 2004.
- M. Bardet, J.-C. Faugère, B. Salvy, B.-Y. Yang, *Asymptotic behaviour of the
  degree of regularity of semi-regular polynomial systems*, MEGA 2005 — the
  `F_2` criterion whose Boolean-ring form §1.1 corrects.
- C. Eder, J.-C. Faugère, *A survey on signature-based algorithms for computing
  Gröbner bases*, J. Symbolic Comput. 80 (2017).
- J.-C. Faugère, *A new efficient algorithm for computing Gröbner bases (F4)*,
  J. Pure Appl. Algebra 139 (1999).
