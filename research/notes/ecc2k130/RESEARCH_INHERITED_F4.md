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
the specialisation of a row its parent has already reduced, and the rows whose
pivot does not contain the assigned variable are *still reduced* after
specialisation.  Building nothing below the root and re-reducing only the
displaced rows, on the tree the solver has always walked, cuts the stage's
counted work by **`6.93×`** across the frozen ladder, deciding all 176 targets
identically.  Two further facts, both measured before they were used, take it
from there.  *Which variable the solver splits on* decides how many rows are
displaced: a pivot is a row's largest monomial, so the historical rule — split
on the lowest-indexed free variable, the largest in the order — picks the
variable most pivots contain and displaces half the basis per level deep in
the tree; splitting on the smallest free variable halves the engine's work at
every rung with a tree and turns its one loss (the chained cubic systems,
`0.60×`) into a win.  And *three quarters of the rows a level keeps are never
touched at that level*, so a kept row is not rewritten into the child's layout
at all: it is shared by reference and materialised only when something needs
it, through the composition of every layout step since, charged once.
Against the pre-round default the shipped configuration costs **`20.95×`**
fewer word operations on the ladder (`27.7×` on the rung that dominates it,
`29.9×` on the deep cubic cell) and `5.87×` on the holdout, with the same
targets decomposed and every decomposition verified in the group, on a tree of
the same size but a different shape; `26.1×` against the from-scratch engine
on the very same tree.  Class: **engineering**, per §3 of `AGENTS.md` — the
floor below has no term this work can move.

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

**Falsification target**, inherited from the predecessor note: a variant is an
improvement only if (a) `compare.py` accepts it — identical verdict digest,
reductions, infeasibility certificates, propagations, splits, F4 calls and
oversize count on every rung — and (b) the counted word operations fall on
every rung, not merely in total.  A rung with no splitting tree cannot fall and
is reported flat.  Inadmissible: changing the ladder, targets, node budget,
degree or size caps; dropping the specialisation from the count; reporting the
stage as an attack cost.

**One admitted departure, declared here.**  §3.5 changes the variable the
solver splits on.  That keeps the algebra, the targets and the answers but not
the tree, so clause (a) cannot hold across the change: when several
decompositions exist a different one is found first, and the node counts
move.  Two comparisons are therefore reported for it and labelled: a
**same-tree** one, reference and candidate both under the new rule, which
`compare.py` accepts and to which clauses (a) and (b) apply; and a
**cross-tree** one against the pre-round default, checked by
`research/inherited_f4_20260922/compare_cross_tree.py` for the same
instances, the same number of targets decomposed, no exhausted budget and no
skipped matrix on every rung — every decomposition the harness counts having
been lifted and verified against the group identity by `groebner_decompose`.
The cross-tree ratio is the whole-method `baseline / candidate` of `AGENTS.md`
§8 for this stage; it is not a same-tree identity and is never presented as
one.

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

**Lazy materialisation.**  Instrumented on `K_1/2^23`, `148` of the `194`
rows a specialisation keeps — `76%` — are neither displaced nor hit as a
pivot while a displaced row reduces, nor read for the tail; the eager
implementation rewrote every one of them into the child's layout and charged
one word read and one written per word for it.  Now a kept row is shared by
reference between the parent and both children, tagged with the layout epoch
it was last written in; each basis records its layout steps, and a row is
materialised only when its content is needed, through the composition of
every step since it was last written, in one pass charged once.  Its pivot
column is tracked through the steps meanwhile, which is all `specialise`
needs to decide whether it stays a pivot.  The composed maps are memoised on
the shared steps, so a child composes one step on top of what its parent
already built.  Nothing about the row space changes; the same rows are
reduced by the same pivots, later.

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
Ratios are reference / candidate.  The reference is the pre-round default,
`MatrixF4` splitting on the lowest free variable.  The candidate columns are:
the inherited engine on the **same tree** (`compare.py` accepts it, clauses (a)
and (b) apply); the from-scratch engine on the **new tree** of §3.5; the
inherited engine on the new tree with eager materialisation, the default of the
round before this one, kept as the before-mark; and the inherited engine on the
new tree with lazy materialisation (§2), which is the configuration shipped.
The last three are cross-tree against the reference and same-tree against each
other; the `nodes` column shows the size of each tree.

| rung | `m` | `ℓ` | targets | decomposed | nodes L → H | reference (F4, lowest) | F5 criterion | inherited, lowest (same tree) | F4, highest | inherited, highest, eager (before) | **inherited, highest, lazy (default)** | ratio to reference | ratio to floor | wall | class | correct |
|:--|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|:--|--:|:--|:--|
| `K_0/2^9`  | 2 | 6 | 40 | 39 | 386 → 338 | 1,352,350 | 1,348,747 | 606,611 (2.23×) | 1,257,665 | 475,686 | **391,693** | **3.45×** | flat | 2.15× | engineering | ✓ verified |
| `K_0/2^9`  | 3 | 6 | 16 | 15 | 285 → 382 | 4,164,666 | 4,149,441 | 2,760,501 (1.51×) | 2,026,237 | 2,524,783 | **1,783,018** | **2.34×** | flat | 1.51× | engineering | ✓ verified |
| `K_0/2^13` | 2 | 12 | 40 | 40 | 599 → 520 | 20,668,245 | 19,776,482 | 10,777,601 (1.92×) | 25,336,199 | 8,302,715 | **5,588,277** | **3.70×** | flat | 2.02× | engineering | ✓ verified |
| `K_1/2^15` | 2 | 4 | 32 | 0 | 32 → 32 | 107,470 | 107,470 | 107,470 (1.00×) | 107,470 | 107,470 | 107,470 | 1.00× | flat | 0.90× | — (no tree) | ✓ identical |
| `K_1/2^17` | 2 | 8 | 32 | 7 | 1,763 → 1,769 | 55,594,932 | 55,594,932 | 8,733,728 (6.37×) | 62,049,801 | 6,535,588 | **4,205,307** | **13.22×** | flat | 6.85× | engineering | ✓ verified |
| `K_1/2^23` | 2 | 11 | 16 | 9 | 4,934 → 4,876 | 699,511,508 | 699,511,508 | 89,762,781 (7.79×) | 883,061,010 | 44,225,148 | **25,226,634** | **27.73×** | flat | 12.33× | engineering | ✓ verified |
| **total**  |   |   | 176 | 110 | 7,999 → 7,917 | 781,399,171 | 780,488,580 (1.001×) | 112,748,692 (**6.93×**) | 973,838,382 (0.80×) | 62,171,390 (12.57×) | **37,302,399** | **20.95×** | flat | **9.29×** | engineering | ✓ |

*Same-tree, under the new rule:* `F4, highest` against the default is
`26.11×` in total, accepted by `compare.py` rung for rung, and the inherited
engine now wins every rung with a tree — including the shallow cubic
`K_0/2^9`, `m = 3`, where the eager engine had lost `0.80×` on the same tree.
*Same-tree, under the old rule:* `6.93×`, every rung with a tree improved
(`6.29×` before lazy materialisation).  *Correct* reads "verified" rather than
"identical" wherever the tree changed: the same targets decompose and none is
exhausted, and every decomposition returned was lifted and checked against
the group identity; where several exist, a different one may be found first,
so the verdict digest differs.  `K_1/2^15` refutes every target at its root
(`32` reductions, `0` splits): there is no tree to inherit along and no choice
of split variable to make, so the row is flat rather than improved.

### 3.1 The root controls

All under the shipped rule with lazy materialisation, same tree as the default:

| variant | total word ops | ratio to reference | min rung ratio vs default | wall |
|:--|--:|--:|--:|--:|
| root as the from-scratch shape policy (**default**) | 37,302,399 | 20.95× | — | 9.29× |
| echelon-only root everywhere | 39,060,565 | 20.01× | 0.93× (`K_1/2^23`) | 8.83× |
| fully reduced root everywhere | 45,006,527 | 17.36× | 0.43× (`K_0/2^13`) | 9.16× |

The echelon-only root wins the two shallow rungs and loses `8%` on
`K_1/2^23`; the fully reduced root loses `57%` on `K_0/2^13`.  The default
takes each regime's better half and is the only one of the three that never
regresses a rung.

### 3.2 Holdout

Two instances the tuning never saw, from the predecessor note's holdout ladder:

| instance | targets | decomposed | nodes L → H | reference | F4, highest | **default** | ratio to reference | same-tree ratio | wall | correct |
|:--|--:|--:|:--|--:|--:|--:|--:|--:|--:|:--|
| `K_0/2^19`, m=2 (`ℓ = 18`, 36 unknowns, the widest matrices) | 12 | 12 | 228 → 228 | 116,912,104 | 163,178,360 | **22,640,333** | **5.16×** | 7.21× | 3.01× | ✓ identical tree |
| `K_1/2^17`, m=2, divisor 1 | 24 | 10 | 1,062 → 1,198 | 33,578,834 | 42,029,541 | **3,011,390** | **11.15×** | 13.96× | 5.94× | ✓ verified |
| **total** | 36 | 22 | | 150,490,938 | 205,207,901 | **25,651,723** | **5.87×** | **8.00×** | 4.21× | ✓ |

`K_0/2^19` has 18 F4 calls per target and its tree does not change shape under
the new rule at all; the earlier `1.56×` on this rung, the floor of the method
when almost all the work is the root, is `5.16×` once the root's rows are
displaced less often and rewritten only when touched.

### 3.3 Where the inherited engine's operations go

The node-level probe (`research/inherited_f4_20260922/probes/`) captures every
system the production solver reduces, reduces each from scratch, specialises
it on the variable the production rule would split on both ways, and compares
against a from-scratch reduction of each child.  Every child's decisive rows
agree with the from-scratch step's through the solver's view (a refutation is
a refutation; otherwise the same forced variables):

| cell | `m` | nodes | mean rank | children from scratch | children inherited | of which re-reduction | of which specialisation | ratio | mismatches |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `K_1/2^23` | 2 | 1,022 | 338 | 138,248 / node | 3,436 / node | 949 | 2,487 | **40.2×** | 0 of 2,044 |
| `K_1/2^17` | 2 | 442 | 197 | 34,099 / node | 1,169 / node | 265 | 903 | **29.2×** | 0 of 884 |
| `K_0/2^13` | 2 | 96 | 233 | 89,913 / node | 20,556 / node | 18,805 | 1,751 | **4.37×** | 0 of 192 |
| `K_0/2^9`  | 2 | 111 | 70 | 4,934 / node | 345 / node | 133 | 212 | **14.3×** | 0 of 216 |
| `K_0/2^9`  | 3 | 165 | 154 | 13,140 / node | 3,260 / node | 2,009 | 1,251 | **4.03×** | 0 of 318 |
| `K_0/2^15` | 3 | 3,672 | 289 | 38,469 / node | 4,758 / node | 1,594 | 3,164 | **8.09×** | 0 of 7,344 |

One level below a root, the specialisation charge is now what the displaced
rows and the pivots they hit cost to materialise and nothing else.  The whole
ladder run, where a row skipped at one level is materialised through several
steps at once when it is finally needed, splits the dominant rung as follows
(the eager column is the round before, for the before-mark):

| `K_1/2^23`, whole ladder run | eager (before) | share | **lazy (default)** | share |
|:--|--:|--:|--:|--:|
| root reductions (16 targets, degrees 2 and 3) | 6,719,658 | 15% | 6,719,658 | 27% |
| re-reduction of displaced rows | 7,150,501 | 16% | 7,150,501 | 28% |
| specialisation (word reads and writes) | 30,228,073 | 68% | **11,229,559** | 45% |
| tail RREF | 126,916 | 0.3% | 126,916 | 0.5% |
| **inherited F4 total** | **44,225,148** | | **25,226,634** | |

Of the `11.2 M` that remain, `7.9 M` materialise displaced rows and `3.3 M`
the pivots they hit; the tail rows cost nothing, because a row leading in a
low column is rarely disturbed.  Both remaining components are proportional
to the number of displaced rows, as is the re-reduction — which is what §3.6
tried, and failed, to lower further.  The roots are now the largest fixed
cost.

### 3.4 Where it lost, and why: the cubic systems under the old rule

Under the historical `LowestFree` rule the inherited engine lost on the chained
cubic `m = 3` systems, and that loss fixed an intermediate revision's policy of
reducing cubic systems from scratch.  Priced through `ic boundary --oracles`
on the frozen oracle ladder (8 targets per cell, seed `123212651130`, both
engines on the same binary, identical verdicts and solver counters on every
target):

| cell | unknowns | equations | reference (F4, lowest) | inherited, lowest | ratio | class |
|:--|--:|--:|--:|--:|--:|:--|
| `K_0/2^9`, m=3 | 27 | 18 | 2,134,142 | 1,394,941 | 1.53× | engineering |
| `K_0/2^15`, m=3 | 30 | 30 | 1,653,034,947 | 2,670,183,123 | **0.62×** | **regression, superseded by §3.5** |

The `K_0/2^15` cell is 8,920 reductions deep and its degree-3 matrices are
`≈ 1,020` rows of rank `≈ 600`: some `40%` rank-deficient, against `< 5%` on
the quadratic cells.  The reduced rows of a rank-deficient matrix are dense
combinations, so every displaced row re-reduces through a hundred-odd dense
pivots — and under `LowestFree` half the rows are displaced per level.  Two
per-node rescue policies were tried and rejected on this cell — dropping a
basis whose specialisation cost more than a fixed multiple of the root's
per-unit reduction cost (`0.48×`) and switching the subtree to the
from-scratch engine at the same trigger (`0.45×`): by the time a per-node
estimate fires, most of the subtree has paid.  A fully reduced root did not
help either (`0.67×`).  The evidence for all of them is retained under
`oracle_ladder/`.  What did help was not a per-node quantity at all (§3.5).

### 3.5 The split variable

The tree is built by splitting on a free variable; the historical rule takes
the lowest-indexed one.  The Macaulay columns are ordered DegRevLex with
`v_0 > v_1 > …`, a reduced row's pivot is its *largest* monomial, and the
largest monomials are made of the largest variables — so `LowestFree` splits
on exactly the variable most pivots contain.  Instrumented on `K_1/2^23`, it
displaces `94` of `180` rows per level on average, and each displaced row
then re-reduces through `40` pivots.  Restoring reduced form after every level
was the first thing tried: it does cut the cascade to `9.8` pivots per
displaced row, but the pass costs `71.9 M` word XORs against the `49 M` it
saves, because a reduced basis is denser than an echelon one and every
specialisation refills some twenty pivot columns per row; it is retained as
the `KIC_F4_INHERIT_RREF` control (`94,952,249` on the ladder against the
default's `37,302,399`).

`SplitRule::HighestFree` splits on the **smallest** free variable that still
occurs.  It is as exhaustive as every other rule — a node's subtree still
covers both values of some free variable — and the algebra at each node is
the same reduction of the same kind of system; only the order of the tree's
levels changes.  Its effect, per engine, on the same targets (all inherited
figures with lazy materialisation):

| configuration | ladder total | ratio to reference | cubic `K_0/2^9` | cubic `K_0/2^15` |
|:--|--:|--:|--:|--:|
| `MatrixF4`, lowest (reference) | 781,399,171 | 1.00× | 2,134,142 | 1,653,034,947 |
| `MatrixF4`, highest | 973,838,382 | 0.80× | 1,026,870 (2.08×) | 305,067,403 (5.42×) |
| `InheritedF4`, lowest | 112,748,692 | 6.93× | 1,394,941 (1.53×) | 2,670,183,123 (0.62×) |
| **`InheritedF4`, highest (shipped)** | **37,302,399** | **20.95×** | **942,041 (2.27×)** | **55,202,947 (29.94×)** |
| *before-mark, 2026-09-24: the same engine on the current main (fresh reference of [`RESEARCH_CHAIN_SPLIT_ORDER.md`](RESEARCH_CHAIN_SPLIT_ORDER.md) §6.4)* | *16,696,846* | | *560,420* | *25,639,588* |
| *`InheritedF4`, interleaved chain order + linear elimination (default since 2026-09-24)* | *16,125,634* | | *354,877* | *7,187,062* |

The from-scratch engine pays `25%` more on the quadratic rungs under the new
rule — same matrix shapes, same tree size within `2%`, `28%` more XORs per
elimination — and `2–5×` less on the cubic cells; the inherited engine gains
everywhere, because the number of displaced rows falls: on `K_1/2^23` its
re-reduction goes from `60.4 M` to `7.15 M` word operations.  On the cubic
cells it now wins on the same tree as `MatrixF4, highest` on both
(`K_0/2^9`, `m = 3`: `1.09×`; `K_0/2^15`: `5.53×`), so the engine inherits on
every degree and the intermediate degree routing is gone.  `SplitRule::Auto`,
the new default of `split_rule_default()`, resolves to `HighestFree` under the
inherited engine and to `LowestFree` otherwise, which is why
`KIC_F4_INHERIT=0` still reproduces the reference to the digit.

### 3.6 Rejected: splitting on the variable in the fewest pivots

`HighestFree` picks the variable *likely* to be in the fewest pivots; the
basis knows the actual count.  A rule that picked, at every node, the free
occurring variable contained in the fewest pivot columns of the current bases
was measured on the ladder against the default and rejected: the per-level
work falls as intended, but the variable that displaces the fewest rows is
also the one that constrains the system least, and the tree grows to
compensate — `4,876 → 8,660` reductions on `K_1/2^23`, `382 → 8,904` on the
cubic `K_0/2^9`, `m = 3`.  Ladder total `73,616,817` against `37,302,399`
(`0.51×`), every rung with a tree worse, the same targets decomposed.  The
split variable must be chosen for the algebra first; ~~`HighestFree` happens to
be good for both~~ *(struck 2026-09-24: for a chain, `m ≥ 3`, it is not — in the
layout order it branches first on the `n` coordinates of the intermediate
point, the least constrained unknowns in the system.  Handing it the
variables in the interleaved order, so that it fixes the last summand first,
and eliminating the linear generators that leaves costs `4.91×` fewer word
operations on a registered holdout; see
[`RESEARCH_CHAIN_SPLIT_ORDER.md`](RESEARCH_CHAIN_SPLIT_ORDER.md) §6)*.

## 4. The phases

Stage wall on the `K_1/2^23` rung, the one that dominates, on the merged tree
(the reference arm carries `main`'s fused flat-matrix packing, #601):

| phase | reference (F4, lowest) | inherited F4 (default) |
|:--|--:|--:|
| build (reference: Macaulay build; inherited: root build, every layout step **and every materialisation**) | 1,981 ms | 333 ms |
| reduce (reference: elimination; inherited: tail RREF of the low block) | 3,164 ms | 10 ms |
| readback | 70 ms | 0 ms |
| **stage wall** | **5,352 ms** | **434 ms** |

The elimination phase, which the unit has always measured and every earlier
round optimised, is `0.3%` of what it was: there is almost nothing left to
eliminate at a child.  What remains is booked under build because that is
what it replaces, and it is now mostly bookkeeping outside the unit:
materialising and inserting the displaced rows (`207 ms`), substituting the
system at every level (`64 ms`, once in the solver and once in the basis),
and the solver's own overhead around the reductions.  A round of
metric-neutral work took the build phase from `609` to `333 ms` — the new
layout is a merge of two sorted lists rather than a comparison sort, the
column index is built only when completion rows need it, the composed maps
are shared with the parent — and none of it moved a word operation.  The
`27%` readback the predecessor note found is gone with the matrices it read
from.

## 5. What this does not establish

Per `AGENTS.md` §8, and because the number is large enough to tempt:

- **This is a stage diagnostic.**  The Gröbner stage is one phase of relation
  collection, which is one phase of the attack.  No crossover, exponent or rho
  comparison follows; none is claimed.  The oracle's cost at `n = 131` is set
  by the counting argument of `RESEARCH_ECC2K130_DECOMPOSITION.md` §5, and
  `21×` off the stage is `2^{4.4}` off `2^{131}`.
- **The algebra is unchanged; the tree is not.**  Same degree, same row space
  at every node (up to the tail-preserving sandwich of §2), same node budget,
  same refutation power at a node.  The split rule changes the order of the
  tree's levels, so a target with several decompositions may return a
  different one and the node counts move by a few percent in either
  direction; that is why the shipped configuration is compared cross-tree
  against the reference (§0) and same-tree against `MatrixF4` under its own
  rule.  Nothing here changes the first fall degree or the solving degree.
- **The from-scratch engine is not helped by the new rule on quadratic
  systems.**  `MatrixF4, highest` costs `0.80×` the reference on the ladder;
  the rule is the inherited engine's, and `SplitRule::Auto` gives it to that
  engine only.
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
  everywhere); both comparison scripts refuse a comparison in which it does.
- **The accounting is asymmetric against the candidate**, deliberately: the
  reference's matrix build is not in the unit, the candidate's specialisation
  is — `11.2 M` of its `25.2 M` operations at `n = 23`.  Dropping that charge
  would read `50×` on that rung; it is not dropped, and the `27.7×` stands.
  Lazy materialisation lowered that charge by doing less of the work, not by
  charging less for it: a row is still charged one read and one write per
  word every time it is rewritten, and index bookkeeping (the layout steps
  and their composition) is uncharged on both sides, as the from-scratch
  build's column indexing always was.
- **The tree is the same size, not the same shape.**  Reductions differ from
  the reference by `−1%` to `+34%` per rung (`+1%` over the ladder), in both
  directions; the ratios above are whole-solve totals that include that.
- **Wall time is a practicality note.**  `9.3×` on the ladder and `12.3×` on
  its dominant rung on this host, single thread, no pinning.  The `ic run`
  end-to-end effect is not measured here.

## 6. Reproducing

```bash
cargo build --release --bin ic --example groebner_stage_bench --example inherited_f4_probe

# reference (the pre-round default: MatrixF4, lowest free variable), the shipped default,
# the same-tree controls, and the F5 criterion
KIC_F4_INHERIT=0 ./target/release/examples/groebner_stage_bench --label reference --out /tmp/gs/baseline/rep1
./target/release/examples/groebner_stage_bench --label default --out /tmp/gs/inherited_f4/rep1
SOLVER_SPLIT_RULE=lowest ./target/release/examples/groebner_stage_bench --out /tmp/gs/inherited_f4_lowest/rep1
SOLVER_SPLIT_RULE=highest KIC_F4_INHERIT=0 ./target/release/examples/groebner_stage_bench --out /tmp/gs/matrix_f4_highest/rep1
KIC_F4_INHERIT=0 KIC_F4_CRITERION=f5 ./target/release/examples/groebner_stage_bench --out /tmp/gs/f5_criterion/rep1

# same-tree comparisons (the frozen suite's script) and the cross-tree one
python3 research/groebner_stage_20260915/compare.py /tmp/gs/baseline /tmp/gs/inherited_f4_lowest
python3 research/groebner_stage_20260915/compare.py /tmp/gs/matrix_f4_highest /tmp/gs/inherited_f4
python3 research/inherited_f4_20260922/compare_cross_tree.py /tmp/gs/baseline /tmp/gs/inherited_f4

# holdout, root and RREF controls, node-level probe (args: a n targets degree m)
# (the rejected fewest-pivots rule of §3.6 was an uncommitted experiment; its numbers are in the note only)
KIC_F4_INHERIT=0 ./target/release/examples/groebner_stage_bench --ladder holdout --out /tmp/gs/holdout_baseline/rep1
./target/release/examples/groebner_stage_bench --ladder holdout --out /tmp/gs/holdout_inherited_f4/rep1
KIC_F4_INHERIT_ROOT=ref  ./target/release/examples/groebner_stage_bench --out /tmp/gs/ref_root/rep1
KIC_F4_INHERIT_ROOT=rref ./target/release/examples/groebner_stage_bench --out /tmp/gs/rref_root/rep1
KIC_F4_INHERIT_RREF=1    ./target/release/examples/groebner_stage_bench --out /tmp/gs/rref_every_level/rep1
./target/release/examples/inherited_f4_probe 1 23 2 3 2
./target/release/examples/inherited_f4_probe 0 15 8 3 3

# the oracle ladder, both engines on one binary (§3.4, §3.5 and the scoreboard rows)
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
