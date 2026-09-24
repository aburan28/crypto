# The chained decomposition system: split order, linear elimination, and the degree-drop rule

**Modules:** `src/cryptanalysis/koblitz_groebner.rs` (`eliminate_linear_generators`,
`InheritPolicy`, `DecompositionSystem::interleaved_order`),
`src/cryptanalysis/koblitz_index_calculus.rs` (`groebner_decompose`),
`src/cryptanalysis/inherited_f4.rs` (`ChildSystem::has_dropped`).
**Harness:** `examples/groebner_stage_bench.rs` (the frozen ladder of
[`RESEARCH_GROEBNER_STAGE.md`](RESEARCH_GROEBNER_STAGE.md) plus the two chained
ladders registered below), `ic boundary --oracles` (the oracle-pricing ladder of
[`RESEARCH_INHERITED_F4.md`](RESEARCH_INHERITED_F4.md) §3.4), `ic run --solver
groebner --summands 3` (whole logarithms).
**Frozen evidence:** `research/chain_split_order_20260924/`.
**Predecessor:** [`RESEARCH_INHERITED_F4.md`](RESEARCH_INHERITED_F4.md), whose engine
this round keeps and whose split rule it keeps; only the order the variables are
handed to that rule, and what the solver does with a linear generator, change.

**Status: registered, then amended (§5), then measured.**  §0–§4 were
committed before any registered run (`2809b498`); §1 lists every run made
before registration.  §5 is a dated amendment made after the registered runs
and before any run of the supplementary holdout it registers.  Results are
appended below §5 and do not edit it; corrections to §0–§4 are struck in place
and point to §5.

## The question

For `m ≥ 3` the decomposition oracle solves the chain

```text
    S₃(x₁, x₂, e₁),  S₃(e₁, x₃, e₂),  …,  S₃(e_{m−2}, x_m, x_R)
```

in `m·ℓ + (m−2)·n` Boolean unknowns: `ℓ` per summand `x_i` (its coordinates in
the factor base's subspace `V`) and `n` per intermediate point `e_j`.  The layout
puts the summands first and the intermediate points last.  The inherited engine
splits on the highest-indexed free variable (`SplitRule::HighestFree`,
[`RESEARCH_INHERITED_F4.md`](RESEARCH_INHERITED_F4.md) §3.5), so for `m ≥ 3` it
branches first on the `n` coordinates of `e_{m−2}` — the least constrained
unknowns in the system, `n > ℓ` of them.

The system is multilinear in its blocks (`DecompositionSystem::blocks`): squaring
is `F_2`-linear, so the last link `S₃(e_{m−2}, x_m, x_R)` is bilinear in
`(e_{m−2}, x_m)`.  Once the `ℓ` coordinates of `x_m` are fixed it is **linear** in
`e_{m−2}`: `n` linear equations that fix the intermediate point up to the kernel
of a quadratic's `F_2`-linearisation (one or two points).  The link before it is
then bilinear in `(e_{m−3}, x_{m−1})`, and so on down to an `m = 2` system in
`x₁, x₂`, the one the quadratic rungs of the frozen ladder already solve.  That
is the hybrid *guess one summand, eliminate the intermediate point, reduce the
rest* — a tree of at most `2^{(m−1)ℓ}` leaves instead of one over
`(m−2)·n + (m−1)·ℓ` split variables.

Two things stand between the engine and that tree:

1. **The order.**  `HighestFree` reaches `x_m` only after `e_{m−2}`.
   `DecompositionSystem::interleaved_order` renames the variables
   `x₁, x₂, e₁, x₃, e₂, …, e_{m−2}, x_m` from the lowest index to the highest, so
   the same rule fixes `x_m` first; `groebner_decompose` solves the renamed system
   and renames each root back before lifting it.  The system's own layout — which
   `blocks()`, `summand_x`, the polynomial-reuse cache and every research example
   read — is untouched.  `m = 2` is the identity.
2. **The linear generators.**  When `x_m` is fixed the last link's `n` generators
   drop from degree 2 to degree 1.  The solver only ever acts on tail rows of the
   form `v` or `v + 1`; a linear generator `v + w + …` is otherwise left in the
   system, where the degree-3 Macaulay matrix multiplies it by every monomial of
   degree ≤ 2, and the inherited engine inserts those products one by one against
   the parent's basis (its completion rows, `ReducedBasis::specialise_shared`).
   **Linear elimination** instead Gauss–Jordans the node's linear generators and
   substitutes each pivot `v := Σ rest + c` into the rest — an affine
   substitution is a ring endomorphism of the Boolean ring, so the roots are
   exactly preserved through the definitions — and rebuilds the node's bases
   from the rewritten system, which has lost the intermediate point's variables.

A third, separable question came out of the exploration: whether a child whose
generators merely *drop in degree* should rebuild its bases from its own system
rather than complete the parent's.  That keeps the tree (the rebuilt basis is the
from-scratch step's row space, one end of the sandwich completion lands in, so
the tail and every decision are identical) and is measured as its own factor.

## 0. The boundaries, stated before anything is measured

**The floor.**  Unchanged from [`RESEARCH_INHERITED_F4.md`](RESEARCH_INHERITED_F4.md)
§0: [`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md) §5
bounds the whole attack at `m·2^131` group operations independently of `ℓ`, a
count of candidate tuples with no term for how a Boolean system is split or
reduced.  **No result below can move it**; its ratio column is flat by
construction and the best class available to this round is **engineering**.

**The reference.**  The shipped default on `main` at `9250d0cc`
(`SolverEngine::InheritedF4 { max_degree: 3 }`, `SplitRule::Auto` →
`HighestFree`, layout order, completion on a degree drop, no linear
elimination), selected on the candidate's own binary with
`KIC_F4_DROP=complete KIC_LINEAR_ELIM=0 KIC_CHAIN_ORDER=layout`.  Checked before
registration: that selection reproduces `main`'s frozen ladder to the digit —
word operations, specialisation operations, reductions, propagations, splits,
refutations, F4 calls, oversize and verdict digest on all six rungs.

On this host the frozen ladder's reference total is ~~`16,734,846`~~
**`16,696,846`** word operations (`K_1/2^23`: `7,982,169`, `4,876` reductions;
*arithmetic correction, §5: the six rungs sum to the second figure*).
[`RESEARCH_INHERITED_F4.md`](RESEARCH_INHERITED_F4.md) records `37,302,399` for the
same configuration on the same tree (`4,876` reductions); ~~the rounds merged since
(#620 and the rest) lowered it~~ *(struck, §5: not bisected; the rounds since that
were checked all report their counts unchanged, so which one lowered it is not
known)*.  Every ratio below is against the figure measured in this round, never
the note's.

**The unit.**  64-bit word operations, as the stage has always counted them:
elimination XORs and the inherited engine's specialisation reads and writes.
Linear elimination is charged in the same unit and conservatively: one word
operation per row operation of its Gauss–Jordan (a row is one 64-bit mask) and
one per monomial its substitution writes.  The solver's ordinary substitution
`v := c`, the renaming of a system into the interleaved order and back, and the
from-scratch root build's column indexing are uncharged on every side, as they
always were.  Wall time is a practicality note (`AGENTS.md` §6).

**Why not the WDSat suite.**  `AGENTS.md` §8 names
`research/index_calculus_baseline_20260914/regression/` as the frozen suite.  It
measures WDSat's conflicts on a fixed CNF-XOR corpus; this change is in the
repository's Gröbner splitting solver, which WDSat never calls, so every counter
that suite records is untouched by construction.  The equivalent matched suite
for this stage, under the parent accounting contract
(`research/index_calculus_baseline_20260914/ec_index_calculus_contract.json`:
encoding `chained_S3`, solver `hybrid_guess_and_F4`, counter `gf2_word_xor`), is
the frozen Gröbner-stage ladder with the two chained ladders registered in §2,
and the oracle-pricing ladder of §3.

## 1. What was run before registration

Every run below was made on `2026-09-24`, on scratch code that is not the code
registered here, and is disclosed so that the cells it touched are read as a
tuning set:

- A probe of the tails the inherited engine reads.  On the frozen ladder, a
  non-unit linear row with no unit row appears on `115/550` reads at
  `K_0/2^9, m = 2`, and on `0.2–5.5%` of reads at `n ≥ 13`; on
  `K_1/2^17, m = 3` (4 targets) on `212` of `66,371`.  Using *tail* rows was
  therefore not pursued; this round uses linear *generators*.
- `koblitz_decompose_bench` (the stage bench's target scalars, node budget
  `20,000`) on `(a, n, targets, m)` = `(1, 17, 4, 3)`, `(0, 15, 8, 3)`,
  `(0, 13, 8, 3)`, `(0, 9, 8, 4)`, with the order given by an environment
  permutation, linear elimination uncharged, and children with a linear
  generator not specialised:

  | cell | default | linear elimination | + interleaved order | reductions default → both |
  |:--|--:|--:|--:|:--|
  | `K_1/2^17`, m=3, 4 targets | 202,667,232 | 204,871,791 | 35,239,803 | 66,371 → 2,324 |
  | `K_0/2^15`, m=3, 8 | 4,342,676 | 4,049,014 | 1,439,390 | 3,672 → 136 |
  | `K_0/2^13`, m=3, 8 | 9,743,787 | 8,408,962 | 6,009,659 | 200 → 200 |
  | `K_0/2^9`, m=4, 8 | 2,868,312 | 1,473,191 | 909,554 | 758 → 156 |

  with the same number of targets decomposed in every column (0, 0, 8 and 7).
  Other orders on `(1, 17, 4, 3)`: splitting the summands `x₁` first or the
  intermediate point first hits the node budget; the from-scratch engine
  (`MatrixF4`, highest) costs `4.88 × 10⁹` in the layout order and
  `2.08 × 10⁹` in the interleaved one.  The inherited engine in the interleaved
  order *with* children completed rather than rebuilt cost `1.19 × 10⁹` on the
  same `2,324`-node tree: that is what item 2 of the question removes.
- The frozen ladder with linear elimination (scratch): identical verdict digests
  on every rung, `K_0/2^9, m = 3` `1,142,867 → 888,008` (`→ 611,213` with the
  order, digest changed, the same 15 of 16 decomposed), the quadratic rungs
  `1.00–1.08×`.
- On the registered code, before registration: the reference selection
  reproduces `main` (above); on `(1, 17, 4, 3)` the candidate costs `38,464,003`
  (same `2,324`-node tree; `≈ 3.2 M` of it the charged substitution), and with
  `KIC_F4_DROP=rebuild` on top it cost more (`100.3 M` against `71.5 M` in a
  build that over-charged the substitution identically on both sides).  That is
  why the registered candidate completes on a drop and rebuilding is a control.

## 2. The variants, and the suites they run on

**Factors.**  `O` — chain order: `KIC_CHAIN_ORDER=layout` against the default
(interleaved; applies only when the resolved split rule is `HighestFree` and
`m ≥ 3`).  `L` — linear elimination: `KIC_LINEAR_ELIM=0|1`.  `D` — degree-drop
rule: `KIC_F4_DROP=complete|rebuild`.  The **reference** is `O = layout, L = 0,
D = complete`; the **registered candidate** is `O = interleaved, L = 1,
D = complete` — the engine's new default.  All eight combinations are run and
tabled; every environment variable is set explicitly on every arm, and each run
records them in its `stage.json` (`policy`).

**Suites**, each run three times per arm (the counters are deterministic; the
repetitions exist to show it and to take a wall median):

1. **Frozen ladder** (`groebner_stage_bench`, unchanged): five quadratic rungs
   and `K_0/2^9, m = 3`.  `O` is the identity on the quadratic rungs.
2. **Chain ladder** (`--ladder chain`, new): `K_0/2^9 m=3` (16 targets),
   `K_0/2^13 m=3` (8), `K_0/2^15 m=3` (8), `K_1/2^17 m=3` (4), `K_0/2^9 m=4` (8)
   — the cells of §1, targets `0…`.  The tuning set.
3. **Chain holdout** (`--ladder chain-holdout`, new): `K_1/2^9 m=3` (16),
   `K_0/2^11 m=3` (8), `K_1/2^11 m=3` (8), `K_1/2^13 m=3` (8), `K_1/2^15 m=3` (8),
   `K_0/2^17 m=3` (4), `K_1/2^9 m=4` (8), `K_0/2^15 m=4` (8), `K_1/2^15 m=4` (8),
   and fresh targets `1000…` on `K_0/2^13 m=3` (8) and `K_1/2^17 m=3` (4).  No arm
   of this change has run on any of them.  A cell the harness skips (no curve,
   no factor base, `m` inadmissible for the cofactor, or too many unknowns for a
   64-bit monomial) is listed as skipped, not dropped.

Same-tree comparisons (`D` against its `complete` twin, every other factor
fixed) go through `research/groebner_stage_20260915/compare.py`, which demands
identical verdict digests, reductions, refutations, propagations, splits, F4
calls and oversize.  Every comparison that changes `O` or `L` changes the tree
and goes through `research/inherited_f4_20260922/compare_cross_tree.py`: the same
instances, the same number of targets decomposed, none exhausted on either side,
the same oversize count, word operations reproducible across repetitions — every
decomposition the harness counts having been lifted and verified against the
group identity by `groebner_decompose`.

## 3. The falsification target

The change **succeeds** only if all of the following hold, and is **abandoned**
(not shipped as the default; the negative result recorded) otherwise:

- **T1, confirmatory.**  On the chain holdout, candidate against reference:
  `compare_cross_tree.py` accepts, the total word-operation ratio is **at least
  `2.0×`**, and the ratio exceeds `1.0` on every rung whose reference run splits
  at least once.  The same is reported for the chain ladder, where it is not
  deciding.  A rung on which the *reference* exhausts its node budget on some
  target cannot be compared (the script refuses unequal exhaustion, and rightly:
  the two sides did not answer the same question).  Such a rung is reported with
  both sides' verdicts and counts and left out of T1's ratio, which is computed
  over the comparable rungs; T1 needs at least six comparable holdout rungs.
- **T2, no regression where the change should not act.**  On the frozen ladder,
  candidate against reference: `compare_cross_tree.py` accepts, and no rung's
  word operations rise by more than `2%`.
- **T3, correctness.**  The unit tests
  (`linear_elimination_keeps_every_root`, `interleaved_chain_solve_has_the_same_roots`,
  the pinned same-tree tests) pass; every harness decomposition is verified in
  the group; on the oracle-pricing ladder
  (`ic boundary --regime koblitz --koblitz-degrees 11 --repeats 1 --seed
  123212651130 --no-fold-max-degree 31 --s4-max-degree 31 --oracles
  --oracle-targets 8`, both arms on one binary) the Gröbner oracle's found/refuted
  verdicts are identical between the arms on every cell and it disagrees with
  enumeration on no target.
- **Abandon** also if the chain-holdout total ratio is below `1.2×`, if any
  root is lost in any test, or if any rung of any suite decomposes a different
  number of targets.

**The degree-drop rule** ships as the default only if, same-tree against its
`complete` twin, `D = rebuild` lowers word operations on every rung with a tree
of all three suites, under both `L = 0` and `L = 1`; otherwise it stays a
retained control.

**Inadmissible:** changing any ladder, target, node budget, Macaulay degree or
size cap after this registration; dropping the specialisation or the linear
elimination from the charge; reporting the stage as an attack cost; choosing
among repetitions.

**Class, stated in advance:** at most **engineering** — the floor has no term
this can move, and nothing here changes the first fall degree or the solving
degree of any system.

## 4. Whole logarithms

`AGENTS.md` §8 asks for the complete pipeline, cold.  `ic run --solver groebner
--summands 3 --random-target --batch 1` runs one: factor base, targets,
Gröbner decomposition of every trial, lifting and verification of every
relation, the relation matrix and the final scalar, checked as `[k]G = Q` against
the planted `k`.  Both arms run on the same binary with the same `--degree`,
`--curve-a`, `--seed` and factor base, so they draw the same trial points; only
the oracle's tree differs.  The cells were sized by a reference-only
feasibility run (seed 1, factor index 0), before registration: of the Koblitz
curves `K_a/F_{2^n}`, `n ∈ {9, 11, 13, 15, 17}`, only `K_0/2^9` (3 columns,
4 trials, `0.01 s`) and `K_0/2^13` (77 columns, 14 trials, `1.6 × 10⁷` word
operations, `0.27 s`) run a logarithm with `m = 3` at all.  At `K_1/2^9`,
`K_1/2^11`, `K_0/2^15`, `K_1/2^15` and `K_1/2^17` the pipeline draws no trial
(the base admits no `m = 3` relation), and `K_0/2^11`, `K_1/2^13` and
`K_0/2^17` have no usable prime-order subgroup in the constructor.  The
registered cells are therefore `K_0/2^13`, seeds `1…10` with holdout seeds
`101…105`, and `K_0/2^9`, seeds `1…5`, as a small-instance control.  Reported per
run: status, recovered and planted `k`, trials, relations, oracle reductions and
word operations, and each phase's wall time.

What this can and cannot say is fixed now.  `ic run` counts the oracle in word
operations and the rest of the pipeline in wall time only; no measured
conversion between them is recorded for it, so the method's `S` stays **null**
(`AGENTS.md` §8: leave a missing conversion null).  If the trial and relation
counts are identical between the arms, every phase but the oracle does identical
work, and `baseline_total / candidate_total` is bounded by the oracle's ratio;
that inference is stated as such, beside the measured whole-process wall ratio,
which is a practicality note.  The per-target oracle price in group-addition
equivalents, and the relation phase it would imply, come from the oracle-pricing
ladder of T3, which converts at its host's measured factor and marks its
projection as an extrapolation.

## 5. Amendment, 2026-09-24: T1 as registered has four comparable rungs, not six

*Written after the registered runs of §2 (committed as `d0b048ac`) and before
any run of the supplementary holdout below.*

**What happened.**  Seven of the eleven holdout cells of §2 were skipped by the
harness.  `K_0/2^11`, `K_1/2^13` and `K_0/2^17` have no curve in the
constructor (no usable prime-order subgroup).  At `K_1/2^9`, `K_1/2^11` and
`K_1/2^15` with `m = 3`, and at `K_1/2^17` with fresh targets, the factor base
admits no `m = 3` relation for its cofactor
(`FrobeniusFactorBase::m_can_decompose`), so the decomposition oracle never
enters the Gröbner stage there.  The list was written without either check.  Four
holdout rungs ran — `K_1/2^9 m=4`, `K_0/2^15 m=4`, `K_1/2^15 m=4` and fresh
targets on `K_0/2^13 m=3` — and T1 needs six.  **T1 is not met as registered**,
on its count clause alone; it is recorded that way and not re-read.

The same check removes two cells from the chain (tuning) ladder: `K_0/2^15`
and `K_1/2^17` with `m = 3` are inadmissible.  Two of §1's exploratory rows —
including its largest ratio, `5.75×` on `K_1/2^17, m = 3` — therefore describe
refutations the pipeline never asks for.  They stay in §1 as what was run; they
are not evidence for this change.

**T1′, the supplementary holdout.**  `examples/chain_ladder_screen.rs` screens
the grid `n ∈ {9, 11, 13, 15, 17, 19, 23}`, `a ∈ {0, 1}`, factor index `0…3`,
`m ∈ {3, 4}` for a curve and a factor base that exist, admissibility for `m`,
and at most 64 unknowns; it decides no target
(`research/chain_split_order_20260924/screen.json`).  Fifteen cells pass.  T1′
is **every one of them on which no arm of this change has run** — nine cells,
no other selection:

| cell | factor index | `ℓ` | base points | unknowns | targets |
|:--|--:|--:|--:|--:|--:|
| `K_1/2^11`, m=4 | 0 | 10 | 991 | 62 | 4 |
| `K_0/2^15`, m=3 | 1 | 4 | 31 | 27 | 8 |
| `K_0/2^15`, m=4 | 1 | 4 | 31 | 46 | 8 |
| `K_0/2^15`, m=3 | 2 | 4 | 21 | 27 | 8 |
| `K_0/2^15`, m=4 | 2 | 4 | 21 | 46 | 8 |
| `K_1/2^15`, m=4 | 1 | 4 | 1 | 46 | 8 |
| `K_1/2^15`, m=4 | 2 | 4 | 11 | 46 | 8 |
| `K_0/2^23`, m=3 | 0 | 11 | 2,025 | 56 | 4 |
| `K_0/2^23`, m=3 | 1 | 11 | 2,071 | 56 | 4 |

(`groebner_stage_bench --ladder chain-holdout-2`, targets `0…`, node budget
`20,000`, four targets where a cell has 56 or more unknowns.)  Arms: the
reference and the registered candidate only, three repetitions each
(`run_t1prime.sh`).  **T1′ holds** if `compare_cross_tree.py` accepts, the total
word-operation ratio is at least `2.0×`, the ratio exceeds `1.0` on every rung
whose reference run splits, and at least six rungs are comparable (a rung whose
reference exhausts its budget is reported and left out, as in T1).  The
candidate ships as the default only if T1′, T2 and T3 hold; otherwise the
engine's default reverts to the reference and the rest stays as retained
controls.
