# Support-local multipliers for the inherited engine

**Modules:** `src/cryptanalysis/inherited_f4.rs` (`ReducedBasis::from_system_with`,
the support-local completion rows), `src/cryptanalysis/koblitz_groebner.rs`
(`build_inherited_macaulay_support_local`, `InheritPolicy::support_local`).
**Harness:** `examples/groebner_stage_bench.rs` (five ladders, one of them new
here), `examples/chain_ladder_screen.rs --wide`, `ic boundary --oracles`,
`ic run --solver groebner --summands 3`.
**Frozen evidence:** `research/support_local_multipliers_20260924/`.
**Predecessor:** [`RESEARCH_CHAIN_SPLIT_ORDER.md`](RESEARCH_CHAIN_SPLIT_ORDER.md),
whose default this round measures against and changes in one place only.

**Status: registered, not yet measured.**  §0–§4 were committed before any
registered run; §1 lists every run made before registration.  Results are
appended below §4 and do not edit it.

## The question

The inherited engine builds a node's Macaulay matrix from scratch at a root
(and, since the previous round, after linear elimination rewrites a node's
system), and it inserts completion rows when a generator's degree drops.  In
both places it multiplies **every** generator by every monomial of the allowed
degree over **every variable occurring in the system** — the from-scratch
step's active-multiplier policy.

For the chained decomposition systems that is mostly wasted.  In the
interleaved order the chain's links are block-multilinear and the tree fixes
the last summand first, so the root's rows `x·L` that multiply the last link `L`
by a variable of the first summands contribute monomials that nothing the tail
reads can use.  **Support-local multipliers** multiply each generator `f` only
by monomials in the variables of its own support: the rows
`{t·f : t ⊆ supp(f), deg(t·f) ≤ D}`.  That is a subset of the occurring-variable
rows, so every consequence is still in the ideal and the splitting solver still
decides every target exactly; the tail it reads at a node can only be weaker,
and the tree can only grow, never lose a root.  Where every generator already
spans every occurring variable — every quadratic `m = 2` root — the rows are
identical and the same code path (with its layout cache) builds them.

## 0. The boundaries, stated before anything is measured

**The floor.**  Unchanged: `m·2^131` group operations
([`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md) §5),
no term for how a Macaulay matrix is built.  Ratio to it flat by construction;
the best available class is **engineering**.

**The reference.**  `main`'s default at `9ca0b991` — the candidate merged in
#690 (interleaved chain order, linear elimination, completion on a drop,
occurring-variable multipliers) — selected on the candidate's own binary with
`KIC_F4_MULTIPLIERS=occurring`.  Checked before registration: that selection
reproduces the round-1 candidate's frozen runs
(`research/chain_split_order_20260924/{frozen,chain,chain-holdout,chain-holdout-2}/candidate/rep1`)
to the digit — word operations, specialisation operations, reductions,
propagations, splits, refutations, eliminations and verdict digest on every
rung of all four ladders.

**The unit.**  64-bit word operations, as in the predecessor: elimination XORs,
specialisation reads and writes, and linear elimination's row operations and
written monomials.  Building a row is uncharged on both sides, as it always
was.  Wall time is a practicality note.

**Why not the WDSat suite.**  As in the predecessor (§0 there): the change is in
the Gröbner splitting solver, which WDSat never calls; the matched suite is the
Gröbner-stage ladders below and the oracle-pricing ladder, under the parent
accounting contract (encoding `chained_S3`, solver `hybrid_guess_and_F4`,
counter `gf2_word_xor`).

## 1. What was run before registration

All on `2026-09-24`, on scratch code (environment switches not in this
commit), disclosed so that the cells they touched are read as tuning data:

- A cost probe of the round-1 default on the four existing ladders.  The roots
  cost `43–88%` of the quadratic rungs and of `K_0/2^13, m = 3` and
  `K_0/2^23, m = 3`; on the heaviest `m = 4` cells re-reduction during
  specialisation dominates (`K_0/2^15 m=4 divisor 2`: `40.7 M` of `64.8 M`).
- The Macaulay degree cap on the T1′ ladder: degree 2 costs about what degree 3
  does (the tree grows `2,376 → 10,493` F4 calls on `K_0/2^23 m=3`), degree 4
  costs `20–200×` more.  Not pursued.
- On `koblitz_decompose_bench 0 23 4 3` (= the T1′ cell `K_0/2^23 m=3`): the
  four initial roots (56 unknowns, rank 1,333, 13,884 columns) cost
  `6.9–7.3 M` each, `73%` of the cell.  Excluding the first two summands'
  variables from the root's multipliers: `38,953,101 → 19,696,631` on the same
  tree (2,380 reductions, identical decompositions); excluding only the first
  summand's: `→ 29,134,985`.
- Support-local multipliers at the roots only, then at the roots and in
  completion rows, on all four ladders (candidate / round-1 default, verdict
  digests identical on every rung of every ladder):

  | ladder | roots only | roots + completion |
  |:--|--:|--:|
  | frozen | 1.02× | 1.02× |
  | chain | 2.25× | 2.25× |
  | chain-holdout | 1.27× | 1.60× |
  | chain-holdout-2 (T1′) | 1.37× | 1.65× |

  The registered candidate is roots + completion.  On the registered code, the
  same four ladders reproduce the right-hand column exactly.

## 2. The variants, and the suites

**Arms.**  Reference: `KIC_F4_MULTIPLIERS=occurring`.  Candidate:
`KIC_F4_MULTIPLIERS=support` (the new default).  Both with
`KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_DROP=complete` set
explicitly; each run records its policy.  Three repetitions each
(`run.sh`).

**Suites.**  The frozen ladder and the three chained ladders of the
predecessor — all tuning data now (§1) — and one new ladder:

- **R2 holdout** (`--ladder r2-holdout`), fixed here before any arm runs on it.
  `chain_ladder_screen --wide` screens `n ∈ {9, 11, …, 31}`, `a ∈ {0, 1}`,
  factor index `0…7`, `m ∈ {3, 4}` for existence, admissibility and at most 64
  unknowns, deciding no target (`screen_wide.json`).  It finds nineteen
  admissible cells: the fifteen of the predecessor's screen, and four no
  ladder has used — `K_0/2^31, m = 3` at factor indices 0, 1 and 5 (`ℓ = 5`,
  63 base points, 46 unknowns) and `K_1/2^29, m = 4` (`ℓ = 1`, one base point,
  62 unknowns).  The R2 holdout is those four cells with targets `0…`, and
  **fresh targets `5000…`** on each of the fifteen: 19 rungs, 8 targets each
  (16 on `K_0/2^9 m=3`, 4 where a cell has 56 or more unknowns).

## 3. The falsification target

The candidate becomes the default only if all of these hold; otherwise the
default stays the reference and the policy is a retained control:

- **T1, confirmatory.**  On the R2 holdout, `compare_cross_tree.py` accepts
  (same instances, same targets decomposed, none exhausted on either side,
  same oversize), the total word-operation ratio reference / candidate is **at
  least `1.3×`**, and no rung's word operations rise by more than `2%`.  At least
  fifteen comparable rungs (a rung on which the reference exhausts its budget
  is reported and left out).
- **T2.**  On the frozen ladder, the cross-tree script accepts and no rung's
  word operations rise by more than `2%`.
- **T3.**  The module tests pass (including
  `support_local_rows_are_occurring_rows` and the root-preservation tests,
  which now run every policy with both multiplier choices); on the
  oracle-pricing ladder both arms give identical found/refuted verdicts for
  matrix-F4 on every cell and no oracle disagrees with another on any target.
- **Abandon** if the R2 total is below `1.1×`, if any test loses a root, or if
  any rung of any suite decomposes a different number of targets.

Inadmissible: changing any ladder, target, budget, degree or cap after this
commit; charging the candidate less than the reference for the same work;
choosing among repetitions.  **Class, stated in advance:** at most
engineering.

## 4. Whole logarithms

`ic run --solver groebner --summands 3 --random-target --batch 1`, both arms on
one binary, on the two cells the predecessor sized (`K_0/2^13`, `K_0/2^9`), with
seeds no earlier run used: `201…210` and holdout `301…305` on `K_0/2^13`,
`201…205` on `K_0/2^9` (`e2e.sh`).  Reported as in the predecessor: status,
recovered and planted `k`, trials, relations, oracle word operations, wall.  The
method's `S` stays null for the reason given there (no measured conversion
between the oracle's unit and the rest of `ic run`), and any statement about
the whole pipeline is limited to the sign of `baseline_total /
candidate_total` when every other counted phase does equal work, with the wall
ratio as a practicality note.
