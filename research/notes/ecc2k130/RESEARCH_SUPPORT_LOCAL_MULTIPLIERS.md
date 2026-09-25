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

**Status: registered, measured (§5), amended (§6): T1, T2 and T3 met; engineering.**
On the registered code the whole logarithms ran slower in wall time
(`0.82×`) while their oracle word operations fell `2.19×`, because the row
packer's hash cost more than the elimination saved.  The unit does not charge
that cost.  §6 fixes the hash with every counter unchanged; the whole
logarithms then run `1.41×` faster.  §0–§4 were committed before any
registered run (`1f0d9751`); §1 lists every run made before registration.
Results are appended below §4 and do not edit it.

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

## 5. Results

*Appended after every registered run; §0–§4 are not edited.  Every number below
is read from `research/support_local_multipliers_20260924/`: `tables.md` for the
stage ladders, `comparisons/` for the acceptance of each pair, `oracle_ladder/`
for T3, and `e2e/` with `e2e_tables.md` for §4.  The registered runs were made
on binaries built from `1f0d9751`.  Word operations are deterministic and
identical across the three repetitions of every arm.*

**Bottom line.**  On the registered holdout, support-local multipliers cost
**`1.50×`** fewer word operations than the reference: 19 comparable rungs of
19, `1.40–3.57×` per rung, none rising, and the same targets decomposed on
every rung.  The other suites:

- the chained ladders, all tuning data: `1.60–2.25×`;
- the frozen ladder: `1.02×`, where only its cubic rung moves (`2.01×`) and two
  quadratic rungs rise by `0.12%` and `0.05%`;
- the oracle ladder: `1.52×` at `n = 15, m = 3` and `2.06×` at `n = 9, m = 3`,
  with identical verdicts.

T1, T2 and T3 are met, so the candidate is the engine's default.  Class:
**engineering**.

The unit does not show the whole of it.  On the registered code the whole
logarithms did `2.19×` less oracle work and ran **slower**: `0.82×` in wall
time on `K_0/2^13`, confirmed with the arms interleaved (§5.4).  The row
building that both arms leave uncharged (§0) had grown in the candidate by
more than the elimination it saved.  In AGENTS.md §3's terms that is the
relabelling pattern, seen in wall time rather than in `S`, which is null here.
§6 finds the cause, a hash that does not suit monomial masks, and removes it
without changing a counter.

### 5.1 The table

Candidate against reference; ratios are reference / candidate; "nodes" are
the reductions each tree took.  The ratio to the floor is flat on every row by
construction (§0).  Correctness: every pair was accepted by
`compare_cross_tree.py`.  That means the same instances and the same targets
decomposed, none exhausted on either side, and the same oversize count, and
every counted decomposition was lifted and verified in the group.

| suite | rungs | targets | decomposed ref / cand | nodes ref → cand | reference | **candidate** | per rung | class |
|:--|--:|--:|:--|:--|--:|--:|:--|:--|
| frozen (5 quadratic rungs, 1 cubic) | 6 | 176 | 110 / 110 | 7,701 → 7,704 | 16,125,634 | **15,816,123 (1.02×)** | 0.9988–2.01× | engineering |
| chain (tuning) | 3 | 32 | 30 / 30 | 571 → 573 | 7,547,415 | **3,351,378 (2.25×)** | 2.01–3.12× | engineering |
| chain-holdout (tuning) | 4 | 32 | 21 / 21 | 8,496 → 8,496 | 52,305,605 | **32,779,266 (1.60×)** | 1.43–3.30× | engineering |
| chain-holdout-2, the predecessor's T1′ (tuning) | 9 | 60 | 25 / 25 | 24,181 → 24,189 | 225,078,474 | **136,355,549 (1.65×)** | 1.46–3.55× | engineering |
| **R2 holdout** (registered here) | 19 | 148 | 65 / 65 | 39,156 → 39,168 | 888,630,176 | **591,938,031 (1.50×)** | 1.40–3.57× | engineering |

The R2 holdout rung by rung is in `tables.md`:

- The four cells no ladder had used: `K_0/2^31, m = 3` at factor indices 0, 1
  and 5 run `1.45×`, `1.43×` and `1.44×` (about `1.9 × 10⁸` word operations
  each in the reference), and `K_1/2^29, m = 4` runs `1.63×`.
- The fresh targets on the fifteen known cells run between `1.40×`
  (`K_1/2^15 m=4`) and `3.57×` (`K_1/2^11 m=4`).

The tree grows by at most ten nodes on any rung, as §0 predicts: the tails are
weaker, and the tree can only grow.

Where the saving comes from: the roots, which §1 found cost `43–88%` of the
chained cells, and the specialisation that inherits their rows.  The
specialisation's reads and writes fall `114.5 M → 68.2 M` over the R2
holdout.  The largest gains, `2.5–3.6×`, are on the `m = 4` rungs with trees
of 124–264 nodes, where the roots are most of the cost.  `K_1/2^29 m=4`, with
16 nodes and 62 unknowns, gains `1.63×`.

### 5.2 Against the registered targets

- **T1: met.**  Nineteen comparable rungs of nineteen (fifteen required).  The
  total is `1.50×` against a threshold of `1.3×`.  No rung rises; the least
  gain is `1.396×`, on `K_1/2^15 m=4`, targets `5000…`.  Nothing is exhausted
  on either side, and the same targets are decomposed on every rung.
- **T2: met.**  The frozen pair is accepted.  Two rungs rise, both by far less
  than the `2%` allowed: `K_1/2^17 m=2` by `0.12%` (1,982,352 → 1,984,783) and
  `K_1/2^23 m=2` by `0.05%` (7,974,980 → 7,979,328).  Both trees have the
  same size, and at `m = 2` only completion rows differ (§0), so the rise is
  theirs.
- **T3: met.**  All 153 tests of the three modules pass, among them
  `support_local_rows_are_occurring_rows` (a strict subset at `m = 3, 4`,
  equal at `m = 2`) and the root-preservation tests over every policy with
  both multiplier choices.  On the oracle-pricing ladder both arms report
  identical found / refuted / inconclusive counts for matrix-F4 on every cell,
  and no two oracles disagree on any target (§5.3).
- **Abandon conditions: none met.**  The R2 total is `1.50×` against an
  abandon floor of `1.1×`.  No test loses a root, and no rung of any suite
  decomposes a different number of targets.

### 5.3 The oracle-pricing ladder (T3)

`ic boundary --regime koblitz --oracles`, seed `123212651130`, 8 targets per
cell, both arms on one binary (`oracle_ladder/{reference,candidate}.json`).
The GAE figures are converted at each run's own host-measured factor, so the
word-XOR ratio is the comparison.  The projected `S` is the ladder's own
extrapolation.

| cell | matrix-F4 found / refuted | word XORs, reference | candidate | ratio | GAE/target ref → cand | projected `S` ref → cand |
|:--|:--|--:|--:|--:|:--|:--|
| `n = 9`, m=3 | 8 / 0 | 354,877 | 172,567 | **2.06×** | 198 → 91 | 88 → 41 |
| `n = 15`, m=3 | 3 / 5 | 7,187,062 | 4,736,391 | **1.52×** | 2,256 → 1,488 | 878 → 579 |
| `n = 9, 11, 13, 15, 17, 23`, m=2 | identical | | | 1.00–1.03×, none rising | | |

The reference's `7,187,062` at `n = 15, m = 3` is, to the word, the
predecessor's candidate figure (its §6.4).  After this round matrix-F4 there
is still about `3×` enumeration (467 GAE per target) and `60×` meet in the
middle (24.7).  The Gröbner oracle remains the slowest route that finishes,
and the page's verdict is unchanged.

### 5.4 Whole logarithms (§4), on the registered code

`ic run --solver groebner --summands 3 --random-target --batch 1`, both arms on
one binary, the same seeds (`e2e/`).  All 40 runs **complete, with the planted
`k` recovered and `[k]G = Q` verified**.  Trials, relations and matrix are equal
on every seed.

| cell | runs | trials = relations | oracle word ops ref → cand | ratio (per run) | wall ref → cand | wall ratio, geometric mean [95% paired bootstrap] |
|:--|--:|--:|:--|:--|:--|:--|
| `K_0/2^13`, seeds 201–210, holdout 301–305 | 15 + 15 | 141 | 104,865,593 → 47,966,249 | **2.19×** (2.15–2.23) | 2.07 → 2.61 s | **0.82×** [0.77, 0.87] |
| the same, rerun with the arms interleaved, 3 repetitions (§6) | 45 + 45 | 141 | the same | the same | 1.89 → 2.33 s | **0.82×** [0.80, 0.84] |
| `K_0/2^9`, seeds 201–205 | 5 + 5 | 17 | 750,778 → 361,523 | **2.08×** (2.02–2.16) | 0.05 → 0.05 s | 0.96× [0.79, 1.22] |

By the registered rule of §4, `baseline_total / candidate_total > 1` on every
run: every counted phase but the oracle does equal work, and the oracle does
less.  The wall contradicts it.  The candidate is `18%` slower on
`K_0/2^13`, and the interleaved rerun rules out host drift.  So some work
outside the count grew by more than the count fell.  The stage ladders show
the same thing more weakly: their wall ratios are `0.90–1.03×` against
`1.02–2.25×` in word operations (`postfix_identity.md`, registered binary).

## 6. Amendment, 2026-09-24: the row packer's hash

*Made after every registered run of §5 and after §5.4's wall result.  The
diagnosis and the fix came first.  Then came one exploratory check on the
fixed binary: `K_0/2^13` seed 209, five interleaved pairs, the candidate
`0.35 → 0.19 s`, the reference unchanged at `0.27 s`, word operations
identical.  The reruns below were then fixed in `run_postfix.sh` before any of
them ran.  The amendment changes no ladder, target, budget, degree, cap or
counter definition.*

**What the profile found.**  I ran callgrind on one whole logarithm
(`K_0/2^13`, seed 210, both arms; `profile/`):

| | reference | candidate, registered | candidate, fixed |
|:--|--:|--:|--:|
| whole run, instructions | 417.3 M | 456.1 M | 312.8 M |
| root build (`ReducedBasis::from_system_with`) | 199.8 M | 253.6 M | 110.2 M |
| elimination within it (`echelon_f2_suffix_counted`) | 87.1 M | 39.4 M | 39.4 M |
| `pack_rows` | — | 151.4 M | 8.1 M |

The candidate halves the elimination, but its roots are built by
`build_inherited_macaulay_support_local`, which packs rows through
`pack_rows`.  That function indexes columns in an `FxMap<u64, usize>`.
`FxHasher` on one word is `x · K`, and a product keeps the trailing zeros of
`x`.  hashbrown takes the bucket from the low bits.  In the interleaved order
the first summand's variables are the lowest bits, and the tail links'
monomials do not contain them, so those monomials share a few buckets and
every lookup probes a long run.  The reference's roots never went that way:
`build_inherited_macaulay` packs through `F4ColumnLayout`, whose index was
already hashed with a splitmix64 finaliser.  The same `FxMap` also served the
inherited engine's completion-row index in both arms.

**The fix** (`d69936e1`).  `F4ColumnLayout`'s splitmix hasher moves into
`fx_hash` as `MaskHasher` / `MaskMap`.  `pack_rows` and the inherited engine's
column index use it.  Every map it touches is lookup-only, so every matrix,
every pivot and every counter is unchanged; only uncharged row-packing time
moves.

**Admissibility, checked** (`run_postfix.sh`, `check_identity.py`,
`postfix_identity.md`).  All five suites were rerun with four arms,
{binary of `1f0d9751`, binary of `d69936e1`} × {reference, candidate}, with the
arm order rotating across three repetitions.  **Every field of every rung of
all 60 runs**, timings aside, equals the registered run of the same arm.  The
same holds for the whole logarithms: trials, relations and oracle word
operations on every seed, both binaries (`e2e_table.py` asserts it).  So
§5.1–§5.3 stand for the shipped code as they are.  The oracle ladder was not
rerun: its verdicts and word XORs cannot depend on a lookup-only hash, and the
same engine's counters were just shown identical on 60 ladder runs.

**Wall time after the fix**, a practicality note: interleaved, medians of
three repetitions, the same host.

| | frozen | chain | chain-holdout | chain-holdout-2 | R2 holdout | whole logs `K_0/2^13` | whole logs `K_0/2^9` |
|:--|--:|--:|--:|--:|--:|--:|--:|
| reference / candidate, registered binary | 0.96× | 0.90× | 1.03× | 0.98× | 1.01× | 0.82× [0.80, 0.84] | 1.05× [1.00, 1.10] |
| reference / candidate, fixed binary | 1.00× | 1.49× | 1.29× | 1.21× | 1.08× | **1.41× [1.37, 1.44]** | **1.24× [1.18, 1.30]** |
| hash fix alone, reference arm | 1.13× | 0.94× | 1.02× | 1.06× | 1.06× | 1.01× [0.99, 1.03] | 1.03× [0.96, 1.11] |
| hash fix alone, candidate arm | 1.17× | 1.56× | 1.28× | 1.31× | 1.12× | 1.73× [1.67, 1.79] | 1.22× [1.11, 1.35] |
| **round as shipped**: main's default → fixed candidate | | | | | | **1.42× [1.38, 1.46]** | **1.28× [1.20, 1.39]** |

Whole-log intervals are 95% paired bootstraps over seeds (`e2e_tables.md`).
The stage-ladder cells are ratios of suite medians.

What this means:

- The fix makes neither arm slower anywhere it matters.  The one place it
  could have is the reference's frozen ladder: its quadratic rungs are where
  Fx's low bits were adequate, and splitmix costs a few more cycles per
  lookup.  That ladder runs `1.13×` faster.  The `0.94×` on the chain ladder
  is a `0.16 → 0.17 s` suite.
- The whole logarithm's `18%` loss becomes a `41%` gain against the same
  reference.
- The stage ladders still gain much less in wall time than in word operations
  (R2: `1.08×` against `1.50×`).  Where the rest of their wall time goes is not
  profiled here; §7 lists it as open.

**Class.**  The support-local change is **engineering** in the unit, as
registered, and stays engineering.  The hash fix is **engineering** too, on a
cost the unit does not see: no counter moves, so it claims nothing in word
operations.  The episode is recorded because it is the failure AGENTS.md §3
describes.  The headline count fell `2.19×` while the uncounted cost rose by
more.  Only the wall time, kept as a practicality note, showed it.

## 7. What this does not establish, and what is next

- **An advance.**  The floor has no term for which rows a root is built from,
  or for a hash.  The ratio to it is flat on every row.  No first fall degree,
  solving degree or fitted exponent moves.
- **An end-to-end speedup in one unit.**  `ic run` counts the oracle in word
  operations and everything else in wall time only, so `S` stays null (§4).
  The whole-log wall ratios above are paired, interleaved and bootstrapped,
  and they remain a practicality note on one shared host.
- **Anything about `m = 2`.**  Quadratic roots are built identically.  The
  quadratic rungs move `0.9988–1.04×` through completion rows alone.
- **That row building is priced.**  It never was, in either arm or in the
  predecessor.  This round is the evidence that it should be: an uncharged
  phase can grow and hide a regression.  Charging it (words written by
  packing, per row) is an accounting change for the next round, to be made
  before its measurements.  It would move every baseline, so it gets its
  own note.
- **Next levers,** all measured in §1 or §6, none claimed:
  - re-reduction during specialisation on the heaviest `m = 4` cells
    (`40.7 M` of `64.8 M` on `K_0/2^15 m=4 divisor 2` in §1's probe);
  - the quadratic roots, which §1 found cost `43–88%` of the `m = 2` rungs;
  - the gap on the R2 holdout between word operations (`1.50×`) and wall
    time (`1.08×`), not yet profiled.

### 7.1 Reproducing

```bash
cargo build --release --example groebner_stage_bench --example chain_ladder_screen --bin ic
./target/release/examples/chain_ladder_screen --wide > research/support_local_multipliers_20260924/screen_wide.json
research/support_local_multipliers_20260924/run.sh           # five suites, both arms, 3 reps
research/support_local_multipliers_20260924/oracle_ladder.sh # T3
research/support_local_multipliers_20260924/e2e.sh           # §4
python3 research/support_local_multipliers_20260924/table.py > research/support_local_multipliers_20260924/tables.md
# §6: the registered commit rebuilt beside the fixed one
git worktree add ../wt-1f0d 1f0d9751 && (cd ../wt-1f0d && CARGO_TARGET_DIR=../target-1f0d cargo build --release --example groebner_stage_bench --bin ic)
OLD=../target-1f0d/release research/support_local_multipliers_20260924/run_postfix.sh
python3 research/support_local_multipliers_20260924/check_identity.py > research/support_local_multipliers_20260924/postfix_identity.md
python3 research/support_local_multipliers_20260924/e2e_table.py > research/support_local_multipliers_20260924/e2e_tables.md
cargo test --release --lib -- cryptanalysis::koblitz_groebner cryptanalysis::inherited_f4 cryptanalysis::fx_hash
```

The retained control is `KIC_F4_MULTIPLIERS=occurring`, which reproduces the
predecessor's default to the digit.
