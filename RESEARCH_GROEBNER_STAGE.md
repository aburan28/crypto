# Optimising the Gröbner stage of the ECC2K-130 decomposition oracle

**Module:** `src/cryptanalysis/koblitz_groebner.rs`
**Harness:** `cargo run --release --example groebner_stage_bench -- --out DIR`
**Frozen evidence:** `research/groebner_stage_20260915/`
**Follow-on:** [`RESEARCH_FACTOR_BASE_SOLVE_COST.md`](RESEARCH_FACTOR_BASE_SOLVE_COST.md)
— the same profiler turned on the factor base, asking which subspace the stage
is cheapest over.
**Related:** [`RESEARCH_KOBLITZ_INDEX_CALCULUS.md`](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)
(the pipeline this stage sits in), [`RESEARCH_GROEBNER_F4.md`](RESEARCH_GROEBNER_F4.md)
(the engine), [`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md)
(what the oracle would have to cost at `n = 131`, and why it cannot).

**The question.**  The decomposition oracle asks *is `R = P_1 + … + P_m` with
every `P_i` in the Frobenius-invariant factor base?* and answers it by
Weil-restricting the Semaev condition to a Boolean system and reducing it with
matrix-F4 plus splitting.  On the Koblitz note's own table that oracle is the
slowest of the three (`206 s` against enumeration's `4.3 s` at `K_1/2^23` before
the projection work, `20.6 s` after).  **How much of that is the algebra, and
how much is the implementation?**

**Bottom line.  Two thirds of it was implementation.**  The stage now costs
`1.60×` fewer counted word XORs and runs `2.84×` faster across the frozen
ladder, deciding every target identically.  That is an **engineering** result in
the sense of `AGENTS.md` §3 and nothing more: the ratio to the floor below is
unchanged, because the floor does not contain a constant this work could move.

## 0. The boundaries, stated before anything is measured

**The floor.**  `RESEARCH_ECC2K130_DECOMPOSITION.md` §5 derives, for an
`m`-summand subspace factor base on the challenge curve, that the whole attack
costs `m·2^131` group operations *independently of `l`*, and §5.3 restates it as
a demand on the oracle: it would have to beat exhaustive search over its own
candidate set by `2^{70.19 + log₂ m}`.  That bound counts candidate tuples.  It
contains no term for how fast one Macaulay matrix reduces, so **no result in
this note can move it**, and the ratio to it is flat by construction.  Stating
that first is the point: the work below is worth doing because the stage is what
every toy-size experiment in this repository waits on, not because it bears on
`n = 131`.

**The reference.**  The frozen baseline measurement in
`research/groebner_stage_20260915/baseline/rep{1,2,3}/stage.json`: the same
harness, the same ladder, the same targets, run against the working tree as it
stood before this round.  `2,475,304,254` word XORs and `11.5 s` of stage time
over 176 targets.

**The unit.**  64-bit word XORs performed by the Macaulay elimination — the unit
`crossbred_bench` and `blocked_macaulay_bench` already use.  Wall time is
reported beside it as a practicality note and never as the metric (§6).

**Falsification target.**  A variant is an improvement only if it (a) decides
every target of the frozen ladder identically — same verdict digest, same
reductions, infeasibility certificates, propagations and splits — and (b)
reduces the counted word XORs on **every** rung, not merely in total, and not
merely in wall time.  Anything that lowers wall time while raising the counted
unit is a relabelling and is recorded as one.  Inadmissible: changing the
ladder, the targets, the node budget or the size caps; comparing medians of
different repetition counts; reporting the stage as an attack cost.

## 1. Where the stage's time actually went

The harness reports the stage in three phases, because `word_ops` only ever
counted the middle one:

| phase | share of baseline stage wall | counted by `word_ops`? |
|:--|--:|:--|
| build the Macaulay matrix | 30.2% | no |
| reduce it | 39.8% | yes |
| read the reduced rows back | 27.3% | no |

**That is the first finding, and it is about the accounting, not the algebra.**
The unit this repository prices the stage in covers two fifths of the stage.  A
round that optimised only what the unit could see would have been capped at
`1.7×` before it started, and the `27.3%` spent turning bit-rows back into
polynomials would never have appeared in any table.

## 2. The table

One unit, one table, every variant a row.  Ratios are baseline / candidate, so
above `1` is an improvement.  Medians of three repetitions; word XORs are
deterministic and identical across all three on both sides.

| instance | `m` | `ℓ` | targets | decomposed | word XORs baseline → candidate | ratio | ratio to floor | wall | correct |
|:--|--:|--:|--:|--:|--:|--:|:--|--:|:--|
| `K_0/2^9`  | 2 | 6 | 40 | 39 | 1,958,318 → 1,625,348 | **1.20×** | flat | 1.82× | ✓ identical |
| `K_0/2^9`  | 3 | 6 | 16 | 15 | 36,193,798 → 32,991,833 | **1.10×** | flat | 2.22× | ✓ identical |
| `K_0/2^13` | 2 | 12 | 40 | 40 | 108,572,993 → 79,728,473 | **1.36×** | flat | 2.56× | ✓ identical |
| `K_1/2^15` | 2 | 4 | 32 | 0 | 135,838 → 88,520 | **1.53×** | flat | 1.24× | ✓ identical |
| `K_1/2^17` | 2 | 8 | 32 | 7 | 139,256,561 → 91,932,686 | **1.51×** | flat | 2.56× | ✓ identical |
| `K_1/2^23` | 2 | 11 | 16 | 9 | 2,189,186,746 → 1,343,398,546 | **1.63×** | flat | 2.95× | ✓ identical |
| **total** | | | 176 | 110 | 2,475,304,254 → 1,549,765,406 | **1.60×** | flat | **2.84×** | ✓ |

**Class: engineering.**  `S` fell; the ratio to the floor did not move, and
could not have.

*Correctness* is not a summary judgement: `compare.py` refuses the comparison
unless the verdict digest over every target and every solver counter agree rung
for rung, and it did not refuse.  `K_1/2^15` decomposes nothing — its subspace
admits no `m = 2` relation — so it prices refutation alone, which is most of an
attack's work and the case the algebra is supposed to be for.

### Holdout

Two instances the tuning never saw, including a different factor-base divisor at
`n = 17` and the widest matrices the oracle reaches here (`K_0/2^19`, a few
hundred rows against a couple of thousand columns).  Only the elimination kernel
is switchable at run time (`F4_F2_RREF=ref`), so this A/B isolates that one
change rather than the whole round:

| instance | word XORs | wall | correct |
|:--|--:|--:|:--|
| `K_0/2^19`, `m = 2` | **1.58×** | 1.70× | ✓ identical |
| `K_1/2^17`, `m = 2`, divisor 1 | **1.51×** | 1.43× | ✓ identical |

## 3. The engineering ledger

Each row is cumulative on the one above, measured on the same ladder; single
repetitions, kept in `research/groebner_stage_20260915/steps/`.  `n = 23` wall
is broken out because it is the rung that dominates.

| step | what changed | total word XORs | ratio | `n=23` wall | build / reduce / read (ms) |
|:--|:--|--:|--:|--:|:--|
| v0 | baseline | 2,475,304,254 | 1.00× | 9.17 s | 2654 / 3913 / 2405 |
| v1 | readback scans set bits instead of testing every column; matrix built into one flat buffer | 2,475,304,254 | 1.00× | 6.74 s | 2162 / 4030 / 371 |
| v2 | Four Russians elimination (`k` from the sweep below) | 1,705,630,354 | 1.45× | 5.01 s | 2163 / 2303 / 375 |
| v3 | its pivot scan made read-only (see §4) | 1,639,596,956 | 1.51× | 4.78 s | 2154 / 2090 / 371 |
| v4 | block width fixed at `k = 4` | 1,549,765,406 | 1.60× | 4.32 s | 2149 / 1628 / 371 |
| v5 | monomials interned on sight instead of sorting every nonzero and binary searching it back | 1,549,765,406 | 1.60× | 3.49 s | 1312 / 1634 / 375 |
| v6 | per-row cancellation by parity flip instead of a sort | 1,549,765,406 | 1.60× | 3.26 s | 1090 / 1621 / 375 |
| **v7** | readback reserves the row's exact weight | **1,549,765,406** | **1.60×** | **3.11 s** | 1100 / 1636 / 207 |

### Rejected, and why

- **Bucketed pivot search** (rows carried their leading column, so a pivot was
  found in `O(1)` amortised instead of by rescanning).  The bookkeeping cost more
  than the scan it replaced at these shapes: `n = 23` reduce went `3913 → 5088 ms`.
  Recorded in `steps/stepB.json`.
- **Packed storage with the reference elimination** (`steps/stepB2.json`) and
  **in-place row slices** (`steps/stepB3.json`) — both real but small, and both
  superseded by v2.
- **Block width `k = 6`.**  Faster wall (`2.12×` total against `k = 4`'s `2.06×`
  at the time of the sweep) but **more counted word XORs on the small rungs**
  (`0.99×` and `0.85×`).  That is the relabelling case in §3 of `AGENTS.md`
  exactly — a headline that moves while the work grows — so `k = 4`, the widest
  block that reduces the unit on every rung, is the default.  The full sweep at
  `k = 3…8` is in `steps/blockwidth_k*.json`.

## 4. What the Four Russians step does here, and the one thing it needed

The reference reducer eliminates one column at a time, so it *touches every row
once per pivot column*.  These Macaulay matrices are a few hundred rows of a
handful of 64-bit words each — at `K_1/2^23`, `276 × 439` on average, seven
words per row — so the per-row touch, not the XOR it guards, is the cost.  Four
Russians eliminates `k` columns per pass: tabulate all `2^k` combinations of the
`k` pivot rows once, then fix every other row with one indexed XOR.

Textbook M4RI (`_mzd_gauss_submatrix_full`) reduces every row its pivot search
walks past, because a row's bit at column `j` is only meaningful after that
reduction.  On a tall, rank-deficient matrix that turns a search into a pass over
the whole matrix, and it showed: on the widest shapes the first implementation
did **six times more** word XORs than the reference.  The block's pivots are an
identity on the block's own columns, so the reduced bit is

```text
    bit(i, j) ⊕ parity( bits(i, c..j) & the pivots' bits at j )
```

— one word test per candidate row and no writes, with only the row that turns
out to be the pivot actually reduced.  That is v3, and it is what made the
method pay on every rung rather than on the square ones only.

Reduced row echelon form is unique once the pivot columns are fixed, and the
pivot columns are a property of the row space, so the new kernel must return the
reference kernel's matrix *exactly*, not an equivalent basis.
`m4ri_matches_the_reference_reducer` pins that over 81 random matrices from
`1×1` to `300×300` and `129×1000` at three densities, and
`m4ri_costs_no_more_than_the_reference_on_the_shapes_it_takes` pins the unit.

## 5. A contract that had to change

`ic_f2_dispatch` calibrates a CUDA or RDMA offload against the CPU reducer and,
until now, required *identical matrix, rank and word count*.  The count
requirement is no longer meaningful: the CPU runs a Four Russians elimination and
`gpu/koblitz_f2/reduce.cu` runs the column-at-a-time one, and two eliminations
that agree on the answer can legitimately disagree on how much work reaching it
took.  Calibration now requires an identical reduced matrix and rank, charges
whichever backend ran for its own count, and records a `count_divergences`
statistic — because a *silent* divergence is the one way an offload could make
the accounting unit mean two things in one run.  The CUDA path is feature-gated
(`ic-cuda`) and its tests are `#[ignore]`-gated on a device; **they were not run
in this round**, and the kernel itself is unchanged.  Nothing selected changes
either way: the dispatch panel on the scoreboard measures the CUDA backend at
`21–28×` the CPU's time on exactly these shapes, so calibration was already
choosing the CPU.

## 6. What this does not establish

Per `AGENTS.md` §8, and stated plainly because the temptation runs the other way:

- **This is a stage diagnostic, not an end-to-end ECDLP speed-up.**  The stage is
  one phase of relation collection, which is one phase of the attack.  No
  crossover, exponent or rho comparison follows from it, and none is claimed.
- The WDSat frozen regression suite
  (`research/index_calculus_baseline_20260914/regression/`) does not apply: it
  replays a SAT solver against a conflict counter, and this is a Gröbner engine
  measured in word XORs.  This round freezes its own matched suite under the
  parent accounting contract instead, with the same discipline — frozen inputs,
  three repetitions, a comparison that refuses a changed answer, and the
  rejected variants kept.
- **Practicality note only:** `ic run --degree 23 --curve-a 1 --solver groebner
  --batch 1` recovers the same logarithm (`53`, verified) with relation
  collection at `15.56 s` (median of three) against `43.36 s` measured on this
  host before the round.  The baseline there is a **single sample** — no binary
  of the pre-change tree was retained — so it is a sanity check that the stage
  ratio survives into a whole run, not a timing claim.
- Nothing here changes what `RESEARCH_ECC2K130_DECOMPOSITION.md` concludes.  At
  `n = 131` the oracle's cost is set by how many candidate tuples must be ruled
  out, and `1.6×` off a Macaulay reduction is `2^0.68` off `2^131`.

## 7. Reproducing

```bash
cargo build --release --example groebner_stage_bench
./target/release/examples/groebner_stage_bench --label candidate --out /tmp/gs-cand
python3 research/groebner_stage_20260915/compare.py \
    research/groebner_stage_20260915/baseline /tmp/gs-cand --output /tmp/gs.json

# the elimination-kernel A/B, on instances the tuning never saw
F4_F2_RREF=ref ./target/release/examples/groebner_stage_bench --ladder holdout --out /tmp/gs-ref
./target/release/examples/groebner_stage_bench --ladder holdout --out /tmp/gs-m4ri
python3 research/groebner_stage_20260915/compare.py /tmp/gs-ref /tmp/gs-m4ri

cargo test --release --lib cryptanalysis::koblitz_groebner
```

`F4_F2_M4RI_K` overrides the block width for the sweep; `F4_SHAPE_LOG=path`
records every matrix shape the solver produces, which is how §4's diagnosis was
made.

## References

- G. Bard, *Accelerating cryptanalysis with the Method of Four Russians*, 2006.
- M. Albrecht, G. Bard, W. Hart, *Efficient multiplication of dense matrices over
  GF(2)*, ACM TOMS 37(1), 2010 — `_mzd_gauss_submatrix_full` and the block
  elimination this follows.
- J.-C. Faugère, *A new efficient algorithm for computing Gröbner bases (F4)*,
  J. Pure Appl. Algebra 139 (1999).
