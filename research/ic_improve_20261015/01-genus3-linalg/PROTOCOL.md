# Protocol: round 2, item 01: the linear-algebra bottleneck of genus-3 IC-2LP after the cubic GHS transfer

Frozen 2026-10-15 before any training or test data. Round rules: [`../RULES.md`](../RULES.md).
Prior work: [`../../ic_improve_20261014/04-genus3-2lp/`](../../ic_improve_20261014/04-genus3-2lp/README.md)
(IC-2LP scales as q^1.385 over 2^10..2^18; its linear algebra as q^1.519 and is 65% of the cost at 2^20).

## Question

On item 04's instances, does replacing the incumbent's linear algebra (and, separately, re-tuning the
factor-base constant kappa) lower the cost of the whole IC-2LP solve, at the exponent level over
2^10..2^20, and does the whole exponent move toward 2 - 2/g = 4/3?

## Instrument

`examples/ic_improve2_01_genus3_linalg.rs`, a copy of `examples/ic_improve_04_genus3_2lp.rs` with the
relation collection, counters and random streams unchanged. `src/` is not edited.

**Shared collection.** One walk collects relations exactly as item 04 does. At every check (the
incumbent's schedule: when the singleton-filtered system has `nrows >= ncols + 5`, then every 10%
more rows), every linear-algebra arm that has not yet produced a verified `k` solves its own system.
Its counters are booked to that arm and removed from the shared counters. An arm's whole cost is
(transfer + relation-phase counters at the check where it succeeded) + (all of its own
linear-algebra counters, over every check it ran). The walk continues until every arm has succeeded
(or made 8 checks: then that arm fails). So every arm sees the same instance, target, walk, relations
and check points; arms differ only in the linear algebra. The incumbent arm reproduces item 04
operation for operation.

## Arms

Incumbent **`inc`**: item 04's IC-2LP. kappa = 2 (factor base `ceil(2 q^{2/3})`), large primes
eliminated by the DFS spanning-tree substitution of item 04, singleton filter, Lanczos on the normal
equations mod ell over all filtered rows.

Linear-algebra candidates (same collection, kappa = 2):

| arm | rows | SGE | solver |
|---|---|---|---|
| `wied` | incumbent rows | no | Wiedemann |
| `sgelz` | incumbent rows | yes | Lanczos (normal equations) |
| `sgewd` | incumbent rows | yes | Wiedemann |
| `bfslz` | BFS cycle basis | no | Lanczos |
| `bfswd` | BFS cycle basis | no | Wiedemann |
| `bfssgelz` | BFS cycle basis | yes | Lanczos |
| `bfssgewd` | BFS cycle basis | yes | Wiedemann |
| `lpsgelz` | large primes kept as columns | yes | Lanczos |
| `lpsgewd` | large primes kept as columns | yes | Wiedemann |

- **BFS cycle basis**: the same consumed relations, re-based on a breadth-first spanning tree of the
  large-prime graph rooted at the "1" vertex; one row per non-tree edge. Same row count and same row
  space as the incumbent's rows; built with additions and sign changes only, like the incumbent's
  substitution (neither is charged products).
- **Large primes as columns**: every consumed relation is a row with its large primes as extra
  columns; singleton filtering and SGE decide which large primes to merge away.
- **SGE** (structured Gaussian elimination, ideas ported from `src/cryptanalysis/koblitz_sparse_la.rs`:
  singleton cascade, Markowitz-order merging of light columns, excess removal of heaviest rows with
  singleton cascade): a column is eliminated, smallest estimated fill first, while its fill is below
  the current average column weight + `cx` (the condition under which `2 nnz ncols + cx ncols^2`
  falls). The `k` column is never eliminated. Each product mod ell is counted (products by +-1 are
  additions and are not counted, as in the incumbent's substitution); each non-unit pivot inversion
  is counted as an inversion mod ell.
- **Wiedemann** (scalar): the system is made square by removing heaviest rows (singleton cascade);
  the sequence `s_j = (A^j b)[k]` feeds Berlekamp-Massey (early termination after 20 terms with
  `L` unchanged past `2L`), and `k = -(1/g_0) sum_{i>=1} g_i s_{i-1}` with `g` the minimal
  polynomial. One matrix-vector product per term (`nnz` products), BM products and inversions counted.

Other arms:
- **A/A** (`aa:inc`): the incumbent on an independent walk stream (item 04's A/A construction), same
  instance and target; also `aa:<best>` for the chosen candidate.
- **kappa re-tuned** (`tuned:inc`, `tuned:<best>`): the frozen kappa*, run as its own collection on
  the same instance, target and walk seed. `tuned:<best>` is the combination arm.
- Rho is not re-run: item 04's rho records on the same instances (same seeds) are the reference column.

## Tuning (training curves 5, 6, 7 only; never used in fits)

1. LA arm and `cx`: all nine candidates with `cx in {0, 2, 8}` at b = 14, 16, 18, kappa = 2.
   Frozen: the (arm, cx) with the lowest geometric-mean **whole** cost at b = 18.
2. kappa: with the frozen LA arm and also with `inc`, kappa in {0.5, 0.71, 1, 1.41, 2} at b = 16, 18, 20.
   Frozen kappa*: lowest geometric-mean whole cost of the frozen LA arm at b = 20.
3. Choices are recorded in an amendment before any test-curve data.

## Sizes, instances, budget

`q` = next prime above `2^b + 1`, seed `20261014*1000 + 10 b + curve` (item 04's instances). Test
curves 0, 1, 2, b = 10..20, and b = 21 (and 22) if a job fits in 60 min and 4 GB. Single-threaded,
`nice`, `prlimit --as=4000000000`, `timeout 3600`.

## Metric

Item 04's unit: `cost = M_q + ceil(log2 q) I_q + c_ell (M_ell)`, `c_ell = 1` primary, `c_ell = 9`
sensitivity, over all phases. `M_ell` counts products mod ell. Inversions mod ell (counter 3)
are not priced, as in item 04; they are reported (each arm makes at most a few thousand). Also reported: linear-algebra cost alone (`la_cost`, `la_cost_cl9`), matrix
dimensions and nonzeros given to the solver. Wall time is recorded but not used.

## Analysis (bootstrap: 2000 resamples of curves within each size, paired across arms)

- Whole and LA exponents: least-squares slope of `ln cost` on `ln q` over b = 10..20 (and 14..20), CI.
- Exponent difference `Delta = slope(cand) - slope(inc)`, paired bootstrap.
- Local exponents between consecutive sizes and their trend (slope per bit) with CIs.
- Paired ratio cand/inc per size (geometric mean of per-curve ratios), CI.

## Decision rule (per level)

- **Exponent**: ADVANCES if the CI of `Delta` (whole cost, b = 10..20) lies below 0; REGRESSES if above;
  MATCHES if it covers 0.
- **Constant**: paired ratio at the largest common size; ADVANCES if its CI lies below 1, REGRESSES if
  above, MATCHES if it covers 1.
- **Primitive**: the same two rules on the linear-algebra cost alone.
- TRADES when levels point in opposite directions.
- **"Moves toward 4/3"**: yes if `Delta` (whole, 10..20) has its CI below 0 **and** the candidate's
  14..20 slope is closer to 4/3 than the incumbent's with the CI of the paired difference below 0.

## Controls and gates

- **Gate 0**: `main:inc` must reproduce item 04's `ic2lp` records on curves 0-2 (cost and k equal at every size).
- **A/A**: `aa:inc / main:inc`: exponent difference CI must cover 0. Per-size ratios read with item 04's
  noise floor (about +-2% with 3 curves).
- **Correctness**: every arm's `k` is checked as `k phi(G) = phi(T)` and against the planted `k`; all
  arms agree with the incumbent, or the arm fails.

## Amendment A1 (2026-10-15, after tuning stage 1, before any test-curve data)

- Stage 1 ([results/tune1/](results/tune1/), training curves 5-7, kappa = 2): lowest geometric-mean
  whole cost at b = 18 is **`lpsgewd` with `cx = 8`** (3.371e9, against 4.995e9 for `inc`; next
  `lpsgelz` cx 8 3.571e9, `bfssgewd` cx 8 3.574e9). Frozen: `cx = 8` for every SGE arm, and
  `lpsgewd` is "the best LA arm" for the A/A and kappa-tuned arms. `cx = 8` is the edge of the grid.
- Stage 2 grid, changed before any stage-2 data for the machine budget (load 9 on 4 cores, three
  items): kappa in {0.71, 1, 1.41, 2, 2.83} at b = 18 and kappa in {1, 1.41, 2, 2.83} at b = 20,
  training curves 5-7, LA arms `inc` and `lpsgewd`. 2.83 is added because the cheaper linear
  algebra may move the optimum above 2; 0.5 is dropped (it can only win if LA dominates, which stage 1
  says it no longer does). b = 16 is dropped. kappa* = lowest geometric-mean whole cost of `lpsgewd` at b = 20.
- Test run: `main` stream with all ten LA arms at b = 10..20; at b = 21 only `inc, wied, bfssgewd,
  lpsgewd` (time cap); `aa` stream with `inc, lpsgewd` at b = 10..21.

## Amendment A2 (2026-10-15, during the test run, before any test result was read)

The host ran at load 8-12 on 4 cores (three items). The first two b = 20 jobs with all ten arms had
used 12 CPU minutes in 41 wall minutes and would have reached the 60-minute cap; one was stopped
(curve 0) and one was left running (curve 1; its result is kept if it finishes, and is identical by
construction to the split jobs). At b >= 19 the arms are split over two collections per curve
(`inc, wied, sgewd` and `bfssgewd, lpsgelz, lpsgewd`). This does not change any arm's cost: an arm's
cost is booked at its own successful check and does not depend on which other arms share the
collection (`main:inc` is in the first group, so gate 0 still applies). The four arms that tuning
stage 1 found dominated (`sgelz`, `bfslz`, `bfswd`, `bfssgelz`) are not run at b >= 19. b = 21 is
dropped: at this load a b = 21 job would not fit in 60 minutes.

## Amendment A3 (2026-10-15, after tuning stage 2 at b = 20, before any analysis of test data)

- Stage 2 at b = 20 (training curves 5-7; [results/tune2/](results/tune2/)), geometric-mean whole cost
  of `lpsgewd`: kappa 1: 2.831e10, 1.41: 2.469e10, **2: 2.446e10**, 2.83: 2.789e10 (curve 5; curves 6, 7
  were still running and cannot change the order). **kappa* = 2**, item 04's value, so the
  kappa-tuned arms `tuned:inc` and `tuned:lpsgewd` are `main:inc` and `main:lpsgewd` and are not run
  separately. The b = 18 part of the stage-2 grid was stopped to free the host.
- For the incumbent's own linear algebra the stage-2 data pick a different kappa: `inc` costs 3.272e10
  (kappa 1), **3.204e10 (1.41)**, 3.695e10 (2), 4.946e10 (2.83, curve 5). So "re-tune kappa alone" is a
  separate candidate. Added arm, exploratory and secondary: stream `tuned` at **kappa = 1.41** with LA
  arms `inc` (= kappa re-tune only) and `lpsgewd` (= kappa re-tune + new LA), same instances, targets
  and walk seed as `main`, b = 10..20 as the budget allows, paired against `main:inc`.
