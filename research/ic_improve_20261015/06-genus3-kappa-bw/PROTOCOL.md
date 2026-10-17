# Protocol: round 2, item 06: a kappa schedule and block Wiedemann for genus-3 IC after the cubic GHS transfer

Frozen 2026-10-15 before any training or test data of this item. Round rules: [`../RULES.md`](../RULES.md)
(this item runs alone: up to 3 single-threaded jobs in parallel, 8 GB and 60 minutes per job).
Prior work: [item 01](../01-genus3-linalg/README.md) and its [PROTOCOL](../01-genus3-linalg/PROTOCOL.md) (A1-A5).

Already seen before this protocol (disclosed, not re-used as evidence):
- item 01's results, including its training records at b = 20 (curves 5-7; `lpsgewd` cx 8 at kappa 1, 1.41, 2, 2.83:
  1.41 and 2 tie) and its exploratory kappa = 1.41 test arm (cheapest at 2^21 on test curves 0-2);
- a smoke test of the new binary on test curve 0 at b = 12, 14, 16 (kappa 2; arms `lpsgewd`, `lpsgebw2/3/4`,
  `lpsgewd@16`): `lpsgewd` reproduced item 01's records exactly and every arm returned the verified `k`. Those
  three records are discarded; the test run regenerates them.

## Question

Against the incumbent IC component, **`lpsgewd` at kappa = 2, cx = 8** (item 01: large primes kept as columns,
structured Gaussian elimination, scalar Wiedemann; q^1.380 over 2^10..2^20), do

1. a **kappa schedule** kappa(q) (with the SGE constant cx allowed to co-vary), fitted on training curves and frozen,
2. **block Wiedemann** mod ell on the post-SGE matrix in place of scalar Wiedemann,
3. their **combination**,

lower the whole-solve cost on held-out curves, at the exponent and constant levels, and does the trend of the local
exponent flatten toward 2 - 2/g = 4/3?

## Instrument

`examples/ic_improve2_06_genus3_kappa_bw.rs`, a copy of item 01's `ic_improve2_01_genus3_linalg.rs` (with its A5
overflow fix). Relation collection, counters, random streams and every existing arm are unchanged; `src/` is not
edited. Added:

- **`<arm>@<cx>`**: a per-arm SGE constant, so arms that differ only in cx share one collection.
- **`lpsgebw<m>`**: `lpsgewd` with the scalar Wiedemann replaced by **block Wiedemann with blocking (m, 1)**:
  m coordinate projections (the `k` coordinate and m-1 pseudo-random others; reading a coordinate costs nothing)
  and one right vector b. The vector sequence X^T A^i b determines the minimal recurrence of b after about
  N (1 + 1/m) terms instead of 2N, so the matrix-vector products fall from ~2N to ~N (1 + 1/m). The generator is
  an order basis of [A(x) | -I_m] (iterative mbasis, shifts (0, 1, ..., 1)) whose work grows like m (m+1) N^2 / 2
  products. `k = -(1/g_0) sum_{i>=1} g_i (A^{i-1} b)[k]` as in the scalar arm. Early termination: stop once the
  number of terms reaches L + ceil(L/m) + 20, with L the least shifted degree of a generator candidate.
  Counted in item 01's unit: one `M_ell` per stored nonzero per matrix-vector product, every product of the
  generator computation, every inversion mod ell as `I_ell` (reported, not priced, as in items 01/04).
- **Symmetric blocking (n = m right vectors)** is not built: single-threaded it needs n (N/m + N/n) = 2N
  single-vector products, the same as scalar Wiedemann, plus a larger generator computation, so its counted cost
  cannot fall; it only adds parallelism. This is stated, not measured.

**Large-prime bound.** The instrument has no large-prime bound (every prime outside the factor base is a large prime,
as in items 04 and 01). Adding one changes the relation collection; it is not done here. Only cx co-varies with kappa.

## Training (curves 5, 6, 7 only; never in any test fit)

Each training job is one collection (`main` stream) at a given (b, curve, kappa) with several LA arms sharing it.

| b | kappa grid | LA arms |
|---|---|---|
| 16, 18 | 1.0, 1.19, 1.41, 1.68, 2.0, 2.38 | `lpsgewd@8`, `lpsgewd@16`, `lpsgewd@32`, `lpsgebw2@8`, `lpsgebw3@8`, `lpsgebw4@8` |
| 20, 21 | 1.19, 1.41, 1.68, 2.0 | `lpsgewd@8`, `lpsgewd@16`, `lpsgebw2@8`, `lpsgebw3@8` |

If a b = 20 job takes more than 20 minutes of wall time, the b = 21 grid drops kappa 1.19 (recorded as an amendment
before any b = 21 result is read).

Frozen rules (applied mechanically, recorded as amendment B1 before any test data):

1. **Block size m\***: the m in {2, 3} with the lowest geometric-mean LA cost of `lpsgebw<m>@8` at kappa = 2 on the
   largest training size with all three curves complete (b = 21, else 20). m = 4 is kept only if it is the lowest
   at b = 18 (then compared at b = 18 instead).
2. **Scalar or block** for the combination: block if `lpsgebw<m*>@8` has a lower geometric-mean whole cost than
   `lpsgewd@8` at kappa = 2 on the largest complete training size; else scalar.
3. **cx(b)**: for each training size, the cx whose `lpsgewd@cx` has the lowest geometric-mean whole cost at that
   size's best kappa. Applied as a step function: test b <= 17 uses the b = 16 choice, 18-19 the b = 18 choice,
   20 the b = 20 choice, >= 21 the b = 21 choice.
4. **kappa(b)**: at each training size, with `lpsgewd@cx(b)`, fit a quadratic in ln kappa to the ln geometric-mean
   whole cost over the kappa grid; kappa\*(b) is its minimiser, clamped to the grid range (grid argmin if the
   quadratic is not convex). The schedule is the least-squares line `ln kappa(b) = alpha + beta (b - 20)` through the
   four points (b, ln kappa\*(b)). It is used, frozen, at every test size, including b = 22 (an extrapolation of 1 bit).

## Test (held-out curves)

Instances: `q` = next prime above 2^b + 1, seed `20261014*1000 + 10 b + curve` (items 04/01). **Test curves 0, 1, 2
at b = 12..22; curves 3 and 4 (never used by item 01 or here before) at b = 18..22**, as the budget allows.
Every job: single-threaded, `nice`, `prlimit --as=8000000000`, `timeout 3600`.

Collections per (b, curve), all with the same instance, target `k`, and walk seed (`main` and `sched` streams start
from the same random state, as item 01's `tuned` stream):

| collection | kappa | arms | role |
|---|---|---|---|
| `main` | 2 | `lpsgewd` (cx 8), `lpsgebw<m*>@8` | **incumbent**; candidate 2 (block Wiedemann) |
| `sched` | kappa(b) | `lpsgewd@cx(b)`, `lpsgebw<m*>@cx(b)` | candidate 1 (schedule); candidate 3 (combination, if rule 2 picks block; else reported) |
| `aa` | 2 | `lpsgewd` (cx 8) on an independent walk stream | A/A null |

Run order (budget priority): `main` and `sched` on curves 0-2 by increasing b; then `aa` on curves 0-2; then curves 3-4.
A size enters the primary fit only if every arm has 3 or more curves there. If a b = 22 job cannot finish in 60 minutes,
b = 22 is dropped for that collection and the fit range ends at 21; this is reported.

## Metric

Item 01's unit: `cost = M_q + ceil(log2 q) I_q + c_ell M_ell`, `c_ell = 1` primary, `c_ell = 9` sensitivity. Phases:
relation (transfer + walk, `rel_cost`) and linear algebra (`la_cost`). For the block arm, LA products are split into
matrix-vector products (`matvecs x nnz`) and generator products (the rest, SGE excluded).

## Analysis (bootstrap: 2000 resamples of curves within each size, paired across arms)

- Exponents: least-squares slope of ln cost on ln q, every solve a point, over **b = 12..top** (top = largest size in
  the fit) and **16..top**; whole, relation phase, LA phase. Paired Delta = slope(cand) - slope(`main:lpsgewd`).
- Local exponents between consecutive sizes; their **trend** = slope of local exponent on b over 12..top, per bit;
  paired trend difference.
- Paired ratio cand / incumbent per size (geometric mean of per-curve ratios) with CI.
- Rho and item 01's old IC-2LP (`main:inc`) are reference columns only (item 01's records, same instances, b <= 21).

## Decision rule (per candidate, per level)

Noise floors from the A/A arm on the same range: eps_e = the larger absolute end of the A/A exponent-Delta CI;
eps_c = the largest |A/A ratio - 1| over the fit sizes.

- **Exponent**: ADVANCES if the Delta CI (whole, 12..top) lies below 0 and its point is below -eps_e; REGRESSES if
  the CI lies above 0 and the point is above +eps_e; else MATCHES.
- **Constant**: the paired ratio at the top size; ADVANCES if its CI lies below 1 and the point is below 1 - eps_c;
  REGRESSES symmetrically; else MATCHES.
- **Primitive** (block Wiedemann only): the same rule on the LA cost ratio at the top size.
- TRADES when levels (or phases) point in opposite directions.
- **Flattens toward 4/3**: yes if the paired trend difference has its CI below 0 **and** the candidate's local exponent
  at the top interval is closer to 4/3 than the incumbent's.

## Controls and gates

- **Gate 0**: `main:lpsgewd` must reproduce item 01's `main:lpsgewd` records on curves 0-2, b = 12..21 (k, cost and
  LA counters equal), and the training `lpsgewd@8` records at b = 20, kappa 1.41 and 2 must equal item 01's
  `tune2` records.
- **A/A**: `aa:lpsgewd / main:lpsgewd`; its exponent Delta CI should cover 0 (item 01 found that 3 curves undercover;
  eps_e and eps_c above absorb that).
- **Correctness**: every solve's `k` is verified (`k phi(G) = phi(T)`) and compared with the planted `k`; a wrong or
  missing `k` fails the arm.

## Amendment B0 (2026-10-15, after the b = 16 training jobs and before any test data)

For the machine budget the `main` and `aa` test collections, which do not depend on the training rules, start while
training runs. Because m\* is not yet known, every `main` test collection carries **both** `lpsgebw2@8` and
`lpsgebw3@8`; the frozen rule 1 then names which of them is the block-Wiedemann candidate, and the other is reported
only. `aa` runs on curves 0-2 only. This changes no arm's cost (arms sharing a collection are booked separately).

## Amendment B2 (2026-10-15 05:25, before the training rules were applied to b = 21 and before any test result was read)

**b = 22 cannot be run by this instrument.** Instance building parses the order of E over F_{q^3} (about q^3 = 2^66)
into a `u64` and panics (`ParseIntError PosOverflow`, [results/b22_overflow/panic.txt](results/b22_overflow/panic.txt));
the Montgomery arithmetic mod ell also requires ell < 2^62. Widening it changes code shared with every arm, so the
top size is **b = 21** (an instrument limit, not a time limit). The freed budget goes to `aa` at b = 21 and to test curves
3 and 4 at b = 18..21 (`main` and `sched`). The schedule still covers b = 12..21; it is then an interpolation inside the
training range 16..21 above b = 16 and an extrapolation below it.

## Amendment B1 (2026-10-15 06:20, training complete, before any `sched` test data and before any test result was read)

Training ([results/train/train.md](results/train/train.md), curves 5-7, all 60 jobs correct (the original text said 66, a miscount, corrected 08:40); the b = 20 `lpsgewd@8`
records at kappa 1.41 and 2 equal item 01's `tune2` records, 6/6). The b = 21 jobs at kappa 2 on curves 6 and 7 were
killed by a session time limit of the orchestration shell and rerun unchanged. The frozen rules give:

1. **m\* = 2** (LA geometric mean at b = 21, kappa 2: `lpsgebw2@8` 2.843e10 against `lpsgebw3@8` 2.989e10; at b = 18
   m = 4 is the most expensive).
2. **Block** for the combination (whole cost at b = 21, kappa 2: `lpsgebw2@8` 6.162e10 against `lpsgewd@8` 6.515e10).
3. **cx(b) = 16 for b <= 17, 8 for b >= 18** (b = 16: cx 16 best; b = 18, 20, 21: cx 8 best).
4. kappa\*(16) = 2.380 (clamped at the grid edge), kappa\*(18) = 2.081, kappa\*(20) = 1.719, kappa\*(21) = 1.607;
   **schedule ln kappa(b) = 0.5533 - 0.0806 (b - 20)**, i.e. kappa(b) = 1.739 x 0.9225^(b-20):
   b = 12: 3.315, 14: 2.821, 16: 2.401, 18: 2.043, 20: 1.739, 21: 1.604.

Candidates on test: (1) `sched:lpsgewd@cx(b)`, (2) `main:lpsgebw2@8`, (3) `sched:lpsgebw2@cx(b)`; all paired against
`main:lpsgewd`. `main:lpsgebw3@8` is reported only.
