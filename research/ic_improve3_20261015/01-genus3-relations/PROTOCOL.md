# Protocol: round 3, item 01: the relation phase of genus-3 IC after the cubic GHS transfer

Frozen 2026-10-15 (about 16:55 host time), before any test-curve data of this item. Round rules:
[`../RULES.md`](../RULES.md). This item shares the host with one other item: up to 2 single-threaded jobs at a
time, `prlimit --as=6000000000`, `timeout 3600`, `nice`, `RAYON_NUM_THREADS=1`, `OMP_NUM_THREADS=1`.
Prior work: [round 2 item 06](../../ic_improve_20261015/06-genus3-kappa-bw/README.md) and its PROTOCOL (B0-B2), and
[round 2 item 01](../../ic_improve_20261015/01-genus3-linalg/README.md).

## Disclosed before freezing (not test evidence)

- Item 06's records (test curves 0-4 and training curves 5-7), including its training kappa grids and frozen schedule.
- A tabulation of item 06's test records, made to choose the cause to attack. It showed:
  - at kappa = 2 the trial count per q^{4/3} falls from 0.28 (2^12) to 0.18 (2^21);
  - the counted cost per trial rises from 494 to 689 units;
  - 6.64 F_q inversions per trial, each priced at ceil(log2 q) multiplications.
- Smoke tests of the new binary on test curve 0 (b = 12 `sched`, b = 16 `main`):
  - the incumbent reproduced item 06's records exactly;
  - every shadow assertion held;
  - R1's relation cost was 0.50 and 0.48 of the incumbent's.

  These two records are discarded; the test run regenerates them.
- Diagnosis runs on training curves 5-7 (relation phase only, below), started before this file was frozen. They are
  training data by design.

## Question

The incumbent is item 06's cheapest frozen arm at large q, the **combination `sched:lpsgebw2`**:
- the kappa schedule kappa(b) = 1.739 x 0.9225^(b-20);
- cx(b) = 16 for b <= 17 and 8 for b >= 18;
- large primes kept as columns, structured Gaussian elimination, then block Wiedemann with blocking (2, 1).

Its relation phase has a local exponent of about 1.43. The question is whether a cheaper relation phase lowers the
whole-solve cost on held-out curves at the exponent and constant levels, and whether the local exponent moves toward
2 - 2/g = 4/3. The relation phase is split into (trials needed) x (cost per trial).

## Instrument

[`examples/ic_improve3_01_genus3_rel.rs`](../../../examples/ic_improve3_01_genus3_rel.rs) is a copy of item 06's
instrument. `src/` is not edited. The incumbent code path (walk, counters, random streams, linear algebra) is
unchanged. The changes are:

1. **Widening (step 2 of the task).**
   - `N = #E(F_{q^3})` is a `u128`.
   - `F_{q^3}` square roots take a `u128` exponent; q^3 overflows `u64` above 2^21.33.
   - For q >= 2^22 the duplicate-divisor key packs 25 bits per coefficient into two words. Below 2^22 item 06's 21-bit key is kept, so those walks are unchanged.
   - Instance selection adds `ell < 2^62`, the existing Montgomery bound. It holds for every candidate cofactor at q <= 2^21 (N < 2^64, c >= 4), so instances there are unchanged. At q >= 2^22 it forces c >= 17.
   - **Gate 0** (below) verifies the widened build against item 06's records, bit for bit.
2. **Per-primitive counters** of the relation phase (diagnosis, step 1). Saved counter deltas are recorded for:
   - the factor-base set-up;
   - the multiplier set-up;
   - walk additions (steps and restarts);
   - the complete-splitting test of u(x);
   - root finding of split u(x).

   The record also carries the trial statistics: trials, splits, full/1LP/2LP/3LP relations, cycles and restarts. Large-prime bookkeeping is arithmetic mod ell on rows and is not charged in the unit (`rel_counters[2] = 0` in every item 06 record). Its counted share is 0 by construction and its wall time is not part of any claim.
3. **Candidate primitives R1** (step 3), aimed at the cause found: the cost per trial, not the trial count.
   - `cadd_fast`: Cantor addition in the generic case (deg u1 = deg u2 = 3, Res(u1, u2) != 0) with **one inversion**.
     - s u1 = r (mod u2) comes from the adjugate of the multiplication-by-u1 matrix, with r = Res(u1, u2).
     - The composed V = Vt / r is kept with its denominator.
     - Only the top five coefficients of r^2 h - Vt^2 are formed for the exact division by U = u1 u2.
     - One joint inversion of r * lc(U') gives both 1/r and 1/lc(U').
     - The second quotient has leading coefficient lc(h), whose inverse is computed once per instance and charged to R1.
     - Any other case, or a degenerate intermediate degree, falls back to `cadd`. Work already done before the fallback stays charged.
   - `split_fast`: the **Stickelberger pre-filter**. A cubic with three distinct roots in F_q has a square discriminant; type (1,2) has a non-square one.
     - Cost: 9 products, plus the 4 constant multiples charged as products, plus a Legendre exponentiation charged at its square-and-multiply count.
     - It then computes w = x^((q-1)/2) mod u and x^q = x w^2.
   - `roots_fast`: reuses w.
     - When the roots' quadratic characters are mixed, monic(w - 1) or monic(w + 1) is the quadratic factor. It is confirmed by one remainder, solved by the incumbent's quadratic formula (same charge), and the third root comes from the trace.
     - Otherwise it falls back to Cantor-Zassenhaus on a private stream that never touches the walk's stream.
   - **Shadow mode** (default):
     - The incumbent primitives drive the walk.
     - The candidate primitives run beside them on every trial, under saved counters.
     - Their outputs are asserted equal: the sum (u, v), the splitting decision and the root set. Any disagreement aborts the job.
   - The candidate arm `<arm>+R1` therefore has exactly the incumbent's walk, relations, matrix and logarithm. Its relation counters are the incumbent's, with the incumbent primitive counts replaced by R1's plus R1's set-up.
   - `--relprim fast` drives the walk with R1 alone (arm `+R1s`). This is a check of the accounting, not a test arm.

The counted unit is unchanged (item 01): `cost = M_q + ceil(log2 q) I_q + c_ell M_ell`, with `c_ell = 1` primary and
`c_ell = 9` as sensitivity.

## Step 1: diagnosis (training curves 5, 6, 7 only)

The relation phase only is run (`--las none --max-checks 1`) with shadow primitives:
- kappa = 2 at b = 12..21, 3 curves each;
- the kappa grid {1.0, 1.19, 1.41, 1.68, 2.38} at b = 16 and 18;
- kappa = 1.19 at b = 20.

Item 06's test records are not re-used for this; its training records are used only for the R2 rule below.
Reported, with exponents fitted by OLS of ln y on ln q over 12..21 (every solve a point, 2000 bootstrap resamples of
curves within size):
- trials T, and T kappa / q^{4/3};
- per-trial cost, split into walk additions, splitting test and root finding;
- the per-trial cost fitted as A + B log2 q;
- the splitting probability, against 1/6;
- the 2LP-graph mean degree 2E/V (with V = (q-1)/2), and cycles over the factor-base size;
- for the relation-phase cost, q^a (log q)^beta against q^c.

The diagnosis is descriptive; it decides no verdict.

## Candidates

- **R1**: the incumbent configuration with R1 primitives. It is the same stream (`sched`), so it is exactly paired.
- **R2**: R1 with kappa re-tuned jointly with the cheaper relation phase. It runs as a new stream `r2` and starts from the same random state as `sched`, so the walk's divisor sequence is identical and only the factor base differs.
  - The R2 kappa schedule comes from a frozen training rule (below).
  - The stream also reports `r2:lpsgebw2` (re-tuned kappa with the incumbent primitives) as a decomposition arm, not a candidate.

**Frozen R2 rule** (applied mechanically to training data and recorded as amendment C1 before any `r2` test job):
1. rho(b, c, kappa) = R1 relation cost / incumbent relation cost, from the relation-only training runs above.
2. For each training size b in {16, 18, 20, 21}, each curve c in {5, 6, 7} and each kappa in item 06's training grid at b:
   - the item-06 training record gives the incumbent relation cost and the LA cost of `lpsgebw2@8`;
   - the predicted R1 whole cost is rho(b, c, kappa') x (incumbent relation cost) + (LA cost);
   - kappa' is the measured kappa nearest in ln kappa;
   - the incumbent relation cost from the item-06 record is used only if it equals this item's relation-only incumbent rerun at the measured kappa (gate), else the size is dropped and reported.
3. At each b: geometric mean over the 3 curves; least-squares quadratic in ln kappa; kappa\*_R1(b) is its minimiser, clamped to the grid range (grid argmin if not convex).
4. The schedule is the least-squares line ln kappa_R2(b) = alpha + beta (b - 20) through the four points. It is used frozen at every test size, with cx(b) as the incumbent.

The linear-algebra cost does not depend on the relation primitives, which is why this prediction is exact up to the
kappa-dependence of rho. That dependence is measured at b = 16, 18 and 20 and reported.

**Considered and not built** (with reasons in the README):
- a large-prime bound;
- 3LP;
- an r-adding walk with more multipliers;
- batching inversions across parallel walks. R1's single inversion per addition already removes most of what batching would remove.

## Test (held-out curves)

Instances: `q` = next prime above 2^b + 1, seed `20261014*1000 + 10 b + curve` (items 04/01/06).

| collection (stream) | kappa | LA arm | arms produced | role |
|---|---|---|---|---|
| `sched` | item 06 kappa(b) | `lpsgebw2@cx(b)` | `sched:lpsgebw2` (incumbent), `sched:lpsgebw2+R1` | incumbent; candidate R1 |
| `r2` | kappa_R2(b) | `lpsgebw2@cx(b)` | `r2:lpsgebw2+R1` (candidate R2), `r2:lpsgebw2` (decomposition) | candidate R2 |
| `aa` | item 06 kappa(b) | `lpsgebw2@cx(b)` | `aa:lpsgebw2`, `aa:lpsgebw2+R1` | A/A: independent walk stream (item 04's `aa` rule) |

Sizes and run order (budget priority):
1. `sched` on curves 0-2, b = 12..21;
2. `sched` at b = 22 on curves 0-2;
3. `r2` on curves 0-2, b = 12..21, then b = 22;
4. `aa` on curves 0-2, b = 12..21;
5. `sched` and `r2` on curves 3 and 4, b = 18..21.

A size enters an arm's fit only if that arm and the incumbent have 3 or more correct curves there. The top size `top`
is 22 if `sched` (and, for R2, `r2`) completed 3 curves there within 60 minutes each, else 21. Any job that exceeds 60
minutes or 6 GB is dropped and reported with where the machine was the limit. b = 23 is attempted only if every
b = 22 job finished in under 25 minutes.

## Metric and analysis

Phases:
- relation phase: the transfer plus the walk, `rel_cost`;
- linear algebra: `la_cost`.

Analysis follows [`examples/ic_improve3_01_analyze.rs`](../../../examples/ic_improve3_01_analyze.rs), with 2000 bootstrap resamples of curves within each size, paired across arms:
- **Exponents**: OLS slope of ln cost on ln q, every solve a point, over 12..top and 16..top, for the whole solve and both phases. Delta = slope(candidate) - slope(incumbent).
- **Local exponents** between consecutive sizes. Their **trend** is the slope of the local exponent on b, per bit, with the paired trend difference.
- **Paired ratio** candidate / incumbent per size: the geometric mean of per-curve ratios, with a CI.
- Rho is not run. Item 06's references are quoted only.

## Decision rule (per candidate, per level)

Noise floors come from the A/A arm over the same range:
- eps_e = the larger absolute end of the A/A exponent-Delta CI (whole, 12..21);
- eps_c = the largest |A/A ratio - 1| over the fit sizes.

Levels:
- **Exponent**: ADVANCES if the Delta CI (whole, 12..top) lies below 0 and its point is below -eps_e. REGRESSES symmetrically. Else MATCHES.
- **Constant**: the paired ratio at the top size. ADVANCES if its CI lies below 1 and its point is below 1 - eps_c. REGRESSES symmetrically. Else MATCHES.
- **Primitive** (R1): the same rule on the relation-phase ratio at the top size.
- **TRADES** when levels or phases point in opposite directions.
- **Moves toward 4/3**: yes if the paired trend difference of the whole-solve local exponent has its CI below 0 **and** the candidate's local exponent on the top interval is closer to 4/3 than the incumbent's.
- Separately, and only descriptively, the same is reported for the relation phase.

## Controls and gates

- **Gate 0 (widened build).** Every `sched:lpsgebw2` test record at b <= 21 must equal item 06's `sched:lpsgebw2` record on the same (b, curve): `k`, cost, relation cost and LA counters.
  - Every relation-only incumbent training record at kappa in item 06's grid must equal item 06's training relation cost on the same (b, curve, kappa).
  - A difference fails gate 0 and stops the item.
- **Shadow equality**: every trial of every job asserts that R1's outputs equal the incumbent's (job aborts otherwise).
- **Accounting check**: on test curve 0 at b = 14 and 16, `--relprim fast` (R1 driving the walk alone) must give the same relation cost as the shadow-accounted `+R1` arm, and the verified `k`.
- **A/A**: `aa:lpsgebw2` against `sched:lpsgebw2`; its exponent Delta CI should cover 0.
- **Correctness**: every solve's `k` is verified (`k phi(G) = phi(T)`) and compared with the planted `k`. A wrong or missing `k` fails the arm.

## Amendment C0 (2026-10-15 17:15, before any test result was read)

The job-chaining scripts waited on `pgrep -f run_diag.sh`, which also matched the waiting shells themselves, so the
test queue did not start when the training queue finished (about 5 minutes lost). The patterns are now anchored
(`pgrep -f "^sh run_diag.sh"`). No record is affected.

## Amendment C1 (2026-10-15 17:16, training complete, before any `r2` test job and before any test result was read)

Training, all on curves 5-7 ([results/train/](results/train/)):
- 63/63 relation-only jobs ran, 0 failed.
- **Training gate**: 45/45 incumbent relation-only reruns equal item 06's training relation cost and step count on the same (b, curve, kappa) ([gatetrain.md](results/train/gatetrain.md)).
- rho does not depend on kappa: 0.4797-0.4802 over kappa 1.0-2.38 at b = 16; 0.4767-0.4769 at b = 18; 0.4739 at both kappa 1.19 and 2 at b = 20.

The frozen rule ([r2rule.md](results/train/r2rule.md)) gives:
- kappa\*_R1(16) = 1.969;
- kappa\*_R1(18) = 1.495;
- kappa\*_R1(20) = 1.289;
- kappa\*_R1(21) = 1.227, inside the grid (1.19-2.0), so not clamped.

The schedule is **ln kappa_R2(b) = 0.2669 - 0.0941 (b - 20)**, i.e. kappa_R2(b) = 1.306 x 0.9102^(b-20):

| b | 12 | 14 | 16 | 18 | 20 | 21 | 22 |
|---|---|---|---|---|---|---|---|
| kappa_R2(b) | 2.771 | 2.296 | 1.902 | 1.576 | 1.306 | 1.189 | 1.082 |

Stored in `r2.env`. cx(b) is the incumbent's.
