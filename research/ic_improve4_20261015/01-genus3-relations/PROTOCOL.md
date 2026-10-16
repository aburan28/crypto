# Protocol: round 4, item 01: the genus-3 relation phase, second pass

Written 2026-10-15 (about 23:55 host time) before any test-curve data of this item; the training-derived
values are filled in by amendment D1 before any test job. Round rules: [`../RULES.md`](../RULES.md): at most 3
single-threaded jobs at a time, `prlimit --as=4500000000`, `timeout 3600`, `nice`, `RAYON_NUM_THREADS=1`,
`OMP_NUM_THREADS=1`. Predecessor: [round 3 item 01](../../ic_improve3_20261015/01-genus3-relations/README.md)
and its [PROTOCOL](../../ic_improve3_20261015/01-genus3-relations/PROTOCOL.md) (amendments C0-C2).

## Disclosed before freezing (not test evidence)

- Round 3's records (test curves 0-4, training curves 5-7) and its analysis.
- Smoke runs of the new binary on test curve 0 at b = 12..16 (`r2` stream, arms `lpsgebw2@cx` and
  `lpsgebw2o20@cx`): the incumbent arm reproduced round 3's `r2:lpsgebw2+R1` cost exactly at 12, 14 and 16; the C1
  relation phase was 0.71-0.86 of R1's. These records are discarded (`results/smoke/`, not in any fit); the test
  run regenerates them.
- Probes on training curve 5 at b = 18 and 20 (oversampling, SGE constant, clique removal, W = 1 against W = 32),
  used to choose the training grid below. Training data by design.

## Question

The incumbent is round 3's **R2**: round 3's R1 relation primitives (one-inversion Cantor addition, the
Stickelberger discriminant pre-filter, root finding reusing `x^((q-1)/2) mod u`) with the re-tuned schedule
ln kappa_R2(b) = 0.2669 - 0.0941 (b - 20), cx(b) = 16 for b <= 17 and 8 above, on `sched:lpsgebw2` (large primes as
columns, SGE, block Wiedemann (2, 1)). Unit: `M_q + ceil(log2 q) I_q + M_ell` (`c_ell = 9` as sensitivity).

Round 3 found that at the optimal kappa the whole-solve exponent is about 4/3 + (2/3)(per-trial log growth) +
(1/3)(LA excess) + trial drift, with the LA excess the largest term. The question: do relation-phase changes that
(a) cut the per-trial cost and its log growth and (b) change what the relation phase hands to linear algebra lower
the whole-solve cost on held-out curves, at the primitive, constant and exponent levels, and move the local exponent
toward 4/3?

## Instrument

[`examples/ic_improve4_01_genus3_rel.rs`](../../../examples/ic_improve4_01_genus3_rel.rs), a copy of round 3's
instrument (`src/` untouched). Changes:

1. **R1 drives the walk.** Round 3 ran the incumbent primitives and accounted R1 beside them. Here R1's primitives
   (`cadd_fast`, `split_fast`, `roots_fast`) are charged directly. Round 3's `roots_split` still runs on every new
   split, with its counters restored, so the root order and the use of the walk's random stream are exactly round
   3's; R1's root set is asserted equal to it. The incumbent arm `r2:lpsgebw2` must therefore equal round 3's
   `r2:lpsgebw2+R1` record bit for bit (**gate 0**).
2. **C1 (field-level splitting test, identical relations).** For a cubic u = x^3 + a x^2 + b x + c over F_q (q prime,
   q > 3; every instance field here is a prime field):
   - as R1: the discriminant (13 products) and its Legendre symbol (square-and-multiply count); a non-square
     discriminant means type (1,2);
   - square discriminant: u is split or irreducible. With the depressed cubic t^3 + P t + Q and Cardano's
     resolvent roots z1, z2 (z1 z2 = -P^3/27), u splits iff y = z2/z1 is a cube in the cyclic group of order
     n = q - (-3|q) (F_q^* when q = 1 mod 3; the norm-1 torus of F_{q^2} when q = 2 mod 3). With P' = 3P = 3b - a^2,
     y + 1/y = 2 + 27 disc / P'^3 lies in F_q, so the test is **V_{n/3}(y + 1/y) = 2** for the Lucas sequence V:
     2 products for P'^3, one inversion, 2 products, then a ladder of 2 products per bit. R1 instead computes
     `x^((q-1)/2) mod u` (12-15 products per bit).
   - root finding (only on splits): when 9 does not divide n, gamma = y^e with e = 3^{-1} mod n/3 is a cube root of
     y, and a root of the depressed cubic is t1 = +-sqrt(-P' (V_e(y + 1/y) + 2)) / 3 (Lucas ladder, one square root
     counted product by product: `(q+1)/4` power or Tonelli-Shanks), the sign by evaluating u (2 products), the
     other two roots from the quotient quadratic by the incumbent's quadratic formula and charge.
   - when 9 divides n (q = 8209 at b = 13, q = 65537 at b = 16) or P' = 0, root finding falls back to R1's route
     (`x^((q-1)/2) mod u`, counted, then `roots_fast`). The splitting test is the same in every residue class.
   - the task's suggestion to share one exponentiation between the Legendre symbol and sqrt(disc) does not apply:
     the Lucas test needs no sqrt(disc). The per-q constants (1/3, Tonelli-Shanks data) are counted and booked to C1.
   - **Shadow**: C1 runs beside R1 on every trial under saved counters; the splitting decision and the root set
     are asserted equal on every trial of every job (the job aborts otherwise). Arm `<arm>+C1` replaces R1's
     splitting and root-finding counts by C1's.
3. **C2 (batched walks).** `--walks W` runs W walks in lockstep (walk 0 is round 3's walk; walks 1..W-1 start
   from further draws of the stream; restarts per walk as round 3). R1's Cantor addition is split at its one
   inversion; the W generic additions of a lockstep step share one inversion by Montgomery's trick (3(W-1)
   products + 1 inversion), and each batched inverse is checked (`r q4 z = 1`), so every sum equals R1's.
   Arms: `<stream>:<la>` is the W-walk collection charged with R1's unbatched additions (an independent walk of the
   incumbent algorithm); `+C2` charges the batched additions actually computed; `+C1C2` both changes.
4. **C3 (oversampling and filtering).** An LA arm `lpsgebw2o<pct>[c]@<cx>` checks first when
   `ceil((1 + pct/100) nfb) + 10` rows are collected (pct = 0: round 3's `nfb + 10`), on the same walk. `c` adds
   clique removal before SGE (components of rows linked through weight-2 columns, largest first, a quarter of the
   surplus per round, singleton cascades, down to an excess of max(20, 1% of columns); bookkeeping only, not
   charged, like the incumbent's singleton filter). The SGE constant cx is part of the arm. The arm's relation
   cost is the collection cost up to its own check.

## Training (curves 5, 6, 7 only)

- **T-W (C2's W).** Relation phase only (`--las none --max-checks 1`) at b = 18, kappa_R2(18), W in {1, 8, 32, 128}.
  Rule: W\* = the smallest W whose batched relation cost (geometric mean over the 3 curves) is within 0.5% of the
  smallest.
- **T-C3/C4.** Stream `w` (W\* walks) at b in {16, 18, 20, 21}, kappa = kappa_R2(b) x {0.8, 1, 1.25, 1.6}, with
  the LA arms of amendment D1's grid (fixed from the probes before any training job). Every job reports R1, C1, C2
  and C1+C2 for every LA arm.
  - **C3 rule**: at the grid point kappa_R2(b), the LA arm (R1 primitives) with the smallest geometric-mean cost
    ratio to `lpsgebw2@cx(b)` over the sizes and curves is C3's arm. If that ratio is not below 1, C3 is dropped and
    reported.
  - **C4 rule**: among the `+C1C2` arms, for each LA arm and size, a least-squares quadratic of ln(geometric-mean
    cost over the curves) in ln kappa gives kappa\*(b) (clamped to the grid; grid argmin if not convex) and its
    predicted minimum. The LA arm with the smallest sum of predicted minima over the sizes is C4's arm; its
    schedule is the least-squares line ln kappa_C4(b) = alpha + beta (b - 20) through its four kappa\*(b).
- The training tables are produced by `examples/ic_improve4_01_analyze.rs` (`wtune`, `tune`, `pertrial`).

## Test (held-out curves)

Instances: q = next prime above 2^b + 1, seed `20261014*1000 + 10 b + curve` (as rounds 2-3).

| stream | walks | kappa | LA arms | arms produced | role |
|---|---|---|---|---|---|
| `r2` | 1 (round 3's walk) | kappa_R2(b) | `lpsgebw2@cx(b)`, C3's arm | `r2:lpsgebw2` (incumbent), `+C1`; C3 arm, its `+C1` | incumbent; C1; C3; C1+C3 |
| `w` | W\* | kappa_R2(b) | `lpsgebw2@cx(b)` | `w:lpsgebw2` (walk control), `+C2`, `+C1C2` | C2 |
| `c4` | W\* | kappa_C4(b) | C4's arm | `c4:<arm>+C1C2` (C4), and its R1/C1/C2 variants (decomposition) | C4 |
| `aa` | 1, independent stream | kappa_R2(b) | `lpsgebw2@cx(b)` | `aa:lpsgebw2` | A/A |

C3's arm uses cx as tuned where its name fixes it; otherwise cx(b).

Sizes and order (budget priority), each job alone under the caps:
1. `r2` on curves 0-2 at b = 12..21, then b = 22;
2. `c4` on curves 0-2 at b = 12..21, then b = 22;
3. `w` on curves 0-2 at b = 12..21, then 22;
4. `aa` on curves 0-2 at b = 12..21;
5. `r2` and `c4` on curves 3-4 at b = 18..21.

A size enters an arm's fit only with 3 or more correct curves for that arm and the incumbent. Top size = 22 where
the arm and the incumbent completed 3 curves there, else 21. A job over 60 minutes or 4.5 GB is dropped and
reported, with where the machine was the limit.

## Metric and analysis

`analyze` of [`examples/ic_improve4_01_analyze.rs`](../../../examples/ic_improve4_01_analyze.rs) (round 3's
analysis): 2000 bootstrap resamples of curves within size, paired across arms on (b, curve, target):

- exponents: OLS of ln cost on ln q, every solve a point, over 12..top and 16..top, for the whole solve and both
  phases; Delta against the incumbent;
- local exponents between consecutive sizes, their trend per bit and the paired trend difference;
- paired ratios per size (geometric mean of per-curve ratios) with CIs;
- per-trial cost by primitive with A + B log2 q fits (`pertrial`), for R1, C1, C2, C1+C2.

## Decision rule

Noise floors from the A/A arm (`aa:lpsgebw2` against `r2:lpsgebw2`, 12..21): eps_e = the larger absolute end of its
exponent-Delta CI; eps_c = the largest |ratio - 1| over its sizes. C1 and C3 run on the incumbent's walk, so their
ratios are exactly paired; C2 and C4 run on other walks and are judged against these floors.

- **Exponent** (whole, 12..top): ADVANCES if the Delta CI lies below 0 and the point is below -eps_e; REGRESSES
  symmetrically; else MATCHES.
- **Constant** (whole-solve ratio at the top size): ADVANCES if the CI lies below 1 and the point is below
  1 - eps_c (for the exactly paired C1 and C3: below 1); REGRESSES symmetrically; else MATCHES.
- **Primitive** (C1: splitting test + root finding per trial; C2: walk addition per trial): ratio of counted cost
  per trial against R1, on the same trials (C1) or the same additions (C2); also the change in B.
- **TRADES** when levels or phases disagree in direction.
- **Toward 4/3**: yes if the paired trend difference of the whole-solve local exponent has its CI below 0 and the
  candidate's local exponent on the top interval is closer to 4/3 than the incumbent's.

## Controls and gates

- **Gate 0**: every `r2:lpsgebw2` test record must equal round 3's `r2:lpsgebw2+R1` record on the same (b, curve):
  k, cost, relation cost, LA counters and steps. A difference stops the item.
- **Shadow equality** (C1): asserted on every trial of every job.
- **Batched inverses** (C2): `r q4 z = 1` asserted for every batched addition.
- **A/A**: `aa:lpsgebw2` against `r2:lpsgebw2`; its exponent-Delta CI should cover 0. The `w:lpsgebw2` walk control
  against the incumbent is a second, W-walk A/A comparison.
- **Correctness**: every solve's k is verified (`k phi(G) = phi(T)`) and equal to the planted k; a wrong or missing
  k fails the arm.

## Amendment D1 (2026-10-15 23:30, before any T-C3/C4 training job and before any test job)

- **Probes** (training curve 5, `r2` stream, b = 18 and 20; `results/smoke/`): oversampling alone lowers the LA
  cost by 5-11% for +2-6% relation cost; raising cx from 8 to 16-32 helps once rows are oversampled; clique removal
  before SGE helps further (b = 20: `lpsgebw2o50c@16` LA 0.76x, relation 1.06x, whole 0.95x of `lpsgebw2@8`), and the
  LA cost was still falling at the largest oversampling tried. W = 1 and W = 32 needed the same trials (b = 18).
- **T-W** ([wtune.md](results/train/wtune.md)): batched relation cost / R1 = 0.9527 (W = 8), 0.9477 (32), 0.9467
  (128); the smallest is 1.1446e9 (W = 32), W = 8 is 0.8% above it, so **W\* = 32**.
- **T-C3/C4 grid**: kappa = kappa_R2(b) x {0.8, 1, 1.3, 1.7} at b = 16, 18, 20, 21, curves 5-7, stream `w` with
  W = 32; LA arms `lpsgebw2@cx(b)`, `lpsgebw2o35c@16`, `lpsgebw2o70c@16`, `lpsgebw2o100c@16`, `lpsgebw2o70c@32`.
  The C3 and C4 rules of the protocol pick among these five.

## Amendment D2 (2026-10-16 01:08, training complete, before any test job)

- **Container restart (00:2x).** Every running job was killed: the three b = 21, kappa x 1.7 training jobs (no
  record written; rerun). The nine b = 21 records written before were complete and are kept. The binary was
  rebuilt from the final source (doc comments, `rustfmt` and headers changed since the training binary
  `d3e1d0c5…`; new hash in [binary_sha256_final.txt](results/repro/binary_sha256_final.txt)); it reproduced two
  earlier records exactly ([results/repro/](results/repro/)). All later training jobs and every test job use it.
- **Scheduler fix (00:03).** `run_slots.sh` counted the `timeout`/`prlimit` wrapper as a job, so only 2 jobs ran at
  a time between about 23:40 and 00:03; the pattern is now anchored on the binary. No record is affected.
- **Training** ([tune.md](results/train/tune.md)): 48/48 jobs, every log verified.
  - **C3**: `lpsgebw2o70c@32` (ratio 0.9331 to `lpsgebw2@cx(b)` at kappa_R2, geometric mean over 16, 18, 20, 21).
  - **C4**: arm `lpsgebw2o100c@16+C1C2`, **ln kappa_C4(b) = 0.3478 - 0.1916 (b - 20)**; kappa\*(16) = 3.235 sits at
    the top of the grid (clamped), the others inside (1.895, 1.382, 1.237 at 18, 20, 21).
  - Stored in `d2.env`. W = 32.
- **Order**: inside each priority block the largest size runs first (slot packing); `w` at b = 22 runs last and only
  if time allows.
