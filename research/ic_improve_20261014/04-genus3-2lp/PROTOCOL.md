# Protocol: item 04, genus-3 double-large-prime IC after the cubic GHS transfer, beyond q = 2^15

Frozen 2026-10-14 before any recorded data. Round rules: [`../RULES.md`](../RULES.md).
Prior work: [`../../ic_next10_20261014/07-cubic-ghs-genus3/`](../../ic_next10_20261014/07-cubic-ghs-genus3/)
(Python tool `ghs3.py`, run `run-20261014`).

## Question

On the cubic-GHS split-arm instances of item 7 (`E0 / F_{q^3}` transferred to the
Jacobian of a genus-3 hyperelliptic curve `C / F_q`, `deg h = 7`), does genus-3
index calculus with double large primes (IC-2LP) improve on the full-factor-base
index calculus (IC-basic) that it replaces, at the exponent level, over a common
size range; and does the IC-2LP exponent settle at the Gaudry-Thome-Theriault-Diem
value `2 - 2/g = 4/3` or keep rising?

## Instrument

`examples/ic_improve_04_genus3_2lp.rs`: a line-by-line Rust port of `ghs3.py`
(instance build, transfer, Cantor arithmetic over `F_q` and over `K = F_{q^3}`,
smoothness test, Cantor-Zassenhaus roots, large-prime graph and cycle handling,
singleton filtering, Lanczos mod `ell`), with the **same operation counters**
and the **same random stream** (a port of CPython's MT19937 `random.Random`,
`randrange`, `sample`), so the same seed gives the same instance, the same walk
and the same counted cost. Group order by PARI `ellcard` (SEA) on `E0 / K`, as in item 7.

Gate 0 (cross-check, before any new arm is run): with `--mode python` the port
must reproduce item 7's recorded rows (`ic_b{8..12}_c{0,1,2}.jsonl`, arms ic2lp
and basic): identical `N`, `ell`, `h`, recovered `k`, walk steps and counter
vectors, hence identical `cost`. Any mismatch is fixed before data, or reported
and the item stops.

Repository modules (`src/prime_hyperelliptic`, `hyperelliptic_index_calculus`)
are BigUint-based with checked, uncounted arithmetic; they are not used because
the counted unit must match item 7 operation by operation. Recorded as a
deviation from "reuse", with this reason.

## Arms (all on the same instance, same `G`, same target `k`, same walk seed)

Paired mode: after the instance, `G` and the item-7 checks, one `k` is drawn;
every arm then starts from a copy of the same generator state, so the walk
set-up draws are identical across arms.

1. **Incumbent: IC-basic** (item 7's full-factor-base control). Factor base =
   all degree-1 primes, no large primes.
2. **Candidate: IC-2LP** as frozen in item 7: factor base = the
   `ceil(2 q^{2/3})` degree-1 primes with smallest `x`, relations with up to two
   large primes, spanning-tree elimination.
3. **Middle point: IC-1LP, tuned.** Same code, relations with at most one large
   prime kept (two-large-prime relations dropped), factor base
   `ceil(kappa1 q^{5/7})` (Theriault's single-large-prime balance for g = 3,
   cost `q^{10/7}`).
4. **Optional: IC-2LP-tuned.** Factor base `ceil(kappa2 q^{2/3})`.
5. **A/A null:** IC-2LP rerun with an independent walk seed (seed + 1), same
   instance and target, at every size where IC-2LP runs.
6. **Reference only: rho** on `E0(K)` (short Weierstrass, affine, same `K`
   arithmetic and counters as item 7), r = 64 independent random multipliers,
   negation map (canonical `y`), fruitless 2-cycle escape by doubling,
   distinguished points. Reported as a column; no verdict depends on it.
   Where a rho run would exceed the job cap, the rho model
   `(cost per step) * sqrt(pi ell / 4)` is reported instead and labelled.

Tuning (arms 3 and 4): training curves `curve = 5, 6, 7` (seeds disjoint from the
test curves) at b = 12, 13, 14; grid `kappa in {0.5, 1, 2, 4}`; the `kappa` with
the lowest geometric-mean cost at b = 14 is frozen and recorded before any test
run of that arm. Training rows are never used in the fits.

## Sizes and instances

`q` = next prime above `2^b + 1`, seed `20261014*1000 + 10 b + curve`
(item 7's formula, so b <= 15, curves 0-2 are item 7's instances).
Test curves 0, 1, 2 (3 per size; curve 3 added where the budget allows).
- IC-2LP, A/A, IC-1LP, IC-2LP-tuned: b = 10 .. 20, or up to the largest b whose
  jobs finish within the 60-minute cap under 3 GB.
- IC-basic: b = 10 upward until a job would exceed 60 min or 3 GB (its linear
  algebra is `~ q^2`); the largest completed size is reported with the reason it stopped.
- Common range = the sizes where every arm has all test curves.

## Metric

Item 7's unit: `cost = M_q + ceil(log2 q) I_q + c_ell M_ell` in
`F_q`-multiplication equivalents, all phases (transfer of `G` and `T`, walk,
smoothness tests, root finding, linear algebra). Primary `c_ell = 1`;
sensitivity `c_ell = 9`. Wall time recorded, not used for the verdict.

Correctness gate: every recovered `k` is checked as `k phi(G) = phi(T)` and
against the planted `k`; a wrong or missing answer is a failure, never dropped.
Any arm with a failure at a size is reported as failing there.

## Analysis (bootstrap: 2000 resamples of curves within each size, paired across arms)

- Per-arm exponent: least-squares slope of `log cost` on `log q` over the common
  range, 95% percentile bootstrap CI.
- Exponent difference `Delta = slope(cand) - slope(inc)`, paired bootstrap.
- Local exponents between consecutive sizes, with bootstrap CIs.
- Paired ratio cand/inc per size: geometric mean of per-curve ratios, bootstrap CI.
- IC-2LP alone over b = 14..20 (and the upper half of its range): slope with CI,
  compared with 4/3; trend in local exponents (slope of local exponent on b, CI).

## Decision rule

Exponent level (primary): ADVANCES if the 95% CI of `Delta` lies below 0;
REGRESSES if above 0; MATCHES if it covers 0. Requires >= 4 common sizes.
Constant level: the paired ratio at each common size; ADVANCES at the largest
common size if its CI lies below 1, REGRESSES above 1, MATCHES if it covers 1.
TRADES if the exponent and constant levels point in opposite directions over
the common range (e.g. lower exponent, ratio CI above 1 at the smallest size).

Settling question: "settles at 4/3" if the upper-range slope CI covers 4/3 and
the local-exponent trend CI covers 0; "keeps rising" if the trend CI lies above
0 or the upper-range slope CI lies above 4/3.

Controls: the A/A arm's exponent difference and per-size ratio CIs must cover
0 and 1; the positive control is item 7's separation of IC-2LP from IC-basic
(q^1.30 vs q^1.65 at b = 8..13), which the port must reproduce exactly (gate 0).

## Budget

Memory cap 3 GB per job (`prlimit --as=3000000000`), single-threaded, `nice`,
60 min per job, about 4 h for the item, shared machine (4 cores, four items).

## Amendment A1 (2026-10-14, after tuning, before any test-curve data)

- Gate 0 passed: the port reproduces item 7's records exactly, 42/42 solves
  (b = 8..15 for IC-2LP, b = 8..13 for IC-basic, curves 0-2): same instance,
  `k`, steps, counter vectors and cost ([results/xcheck/xcheck.md](results/xcheck/xcheck.md)).
- Tuning ([results/tune/tune.md](results/tune/tune.md)), training curves 5-7,
  lowest geometric-mean cost at b = 14: **kappa1 = 2 frozen for IC-1LP**
  (1.59e8 against 2.15e8 for 4 and 3.19e8 for 1); **kappa2 = 2** for the tuned
  2LP arm, which is item 7's value, so IC-2LP-tuned is the IC-2LP arm with the
  same seed and is not run separately.
- A 2^16 timing run (curve 0, not part of the record) was used to size the jobs.
