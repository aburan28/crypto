# Item 03 protocol: trace-zero pre-filter and free relations on the compact-orbit IC

Frozen 2026-10-14, before any arm data. Rules: [`../RULES.md`](../RULES.md).

## Incumbent

The compact-orbit / pair-table pipeline behind the n = 41–61 Koblitz complete
solves: `examples/koblitz_orbit_dlp_fast_online.rs`, run as
`construct:<n>:0:<K> <targets> <rank_seed> <out>` with default environment
(`KIC_PAIR_SYMMETRY` unset = 0, no Bloom, no rank probe cap). This is the
configuration of the retained-setup n = 61 row (`construct:61:0:500`, LEADERBOARD
"13.9x") and of item 09's IC arm. Its phases: base (x-scan, cofactor
projection, signed Frobenius orbits), normal basis, index (every ordered pair
state `(l, r, rel)`, both S3 roots, canonical-rotation root table), guided rank
stage (decompose `[a]G - R_j` for the first pivotless column `j` until rank K),
linear algebra, and one online target interval per target.

K per size is the value the repository already uses: n = 41 K = 255,
n = 53 K = 220, n = 61 K = 500 (all `a = 0`).

## Candidates

All candidate code lives in a new example `examples/ic_improve_03_tracezero_freerel.rs`,
a copy of the incumbent's `u64` path with two switches. With both switches off it
must reproduce the incumbent binary's counts and logs exactly (instrument A/A).

- **F (trace filter).** Before every decomposition attempt (rank stage and
  target stage) compute `Tr(sqrt(b)/X)` of the point to decompose; if it is 1,
  reject the attempt without running the oracle. Charged: one inversion and one
  trace per attempt.
- **R (free relations).** While the root table is built, record every
  `insert_if_absent` key collision whose stored state is not the swap mate of the
  inserting state. Verify each as a signed four-point identity
  (`P1 + P2 + P3 + P4 = O` over the base, all Frobenius shifts tried), keep
  non-zero log rows, and insert them into a second echelon. The rank stage then
  runs the incumbent's attempt stream unchanged (same scalars, same column
  choice computed from the incumbent's own echelon) and stops as soon as the
  combined echelon (free rows plus stage rows) reaches rank K. The logs are
  solved from the combined echelon. This "prefix" rule makes R's attempt
  sequence a prefix of the incumbent's, so the counted-unit saving is exact.
  Charged: the collision bookkeeping inside the index build, the verification
  lifts and the echelon inserts (timed as phase `freerel_harvest`).
- **FR.** Both.

## Theory gate for F (before it is used)

The filter is derived, then checked, before any arm uses it:

1. For any solution of `S3(x1, x2, X) = 0` on `y^2 + xy = x^3 + a x^2 + b`,
   `Tr(x1) + Tr(x2) = Tr(sqrt(b)/X)` (quadratic-solvability of `S3` in
   `p = x1 x2`). Verify on 10^6 random solutions at n = 41.
2. Generalisation to m = 3, 4 and to other subspaces: derive the condition, and
   check it exhaustively on a small field (n = 13, subspace bases inside and
   outside the trace-zero hyperplane, all m-sums for m = 2, 3, 4).
3. Hypothesis check on the incumbent: count base representatives with
   `Tr(x) != Tr(a)` and targets (rank-stage and online) with `Tr(sqrt(b)/X) = 1`.
   Positive control: on uniformly drawn curve abscissae (item 08's target law) the
   filter must reject 50% ± 1% (binomial 99% interval).

If item 3 finds the incumbent base is not inside a trace hyperplane, the cost of
switching to such a base (yield and relation count) is measured as an extra arm.

## Paired design

- Sizes: (n, K) = (41, 255), (53, 220), (61, 500); curve `KoblitzCurve::new(0, n)`.
- Targets: the first 8 public points of item 09's frozen files
  `research/ic_next10_20261014/09-equal-advice/data/targets_n<n>.txt`; every
  recovered log is checked by `[d]G = Q` in-process and against the held-back
  scalar file afterwards.
- Replicates: rank seeds 1..10 at every size. A replicate = (n, seed); the same
  8 targets in every replicate.
- Arms per replicate: `inc` (incumbent binary), `A` (new example, switches off),
  `F`, `R`, `FR`; run in that order per seed, one process each, under `nice`,
  `RAYON_NUM_THREADS=1`, `prlimit --as=3 GiB`.

## Metrics

- **Complete-solve cost per target** (primary): setup (base + basis + index,
  plus harvest for R) + rank stage + linear algebra + that target's online
  interval, in wall ms from the binaries' own phase timers. Per replicate the
  value is the mean over its 8 targets.
- **Counted unit**: S3 root calls (`index_pair_candidates + rank_root_calls +
  target root_calls`), reported per phase beside wall ms. Filter checks,
  harvest collisions and harvest lift additions are reported as separate
  counts.
- Per-phase wall ms: index, harvest, rank, LA, target online.

## Analysis and decision rule

Paired ratio candidate / `A` per replicate; ratio of means over the 10 seeds
with a 95% percentile bootstrap (10,000 resamples of seeds), per size and
pooled over sizes. Level: **constant** (paired ratio); also **primitive**
(rank-stage cost per needed rank-stage relation is unchanged, the number of
relations falls by the free rank).

- **ADVANCES**: whole complete-solve ratio CI entirely below 1.
- **TRADES**: one phase better with CI below 1 and another worse with CI above 1,
  whole-cost CI covering 1.
- **MATCHES**: whole-cost CI covers 1 and no phase separates.
- **REGRESSES**: whole-cost CI entirely above 1.

The verdict is given per lever (F, R, FR) and per size; the pooled verdict uses
all sizes.

## Controls

- **A/A null**: `A` against `inc` must have identical rank attempts, root calls,
  relations, logs and recovered scalars on every replicate; its wall ratio CI
  must cover 1 (it measures host noise).
- **Positive control, R**: the harvest must reproduce item 09's census rank at
  n = 61 (K = 256: 3, K = 360: 5, K = 600: 49) within the census's own
  definition (rank of verified non-zero rows), and the R arm's rank-stage
  root-call reduction must equal the incumbent's tail cost over the last
  `F_rank` attempts exactly (it is a prefix).
- **Positive control, F**: the 50% rejection on uniform curve abscissae (theory
  gate 3).

## Correctness gate

Every arm: every target's recovered scalar verified `[d]G = Q` and equal to the
held-back scalar; the candidate's base logs equal the incumbent's logs on every
column; every free relation verified as a group identity. Any disagreement
fails the arm.

## Budget

3 GB address space per job, one thread, < 60 min per job, about 4 h total.

## Amendment 1 (2026-10-14, after two smoke runs, before any arm data)

Smoke runs (n = 41, 53, 61; 2 targets; seed 7; not part of the data) showed:

1. **The prefix rule wastes the free rows.** At n = 41, K = 255 the harvest has
   free rank 161, but the prefix rule still needed 208 of the incumbent's 255
   stage relations: the incumbent's guided rows target columns the free rows
   already pin (114 rows gained nothing in the combined echelon). Lever R is
   therefore run with the **natural rule**: the guided column choice reads the
   combined echelon (first column with no pivot there), so every stage relation
   raises the combined rank and exactly `K - free_rank` are needed. The prefix
   rule stays available (`IC03_PREFIX_RULE=1`) and is reported from the smoke run
   only. The attempt streams of `A` and `R` share the scalar sequence but target
   different columns, so the rank-stage pairing is statistical, not exact; the
   expected rank-stage ratio is `(K - free_rank)/K`, checked against the measured one.
2. **Wall clock on this host is too noisy for the effect size.** Two identical
   n = 61 index builds took 16.6 s and 12.0 s. The verdict is therefore decided on
   the **counted unit** (the RULES' unit), and wall ms is reported beside it:
   `C = index_pair_candidates + rank_root_calls + target_root_calls +
   filter_checks + harvest_additions` (one S3 root solve, one filter inversion and
   one harvest point addition each count 1; an addition is cheaper than a root
   solve, so this over-charges the candidates). Same decision rule, applied to
   the C ratio; the wall-ms ratio must not contradict it (its CI must not lie
   entirely on the other side of 1).
3. Replicates raised to 12 seeds (1..12) per size, since the rank-stage pairing
   is now statistical.
