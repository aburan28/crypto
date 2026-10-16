# Item 02 protocol: free relations in the production compact-orbit IC, K tuned per arm

Frozen 2026-10-15, before any training or test data. Rules: [`../RULES.md`](../RULES.md).
Prior work: [`../../ic_improve_20261014/03-tracezero-freerel/`](../../ic_improve_20261014/03-tracezero-freerel/README.md).

## Incumbent and candidate

- **Incumbent (`inc`).** `examples/koblitz_orbit_dlp_fast_online.rs` as committed at
  `45946b054` (blob `ca551bb`), run as `construct:<n>:0:<K> <targets> <seed> <out>`
  with the default environment. The binary is built from the committed file before
  any edit and frozen (sha256 in `data/binaries.sha256`).
- **Candidate (`free`).** The same file with the free-relation harvest added behind
  `KIC_FREE_RELATIONS=1` (default off). Root-table key collisions between pair states
  that are not swap mates are recorded during the insertion the incumbent already
  does, verified as signed four-point identities by point additions, and inserted as
  homogeneous rows into the rank echelon before the guided rank stage. The guided rule
  is unchanged (first column with no pivot), so it fills only the columns the free
  rows leave open (the "natural rule" of the prior item's Amendment 1). Binary frozen
  after the A/A check.
- K is a free parameter of both arms. The repository's configured values are
  `K_def` = 255 (n = 41), 220 (n = 53), 500 (n = 61).

## Metric

Per run (one process, 8 targets, the first 8 public points of item 09's frozen
`targets_n<n>.txt`), the **complete-solve cost per target**:

- **Counted unit C** (verdict): `index_pair_candidates + rank_root_calls_total +
  mean over targets of root_calls + harvest_additions`. One S3 root solve and one
  harvest point addition count 1 each (over-charges the candidate).
- **Wall ms** (reported beside, must not contradict): `setup_total + rank_stage +
  linear_algebra + mean target online_ms`, from the binary's own timers.
- Per phase: index, rank, target, harvest; rank attempts; free rank; peak RSS.

## Stage 1: A/A (default unchanged)

The new binary with `KIC_FREE_RELATIONS` unset against the frozen incumbent binary,
both with `KIC_DUMP_RANK`, on seed 101 at every training (n, K) and on seeds 1, 2 at
each `K_def`. Pass condition: after deleting wall-clock, RSS and hardware-counter
fields only, the summary line, every target record and the full rank trace are
byte-identical. Any difference stops the item.

## Stage 2: training sweep (K selection)

- Grids: n = 41: K ∈ {64, 80, 100, 125, 160, 200, 255};
  n = 53: K ∈ {160, 200, 240, 280, 330, 400}; n = 61: K ∈ {400, 450, 500, 550, 600}
  (600 is the last grid point under the 3 GB cap; K ≈ 620 is the wall).
- Training seeds 101–104, both arms, interleaved (inc then free per seed and K).
- Selection: per (n, arm), `K*` = argmin over the grid of the mean C over the 4
  training seeds. A K with any unsolved or unverified target, or any run killed by
  the memory cap or the time limit, is ineligible. `K*` is written to
  `data/frozen_K.json` before any test run; test seeds are never used to choose K.

## Stage 3: held-out test

- Seeds 1–16 at n = 41 and n = 53, seeds 1–12 at n = 61 (disjoint from training).
- Per seed, interleaved, one process each: `inc@K_inc*`, `free@K_free*`, and for the
  same-K comparisons `inc@K_free*`, `free@K_inc*`, `inc@K_def`, `free@K_def`
  (each distinct (arm, K) once).

## Comparisons (paired by seed; same targets everywhere)

- **(a) Same K.** `free@K / inc@K` for each K in {K_inc*, K_free*, K_def}.
- **(b) Tuned vs tuned (primary, deployment).** `free@K_free* / inc@K_inc*`.
- **(c) Against today's configuration.** `free@K_free* / inc@K_def` (reported, no verdict).

Statistic: ratio of means over seeds, 95% percentile bootstrap over seeds (10,000
resamples), per size and pooled over sizes (pooling = resampling seeds within each
size, ratio of summed means).

## Decision rule (on C; level = constant)

- **ADVANCES**: CI entirely below 1. **REGRESSES**: CI entirely above 1.
- **TRADES**: CI covers 1 while one phase (index or rank) separates in each direction.
- **MATCHES**: CI covers 1 and no phase separates.
- The wall-ms ratio must not have its CI entirely on the other side of 1; if it
  does, the verdict is reported as contradicted.
- Verdict per size and pooled, for (b) and for each (a).

## Correctness gate

Every run: all 8 targets recovered, verified `[d]G = Q` in process, and equal to
item 09's held-back scalars. Candidate runs: `free_relations.failed = 0`, rank K
reached, and `logs_hash` equal to the blake3 hash of the incumbent's base logs from
the A/A rank trace at the same (n, K) (the logs are unique). Any failure fails the arm.

## Controls

- Null: the A/A of stage 1 (exact identity, not a ratio).
- Positive: at `K_def` the candidate must reproduce the prior item's free ranks
  (161, 4, 25 at n = 41, 53, 61) and its rank-stage attempts `K − free_rank`.

## Optional extension (measure only)

If time allows: a census, outside the production binary, of root-table matches under
the 2-torsion translation `X ↦ sqrt(b)/X` (the abscissa of `Q + T2`), reported as
counts of matches and of verified relations.

## Budget

One thread per job under `nice`, at most 2 jobs at once, `prlimit --as=3 GiB`,
`timeout 3000` s per job, about 4 h in total.
