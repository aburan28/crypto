# N83 precomputed-table rho control: protocol (frozen 2026-10-16, before timing)

## Question

The N83 compact-orbit stream reports index calculus (IC) beating Pollard rho
online on the a=1 Koblitz curve over `F_2^83` (prime subgroup
`r = 8569786107849059`, cofactor 1128547018, curve
`EC1N83Ckb1a1haacb69ee9f37`):

- **12.29×** (median, 95% CI 3.26–32.21): hosted run `gha-38597929540-1`,
  K600 IC (`koblitz_orbit_dlp_fast_online.rs` at `eaa75494`) against strong
  rho rung 3, 16 public hash-to-curve targets
  ([evidence](../sat_factor_base_review_20260908/autolab/evidence/20261014-koblitz-n71-n83-hosted-single-target/README.md)).
- **113.011×** (median, 95% CI 40.2–148.8): Tier 1 run 38660150265, K3300
  root11 IC (`koblitz_orbit_dlp_fast.rs` at `d0e4aa13`) against the fold-backend
  rho worker, five reserve targets
  ([result](../n83_tier1_refine_reserve_results_20261015/RESULT.md)).

In both, IC's online interval starts after 600–760 s of reusable,
target-independent setup holding a multi-GiB state table (K3300: 628 s setup +
50.6 s 4-thread rank stage, peak RSS 10.79 GB; K600: 758 s, RSS 5.8 GB), while
every rho reference starts each target with an empty distinguished-point
table. IC's online phase is a table lookup over that precomputation. The
matched reference for a setup-excluded online interval is rho with
precomputation (Bernstein–Lange, "Computing small discrete logarithms faster",
INDOCRYPT 2012): the same reusable budget spent on a table of `G`-only
distinguished points. AGENTS.md ("Distance to rho") requires this reference
where setup is excluded. It had not been run.

## Arms (all on one hosted runner per job, same 16 targets)

Targets: the 16 public evaluation points of `gha-38597929540-1`
(`inputs/eval_w000.target` … `eval_w015.target`, in that order),
`targets.jsonl` SHA-256
`7c806b865814a21ece7591a111e2f09ce8eaffaeafab847515e24f33ddc457ce`.
No scalar enters any producer. `reference_scalars.txt` (the hosted run's
recovered scalars) is read only by the analyzer after all arms finish.

| arm | producer | precomputation | online |
|---|---|---|---|
| `ic` (job k3300) | frozen `koblitz_orbit_dlp_fast.rs` @ `d0e4aa13` (SHA-256 `5e9fd826…`), `construct:83:1:3300`, Tier 1 root11 environment, rank seed 29, `KIC_RANK_THREADS=4`, `KIC_TARGET_THREADS=1`, CPUs 0–3 | in process, once | per target, its own clock |
| `ic` (job k600) | frozen `koblitz_orbit_dlp_fast_online.rs` @ `eaa75494` (SHA-256 `4729f3b9…`), `construct:83:1:600`, rank seed 7, `RAYON_NUM_THREADS=1`, CPU 3 | in process, once | per target |
| `rho_plain_fold` | frozen `n83_rho_optimization_eval.rs` @ `d0e4aa13` (SHA-256 `fb14072f…`), `research/n83_dense_root_20261013/rho_config.json` (fold, 64 lanes, dp 6), seed `2026101600+i`, CPU 3, one process per target | none (empty table) | per target |
| `bl_a_time_matched` | `examples/n83_precomputed_rho.rs` (this PR), `configs/bl_a_time_matched.json` | 1 thread on CPU 3; table walks stop **650 s** after setup start (≤ either IC's target-independent wall: K3300 ≈ 680 s with 4 rank threads, K600 ≈ 758 s); dp_bits 7 | 32 lanes, per target |
| `bl_b_memory_matched` | same producer, `configs/bl_b_memory_matched.template.json` | 4 threads on CPUs 0–3; table bytes (16 B/entry) = **this job's IC process maximum RSS** (GNU `time -v`), an upper bound on IC's state-table memory; safety cap 12,600 s of walks; dp_bits 6 | 4 lanes, per target |

Both precomputed arms use the repository's fastest rho step: the fold field
backend (`Gf2_128`, two-CLMUL tail folds) of `n83_rho_optimization_eval.rs`
with the production `koblitz_strong_rho` jumps (32 `G`-only), partition,
distinguished-point hash, signed-Frobenius canonicalization (166-element
quotient) and batched slope inversion. The new driver
`StrongRhoG::dp_walks` differs from `solve` only in (i) supplying walk starts
from a closure (`c·G` for the table, `c·G + Q` online) and (ii) escaping a
fruitless cycle to `[2]·min(current, next)` so two walks caught in the same
2-cycle leave it alike, which a shared table needs.

### Online interval (precomputed arms)

From before the target's on-curve and subgroup checks to after the fast
`[d]G = Q` recovery check: target validation, `Q`-dependent walk starts,
walks to a distinguished point, table lookups (binary search over the sorted
fingerprint table, plus a per-target table of the target's own distinguished
points), and recovery. This includes a subgroup check IC's interval does not
include. Excluded, as IC excludes its setup: the table build (reported in the
setup record with wall, CPU, steps and bytes) and the generic-library replay.

### IC change from the original runs (declared)

Each original run started one IC process per target. Here one IC process
receives all 16 targets: setup once, then each target's own online clock
(the producers' multi-target mode, same code path). This changes nothing
inside the online interval; it is what makes IC's setup reusable, the premise
being tested.

## Verification

`examples/n83_precomputed_rho_analyze.rs` replays every recovered scalar of
every arm with the generic library curve (`KoblitzCurve::mul`, independent of
all producers), requires all four arms to agree target by target, and checks
the hosted run's recorded scalar. Any failure is `FAIL`; failed or capped
targets stay in the records.

## Statistics

Per arm: mean and median online ms with 95% percentile bootstrap intervals
(10,000 target resamples, fixed seeds). Pairs (IC/arm): median of per-target
ratios, ratio of means and ratio of medians, each with a target-resampled 95%
interval. Operation counts priced at 0, 0.1 and 1 GAE per unpriced native
operation (AGENTS.md §5a) for every arm that emits them.

## Decision rule (verdict, fixed now)

For each IC configuration (K600, K3300) and each precomputed arm:

- **Claim survives** against that arm if the 95% interval of the median paired
  ratio IC/arm lies entirely **below 1** (IC faster online).
- **Claim does not survive** if that interval lies entirely **above 1**
  (precomputed rho faster online under matched precomputation).
- Otherwise **unresolved**.

The 12.29× / 113× "IC beats rho" readings survive only if they survive the
time-matched arm (a); arm (b) is reported beside it. The plain-rho arm is the
positive control: it must reproduce IC ≪ plain rho on this runner, or the job
is not comparable to the claims.

## Prediction (recorded before the run; a derivation, not a measurement)

`N = r/166 ≈ 5.16·10^13` orbit classes. A table built with `P` walk steps of
mean length `W = 2^d` covers about `P` classes, so an online solve needs about
`N/P + W` steps plus lane overshoot. At an assumed 3–5·10^6 fold steps/s on
one EPYC 7763 thread, arm (a) has `P ≈ 2–3·10^9`, predicting about
`2–3·10^4` online steps (≈ 5–15 ms). Arm (b) at ≈ 6.7·10^8 (K3300) or
≈ 3.6·10^8 (K600) entries with `W = 64` predicts about `1.3–2.3·10^3` steps
(≲ 1 ms). IC online medians in the claims: K600 557.5 ms, K3300 21.1 ms. If
the prediction holds, precomputed rho is faster online than both IC
configurations.

## Limits

Public synthetic points only; no key recovery, asymptotic, or deployed-curve
claim. Wall time on hosted `ubuntu-24.04` runners (AMD EPYC 7763 class, 4
vCPU, shared tenancy); `isolated_bench.py` reservation is not used, so the
isolation tier is L0/cloud-shared, recorded with `lscpu`, `meminfo` and the
cgroup CPU limit. Both jobs run the plain-rho control so each job's ratios
are internal to one runner. Precomputation is target independent and reused
across all 16 targets in both IC and rho; this compares online cost at
matched reusable cost, not cold cost.
