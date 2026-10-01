# n=61 compact orbit vs strong rho, one target per workload (2026-10-01)

Status **`PENDING_INDEPENDENT_VALIDATION`**. No ledger row changes. The replay
ran on the same host, the walls come from a shared, heavily loaded macOS host
(not AGENTS.md §10 evidence), and there is no independent-host validation.

This implements option A of
[`../20261001-koblitz-n61-compact-orbit-panel/VS_RHO_SINGLE_TARGET_OPTIONS.md`](../20261001-koblitz-n61-compact-orbit-panel/VS_RHO_SINGLE_TARGET_OPTIONS.md)
with the option C plumbing.

- Each workload is one previously unseen public hash-to-curve point on
  `EC1N61Ce0hd8c15448081e` (a=0, n=61, 47-bit prime subgroup).
- Each workload uses a fresh process per arm.
- The IC's reusable setup is outside its online interval.
- The comparator is strong rho rung 3 (`rho_reference_minimum`), not KS v2.

## Run

`20261001T200614Z-d687d78737`, from
`boundary_autolab.py launch-single --beat koblitz.compact_orbit.n61_single_target`
at `34f70a32b`. It covered 64 eval workloads and 4 disjoint tune workloads per
K. The run directory is in `autolab_runs/20261001T200614Z-d687d78737/`.

| Producer | Source | What changed vs the frozen original |
|---|---|---|
| IC | `examples/koblitz_orbit_dlp_fast_online.rs` | Split timer in `extract()`. Exclusive query / PDP / relation-check / descent / recovery-check phases that sum to `online_ms`, plus explicit `target_query_begin` and `recovery_check_end` events with ns offsets. Adds a `relation_checks` counter. |
| Rho | `examples/koblitz_rho_batch_ks_strong_online.rs` | The online timer starts after Q is built. Q is either a public point (`KIC_RHO_TARGET_POINT`) or the derived fixture. Walk + collision and the `[d]G = Q` check are reported separately, with `after_target_built` and `recovery_check_true` events. |

The frozen `koblitz_orbit_dlp_fast.rs` and `koblitz_rho_batch_ks_strong.rs`
are unchanged, so every existing pin still verifies. The run's
`producer_identity` step ran each online producer and its frozen original. The
IC ran at K=400 on scalar 123456789012345; rho ran with seed 610000. Their
records are **identical** once timers, `relation_checks` and
`producer_version` are removed (`artifacts/producer_identity.json`).

## K tune (disjoint one-target workloads)

The four tune points come from the same law with role `tune`. They are
disjoint from the 64 eval points. K is picked by the lowest median
whole-process wall, as registered in `protocol.json`.

| K | Median process wall s | Median user s | Median instructions | Median online ms | Est. / measured RSS GiB |
|---:|---:|---:|---:|---:|---:|
| 300 | 89.8 | 32.4 | 221.6 G | 160.0 | 0.69 / 0.63 |
| **400** | **30.0** | 23.0 | 167.1 G | 32.1 | 1.28 / 1.23 |
| 500 | 35.5 | 21.4 | 152.1 G | 29.1 | 1.40 / 1.35 |
| 600 | 32.1 | 21.2 | 164.5 G | 6.6 | 2.55 / 2.51 |
| 700 | 64.1 | 25.6 | 183.4 G | 12.3 | 2.73 / 2.69 |
| 800 | 56.7 | 28.4 | 214.3 G | 19.8 | 4.94 / 4.89 |

The wall rule picked K=400. On user CPU and instructions, K=500–600 would edge
it by about 7–9%, so the pick rests on contended walls. The online column comes
from four targets per K and is dominated by luck, since IC probes are roughly
geometric. The new RSS model (`ic_rss_model`) over-predicts the measured peaks
by 1–8% at every K (0.9% at K=800, 8.4% at K=300).

## Result (64 workloads, K=400)

The online speedup is `rho_online_wall_ms / ic_online_wall_ms` on the same
point:

| Statistic | IC online ms | Rho online ms | Speedup (rho / IC) |
|---|---:|---:|---:|
| Median (95% bootstrap CI) | 38.6 (22.1–49.3) | 265.5 (230.9–362.2) | **9.93 (4.57–14.37)** per workload |
| Mean | 99.0 | 497.9 | 5.03 (ratio of means) |
| Ratio of medians | | | 6.88 |
| Min / max | 0.12 / 2,220.7 | 19.5 / 4,433.8 | 0.17 / 1,374 |

- **IC wins the online interval.** It is faster on 60 of 64 workloads.
- Both distributions are heavy-tailed: the median IC solve takes 95k PDP
  probes, the median rho walk 1.19 M steps. The per-workload extremes are
  lucky IC targets (161 probes) or long rho walks, so the ratio of means
  (5.0×) is the expected-cost reading.
- PDP probing takes essentially all of IC online time. Median phases: PDP 38.49
  ms, recovery check 0.046 ms, relation check 0.017 ms, descent 0.001 ms, query
  < 0.001 ms. Every workload needed exactly one relation check.

Cold, from process start to verified recovery, including the IC's setup:

| | IC | Rho |
|---|---:|---:|
| Median process wall s | 34.7 | 0.295 |
| Median IC setup s (base, rank, index, LA) | 34.6 | — |
| Median instructions per process | 167.1 G | 2.90 G |
| Peak RSS | 1.23 GiB | 14.4 MiB |

The cold ratio (IC / rho, in-process to solved) has median **117×**
(95% CI 96–154×). **IC loses the cold one-target comparison by two orders of
magnitude**, as expected: at L=1 it pays the full K²·n index build plus the
rank stage in every process. Using the means, the setup is amortized after
roughly 34.6 s / (0.498 − 0.099) s ≈ 87 targets on this host. That figure is a
rough break-even, not a measured crossover.

## Verification and claim-check

Every workload passes all of the following (`single_target_summary.json` →
`verification`):

- Both producers exit 0, and IC `group_verified` holds.
- The two recovered scalars agree.
- The IC base hash matches the dumped base.
- The IC's setup completes before its online start.
- Both out-of-process `oracle.py` replay certificates hold.
  - Both arms: digest, target equals the workload point, subgroup membership,
    and `[d]G = Q`.
  - IC also: the relation points are on the curve, the `x_codes` match, and
    the four points sum to the target.

Claim-check `--stage vs_rho`:

- Per workload (`artifacts/claims/wNNN.claim_check.json`): **64 / 64 PASS**,
  0 missing fields, 0 validation errors.
- `claim_draft.json` (workload 0) through the CLI: **PASS**
  (`artifacts/claim_check_cli_vs_rho.json`).

Each claim carries the IC1 candidate, workload and run IDs with manifest
hashes, an identical resource envelope for both arms, the replay certificate
SHA-256 digests, rho's policy (rung 3, 32 lanes, 4 distinguished bits, empty
table) and the cold metrics. `independent_validation` is true in the
same-host, out-of-process sense used by `research/ic_single_target_20260930`.
The `independent_validation_scope` field says so. `verify` re-hashes all
1,000 manifest files with 0 mismatches (`artifacts/verify.json`).

PASS means the schema is complete and the replays hold. It promotes nothing:
the run status stays `PENDING_INDEPENDENT_VALIDATION`.

## Host and timing conditions

Apple M4 Pro, 14 logical CPUs, 48 GiB, macOS 26.6, rustc 1.93.1, one thread
per arm (`RAYON_NUM_THREADS=1`), unpinned. Other agents' jobs ran throughout.
The 1-minute load before each eval workload was 13–337 (median 35), with swap
4.7–5.1 GiB used. Every timed process reported 0 swaps. IC process wall
exceeds user + sys by a median of 10.1 s, rho by 0.06 s, so the IC walls carry
heavy contention. The online intervals are short, but they are wall time on
the same host. Read them as diagnostics until the hosted isolated run
(`.github/workflows/compact-orbit-n61-single-target-isolated.yml`, dispatch
after merge) reproduces them.

## What this is not

- No ledger promotion and no key recovery. Public synthetic points only, with
  no known scalar (`target_scalar_constructed: false`).
- No asymptotic claim. One curve size, one host.
- The online win depends on reusable setup sitting outside the interval, as
  the measurement rules require. The cold ratio beside it is the price.
- Rho has no precomputed distinguished-point table. Bernstein–Lange is a
  different reference and is not measured.
- The claims' `independent_replay_pointer` names the original
  `runs/<run-id>/` path (gitignored). The same files are here under
  `autolab_runs/<run-id>/artifacts/claims/`. Editing the pointer would break
  the manifest hashes.

## Files

`autolab_runs/20261001T200614Z-d687d78737/` holds the following:

- `state.json`
- `inputs/`: protocol snapshot, `ledger_pin.json`, `targets.json`
- `artifacts/`: `k_tune.json`, `producer_identity.json`, `workloads.jsonl`,
  `certificates/`, `claims/`, `manifests/`, `single_target_summary.json`,
  `claim_draft.json`, `claim_check.json`, `claim_check_all.json`,
  `candidate.json`, `review_manifest.json`, `verify.json`
- `receipts/`
- `logs/`

Left out, as the autolab README asks: the 2.9 MB base dump
`logs/base_n61_K400.jsonl` and the ledger copy
`inputs/boundary_targets.json` (pinned in `ledger_pin.json`). The candidate
manifest keeps the base representatives, and rerunning the command regenerates
the dump.
