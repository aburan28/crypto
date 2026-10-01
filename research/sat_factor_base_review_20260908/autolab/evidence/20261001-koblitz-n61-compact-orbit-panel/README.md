# n=61 compact orbit vs KS v2 batched rho through the autolab (2026-10-01)

Status **`PENDING_INDEPENDENT_VALIDATION`**. These are multi-target batch
diagnostics against a superseded comparator, measured on a shared, unisolated
macOS host. They promote no ledger row. Continues PR #1134 (L=65,536), which
ran outside the autolab and tuned K only at L=65,536.

Both panels ran through
`boundary_autolab.py launch-panel --beat koblitz.compact_orbit.n61_panel`
(new batch-panel beat in `protocol.json`; `launch` rejects it). The run
directories are promoted under `autolab_runs/<run-id>/`.

| Run ID | L | K (tuned) | Wall ratio median (range) | User-CPU ratio | Instructions ratio | end_to_end_dlp | vs_rho |
|---|---:|---:|---|---:|---:|---|---|
| `20261001T155442Z-9990f40315` | 1,024 | 700 | 0.282 (0.276–0.306) | 0.289 | 0.132 | PASS | FAIL (closed) |
| `20261001T160852Z-1f8383bd42` | 4,096 | 1,000 | 0.262 (0.247–0.269) | 0.256 | 0.117 | PASS | FAIL (closed) |

All ratios are compact / rho on the same block pair.

## What this is and is not

- The comparator is the frozen KS v2 batched rho (`koblitz_rho_batch_ks_v2_n61`,
  `KIC_RHO_DP_BITS=4`). It is **not** on
  `measurement_schema.vs_rho.rho_reference_minimum`. Against strong rho R3 the
  n=61 L=1,024 instruction ratio is IC/rho = 2.297
  (`research/notes/index-calculus/RESEARCH_STRONG_RHO_SWEEP_PROTOCOL_20260929.md`).
  Nothing here changes that verdict.
- Each run solves one batch of L targets. That is a historical multi-target
  diagnostic, not the one-target primary comparison. See
  [`VS_RHO_SINGLE_TARGET_OPTIONS.md`](VS_RHO_SINGLE_TARGET_OPTIONS.md).
- No key recovery. No asymptotic claim: at optimal K both arms scale as
  √(L·r/n), so any lead is a constant factor.

## Host and timing conditions

Apple M4 Pro, 14 logical CPUs, 48 GiB, macOS 26.6 arm64, rustc 1.93.1
(Homebrew), single thread (`RAYON_NUM_THREADS=1`), source `09a92b266`.
`tools/isolated_bench.py` needs Linux (`sched_getaffinity`, `/proc`, PSI), so
these walls are **not** AGENTS.md §10 evidence. Each run records its load
average, available memory and swap before and after. Other agents' jobs were
running: two single-core Python jobs throughout, and the 1-minute load was
5.9–14.3. Swap was 7.3 of 8 GiB used but stable, and every timed process
reported **0 swaps**. Wall exceeds user CPU by up to 4.2 s (IC) and 17.5 s
(rho, L=1,024 block 1), which is contention. The user-CPU and retired-
instruction ratios are the less contaminated readings. Instruction counts come
from `/usr/bin/time -l` and include kernel work.

## 1,024-target K tune (gap 1)

Disjoint tune corpus `n61-ks-growing-tune-1024-v1` (sha256 `8ec0c043…`). The
eval corpus is `n61-ks-growing-1024-v1`. Both regenerate byte-for-byte to PR
#830's committed scalar files and share no scalar. K is picked by lowest
whole-process wall. All six candidates fit in memory.

| K | Wall s | User s | Sys s | Instructions | Est. RSS GiB | Measured RSS GiB | Load |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 400 | 55.6 | 50.3 | 0.7 | 470.3 G | 0.9 | 1.23 | 7.3 |
| 500 | 37.2 | 36.2 | 0.4 | 342.3 G | 1.3 | 1.35 | 11.9 |
| 600 | 40.9 | 33.6 | 0.8 | 296.5 G | 1.9 | 2.51 | 9.3 |
| **700** | **29.9** | 29.2 | 0.5 | **277.7 G** | 2.6 | 2.68 | 12.5 |
| 800 | 32.2 | 31.4 | 0.6 | 289.3 G | 3.4 | 4.89 | 9.7 |
| 1,000 | 39.4 | 38.3 | 0.9 | 334.1 G | 5.3 | 5.39 | 8.6 |

K=700 is the lowest wall and also the lowest instruction count, so the choice
does not hang on wall noise. K=600's wall is 7.3 s above its user CPU
(contention), and its user CPU is still above K=700's. PR #830 chose K=600 on
the same corpus under load ~35 (walls 158–292 s), on a grid without 500 or
700.

L=4,096 tune (`n61-ks-growing-tune-4096-v1`, grid from PR #830):

| K | Wall s | User s | Instructions | Est. / measured RSS GiB |
|---:|---:|---:|---:|---:|
| 600 | 77.3 | 71.9 | 666.8 G | 1.9 / 2.51 |
| 800 | 60.2 | 56.5 | 510.9 G | 3.4 / 4.89 |
| **1,000** | **54.7** | 52.1 | **475.5 G** | 5.3 / 5.39 |
| 1,200 | 58.7 | 57.0 | 520.0 G | 7.7 / 10.00 |

This is an interior optimum on both wall and instructions. The memory-skip
estimate (94 B per regular state, from PR #1134) under-predicts at K≥800. It
skipped nothing here, but it is optimistic for larger K.

## Paired blocks

Alternating order, the same eval targets in both arms. Each cell shows wall s,
with user s in parentheses.

L=1,024, K=700:

| Block | Order | IC | Rho | Wall | User | Instr. | Peak RSS IC / rho GiB |
|---:|---|---:|---:|---:|---:|---:|---|
| 0 | IC→rho | 29.6 (28.9) | 107.1 (101.9) | 0.276 | 0.283 | 0.132 | 2.71 / 0.36 |
| 1 | rho→IC | 34.8 (30.6) | 123.4 (105.9) | 0.282 | 0.289 | 0.132 | 2.71 / 0.36 |
| 2 | IC→rho | 33.3 (30.4) | 108.9 (103.1) | 0.306 | 0.294 | 0.132 | 2.71 / 0.36 |

L=4,096, K=1,000:

| Block | Order | IC | Rho | Wall | User | Instr. | Peak RSS IC / rho GiB |
|---:|---|---:|---:|---:|---:|---:|---|
| 0 | IC→rho | 53.7 (51.9) | 217.8 (204.6) | 0.247 | 0.254 | 0.117 | 5.42 / 0.75 |
| 1 | rho→IC | 54.8 (52.4) | 203.7 (201.8) | 0.269 | 0.259 | 0.117 | 5.42 / 0.75 |
| 2 | IC→rho | 53.3 (51.6) | 203.5 (201.5) | 0.262 | 0.256 | 0.117 | 5.42 / 0.75 |

The IC arm uses 7.2–7.5× rho's peak memory. Compact is faster on wall in all
six blocks against this comparator.

## Verification

In every block of both runs, both arms solved every target. IC
`recovered_matches_published` and `group_verified` hold everywhere, as does
rho's recovered == published. Both arms walk the same published points, and the
scalars match the eval corpus. Records with their timers removed hash equal
across blocks in each arm. The pure-Python replays
(`autolab_orbit_extract_20260924/independent_replay.py` for IC against the
dumped base, and `growing_n_n61_L65536_20260930_rho_replay.py` for rho) pass
**3,072 / 3,072** (IC) and **3,072 / 3,072** (rho) at L=1,024, and
**12,288 / 12,288** in each arm at L=4,096. These are same-host,
same-campaign replays, so independent-host validation is still outstanding.
Receipts are in `artifacts/replay_*.json`.

## Claim-check

- `end_to_end_dlp` (`artifacts/claim_draft.json` → `claim_check.json`):
  **PASS** in both runs. This checks schema completeness only.
- `vs_rho` (`claim_draft_vs_rho.json` → `claim_check_vs_rho.json`):
  **FAIL, closed by design**. There are 18 missing stage fields plus
  `target_count must equal 1`. No single-target value was filled in.

## Files

`autolab_runs/<run-id>/` holds `state.json`, `artifacts/`, `receipts/` and
`inputs/` (protocol snapshot, `ledger_pin.json`, both scalar corpora). Two
things are left behind, as the autolab README asks: `logs/` (raw per-target
records, the base dump, producer stderr; 4–60 MB per run) and the 930 KB
`inputs/boundary_targets.json` copy. The ledger copy is pinned by SHA-256 in
`ledger_pin.json` and recoverable from git at `09a92b266`. Because of this,
`artifacts/review_manifest.json` lists files that are not in the bundle.
Rerunning the `launch-panel` command regenerates everything from the public
corpus names and seeds.

## Isolated timings (gap 3)

`.github/workflows/compact-orbit-n61-isolated.yml` runs the same
`launch-panel` command on a hosted `ubuntu-24.04` runner under
`tools/isolated_bench.py reserve`, with every producer pinned to the reserved
CPU (`--cpu`). It follows the `compact-ir-cold-gap.yml` pattern. On pull
requests it validates only. The measurement job is `workflow_dispatch`, so it
can only run once the workflow is on the default branch. Hosted runners are
VMs, so expect the ±5–10% residual wall noise §10 describes. The A/A spread
is the per-arm spread across the alternating blocks (default 5).
