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
| `20261001T212647Z-c0875d45ed` | 65,536 | 1,480 | 0.332 (0.139–0.491) | 0.291 | 0.132 | PASS | FAIL (closed) |

All ratios are compact / rho on the same block pair. The L=65,536 run is
described in its own section below; its host was far more loaded.

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

## L=65,536 run `20261001T212647Z-c0875d45ed`

Same command with `--targets 65536`, at source `777e767d4`, on the same
Apple M4 Pro host. The status is **`PENDING_INDEPENDENT_VALIDATION`**. The
host was heavily loaded by other agents' jobs: the 1-minute load was 72–213
before the tune runs and 12–55 before the timed block processes. Swap was
8.6–11.3 GiB used (of 10–12 GiB), but every timed process reported **0 swaps**.
These walls are **not** AGENTS.md §10 evidence. Read the instruction ratio
first, then user CPU, then wall.

The launching agent lost its session during the replay phase, but the
`launch-panel` process kept running. It finished on its own and released its
lock (exit 0), so `--resume` was not needed.

### K tune

The tune used the disjoint corpus `n61-ks-growing-tune-65536-v1` (sha256
`aba453fc…`), which shares no scalar with the eval corpus
`n61-ks-growing-65536-v1` (`162ddce8…`). Its grid comes from `protocol.json`.
K is picked by lowest whole-process wall. The memory gate uses the
`compact_orbit_rss` model with 2 GiB headroom, and every candidate fit
(K=1,480: 11.05 + 2 GiB against 19.1 GiB available), so none was skipped.

| K | Wall s | User s | Sys s | Instructions | Est. / max RSS / peak footprint GiB | Avail. GiB | Load |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 653.9 | 444.0 | 19.2 | 3,340.9 G | 5.43 / 4.55 / 5.39 | 13.1 | 199 |
| 1,200 | 1,420.2 | 399.2 | 291.2 | 3,148.7 G | 10.03 / 8.71 / 10.00 | 15.8 | 124 |
| 1,400 | 891.3 | 331.2 | 194.7 | 2,550.9 G | 10.73 / 8.53 / 10.72 | 17.7 | 213 |
| **1,480** | **281.1** | 255.1 | 17.3 | **2,041.5 G** | 11.05 / 10.48 / 11.03 | 19.1 | 72 |

- K=1,480 has the lowest wall, user CPU and instruction count, so the choice
  does not hang on wall noise.
- The K=1,000–1,400 walls are badly contention-inflated. Their wall exceeds
  user + sys by 190–730 s, and K=1,200 and 1,400 spent 195–291 s in the kernel.
- K=1,480 is the top of the grid, not an interior optimum. The root table
  doubles above K=1,482 (19.1 GiB estimated at K=1,500). With the 2 GiB
  headroom, that exceeds the 19–21 GiB this host had free.
- The `compact_orbit_rss` estimate matches the measured peak footprint within
  0.1–0.8% at every K.

### Paired blocks

Alternating order, with the same 65,536 eval targets in both arms. Each cell
shows wall s, with user s in parentheses.

| Block | Order | IC | Rho (KS v2) | Wall | User | Instr. | Peak footprint IC / rho GiB |
|---:|---|---:|---:|---:|---:|---:|---|
| 0 | IC→rho | 277.6 (238.3) | 835.2 (819.6) | 0.332 | 0.291 | 0.132 | 11.08 / 2.97 |
| 1 | rho→IC | 503.2 (265.4) | 1,023.8 (843.7) | 0.491 | 0.315 | 0.149 | 11.08 / 2.68 |
| 2 | IC→rho | 289.0 (240.2) | 2,072.2 (1,555.4) | 0.139 | 0.154 | 0.132 | 11.08 / 2.30 |

- The wall ratio median is **0.332**, with a range of 0.139–0.491 over 3
  blocks. With three blocks no confidence interval is reported.
- The user-CPU ratio median is 0.291 (0.154–0.315). The instruction ratio
  median is **0.132** (0.132–0.149).
- Block 1's IC wall is 238 s above its user CPU, with 147 s in the kernel.
- Block 2's rho used 1,555 s of user CPU against about 830 s in blocks 0 and
  1, for the same instruction count (16.9 T). That is host contention, so
  block 2's wall and user ratios are contaminated. The instruction ratio is
  the stable reading.
- The IC arm uses 3.7–4.8× rho's peak memory.

### Verification and claim-check

- In every block, both arms solved all 65,536 targets.
- `ic_all_verified`, `rho_all_verified`, scalars == eval corpus and same
  target points hold in every block. Records with their timers removed hash
  equal across blocks in each arm.
- The pure-Python replays pass **196,608 / 196,608** (IC, against the dumped
  base) and **196,608 / 196,608** (rho), with 0 failures. They are
  same-host, same-campaign replays.
- `end_to_end_dlp` passes. Re-running `claim-check` on
  `artifacts/claim_draft.json` gives a result identical to
  `artifacts/claim_check.json`.
- `vs_rho` fails closed (exit 4), as for L=1,024 and 4,096: 18 missing stage
  fields plus `target_count must equal 1`.
- `verify` re-hashes all 54 manifest files with 0 mismatches. `preflight`
  is ok.
- No ledger row changes.

This run additionally keeps the small producer logs: `logs/*.stderr.txt`
(`/usr/bin/time -l`), `logs/*.summary.json` and `logs/replay_*.log`. The
per-target `.jsonl` records and the base dump (about 300 MB) are left out.

## Isolated timings (gap 3)

`.github/workflows/compact-orbit-n61-isolated.yml` runs the same
`launch-panel` command on a hosted `ubuntu-24.04` runner under
`tools/isolated_bench.py reserve`, with every producer pinned to the reserved
CPU (`--cpu`). It follows the `compact-ir-cold-gap.yml` pattern. On pull
requests it validates only. The measurement job is `workflow_dispatch`, so it
can only run once the workflow is on the default branch. Hosted runners are
VMs, so expect the ±5–10% residual wall noise §10 describes. The A/A spread
is the per-arm spread across the alternating blocks (default 5).

### Hosted run `gha-36912327205-1` (Actions run 36912327205, L=1,024)

[Actions run 36912327205](https://github.com/aburan28/crypto/actions/runs/36912327205)
dispatched the workflow at `64e9ee32a` with `--targets 1024 --blocks 5`. Host:
Intel Xeon Platinum 8370C @ 2.80 GHz, 4 logical CPUs, Linux 6.17 (Azure),
`isolated_bench.py reserve` on CPUs {2, 3}, producers pinned to CPU 3.
`isolation.jsonl` reports **0 contended samples** and exit status 0. The
`verify` re-hash of all 76 manifest files passes. The run is archived under
`autolab_runs/gha-36912327205-1/`, with the artifact's `isolation.jsonl`,
`verify.json`, `preflight_1.log` and `lscpu.txt` in `hosted/`. `logs/` and
`inputs/boundary_targets.json` are left out as above.

The tune ran on the disjoint 1,024-target corpus and picked K=700, as on macOS:

| K | Process wall s | Setup s | Rank s | Targets s | Peak RSS GiB |
|---:|---:|---:|---:|---:|---:|
| 400 | 71.1 | 5.9 | 19.6 | 45.5 | 1.23 |
| 500 | 52.8 | 9.3 | 14.4 | 29.2 | 1.35 |
| 600 | 45.3 | 13.4 | 12.0 | 19.9 | 2.50 |
| **700** | **43.1** | 18.2 | 10.4 | 14.4 | 2.68 |
| 800 | 45.0 | 24.1 | 9.2 | 11.6 | 4.89 |
| 1,000 | 52.5 | 37.5 | 7.4 | 7.5 | 5.39 |

The producer's `setup_total` excludes `rank_stage`, so process wall ≈ setup +
rank + targets.

The five alternating blocks, with the same eval targets in both arms, give
wall in seconds:

| Block | Order | IC | Rho (KS v2) | Wall ratio | User ratio |
|---:|---|---:|---:|---:|---:|
| 0 | IC→rho | 42.66 | 183.91 | 0.2320 | 0.2295 |
| 1 | rho→IC | 42.16 | 183.94 | 0.2292 | 0.2268 |
| 2 | IC→rho | 42.12 | 183.79 | 0.2292 | 0.2268 |
| 3 | rho→IC | 42.24 | 186.20 | 0.2269 | 0.2244 |
| 4 | IC→rho | 42.36 | 183.90 | 0.2303 | 0.2280 |

- The wall ratio median is **0.229** (range 0.227–0.232), and the user-CPU
  ratio median is 0.227.
- The A/A spread per arm is 1.3% for IC (42.12–42.66 s) and 1.3% for rho
  (183.79–186.20 s), inside the §10 hosted-VM noise band.
- Instruction counts are not collected on Linux.
- Verification is the same as on macOS: every block solves all 1,024 targets
  in both arms, and both pure-Python replays pass 5,120 / 5,120.
- `end_to_end_dlp` passes. `vs_rho` fails closed for the same reason as
  above.

The isolated hosted ratio (0.229) sits below the contended macOS ratio
(0.282). It is still a multi-target batch against KS v2, so it is a diagnostic
and not the one-target primary comparison. For that comparison, see
`../20261001-koblitz-n61-single-target/`.

## L=65,536 through the autolab (`20261001T212647Z-c0875d45ed`)

This reruns PR #1134's L=65,536 point as a proper autolab run:
`launch-panel --beat koblitz.compact_orbit.n61_panel --targets 65536` at
`777e767d4` (3 blocks).

- K comes from a tune on the disjoint corpus `n61-ks-growing-tune-65536-v1`.
  The eval corpus is `n61-ks-growing-65536-v1`.
- The memory skip uses the new `ic_rss_model`, and the grid is
  [1000, 1200, 1400, 1480].
- K=1,480 is the last K before the IC root table doubles (K ≈ 1,482). Above
  it, peak memory jumps to about 19 GiB, which this shared host cannot hold
  next to other agents' jobs.
- The comparator is still the frozen KS v2 batched rho. It is not
  `rho_reference_minimum`, and this is a multi-target diagnostic.

Host: the same M4 Pro, unpinned. The 1-minute load before and after each timed
process was 12–213, and system swap was 4.7–11.3 GiB used by other jobs. Every timed
process reported **0 swaps**.

| K | Wall s | User s | Sys s | Instructions | Est. / measured peak footprint GiB | Load at start |
|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 653.9 | 444.0 | 19.2 | 3,340.9 G | 5.43 / 5.39 | 199 |
| 1,200 | 1,420.2 | 399.2 | 291.2 | 3,148.7 G | 10.03 / 10.00 | 124 |
| 1,400 | 891.3 | 331.2 | 194.7 | 2,550.9 G | 10.73 / 10.72 | 213 |
| **1,480** | **281.1** | **255.1** | 17.3 | **2,041.5 G** | 11.05 / 11.03 | 72 |

The registered wall rule picks **K=1,480**, which is also the lowest user CPU
and the lowest instruction count.

- At K=1,200 and K=1,400 the large sys times are memory-pressure page work.
  The macOS instruction counter includes kernel instructions, so those two
  counts are inflated too. User CPU is the cleanest column, and it falls
  monotonically to the grid edge.
- The optimum therefore sits at the memory cap, not at an interior K.
- The RSS model predicts the measured peak footprint within 0.1–0.8% at all
  four K. macOS `max_rss` reads lower under memory compression.

Paired blocks, with the same eval targets in both arms (wall s, user s in
parentheses):

| Block | Order | IC | Rho (KS v2) | Wall | User | Instr. | Peak footprint IC / rho GiB |
|---:|---|---:|---:|---:|---:|---:|---|
| 0 | IC→rho | 277.6 (238.3) | 835.2 (819.6) | 0.332 | 0.291 | 0.132 | 11.08 / 2.97 |
| 1 | rho→IC | 503.2 (265.4) | 1,023.8 (843.7) | 0.492 | 0.315 | 0.149 | 11.08 / 2.68 |
| 2 | IC→rho | 289.0 (240.2) | 2,072.2 (1,555.4) | 0.139 | 0.154 | 0.132 | 11.08 / 2.30 |

- The wall ratio median is **0.332** (range 0.139–0.492). The user-CPU ratio
  median is 0.291, and the instruction ratio median is **0.132**.
- Contention dominates the walls. Block 1's IC wall is 238 s above its user
  CPU. Block 2's rho used 1,555 s of user CPU against 820–844 s in blocks
  0–1, which is consistent with efficiency-core scheduling under load ~200.
- The instruction ratio is the stable reading. Blocks 0 and 2 agree to
  0.0003. Block 1 is higher because of kernel page work.
- IC setup (index and rank) took 62–113 s per block, and the targets took
  209–385 s. The IC uses 3.7–4.8× rho's peak memory.

PR #1134 measured the same point outside the autolab at K=1,400, with a wall
ratio median of 0.352. This run's instruction ratio is in line with the
L=1,024 (0.132) and L=4,096 (0.117) runs above. The wall ratio is not
comparable, because of the load.

Verification: in every block, both arms solved all 65,536 targets.

- IC `recovered_matches_published` and `group_verified` hold everywhere.
- Rho's recovered scalar equals the published one.
- Both arms walk the same points, and the scalars match the eval corpus.
- Records with timers removed hash equal across blocks in each arm.

The pure-Python replays pass **196,608 / 196,608** for IC
(`independent_replay.py` against the dumped base) and **65,536 / 65,536** in
each rho block (3 × 65,536). They are same-host and same-campaign, and the IC
replay alone took about 4 h of CPU.

Claim-check:

- `end_to_end_dlp`: **PASS**.
- `vs_rho`: **FAIL** (closed by design; this is a batch run).
- `verify`: PASS, 54 manifest files with 0 mismatches (`artifacts/verify.json`).

The run is archived under `autolab_runs/20261001T212647Z-c0875d45ed/`. Left
out, as above: `logs/` (per-target records and the base dump)
and `inputs/boundary_targets.json`. The two 1 MB scalar corpora are kept.
