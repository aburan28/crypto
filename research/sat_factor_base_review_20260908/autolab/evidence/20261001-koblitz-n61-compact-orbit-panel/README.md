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
| `20261002T150549Z-89bc0af28f` | 1,024 | 700 | 0.285 (0.258–0.295) | 0.277 | 0.132 | PASS | FAIL (closed) |
| `20261002T151923Z-95d6146694` | 4,096 | 1,000 | 0.255 (0.240–0.259) | 0.253 | 0.117 | PASS | FAIL (closed) |
| `20261002T154358Z-c3faa1f97b` | 16,384 | 1,400 | 0.261 (0.255–0.295) | 0.245 | 0.116 | PASS | FAIL (closed) |
| `20261002T210608Z-7c5a54a1ed` | 4,096 | 700 (1,024-target tune) | 0.319 (0.299–0.320) | 0.311 | 0.146 | PASS | FAIL (closed) |
| `20261002T213417Z-6e1f87218b` | 16,384 | 1,000 (1,024-target tune) | 0.414 (0.308–0.443) | 0.321 | 0.139 | PASS | FAIL (closed) |
| `20261002T230454Z-8be7370aae` | 65,536 | 1,000 (1,024-target tune) | 0.604 (0.579–0.722) | 0.491 | 0.224 | PASS | FAIL (closed) |

All ratios are compact / rho on the same block pair. The `20261002T15…`
rows are reruns with producer logs archived; see
[Reruns with producer logs](#reruns-with-producer-logs-and-l16384-2026-10-02).
The rows marked "1,024-target tune" follow the goal's rule that K is tuned on
a disjoint 1,024-target corpus; see
[K tuned on the disjoint 1,024-target corpus](#k-tuned-on-the-disjoint-1024-target-corpus-2026-10-02).

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

The left-out `logs/` of `20261001T155442Z-9990f40315` and
`20261001T160852Z-1f8383bd42` no longer exist. A search of every worktree and
temporary directory on the host on 2026-10-02 found only checked-out copies
of this bundle, none with `logs/`. Both panels were rerun through the
autolab instead (below), and those runs keep their small producer logs. The
two original runs stay as they are.

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

Added 2026-10-02: the small producer logs from the Actions artifact
(`compact-orbit-n61-isolated-L1024`, artifact id 11190340275, digest
`sha256:12a40f31…`, expires 2026-12-30) are now in `logs/`: 33 files, the
`*.stderr.txt`, `*.summary.json` and `replay_*.log`. All 76 entries of the
run's committed `artifacts/review_manifest.json`, including the 44 `logs/`
entries, match the downloaded artifact by SHA-256. The per-target
`*.jsonl` records (13 MB) and `base_n61_K700.jsonl` (5.3 MB) are still left
out. Their hashes are in `review_manifest.json`, and the artifact holds them
until it expires.

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
| 1 | rho→IC | 503.2 (265.4) | 1,023.8 (843.7) | 0.491 | 0.315 | 0.149 | 11.08 / 2.68 |
| 2 | IC→rho | 289.0 (240.2) | 2,072.2 (1,555.4) | 0.139 | 0.154 | 0.132 | 11.08 / 2.30 |

- The wall ratio median is **0.332** (range 0.139–0.491). The user-CPU ratio
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
out, as above: the per-target `logs/*.jsonl` records, the base dump
and `inputs/boundary_targets.json`. The two 1 MB scalar corpora are kept, as
are the small producer logs: `logs/*.stderr.txt` (`/usr/bin/time -l`),
`logs/*.summary.json` and `logs/replay_*.log`.

## Reruns with producer logs, and L=16,384 (2026-10-02)

Three new runs of
`launch-panel --beat koblitz.compact_orbit.n61_panel --targets L` (3 blocks
each), run one after another at `60fb2ae47` (`main`; `preflight` ok, 8/8
checks, 0 dirty paths). L=1,024 and L=4,096 replace the two runs whose logs
are gone. L=16,384 is the first autolab run at that L. PR #830 measured it
outside the autolab (`autolab_orbit_extract_20260924/growing_n_n61_L16384/`,
K=1,400, wall ratio median 0.205).

Host: the same M4 Pro, unpinned, rustc 1.93.1. Other agents' jobs were running
throughout. The 1-minute load before and after each timed process was
7.8–19.2 (L=1,024), 6.9–10.4 (L=4,096) and 6.9–26.2 (L=16,384). System swap
was 11.3–16.4 GiB used by other jobs, and every timed process reported
**0 swaps**. These walls are **not** AGENTS.md §10 evidence. The user-CPU
and instruction ratios are the less contaminated readings.

### K tunes

Each tune ran on the beat's disjoint corpus `n61-ks-growing-tune-{L}-v1`.
That corpus has L targets, so only the L=1,024 tune uses a 1,024-target
corpus (see [Tune corpora](#tune-corpora)). The memory gate
(`ic_rss_model` + 2 GiB headroom against available memory at launch)
skipped no candidate.

| L | K | Wall s | User s | Sys s | Instructions | Est. / measured peak GiB | Available GiB | Load |
|---:|---:|---:|---:|---:|---:|---|---:|---:|
| 1,024 | 400 | 55.3 | 52.2 | 0.8 | 471.0 G | 1.28 / 1.23 | 10.7 | 19.2 |
| 1,024 | 500 | 36.0 | 35.5 | 0.4 | 342.5 G | 1.40 / 1.35 | 10.9 | 16.9 |
| 1,024 | 600 | 31.4 | 30.6 | 0.6 | 296.5 G | 2.55 / 2.51 | 11.6 | 12.4 |
| 1,024 | **700** | **28.9** | 28.3 | 0.5 | **277.7 G** | 2.73 / 2.69 | 11.9 | 12.0 |
| 1,024 | 800 | 33.1 | 31.0 | 1.3 | 291.0 G | 4.94 / 5.08 | 11.9 | 10.8 |
| 1,024 | 1,000 | 45.1 | 39.5 | 2.6 | 338.5 G | 5.43 / 5.39 | 13.8 | 10.1 |
| 4,096 | 600 | 72.5 | 70.4 | 0.9 | 667.0 G | 2.55 / 2.51 | 12.2 | 7.5 |
| 4,096 | 800 | 57.0 | 54.6 | 1.7 | 513.7 G | 4.94 / 4.89 | 12.3 | 6.9 |
| 4,096 | **1,000** | **51.6** | 50.1 | 1.2 | **476.0 G** | 5.43 / 5.39 | 15.7 | 7.3 |
| 4,096 | 1,200 | 65.8 | 56.2 | 7.6 | 540.4 G | 10.03 / 9.99 | 15.8 | 7.4 |
| 16,384 | 800 | 161.7 | 157.1 | 2.8 | 1,384.9 G | 4.94 / 4.89 | 12.5 | 9.8 |
| 16,384 | 1,000 | 116.8 | 114.8 | 1.8 | 1,049.7 G | 5.43 / 5.39 | 14.3 | 10.5 |
| 16,384 | 1,200 | 123.1 | 105.4 | 15.1 | 970.8 G | 10.03 / 10.00 | 13.9 | 9.0 |
| 16,384 | **1,400** | **106.1** | 100.5 | 4.3 | **919.1 G** | 10.73 / 10.72 | 16.7 | 10.4 |

- The wall rule picks K=700, 1,000 and 1,400, the same K as the earlier
  runs (and PR #830 at L=16,384). Each pick is also the lowest instruction
  count.
- At L=16,384 the K=1,200 wall sits above K=1,000 because of 15 s of sys
  (memory-pressure page work). Its user CPU and instructions are lower, and
  K=1,400 is lowest on all three, so the pick does not hang on that noise.
- At L=16,384, K=1,400 is the top of the registered grid. The optimum may
  lie above it. K=1,480 is the last K below the root-table doubling.

### Paired blocks

Alternating order, the same eval targets in both arms. Wall s, with user s
in parentheses.

| L | Block | Order | IC | Rho (KS v2) | Wall | User | Instr. | Peak footprint IC / rho GiB |
|---:|---:|---|---:|---:|---:|---:|---:|---|
| 1,024 | 0 | IC→rho | 31.5 (30.0) | 106.7 (104.7) | 0.295 | 0.287 | 0.132 | 2.71 / 0.29 |
| 1,024 | 1 | rho→IC | 28.6 (28.0) | 110.7 (107.1) | 0.258 | 0.261 | 0.132 | 2.71 / 0.29 |
| 1,024 | 2 | IC→rho | 29.2 (28.2) | 102.5 (101.8) | 0.285 | 0.277 | 0.132 | 2.71 / 0.34 |
| 4,096 | 0 | IC→rho | 52.4 (50.8) | 217.9 (204.8) | 0.240 | 0.248 | 0.116 | 5.61 / 0.58 |
| 4,096 | 1 | rho→IC | 55.4 (51.5) | 217.4 (203.6) | 0.255 | 0.253 | 0.118 | 5.42 / 0.67 |
| 4,096 | 2 | IC→rho | 53.4 (51.6) | 205.8 (203.6) | 0.259 | 0.254 | 0.117 | 5.42 / 0.67 |
| 16,384 | 0 | IC→rho | 106.8 (99.3) | 419.1 (410.4) | 0.255 | 0.242 | 0.115 | 10.76 / 1.34 |
| 16,384 | 1 | rho→IC | 121.4 (101.7) | 411.3 (406.3) | 0.295 | 0.250 | 0.118 | 10.76 / 1.34 |
| 16,384 | 2 | IC→rho | 118.3 (103.4) | 453.9 (421.7) | 0.261 | 0.245 | 0.116 | 10.76 / 1.15 |

- The L=1,024 and L=4,096 reruns reproduce the earlier runs' instruction
  ratios (0.132 and 0.117 medians) and stay within their wall ranges.
- At L=16,384 the instruction ratio median is 0.116, in line with L=4,096.
  The wall ratio median (0.261) is above PR #830's 0.205, which ran under a
  different load and is not comparable.
- The IC arm's peak footprint is 8.0–9.7× rho's.
- These are multi-target batch diagnostics against the superseded KS v2
  comparator. Both arms scale as √(L·r/n) at optimal K, so the lead is a
  constant factor, not an asymptotic one.

### Verification and claim-check

In every block of all three runs, both arms solved all L targets. IC
`recovered_matches_published` and `group_verified` hold everywhere, rho's
recovered scalar equals the published one, both arms walk the same points,
and the scalars match the eval corpus. The pure-Python replays pass:

| Run | IC replay | Rho replay | `verify` | `end_to_end_dlp` | `vs_rho` | Status |
|---|---:|---:|---|---|---|---|
| `20261002T150549Z-89bc0af28f` | 3,072 / 3,072 | 3,072 / 3,072 | PASS (58 files) | PASS | FAIL (closed) | `PENDING_INDEPENDENT_VALIDATION` |
| `20261002T151923Z-95d6146694` | 12,288 / 12,288 | 12,288 / 12,288 | PASS (54 files) | PASS | FAIL (closed) | `PENDING_INDEPENDENT_VALIDATION` |
| `20261002T154358Z-c3faa1f97b` | 49,152 / 49,152 | 49,152 / 49,152 | PASS (54 files) | PASS | FAIL (closed) | `PENDING_INDEPENDENT_VALIDATION` |

The replays are same-host and same-campaign, so independent-host validation
is still outstanding. `end_to_end_dlp` was also re-run standalone with
`claim-check --stage end_to_end_dlp` on each `claim_draft.json` (PASS). It
checks schema completeness only. `vs_rho` fails closed because these are
batch runs (`target_count must equal 1`). No ledger row is promoted.

The runs are archived under `autolab_runs/<run-id>/` in the same layout as
`20261001T212647Z-c0875d45ed`: `state.json`, `artifacts/` (with
`verify.json`), `receipts/`, `inputs/` (without `boundary_targets.json`) and
the small producer logs in `logs/`. Left out: the per-target `logs/*.jsonl`
records and the base dump. Their hashes are in each run's
`review_manifest.json`.

### Tune corpora

Every panel's tune corpus is disjoint from its eval corpus. That is recorded
in `panel_summary.json` and re-checked from the archived scalar files
(0 shared scalars in each). The tune corpus is the beat's
`n61-ks-growing-tune-{L}-v1`, which has **L** targets. Only the L=1,024
panels tune on a 1,024-target corpus.

| Run | L | Tune corpus | Tune targets | Disjoint from eval |
|---|---:|---|---:|---|
| `20261001T155442Z-9990f40315` | 1,024 | `n61-ks-growing-tune-1024-v1` | 1,024 | yes |
| `20261001T160852Z-1f8383bd42` | 4,096 | `n61-ks-growing-tune-4096-v1` | 4,096 | yes |
| `20261001T212647Z-c0875d45ed` | 65,536 | `n61-ks-growing-tune-65536-v1` | **65,536** | yes |
| `gha-36912327205-1` | 1,024 | `n61-ks-growing-tune-1024-v1` | 1,024 | yes |
| `20261002T150549Z-89bc0af28f` | 1,024 | `n61-ks-growing-tune-1024-v1` | 1,024 | yes |
| `20261002T151923Z-95d6146694` | 4,096 | `n61-ks-growing-tune-4096-v1` | 4,096 | yes |
| `20261002T154358Z-c3faa1f97b` | 16,384 | `n61-ks-growing-tune-16384-v1` | 16,384 | yes |

This matches `protocol.json` (`k_selection`: "disjoint tune corpus of the
same L") and PR #830's scripts. The autolab skill's wording, "K is chosen by
lowest wall on a disjoint 1024-target corpus", holds only at L=1,024. The
optimal K grows with L (700 → 1,000 → 1,400 → 1,480 here), so a 1,024-target
tune would pick an L=1,024 K for every panel. Choosing between the two rules
is left open. It would change `protocol.json` and is not done here.

Resolved 2026-10-02: the user chose the 1,024-target rule. See the next
section.

## K tuned on the disjoint 1,024-target corpus (2026-10-02)

The goal's rule is that K is chosen by lowest wall on a disjoint
**1,024-target** corpus. `launch-panel` now has `--tune-targets T`: the tune
corpus becomes `corpora.tune` at L=T (here `n61-ks-growing-tune-1024-v1`,
1,024 targets), whatever the panel's L. The default is unchanged (T = the
panel's L), so the earlier runs and their manifests are untouched. The option
records `tune_targets` in `state.json`, `k_tune.json` and `candidate.json`,
and the run fails if the tune and eval corpora share a scalar. The grids are
`protocol.json`'s new `k_candidates_by_tune_targets["1024"]`. Every K in them
fits under the `ic_rss_model` on this 48 GiB host. `--tune-only` stops after
the tune (status `TUNE_ONLY`); `--resume` continues the same run into the
panel. Each panel below was launched as a tune-only run, and because the
tune picked a different K from the earlier run, it was resumed into the full
3-block panel.

Expected and observed: a 1,024-target tune favours a smaller K than the
panel's L would, because the setup cost (index and rank, which grow with K)
is amortised over 1,024 targets rather than L. That costs wall at large L,
as the ratios below show.

All three runs used the same tune corpus file. Its SHA-256 is
`8ec0c0438ba46a1b…`, byte-identical to the L=1,024 runs' tune corpus.
Re-checked from the archived scalar files, it shares **0** scalars with
each eval corpus (`n61-ks-growing-{4096,16384,65536}-v1`).

### 1,024-target tunes

Wall s, user s, sys s, instructions and est. / measured peak footprint GiB.
The bold row is the lowest wall, which is the pick.

| Panel L | K | Wall s | User s | Sys s | Instructions | Est. / measured GiB |
|---:|---:|---:|---:|---:|---:|---|
| 4,096 | **700** | **30.1** | 29.2 | 0.5 | **293.5 G** | 2.73 / 2.68 |
| 4,096 | 800 | 32.1 | 31.3 | 0.6 | 305.7 G | 4.94 / 4.89 |
| 4,096 | 1,000 | 40.4 | 37.9 | 1.2 | 354.7 G | 5.43 / 5.39 |
| 4,096 | 1,200 | 57.2 | 49.0 | 4.2 | 452.2 G | 10.03 / 9.99 |
| 16,384 | **1,000** | **40.8** | 37.4 | 1.9 | **355.8 G** | 5.43 / 5.39 |
| 16,384 | 1,200 | 57.7 | 49.2 | 5.3 | 454.8 G | 10.03 / 10.00 |
| 16,384 | 1,400 | 73.1 | 63.5 | 6.1 | 580.0 G | 10.73 / 10.71 |
| 65,536 | **1,000** | **100.6** | 64.2 | 5.5 | **361.3 G** | 5.43 / 5.39 |
| 65,536 | 1,200 | 697.3 | 94.4 | 201.0 | 857.2 G | 10.03 / 10.18 |
| 65,536 | 1,400 | 179.6 | 99.1 | 20.4 | 607.7 G | 10.73 / 10.71 |
| 65,536 | 1,480 | 305.9 | 102.9 | 89.5 | 830.4 G | 11.05 / 11.03 |

- Each pick is also the lowest user CPU and the lowest instruction count, so
  it does not hang on wall noise.
- At L=4,096 the pick is K=700, against K=1,000 in `20261002T151923Z-95d6146694`.
- At L=16,384 the pick is K=1,000, against K=1,400 in `20261002T154358Z-c3faa1f97b`.
- At L=65,536 the pick is K=1,000, against K=1,480 in `20261001T212647Z-c0875d45ed`.
  This tune ran under heavy load (1-minute load up to 68). K=1,200 and
  K=1,480 spent 201 s and 90 s in sys (memory-pressure page work), which
  also inflates their kernel-inclusive instruction counts. User CPU, the
  cleanest column, still rises monotonically with K (64 → 94 → 99 → 103 s),
  so K=1,000 is the pick on every column. The launch waited until 13.5 GiB
  was available, so the RSS gate (13.4–19.7 GiB available at each K's start)
  skipped no candidate.
- The L=16,384 and L=65,536 grids start at 1,000 as specified, so 700 was
  not offered there. On the same corpus at L=4,096 (and in the L=1,024 runs)
  K=700 beat K=1,000 by 10 s, so a grid that included it would very likely
  pick K=700 at every L.
- The 1,024-target K=700 instruction count is 293.5 G here, against 277.7 G
  in the L=1,024 runs on the same corpus. The producer sources are
  byte-identical (`koblitz_orbit_dlp_fast.rs` `8dea682a…`,
  `koblitz_rho_batch_ks_v2_n61.rs` `019b2e2d…`), but the binaries differ,
  because library code under `src/` changed on `main` between `60fb2ae47`
  and these commits. Compare ratios within a run, not raw counts across
  runs.

### Paired blocks

Alternating order, the same eval targets in both arms. Wall s, with user s
in parentheses. The IC peak column is the producer's reported peak RSS.

| L | K | Block | Order | IC | Rho (KS v2) | Wall | User | Instr. | IC peak GiB |
|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| 4,096 | 700 | 0 | IC→rho | 63.7 (61.6) | 213.3 (206.1) | 0.299 | 0.299 | 0.146 | 2.68 |
| 4,096 | 700 | 1 | rho→IC | 67.1 (64.2) | 210.2 (206.1) | 0.319 | 0.311 | 0.146 | 2.69 |
| 4,096 | 700 | 2 | IC→rho | 67.9 (65.2) | 211.8 (205.5) | 0.320 | 0.317 | 0.147 | 2.68 |
| 16,384 | 1,000 | 0 | IC→rho | 205.4 (140.5) | 463.9 (419.2) | 0.443 | 0.335 | 0.138 | 4.02 |
| 16,384 | 1,000 | 1 | rho→IC | 167.4 (136.6) | 542.7 (425.9) | 0.308 | 0.321 | 0.139 | 4.03 |
| 16,384 | 1,000 | 2 | IC→rho | 195.0 (136.6) | 470.7 (427.3) | 0.414 | 0.320 | 0.141 | 4.08 |
| 65,536 | 1,000 | 0 | IC→rho | 831.0 (623.0) | 1,376.1 (1,168.3) | 0.604 | 0.533 | 0.220 | 4.16 |
| 65,536 | 1,000 | 1 | rho→IC | 841.9 (497.1) | 1,166.6 (1,021.4) | 0.722 | 0.487 | 0.231 | 3.95 |
| 65,536 | 1,000 | 2 | IC→rho | 612.6 (466.3) | 1,058.7 (950.5) | 0.579 | 0.491 | 0.224 | 4.19 |

Medians (range) against the earlier same-L tuned runs:

| L | Run | K (tune corpus) | Wall ratio | User-CPU ratio | Instructions ratio |
|---:|---|---|---|---:|---:|
| 4,096 | `20261002T151923Z-95d6146694` | 1,000 (4,096 targets) | 0.255 (0.240–0.259) | 0.253 | 0.117 |
| 4,096 | `20261002T210608Z-7c5a54a1ed` | **700 (1,024 targets)** | 0.319 (0.299–0.320) | 0.311 | 0.146 |
| 16,384 | `20261002T154358Z-c3faa1f97b` | 1,400 (16,384 targets) | 0.261 (0.255–0.295) | 0.245 | 0.116 |
| 16,384 | `20261002T213417Z-6e1f87218b` | **1,000 (1,024 targets)** | 0.414 (0.308–0.443) | 0.321 | 0.139 |
| 65,536 | `20261001T212647Z-c0875d45ed` | 1,480 (65,536 targets) | 0.332 (0.139–0.491) | 0.291 | 0.132 |
| 65,536 | `20261002T230454Z-8be7370aae` | **1,000 (1,024 targets)** | 0.604 (0.579–0.722) | 0.491 | 0.224 |

- The 1,024-target rule raises every ratio at L≥4,096, and the gap grows
  with L. The instruction ratio goes from 0.117 to 0.146 at L=4,096, from
  0.116 to 0.139 at L=16,384, and from 0.132 to 0.224 at L=65,536. The
  compact orbit still beats this comparator in every block.
- At L=65,536 and K=1,000, IC setup is short (38 s index, 11 s rank in
  block 0), but the per-target phase is long (782 s in block 0, median
  7.1 ms per target). The rule trades a smaller setup for a slower descent,
  which is the wrong trade at this L.
- The L=65,536 run was heavily contended. The 1-minute load at process
  boundaries was 11.8–68.2. IC wall exceeds IC user CPU by 146–345 s, and
  rho's user CPU (950–1,168 s per block) is above the earlier run's
  820–844 s, consistent with efficiency-core scheduling. User-CPU ratios
  are therefore contaminated too. The instruction ratio is the steadiest
  reading (0.220–0.231).
- The L=16,384 walls are contended. The 1-minute load at process boundaries
  was 11.3–31.0, reaching about 49 between them, and IC wall exceeds IC user
  CPU by up to 65 s. The user-CPU and instruction ratios are the readings to
  use.
- L=4,096 ran under load 9.7–13.5. Every timed process in all three runs
  reported **0 swaps**. These walls are **not** AGENTS.md §10 evidence.

### Verification, claim-check and archive

| Run | L | IC replay | Rho replay | `verify` | `end_to_end_dlp` | `vs_rho` | `preflight` | Status |
|---|---:|---:|---:|---|---|---|---|---|
| `20261002T210608Z-7c5a54a1ed` | 4,096 | 12,288 / 12,288 | 12,288 / 12,288 | PASS (56 files) | PASS | FAIL (closed) | ok | `PENDING_INDEPENDENT_VALIDATION` |
| `20261002T213417Z-6e1f87218b` | 16,384 | 49,152 / 49,152 | 49,152 / 49,152 | PASS (54 files) | PASS | FAIL (closed) | ok | `PENDING_INDEPENDENT_VALIDATION` |
| `20261002T230454Z-8be7370aae` | 65,536 | 196,608 / 196,608 | 196,608 / 196,608 | PASS (56 files) | PASS | FAIL (closed) | ok | `PENDING_INDEPENDENT_VALIDATION` |

In every block both arms solved all L targets, both arms walk the same
points, the scalars match the eval corpus, and timer-stripped records hash
equal across blocks in each arm. `end_to_end_dlp` passes in the run and again
standalone (`claim-check --stage end_to_end_dlp` on `claim_draft.json`); it
checks schema completeness only. `vs_rho` fails closed because these are
batch runs (`target_count must equal 1`). The replays are same-host and
same-campaign, so independent-host validation is still outstanding. No
ledger row is promoted.

Sources: L=4,096 ran at `943268e9f`, which added `--tune-targets`. L=16,384
and L=65,536 ran at `891b55068`, which fixed the manifest self-entry bug
below. Both commits are on this PR's branch.

Manifest repair for `20261002T210608Z-7c5a54a1ed`: its tune-only stage wrote
`review_manifest.json`, and the resumed run's final manifest then hashed that
stale file as an entry for itself before overwriting it. `verify` flagged
exactly that one entry and nothing else. The as-written manifest is kept as
`artifacts/review_manifest_as_written.json`, together with its failing receipt
`artifacts/verification_before_manifest_repair.json`. `review_manifest.json`
is the same file minus that self-entry, and `verify` then passes on all 56
files. No measurement file changed. The code now never lists the manifest in
itself and, on resume, keeps the tune-only manifest as
`artifacts/review_manifest_tune_only.json` (the later runs carry it).

The runs are archived under `autolab_runs/<run-id>/` in the same layout as
the other 2026-10-02 reruns, with the small producer logs (tune logs
included) and without the per-target `logs/*.jsonl`, the base dump and
`inputs/boundary_targets.json`. Their hashes are in each run's
`review_manifest.json`.

### Which runs satisfy the 1,024-target rule

| L | Run | Tune corpus targets | Satisfies the rule |
|---:|---|---:|---|
| 1,024 | `20261001T155442Z-9990f40315`, `20261002T150549Z-89bc0af28f`, `gha-36912327205-1` | 1,024 | yes |
| 4,096 | `20261002T210608Z-7c5a54a1ed` | 1,024 | **yes** |
| 4,096 | `20261001T160852Z-1f8383bd42`, `20261002T151923Z-95d6146694` | 4,096 | no (same-L tune; kept as history) |
| 16,384 | `20261002T213417Z-6e1f87218b` | 1,024 | **yes** |
| 16,384 | `20261002T154358Z-c3faa1f97b` | 16,384 | no (same-L tune; kept as history) |
| 65,536 | `20261002T230454Z-8be7370aae` | 1,024 | **yes** |
| 65,536 | `20261001T212647Z-c0875d45ed` | 65,536 | no (same-L tune; kept as history) |

Every L from 1,024 to 65,536 now has at least one run whose K was chosen
on a disjoint 1,024-target corpus. The same-L-tuned runs stay in the bundle
as history, and they show what the rule costs: at L≥4,096 a same-L tune
picks a larger K and gives lower ratios (0.116–0.132 instructions against
0.139–0.224 here).
