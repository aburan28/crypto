# n=61 compact orbit: hosted isolated reruns at larger sizes (2026-10-03)

Status **`PENDING_INDEPENDENT_VALIDATION`**. These are reruns of the two
n=61 autolab beats on hosted GitHub runners, a different host class from the
macOS campaign, under `tools/isolated_bench.py reserve` with every producer
pinned to the reserved CPU. Each run replays its own records on the same
runner. That is an independent-host rerun of the measurement, not an
independent-author replay, and it promotes no ledger row.

Earlier hosted runs (L=1,024 batch, W=64 single-target) are in
[`../20261001-koblitz-n61-compact-orbit-panel/`](../20261001-koblitz-n61-compact-orbit-panel/README.md#isolated-timings-gap-3)
and
[`../20261001-koblitz-n61-single-target/`](../20261001-koblitz-n61-single-target/README.md#isolated-hosted-run-gha-36933523033-1).

| Actions run | Workflow | Size | K | Headline (median, range or 95% CI) | Replays | `verify` | Claim-check |
|---|---|---|---:|---|---|---|---|
| [37156127161](https://github.com/aburan28/crypto/actions/runs/37156127161) | batch panel | L=4,096, 5 blocks | 700 (fixed) | wall ratio IC/rho **0.213** (0.2125–0.2132) | IC 20,480 / 20,480, rho 20,480 / 20,480 | PASS (64 files) | `end_to_end_dlp` PASS, `vs_rho` FAIL (closed) |
| [37159733133](https://github.com/aburan28/crypto/actions/runs/37159733133) | batch panel | L=16,384, 3 blocks | 700 (fixed) | wall ratio IC/rho **0.323** (0.3225–0.3230) | IC 49,152 / 49,152, rho 49,152 / 49,152 | PASS (46 files) | `end_to_end_dlp` PASS, `vs_rho` FAIL (closed) |
| [38302195019](https://github.com/aburan28/crypto/actions/runs/38302195019) | batch panel | L=65,536, 3 blocks | 1,480 (same-L tune) | wall ratio IC/rho **0.239** (0.2374–0.2426) | IC 196,608 / 196,608, rho 196,608 / 196,608 | PASS (54 files, full artifact) | `end_to_end_dlp` PASS, `vs_rho` FAIL (closed) |
| [38822588626](https://github.com/aburan28/crypto/actions/runs/38822588626) | batch panel | L=65,536, 3 blocks | 700 (fixed) | wall ratio IC/rho **0.805** (0.8019–0.8083) | IC 196,608 / 196,608, rho 196,608 / 196,608 | PASS (46 files, full artifact) | `end_to_end_dlp` PASS, `vs_rho` FAIL (closed) |
| [37156392134](https://github.com/aburan28/crypto/actions/runs/37156392134) | single target | W=256 | 500 (fixed) | online speedup rho/IC **14.67** (12.40–17.21); cold IC/rho 74.8 (70.0–80.0) | 256 / 256 workloads, both arms | PASS (3,616 files) | `vs_rho` PASS 256 / 256 |

The L=4,096 and W=256 runs were dispatched with `workflow_dispatch` on `main`
at `b5077d433`, L=16,384 at `87055a99f` (AMD EPYC 7763). The L=65,536
K=1,480 tune run was dispatched on `main` at `966a0dbc9` (AMD EPYC 9V74);
the L=65,536 fixed-K=700 apples-to-apples rerun on `main` at `700f7d311`
(AMD EPYC 9V45). All: 4 vCPUs, Linux 6.17 (Azure), `ubuntu-24.04`,
reservation on CPUs {2, 3}, producers pinned to CPU 3,
`RAYON_NUM_THREADS=1`.

## Choice of K

The batch workflow takes a fixed K or, when `k` is empty, tunes on the
disjoint same-L corpus (it has no `--tune-targets` input). The L=4,096 and
L=16,384 runs fix K=700, the pick of the 1,024-target tune for L=4,096
(`20261002T210608Z-7c5a54a1ed`) and of the K=700 retunes for L=16,384 /
L=65,536. The first L=65,536 hosted run left `k` empty and tuned on
`n61-ks-growing-tune-65536-v1`, picking **K=1,480** by lowest wall on the
grid {1000, 1200, 1400, 1480}. A later L=65,536 hosted run fixes **K=700**
to match the macOS 1,024-target retune (`20261003T233040Z-454badefbb`).
The single-target run fixes K=500, the pick of the previous hosted
single-target tune (`gha-36933523033-1`), to fit 256 workloads in the job
limit.

## Batch panel, L=4,096 (`gha-37156127161-1`)

`launch-panel --beat koblitz.compact_orbit.n61_panel --targets 4096
--blocks 5 --k 700 --cpu 3`. `isolation.jsonl`: **0 contended samples**,
exit status 0. The measure job took 59 min.

| Block | Order | IC wall s | Rho (KS v2) wall s | Wall ratio | User ratio |
|---:|---|---:|---:|---:|---:|
| 0 | IC→rho | 82.78 | 388.33 | 0.2132 | 0.2110 |
| 1 | rho→IC | 82.82 | 388.44 | 0.2132 | 0.2113 |
| 2 | IC→rho | 82.87 | 389.24 | 0.2129 | 0.2109 |
| 3 | rho→IC | 82.79 | 388.52 | 0.2131 | 0.2111 |
| 4 | IC→rho | 82.68 | 389.17 | 0.2125 | 0.2105 |

- The A/A spread per arm is 0.2% for IC (82.68–82.87 s) and 0.2% for rho
  (388.33–389.24 s).
- IC phases per block: index 16.4–16.5 s, rank 9.8 s, targets 56.4–56.6 s.
  Peak RSS IC 2.71 GiB, rho 0.58 GiB.
- Instruction counts are not collected on Linux.
- In every block both arms solved all 4,096 targets and walk the same
  points; the pure-Python replays pass 20,480 / 20,480 (IC) and
  5 × 4,096 (rho).
- The contended macOS run at the same L and K
  (`20261002T210608Z-7c5a54a1ed`) read wall 0.319, user 0.311 and
  instructions 0.146. The isolated hosted wall ratio is lower, as it was at
  L=1,024 (0.229 hosted against 0.282 on macOS).

## Batch panel, L=16,384 (`gha-37159733133-1`)

`launch-panel --beat koblitz.compact_orbit.n61_panel --targets 16384
--blocks 3 --k 700 --cpu 3`. Three blocks rather than five keep the job,
with its 49,152-record IC replay, inside the 300-minute limit.
`isolation.jsonl`: **0 contended samples out of 1,060**, exit status 0. The
measure job took 92 min.

| Block | Order | IC wall s | Rho (KS v2) wall s | Wall ratio | User ratio |
|---:|---|---:|---:|---:|---:|
| 0 | IC→rho | 249.97 | 774.97 | 0.3226 | 0.3193 |
| 1 | rho→IC | 249.35 | 773.17 | 0.3225 | 0.3193 |
| 2 | IC→rho | 249.58 | 772.71 | 0.3230 | 0.3197 |

- The A/A spread per arm is 0.2% for IC (249.35–249.97 s) and 0.3% for rho
  (772.71–774.97 s).
- IC phases per block: index 16.4–16.5 s, rank 9.7–9.8 s, targets
  223.2–223.6 s. Peak RSS IC 2.71 GiB, rho 1.15 GiB.
- In every block both arms solved all 16,384 targets and walk the same
  points; the replays pass 49,152 / 49,152 (IC) and 3 × 16,384 (rho).
- The timer-stripped records hash to the same values as the macOS run at
  the same L and K (`20261003T211804Z-e9e3827d82`, in the companion PR):
  IC `892a95e5…`, rho `e533b6e3…`. Both hosts produced identical outputs;
  only the timers differ.
- The macOS run read wall 0.481, user 0.496 and instructions 0.225 under
  heavy load. The isolated hosted wall ratio, 0.323, is lower, as at
  L=1,024 and L=4,096.
- At L=16,384 the ratio is higher than at L=4,096 (0.213) because K=700,
  chosen on a 1,024-target corpus, makes the per-target descent the
  dominant IC cost at this L (223 s of 250 s).

It is still a multi-target batch against the frozen KS v2 batched rho, which
is superseded by strong rho R3 and not on
`measurement_schema.vs_rho.rho_reference_minimum`. Against strong rho R3 the
n=61 L=1,024 instruction ratio is IC/rho = 2.297. This run does not change
that verdict.

## Single-target panel, W=256 (`gha-37156392134-1`)

`launch-single --beat koblitz.compact_orbit.n61_single_target --workloads
256 --tune-workloads 4 --k 500 --cpu 3`. Every workload is one public
hash-to-curve target with no known scalar, solved by a fresh process in each
arm against strong rho rung 3. `isolation.jsonl`: **0 contended samples out
of 1,163**, exit status 0. The measure job took 100 min.

| Statistic | IC online ms | Rho online ms | Online speedup (rho / IC) | Cold ratio (IC / rho) |
|---|---:|---:|---:|---:|
| Median (95% bootstrap CI) | 19.97 (16.65–23.47) | 294.8 (273.6–317.8) | **14.67 (12.40–17.21)** | **74.8 (70.0–80.0)** |
| Min / max | 0.19 / 170.6 | 17.5 / 882.4 | 0.46 / 1,601 | 25.6 / 1,145 |

- IC is faster online on **255 of 256** workloads.
- Read the other way round, the median of the 256 per-workload IC/rho
  online ratios (`ic_online_ms / rho_online_ms` from the claim rows in
  `artifacts/single_target_summary.json`) is **0.0682** (min 0.0006,
  max 2.156). This is the figure the scoreboard's online series plots.
- Median IC online phases: PDP 19.94 ms, recovery check 0.025 ms, relation
  check 0.008 ms, descent 0.0004 ms.
- Cold, the IC process takes a median 22.06 s (setup 22.02 s) against rho's
  0.300 s, so IC loses the cold one-target comparison by about 75×.
- Peak RSS: IC median 1.35 GiB, rho median 33 MiB.
- Producer identity holds at K=500: the online producers' untimed records
  equal the frozen originals' on the known-answer target.
- Every workload passes verification and both replay certificates.
  Claim-check `vs_rho` passes **256 / 256**; status stays
  `PENDING_INDEPENDENT_VALIDATION`.

The W=64 hosted run (`gha-36933523033-1`, K=500 from its own tune) read an
online speedup of 12.68 (8.79–18.14) and a cold ratio of 98 (84–115) on an
EPYC 7763 runner. This W=256 run agrees on direction and size: IC wins the
online interval by about an order of magnitude and loses cold by about two.
The setup is shorter here (22.0 s against 33.5 s), on the same CPU model, so
the cold ratio is lower; both are single runner VMs.

## Files

Each `autolab_runs/<run-id>/` holds the run directory from the Actions
artifact: `state.json`, `artifacts/`, `receipts/`, `inputs/` and the small
producer logs in `logs/`. `hosted/` holds the artifact's `isolation.jsonl`,
`verify.json`, `preflight_1.log` and `lscpu.txt`. Left out, as the autolab
README asks: `inputs/boundary_targets.json` (pinned by SHA-256 in
`ledger_pin.json`), the base dumps `logs/base_n61_K*.jsonl`, and for the
batch runs the per-target `logs/*.jsonl` records. Their hashes are in each
run's `artifacts/review_manifest.json`, and the artifacts hold them until
they expire.

| Run | Artifact | Artifact id | Digest | Expires |
|---|---|---:|---|---|
| `gha-37156127161-1` | `compact-orbit-n61-isolated-L4096` | 11287237024 | `sha256:83dd37be3ccbbd60…` | 2027-01-01 |
| `gha-37159733133-1` | `compact-orbit-n61-isolated-L16384` | 11289103797 | `sha256:a42d32a3c2cf8505…` | 2027-01-01 |
| `gha-37156392134-1` | `compact-orbit-n61-single-target-W256` | 11287912438 | `sha256:8cd7bc1b6a2b137f…` | 2027-01-01 |
| `gha-38302195019-1` | `compact-orbit-n61-isolated-L65536` | 11780627009 | `sha256:fd4b739102fbeb94…` | 2027-01-10 |
| `gha-38822588626-1` | `compact-orbit-n61-isolated-L65536` | 12020154541 | `sha256:7bb37fc5ece79f6c…` | 2027-01-13 |

## Batch panel, L=65,536 (`gha-38302195019-1`)

Timeout bump from PR #1736 (`timeout-minutes` 480). Earlier same-day
dispatches with fixed `k=700` ([38288893010](https://github.com/aburan28/crypto/actions/runs/38288893010)
and siblings) sat in a saturated queue and were superseded; this run
completed measure on 2026-10-13.

`launch-panel --beat koblitz.compact_orbit.n61_panel --targets 65536
--blocks 3 --cpu 3` (empty `k` → tune on disjoint same-L corpus). Grid
walls on `n61-ks-growing-tune-65536-v1`: K=1000 487.4 s, K=1200 379.2 s,
K=1400 330.6 s, **K=1480 320.2 s** (chosen). `isolation.jsonl`: **0
contended samples out of 2,736**, exit status 0. Measure job ~3 h 52 min
after validate. Host: AMD EPYC 9V74, Linux 6.17.0-1022-azure.

| Block | Order | IC wall s | Rho (KS v2) wall s | Wall ratio | User ratio |
|---:|---|---:|---:|---:|---:|
| 0 | IC→rho | 319.26 | 1336.73 | 0.2388 | 0.2316 |
| 1 | rho→IC | 323.91 | 1335.17 | 0.2426 | 0.2354 |
| 2 | IC→rho | 316.88 | 1334.70 | 0.2374 | 0.2301 |

- Median wall ratio **0.2388** (0.2374–0.2426); user median **0.2316**.
- IC phases (block 0): index ~95.4 s, rank ~6.0 s, targets ~217.6 s.
  Peak RSS IC ~11.0 GiB, rho ~2.3 GiB.
- Both arms solve all 65,536 targets every block and walk the same points.
  Pure-Python replays: IC 196,608 / 196,608; rho 3 × 65,536.
- Rho untimed transcript hash `b0294263…` matches the macOS K=700 L=65,536
  panel (`20261003T233040Z-454badefbb`); IC hash differs because K=1,480
  here vs K=700 there.
- Full-artifact `verify` PASS (54 files). Trimmed archive omits per-target
  `logs/*.jsonl`, `base_n61_K*.jsonl`, and `inputs/boundary_targets.json`
  (hashes retained in `review_manifest.json`).
- Claim-check: `end_to_end_dlp` PASS; `vs_rho` FAIL-closed (batch
  `target_count ≠ 1`). Status
  **`PENDING_INDEPENDENT_VALIDATION`**. No ledger promotion.

This is still multi-target batch against frozen KS v2 batched rho, not the
strong-rho / single-target primary ladder.


## Batch panel, L=65,536 K=700 fixed (`gha-38822588626-1`)

Apples-to-apples independent-host rerun of the macOS K=700 L=65,536 panel
(`20261003T233040Z-454badefbb`). Prior hosted L=65,536 archive
(`gha-38302195019-1`) used same-L tune K=1,480.

`launch-panel --beat koblitz.compact_orbit.n61_panel --targets 65536
--blocks 3 --k 700 --cpu 3`. `isolation.jsonl`: **1 contended sample**,
exit status 0. Measure job ~2 h 55 min after validate (panel step
23:24–02:18 UTC). Host: AMD EPYC 9V45, Linux 6.17.0-1022-azure,
`main` @ `700f7d311`.

| Block | Order | IC wall s | Rho (KS v2) wall s | Wall ratio | User ratio |
|---:|---|---:|---:|---:|---:|
| 0 | IC→rho | 805.29 | 1000.00 | 0.8053 | 0.7977 |
| 1 | rho→IC | 798.13 | 995.30 | 0.8019 | 0.7944 |
| 2 | IC→rho | 794.19 | 982.54 | 0.8083 | 0.8008 |

- Median wall ratio **0.8053** (0.8019–0.8083); user median **0.7977**.
- Peak RSS IC ~2.72 GiB, rho ~2.30 GiB (K=700; much smaller than K=1,480's
  ~11 GiB IC).
- Both arms solve all 65,536 targets every block and walk the same points.
  Pure-Python replays: IC 196,608 / 196,608; rho 196,608 / 196,608.
- Rho untimed transcript hash `b0294263…` matches the macOS K=700 panel
  and the prior hosted K=1,480 run; IC hash differs with K.
- Full-artifact `verify` PASS (46 files). Trimmed archive omits per-target
  `logs/*.jsonl`, `base_n61_K700.jsonl`, and `inputs/boundary_targets.json`
  (hashes retained in `review_manifest.json`).
- Claim-check: `end_to_end_dlp` PASS; `vs_rho` FAIL-closed (batch
  `target_count ≠ 1`). Status
  **`PENDING_INDEPENDENT_VALIDATION`**. No ledger promotion.
- Against the same KS v2 comparator, fixed K=700 is slower than the
  same-L-tuned K=1,480 hosted run (wall median 0.805 vs 0.239). The
  macOS K=700 wall median was 1.033 under load (not §10 evidence).

This is still multi-target batch against frozen KS v2 batched rho, not the
strong-rho / single-target primary ladder.
