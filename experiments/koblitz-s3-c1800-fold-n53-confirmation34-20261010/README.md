# N53 C1800 plus PMULL-fold fresh-target diagnostic

The fresh N53 comparison exposed a large preparation-time stability problem on the local Mac: verified C1800 runs took 34.1–160.2 seconds of reusable setup, with a 45% coefficient of variation, and T213 expanded reached the frozen 180-second process cap before emitting a report. The analyzer correctly leaves the 34-target online ratio **unknown**. The first attempt preserved 27 verified rows and one truncated-report failure; the explicit recovery preserved those rows, verified the failed point with a new `R2` run ID, advanced to 34 verified rows, and then stopped at the T213 timeout. Checked Sage independently replayed all 34 verified recovery scalars across 17 completed target pairs. Both partial ledgers and their failures are retained.

This study uses [34 previously unrun N53 public points](../koblitz-s3-c1800-confirmation34-20261010/README.md), T213–T246, with unique workload IDs. The same Mach-O executable and source are used for C1200/table (`GF2_ARM_PMULL_FOLD=0`) and C1800/fold (`=1`), each in a fresh process. [FREEZE.json](FREEZE.json) binds the SHA-256 identities, one-process-per-target online timing boundary, seeded target shuffle, alternating arm order, 14 preparation workers, one target worker, 180-second cap, 32 GiB memory envelope, and the frozen 2× gate. The candidate IDs are `IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0h5172f19b88c7` and `IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h5fca79d17564`.

## Recorded execution

| Attempt | Verified rows | Stopping row | Status |
| --- | ---: | --- | --- |
| Original | 27 / 68 | T220 C1200 `R1` | Exit 101; report truncated at exactly 122,880 bytes; stdout and stderr empty |
| [Explicit recovery](recovery/FREEZE.json) | 34 / 68 in the combined ledger | T213 C1800 `R1` | Process killed at 181.2 seconds by the frozen 180-second timeout; no report |

The original failure happened while `/private/tmp` had low free space. That is consistent with an output write failure, but the empty stderr does not prove its cause. The recovery moved raw output to the workspace volume, preserved the original failed row in its archive, assigned T220 C1200 retry `R2`, and kept all other completed and future run IDs. T220 `R2` verified. The T213 timeout is a separate execution failure and remains in the recovery archive. The analyzer returned `complete: false`, with `online_ratio: null`, for both attempts. Its missing-ledger check also returned zero verified rows and no ratio.

| Verified-run preparation diagnostic | C1200/table, 17 rows | C1800/fold, 17 rows |
| --- | ---: | ---: |
| Reusable setup, seconds: min / median / max | 15.96 / 28.02 / 39.15 | 34.13 / 81.71 / 160.16 |
| Setup coefficient of variation | 0.273 | 0.450 |
| Involuntary context switches: min / median / max | 116,117 / 169,559 / 227,462 | 284,267 / 489,650 / 764,569 |

The extreme setup spread and context switching make this host unsuitable for promoting an N53 CPU wall-time ratio. Since T213 produced no report, its exact stopping phase is unknown; the completed runs show that preparation alone can approach this panel's process cap. The next N53 experiment should use a checked host-level isolation receipt, sufficient output space, and a process cap chosen from a fresh preparation pilot before freezing its target panel. Modal guest CPU/memory diagnostics can establish exploratory fidelity, but a CPU reservation alone does not satisfy the repository's host-exclusive isolation gate.

## Raw receipts

- Original [partial analysis](local/partial_analysis.json), [raw ledger and process files](local/partial_raw.tar.zst), and [checked Sage runtime](local/sage_runtime_info.json). Raw archive SHA-256: `a9e587f329bce6c1fc611fcaf6df444f95020f98b4db09a3dc5eae259bd68166`.
- Recovery [partial analysis](recovery/local/partial_analysis.json), [independent Sage replay of 34 verified rows](recovery/local/partial_sage_replay.json), [raw ledger and process files](recovery/local/partial_raw.tar.zst), and [checked Sage runtime](recovery/local/sage_runtime_info.json). Raw archive SHA-256: `4003a262cf260b2d6f943ce9abff2e4e610d176422bc52d30a5d0eb0c318df9d`.
- The archived original and recovery ledgers hash to `0826850be92aa74eed9a0519458dac07ae077cd47f0bf9b3fc25633d30cc1a76` and `ae3a53ed6b7592a246d81d5755bc9e4c238b2740a17c0a3d9023d16d9da232c5`, matching the respective analysis receipts.
