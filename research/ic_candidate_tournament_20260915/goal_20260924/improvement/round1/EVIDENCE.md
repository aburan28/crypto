# Round one: provenance and interpretation

The registered tournament retained the qualified `pairinv` incumbent. The selected `stop6` challenger did not meet the complete cold instruction/time thresholds, and its single-target online time regressed. This is a completed negative improvement round on synthetic toy instances. It is not completion of the three-round goal, a production-curve result, or a global-optimum claim.

[The generated tables](README.md) contain all variants, separately for online native time, cold instructions and cold native time. [RESULTS.json](RESULTS.json) retains full precision, phase diagnostics, per-cell counts and the exact decision. [RUNS.csv](RUNS.csv) retains all 3,243 canonical run keys and their exclusive ledgers. No scheduled pair failed or timed out in this round; the exporter still preserves those statuses when present.

## Frozen execution

[PR 772](https://github.com/aburan28/crypto/pull/772) merged the [protocol](../PROTOCOL.md), [16-arm registry](../round1.json), controls and runner at `5915da9d56758f81ceabbf978e796cd6be9740c3`. The measured reviewed head was `a34c8fb3b43200b1243dbad2d77bd7b621c2569d`; its 19 reviewed implementation files were byte-identical at the merge. [Workflow 36140265516](https://github.com/aburan28/crypto/actions/runs/36140265516) was dispatched once at 2026-09-25 13:19:52 UTC and completed successfully before 14:17 UTC. No measured process, source, fixture or promotion rule was retried or changed.

The measurement environment was Linux amd64, rustc 1.94.1, the musl target and Valgrind 3.22.0. Workers used one thread pinned to CPU index 3, an 8 GiB address-space cap and a 180-second child timeout. The maximum was 3,500 native/profile pairs within the declared workflow time budget. The recorded host/kernel is in `round/tournament/contract.json`; an exact CPU model was not captured, so it is unknown.

Five development cells were `n17a1`, `n19a0`, `n23a0`, `n23a1`, `n31a0`; the final panel added `n29a1`. The `n` label is the binary-field degree, not subgroup bits. Confirmation used 12 fresh points per cell, 72 total, and three processes per point. Replay repeated those same 72 points in fresh processes. The frozen historical census excluded 5,133 prior points, and the round preserves every generated fixture and admission, including candidates that were not measured in the final stages.

| Stage | Native/profile pairs | Verified |
|---|---:|---:|
| A/A | 30 | 30 |
| Smoke | 270 | 270 |
| Development | 810 | 810 |
| Selection | 405 | 405 |
| Confirmation | 864 | 864 |
| Replay | 864 | 864 |
| Total | 3,243 | 3,243 |

The preceding fixed-vector controls independently verified 30/30 pairs. Those controls do not supply natural-query yield or performance selection.

## Decision and interpretation

The primary interval starts with the supplied point after reusable preparation and ends after scalar replay. Fixture generation and reusable preparation are excluded. The complete cold costs remain additional acceptance gates; three process repetitions never become three independent targets.

`stop6` confirmation ratios to the incumbent were 1.064757 for online time, 0.956822 for cold Ir and 0.964563 for cold native time. Replay ratios were 1.071483, 0.956821 and 0.967258 respectively. The predeclared upper bounds for cold Ir exceeded one in both stages (1.008119 and 1.008117); online upper bounds and maximum cell ratios also failed. Neither cold point estimate reached 0.8. The frozen checker therefore retained the incumbent without promotion.

The incumbent's measured confirmation online time was 0.0347468 ms against 0.292392 ms for the separately qualified online rho reference: an observed ratio of 8.41491. This describes the registered one-point online interval on the toy panel. Cold native time did not establish a robust rho win; `beats_rho_strict` is false. The older decision field `beats_rho` describes only the cold-instruction point estimate and must not be promoted into a broader timing or crossover claim.

Development illustrates why stage winners are insufficient. `stop6_full` led online time with ratio 0.961748 but used 1.404833 times the cold instructions. `row_word` led cold instructions at 0.998918, a small whole-pipeline change. `stop4` led cold native time at 0.935889 but had online ratio 1.380165. The retained six were those three leaders, `row_bounded`, `row_suffix`, and the predeclared exploration arm `stop6`. Selection on separate points chose `stop6`. None of these development observations establishes a promotion. Any next-round design must use development evidence only and exclude all prior generated points; confirmation and replay are not tuning inputs.

The K-instruction floor is a weak implementation-specific full-rank collector floor, not a generic IC bound. Changing the factor-base policy is explicit; actual B/K values are retained. Base-only memory was not instrumented and remains unknown; peak process RSS is separate. Stopped-walk yield intervals resample whole target/walk clusters and are only descriptive. F4/F5/SAT and large sparse LA were not executed in this round. The tested mechanisms are pair-table PDP policies, Frobenius-orbit base construction and scalar-field Gaussian row kernels.

## Durable evidence and independent replay

The [archive](../../../evidence/ic-improvement-round1-20260925.tar.zst) is committed with [its manifest](../../../evidence/manifest.json):

- SHA-256: `f2b0f9b5c226ddc9331449d9e37ace1c2d7a94faa5ef74396d9a0a1231b13593`.
- Compressed bytes: 30,030,589.
- Restored regular files: ic-improvement-round1-20260925.tar.zstS; restored file bytes: 443,365,491.

It retains the original sources, worker binaries, fixtures, admissions, native/profile output, mathematical certificates, all raw profiles, stage summaries, decision, command logs, workflow metadata and controls. Build caches, Python bytecode and operation locks are excluded. Expanded profiles in the zstd container reconstruct to their original gzip bytes, which the restorer hashes. The archive also includes the initial dispatch checkpoint; its historical running state is superseded by `workflow-completion.json`.

[Fresh archive extraction](ARCHIVE_RESTORE.json) passed the archive hash, reconstructed-profile hashes, file count and byte count. The measured Linux checker verified all 3,243 records and 1,362 files across the three source snapshots. [Independent macOS transport replay](TRANSPORT_MACOS.json) reproduced the full frozen audit, including every summary and the decision, and rechecked all 30 controls without executing benchmark workers. A fresh Linux/Python 3.12 restore, full audit, byte-exact table export and scoreboard check are required by the evidence PR's CI. Successful verification does not turn the rejected challenger into a winner.

From the repository root, with Python 3.12 and zstd available, choose a new output directory and run:

```sh
python3.12 research/ic_candidate_tournament_20260915/evidence/restore.py \
  --archive ic-improvement-round1-20260925 --out /tmp/ic-round1-replay
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/improvement/audit_archive.py \
  --bundle /tmp/ic-round1-replay/ic-improvement-round1 \
  --out /tmp/ic-round1-audit.json
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/improvement/export_report.py \
  --round /tmp/ic-round1-replay/ic-improvement-round1/round/tournament \
  --out research/ic_candidate_tournament_20260915/goal_20260924/improvement/round1 --verify
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/improvement/render_report.py \
  --results research/ic_candidate_tournament_20260915/goal_20260924/improvement/round1/RESULTS.json \
  --markdown research/ic_candidate_tournament_20260915/goal_20260924/improvement/round1/README.md \
  --scoreboard docs/index-calculus-scoreboard.html --verify
```

These are read-only evidence replays, not new measurements. Do not dispatch the first round again. The repository skill and goal status point to this result and keep the remaining two-round budget open.
