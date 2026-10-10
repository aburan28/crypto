# Immediate S3 target lookup with projective target checks

Immediate root-order lookup reduced the geometric target S3-call count by
**1.629×** on 12 frozen N41 C1200 public points, while preserving every
recovered scalar and the 98,400-point factor base. Its paired single-target
online ratio was **0.896× prefetch/direct** with a target-bootstrap 95%
interval of **[0.187×, 4.856×]**. The predeclared 2× development gate did not
pass. The 24-process panel completed, and checked Sage independently replayed
all 24 scalars. The large timing spread and absent host-exclusive CPU receipt
leave `controlled_speedup` null.

## Frozen C1200 screen

[`FREEZE.json`](FREEZE.json) binds the source and binary hashes, candidate
manifests, 12 existing public points T210–T221, one-target workloads, alternating
arm order, the checked Sage runtime, 14 preparation workers, one target worker,
a 360-second process cap, and a 32 GiB observed-RSS envelope. The comparison
changes only immediate root-order lookup versus eight-state speculative
prefetch with the same ARM64 PMULL fold. The direct candidate is
`IC1N41Ckb1fb98400PDP4rootRCsampleLAgaussTDdirectISO0hb60dc51124bb`;
the prefetch candidate is
`IC1N41Ckb1fb98400PDP4rootRCsampleLAgaussTDdirectISO0hc4ef08037c5a`.
Both used full prepared rank 1,200 and the same base digest. Each online clock
starts after reusable preparation and ends after scalar recovery check; all
target-dependent attempts are charged.

| Frozen measure | Prefetch | Direct lookup |
| --- | ---: | ---: |
| Verified one-target process rows | 12 / 12 | 12 / 12 |
| Median online interval, ms | 0.895 | 1.684 |
| Median reusable setup, seconds | 47.675 | 49.436 |
| Geometric paired prefetch/direct online ratio | — | 0.896× [0.187×, 4.856×] |
| Geometric paired prefetch/direct S3-call ratio | — | 1.629× |

The online intervals ranged from 0.061 to 54.325 ms across the 24 runs. In
the frozen order, the second process took 3.410× as long online as the first
process geometrically, and 9 of 12 second processes were slower. The process
counter recorded 616 major page faults in 14 rows; it cannot locate those
faults within the target interval. Reusable setup ranged from 27.824 to
82.011 seconds. During the run, the Mac's sampled load average was about
86–95 on 14 logical CPUs, and the system volume had under 1 GiB free near
completion. The alternating order balances arm assignment but does not remove
this noise. The exact per-target times, calls, phase accounting, failures,
resource counters, and scalar certificates are retained in the raw reports.

The [analysis](local/analysis.json) binds the complete ledger SHA-256
`bdf7d50f35de93b3155ea215d75c694e5ed8af4f17bf1cd4958c4193f72a07ea`;
the [independent Sage receipt](local/sage_replay.json) checks every scalar.
The [raw archive](local/raw_screen12.tar.zst) has SHA-256
`e52d7ebaf1a5ab3ec10f328e5fe0c8929b97846be6f6b8a21d6ae4d95ce23723`;
extracting its `raw.jsonl` reproduces the analysis ledger hash. These reused
points are a development screen, so a variant clearing the development gate
would need fresh points before confirmation. The next measurement should
prioritize the C1800 plus PMULL candidate, whose separate 48-point panel
cleared the exploratory 2× gate, on a host with an auditable isolation receipt.

The new [direct-lookup source snapshot](../../examples/s3_partial_prepared_fullrank_projective_directlookup_targetspan.rs) checks each S3 partner against the root table as soon as its roots are available. It retains the eight-state shared inversion, the target search order, projective subgroup validation and scalar replay, and the opt-in ARM64 PMULL fold. A verified relation can now stop the remaining roots in its batch. On the frozen **C1200** public target T220, this reduced target S3 calls from **16 to 2** while preserving the recovered scalar and factor-base digest. The release build is bound in [build_receipt.json](local/build_receipt.json).

The bounded C1200 direct-lookup retry `R2` recovered scalar `254961180376` with 98,400 distinct subgroup-usable factor-base points, full prepared rank 1,200, and two target S3 calls. Its reusable preparation took **160.736 seconds**; its target online and charged PDP intervals were **0.114 ms** and **0.069 ms**. The earlier prefetch run of the same target and base made 16 S3 calls. It ran at a different time under substantial host contention, so the two online times are not a paired speedup estimate. Checked Sage independently replayed both scalars; see the [retry receipt](local/c1200_retry_receipt.json), [C1200 Sage receipt](local/c1200_smoke_replay.json), [direct report](local/direct_c1200_T220_R2.json), and [prior prefetch report](local/prefetch_c1200_T220_prior.json).

A subsequent prefetch solve of the same point completed under the same 360-second cap, again with 16 S3 calls. Against that run, direct lookup's charged PDP interval was **0.069 versus 0.098 ms** (1.414× prefetch/direct), and online time was **0.114 versus 0.138 ms** (1.208×). Both scalars passed checked Sage replay. The [one-pair diagnostic](local/c1200_pair_r2_diagnostic.json), [prefetch report](local/prefetch_c1200_T220_R2.json), and [Sage receipt](local/c1200_pair_r2_sage_replay.json) retain the measurements. Reusable setup differed by more than fourfold across these consecutive processes, so the ratios describe this point under this host state and are not an aggregate result.

On public N41 target T220 with a **244-column diagnostic index**, both the direct-lookup and prefetch binaries recovered scalar `254961180376` using 20,008 actual subgroup-usable factor-base points and full prepared rank 244. Checked Sage independently replayed both scalars on the frozen public point. The direct path made **1,170 S3 calls** versus **1,184** for prefetch. This single pair had online intervals of 0.719 ms and 0.403 ms, and charged target-PDP intervals of 0.654 ms and 0.385 ms, respectively. These timing values are diagnostics from the contended Mac, not an aggregate comparison. See the [Sage receipt](local/smoke_replay.json) and the [direct](local/direct_c244_T220.json) and [prefetch](local/prefetch_c244_T220.json) native reports.

The first C1200 T220 smoke reached its 180-second process cap with an empty report. The exact stopping phase is unknown. Its [failure receipt](local/c1200_timeout.json) retains the source and binary hashes, public-point hash, command shape, and cap. This host's earlier C1200 setup times varied sharply, and workspace free space fell from roughly 47 GiB to 10 GiB while the release build and smoke were running. The completed retry's 160.736-second reusable setup remains much slower than the prior panel's 31.722-second median for the prefetch arm, reinforcing the need for a fresh host preflight before a timed panel.

The [preparer](prepare.py), [runner](run.py), [analyzer](analyze.py), and [checked Sage replay](replay_sage.py) produced the complete frozen screen above. The candidate ID is tied to the commit that introduced the source and remains stable across unrelated commits. The local panel confirms the call-count reduction and directs the next experiment toward the stronger C1800 plus PMULL configuration under the CPU isolation gate.
