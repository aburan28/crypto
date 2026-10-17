# N53 fixed-useful-base S3 root batching: nine verified single-target pairs

**Status:** engineering result in source-level root operations; shared-host
wall measurements are L0 diagnostics. The frozen
[protocol](PROTOCOL.md), [inputs and hashes](FROZEN.json),
[native pair analysis](ANALYSIS.json), and
[native aggregate](RESULT_DATA.json) are the authoritative records.
The [run report](reports/run-20261012T1737Z/REPORT.md) and its
[paired comparison](reports/run-20261012T1737Z/paired_ratios.svg) and
[rho comparison](reports/run-20261012T1737Z/rho_reference.svg) show the
measurements. The earlier [sandbox-wrapper failure](reports/run-20261012T1735Z/REPORT.md)
is retained separately.

## Decision against the frozen gates

The control and candidate used the same constructed K220 base of 23,320
points and 220 signed-Frobenius columns on
`EC1N53Ce0hb097de99be9a`, subgroup order 21,044,858,204,113. The three
public target seeds, three rank seeds per target, rho seeds, and cell order
were committed before timing. All 27 fresh processes completed within the
90-second and 16-GiB caps. The nine control/candidate pairs each reached
rank 220/220 with 220 accepted independent relations. Native replay verified
all 18 IC results and the nine matched strong-rho results, every accepted
relation and solved factor log, the target scalar, and the exclusive phase
sums. The base dumps, root-table sizes, consumed-root counts, rank probes,
and selected rank traces match within every pair.

The exact-output gate and every pair's at-least-50%-fewer-inversions gate
passed. The descriptive progression gate also passed: median candidate over
control is 0.738 cold (required at most 0.90) and 0.612 target-online
(required at most 1.10). These ratios include the complete IC process or
the explicitly separated target phase, respectively. They do not price
common-operation `S` against rho and do not qualify a controlled native
speed comparison on this shared host.

| Variant / reference | Cold process, median paired ratio | Target-online, median paired ratio | Correctness | Class |
|:--|--:|--:|:--|:--|
| control / strong rho, same Q | 78.146 [44.220, 92.602] | 0.490 [0.168, 0.763] | 9/9 replayed | baseline diagnostic |
| batch / strong rho, same Q | 52.586 [33.179, 64.668] | 0.378 [0.103, 0.442] | 9/9 replayed | engineering, L0 diagnostic |
| batch / control, same Q and rank seed | 0.738 [0.714, 0.750] | 0.612 [0.569, 0.903] | 9/9 exact pairs | engineering |

Ratios are unitless, from in-process milliseconds. Brackets are descriptive
95% target-block bootstrap intervals: 10,000 deterministic resamples of the
three public target blocks, retaining all three rank seeds within each block.
The complete nine-value arrays and exact medians are in
[`RESULT_DATA.json`](RESULT_DATA.json). A separate wall-clock wrapper and
peak RSS record is retained for each process. The cold ratio against rho
remains above 1 in every paired cell (batch: 21.864–426.310), so this
implementation does not meet the repository's one-target cold reference.

| Rank-stage root primitive, sum over nine IC cells | Control | Batch | Batch / control |
|:--|--:|--:|--:|
| inversion primitives | 80,790,012 | 763,167 | 0.009446 |
| field multiplication calls | 1,211,850,180 | 902,063,394 | 0.744369 |
| field squaring calls | 242,370,036 | 39,684,684 | 0.163736 |

These source-level counts charge batch preparation and roots prefetched
before early exit. The Montgomery batch uses one inversion for each index
state's regular-root group; all generated roots consumed by the original
scan retain their order. The field-call counts come from the instrumented
native source, not an instruction or group-addition-equivalent calibration.
The first seed, for example, makes 9,785,496 candidate root calls (including
11,382 prefetched calls) versus 9,774,114 control calls, while both consume
9,774,114 roots. This is an implementation gain with the same rank path,
not evidence of a new relation-yield exponent or generic-floor crossing.

## Boundary and cost scope

The declared generic-group reference floor is
`S_floor = sqrt(pi/(2*106))` for the 106-element signed-Frobenius orbit,
where `S = total group-addition equivalents / sqrt(r)`. The complete
compact-orbit example has no calibrated conversion from its field operations,
normal-basis transforms, memory accesses, and rho operations into one `S`
unit. Consequently `S`, the ratio to the floor, and any full-method
operation-count speed ratio remain unset. The table reports the same-point
strong-rho timing reference as a diagnostic only. Setup (base construction,
normal-basis checks, index, rank, matrix and LA), target PDP/descent/recovery,
and rho's own setup and recovery were included in their respective cold
intervals. No amortized or multi-target cost is substituted.

The host was an Apple M4 Pro (`Mac16,7`, 48 GiB), macOS 26.6,
`rustc 1.93.1`, release profile. No core reservation or contention-free
receipt was obtained, so this is L0. Wrapper peak RSS was 341–351 MB for
IC and 8–20 MB for rho, all below the frozen 16-GiB cap. The independent
replay was on the same host class; a second-host replay and a calibrated
common-`S` comparison remain open. A source-level operation reduction
does not itself measure hardware throughput.

## Reproduction and custody

The executable source commit is `00a8da2d16be88ef9be799dbd51da3103a81e5bb`;
the three public inputs and binary hashes were pushed in `576f7087b`.
The valid execution began after the failed wrapper attempt was preserved in
`e70a0fc4e`; neither the frozen experiment binary nor the inputs changed.
The 27 intents, raw outputs, stderr, statuses, base dumps, rank traces, and
replay receipts are under [`runs/`](runs/). The analysis command was:

```sh
target/release/examples/n53_s3_root_batch_run analyze \
  target/release/examples/koblitz_s3_root_batch_experiment \
  target/release/examples/koblitz_rho_batch_ks_strong_online \
  target/release/examples/koblitz_rho_fixture \
  target/release/examples/koblitz_one_target_replay
```

The native aggregate was produced by
`examples/n53_s3_root_batch_report.rs` from `ANALYSIS.json` with
`target/release/examples/n53_s3_root_batch_report
research/notes/ecc2k130/n53_s3_root_batch_20261012/ANALYSIS.json
research/notes/ecc2k130/n53_s3_root_batch_20261012/RESULT_DATA.json`.
The
[rendering source](reports/render_figures.py) only draws those recorded
native values. The corresponding [PDF](RESULT.pdf) includes both graphs.
The benchmark copies are under
[`benchmarks/20261012-n53-s3-root-batch/`](../../../../benchmarks/20261012-n53-s3-root-batch/).

The canonical scoreboard and timeline gain the two L0 candidate/rho
diagnostic points and retain the prior n53 points as history. The ecbench
every-candidate panel is unchanged because this is a native example run,
not an ecbench session. Curve-identity and algebraic-transfer diagrams are
unchanged: this experiment constructs no new curve or transfer map.
`docs/ic/boundary_targets.json` is unchanged because the frozen protocol
does not admit a calibrated `single_target_online` row.

## Validation

`cargo test --release --lib`, the touched release example tests,
and the four targeted site/scoreboard checks were run on this branch.
The full site test suite has one unchanged baseline error:
`test_lab_browser_index_is_current_and_self_consistent` rejects
`docs/curves/covers.json` because its recorded registry digest
`89b580c5…` differs from the current registry digest `5520bc7e…`.
Both files are untouched here and have the same mismatch on
`origin/main`. This report does not relabel that preverified cover
record without re-running its native mathematical checks.
