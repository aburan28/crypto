# C1800 32-state target batch screen

The matched eight-point N53 screen completed **16/16 verified one-target
solves**, and checked Sage independently replayed every scalar. The geometric
batch-8/batch-32 online wall ratio was **0.868×**, with a target-bootstrap 95%
interval of **0.729×–1.007×**. This falls below the frozen **1.10×**
advancement gate. Batch-32 retains the same first relation witness on all
eight points, while larger batch rounding changes the number of states
examined on five of them.

| Point | Batch-8 online, ms | Batch-32 online, ms | Batch-8 S3 calls | Batch-32 S3 calls |
| --- | ---: | ---: | ---: | ---: |
| T205 | 0.236 | 0.449 | 128 | 128 |
| T206 | 0.095 | 0.118 | 128 | 128 |
| T207 | 0.103 | 0.132 | 128 | 128 |
| T208 | 0.118 | 0.106 | 160 | 192 |
| T209 | 0.240 | 0.234 | 560 | 576 |
| T210 | 0.129 | 0.124 | 224 | 256 |
| T211 | 0.261 | 0.338 | 624 | 640 |
| T212 | 0.261 | 0.247 | 672 | 704 |

The target query plus PDP phase had a geometric batch-8/batch-32 ratio of
**1.165×**. The charged descent residual had a ratio of **0.281×**. The
extractor computes an 1,800-column target relation row even when its prepared
full-rank logs allow direct recovery from four labels. Removing that unused
online row is the next separately measured implementation change; these
phase data establish the bottleneck to examine, not its cause.

[`PAIR_FREEZE.json`](PAIR_FREEZE.json) binds both candidate manifests,
source and binary hashes, T205–T212 point files, paired order, resource
envelope, and runners. The [`analysis`](analysis.json),
[`run ledger`](local/paired_screen/ledger.jsonl),
[`phase summary`](phase_summary.json), compressed per-run reports, and
[`Sage replay`](sage_replay.json) retain the full evidence. The
[`Sage runtime receipt`](sage_runtime_info.json) was saved before replay.
These local CPU timings have no host-wide isolation receipt and remain
exploratory.

## Initial control

The C1800 extractor now groups up to **32 ordinary target S3 states** into
one normal-basis inversion and one root-table prefetch group. Its frozen N53
T025 control reached prepared rank 1,800, recovered and checked the public
target scalar, used **512 target S3 calls**, and found the same relation
indices `[55035, 67241, 57152, 189185]` as the eight-state C1800 control.
The control stayed inside the 32 GiB observed-RSS envelope and used no
online rank fallback.

The 32-state control's target online diagnostic was 0.367 ms, compared with
0.250 ms for the archived eight-state control on separate local executions.
Those timings motivate a paired screen; they do not establish a wall-time
ratio. The new source preserves the C1200 embedded prefix, C1800 base, exact
root table, projective scalar checks, and 14 preparation workers with one
target worker. Its batching and prefetch width are the algorithmic variable.

[`CONTROL_FREEZE.json`](CONTROL_FREEZE.json) binds the source, binary,
checked public point, parent and reference receipts, resource envelope, and
runner. The [`control receipt`](local/control/receipt.json) and
[`compressed raw report`](local/control/report.json.gz) retain the full
result. All 16 Rust example tests passed, including all batch lengths from
one through 32 over three field degrees.
