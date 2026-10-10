# C1800 full-rank target row elision

The C1800 extractor now recovers a target directly from its four verified
factor-base labels when reusable full-rank logs are ready. It skips the
1,800-column relation row previously built and discarded on that path. The
partial-rank fallback still builds the row for its span checks.

The matched T205–T212 screen completed **16/16 verified one-target solves**,
and checked Sage independently replayed each scalar. The geometric C1800
reference/row-elision online ratio was **1.131×**, above the frozen **1.10×**
advancement gate. Its target-bootstrap 95% interval was **0.557×–2.134×**.
The charged descent phase ratio was **4.470×**; query plus PDP was **1.141×**.
All eight pairs used the same first relation witness and S3 call count. The
local host has no isolation receipt, and whole online timings include large
query/PDP outliers in both arms. Repeated paired confirmation is the next
measurement step.

| Point | Reference online, ms | Row-elision online, ms | S3 calls per arm |
| --- | ---: | ---: | ---: |
| T205 | 0.161 | 0.111 | 128 |
| T206 | 0.497 | 0.097 | 128 |
| T207 | 0.095 | 0.623 | 128 |
| T208 | 0.111 | 0.174 | 160 |
| T209 | 0.222 | 0.219 | 560 |
| T210 | 0.168 | 0.118 | 224 |
| T211 | 0.233 | 0.251 | 624 |
| T212 | 0.832 | 0.297 | 672 |

[`PAIR_FREEZE.json`](PAIR_FREEZE.json) binds both candidate manifests,
source and binary hashes, paired points and order, workers, and resource
limits. The [`analysis`](analysis.json),
[`phase summary`](phase_summary.json),
[`run ledger`](local/paired_screen/ledger.jsonl), compressed raw reports,
[`Sage replay`](sage_replay.json), and
[`runtime receipt`](sage_runtime_info.json) retain the evidence.

## Initial control

The frozen N53 T025 control reached prepared rank 1,800, verified the target
scalar, used the same **512 S3 calls** and first relation indices as C1800,
and stayed inside the 32 GiB observed-RSS envelope. Its charged target
descent residual was **0.0011 ms**, versus **0.0106 ms** in the archived
C1800 control. The whole online interval was 0.269 ms versus 0.250 ms on
separate local executions, so a paired screen is needed for the full-cost
comparison.

[`CONTROL_FREEZE.json`](CONTROL_FREEZE.json) binds the exact source and binary
hashes, T025 point, reference receipt, workers, and limits. The
[`control receipt`](local/control/receipt.json) and
[`raw report`](local/control/report.json.gz) retain the result. All 16 Rust
example tests passed.
