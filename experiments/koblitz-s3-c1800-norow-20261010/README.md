# C1800 full-rank target row elision

The C1800 extractor now recovers a target directly from its four verified
factor-base labels when reusable full-rank logs are ready. It skips the
1,800-column relation row previously built and discarded on that path. The
partial-rank fallback still builds the row for its span checks.

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
example tests passed. The next stage pairs this candidate with C1800 on
T205–T212 in separate fresh-process one-target runs.
