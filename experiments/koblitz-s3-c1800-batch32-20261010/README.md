# C1800 32-state target batch control

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
one through 32 over three field degrees. The next comparison pairs this
candidate against the eight-state C1800 extractor on the frozen T205–T212
single-target points, with each arm run in a fresh process.
