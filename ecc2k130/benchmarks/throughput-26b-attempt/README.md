# 26 B/s attempt

Canonical verdict: [`../../THROUGHPUT-26B.md`](../../THROUGHPUT-26B.md).

The exact DP4A phase selector passed 300-report device replay and 4,096-point
host differential testing, but measured 19.515 B/s against the 20.134 B/s
reference. The 26 B/s target remains above the current walk's derived 22--25
B/s one-addition boundary.

- `build-dp4a.log`: CUDA compiler and resource report.
- `verify-dp4a.log`: 300/300 scalar replay, zero dropped.
- `bench-dp4a.log`: three complete throughput samples.
- `*-fused-weight-phase.log`: exact fused selector build, replay and timing.
- `result.json`: frozen table values and boundary ratios.
