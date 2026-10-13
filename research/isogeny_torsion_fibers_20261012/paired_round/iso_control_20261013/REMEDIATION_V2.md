# Version 2: six-class distinguished-point control

Frozen before v2 execution. The original protocol and v1 binary used the
negation walk's default distinguished-point bits (`floor(bits(r)/4)`). The
[v1 raw row](raw-v1.tsv) for `p=13`, target 4 reports `rho_failed` after
66,667 walk steps, 33,333 distinguished points, 52,778 starts, and the
100,000 charged-walk-operation cap. The other four arms on that target
verified scalar 10. This is a walk-parameter/collision failure in that
bounded run, not evidence against the six-class group action. The v1
source and raw output remain unchanged.

V2 changes only the six-class arm's `RhoWalk.dp_bits` to `Some(0)`, so every
nonidentity class representative reached after a step is distinguished.
All other walk tuning, exact targets, seeds, curves, relation parameters,
map accounting, independent checks, and 100,000-operation cap remain as in
[PROTOCOL.md](PROTOCOL.md). The runner will emit all 80 arm rows even if a
target fails, then exit nonzero; this preserves complete failure coverage.
V2 raw and replay outputs use separate immutable names. Success still
requires all 80 rows to verify, byte-identical replay, and independent
arithmetic validation. The v1 failure is included in the final analysis
regardless of v2 outcome. Neither rho variant establishes a best tuned
six-class rho reference at larger orders.
