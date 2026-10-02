# Sibling Boolean pivot-trace feasibility screen

`PROTOCOL.md` and `protocol.json` freeze a structural question before any
candidate code or results. The inputs are generated public quadratic Boolean
systems. The study will compare exact degree-3 Macaulay matrices after the two
values of one branch variable, measuring the rank of their difference and the
survival of a deterministic pivot trace. The Boolean derivative is affine, but
that does not imply the matrix difference is low rank.

Status: **protocol only**. No discovery or holdout fixture has been generated,
no correctness result has passed, and no speedup has been measured. A native
Rust producer and verifier, immutable discovery evidence, and a gate decision
are required in a follow-on PR. The unused holdout seeds remain sealed in the
protocol. This work does not accept curve or key inputs.
