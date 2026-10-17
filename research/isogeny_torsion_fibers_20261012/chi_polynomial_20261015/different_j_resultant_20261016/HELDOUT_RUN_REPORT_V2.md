# Run receipt: four held-out order-229 targets

Run ID: `different_j_resultant_20261016/heldout-five-arm-v2`.
The V2 spec freezes the source and binary hashes, four public target
positions, five arms, three rotated rounds, trial caps, and 300-second
child cap. The reference is same-target signed rho with an A/A arm.
All 60 children exited zero and passed separate Rust scalar and
relation-witness replay. The raw files and every exit and meter
receipt are in `heldout_five_arm_v2/`; their SHA-256 values are in
`HELDOUT_FIVE_ARM_V2.artifacts.sha256`.

Measured mean cold `S_group` across 12 children per arm was 140.672
for χ resultant, 88.908 for direct pair, 7.979 for signed rho and its
A/A repeat, and 4.361 for source BSGS. χ and direct pair had the
same 33.58 mean trials, 13.83 mean relations, rank nine and verified
scalar in every paired row. Mean whole-child CPU was 14.304, 6.538,
6.428, 6.405 and 5.899 ms in that order on a shared Darwin arm64
host. Solver-thread CPU and full meter values are listed separately
in the derived summaries. No wall ratio is admitted.

The chart is [`HELDOUT_GRAPH_V2.svg`](HELDOUT_GRAPH_V2.svg).
V1 is retained with its setup-charging defect; V2 is the corrected
run. Three planned relation arms and complete portable pricing
remain pending.
