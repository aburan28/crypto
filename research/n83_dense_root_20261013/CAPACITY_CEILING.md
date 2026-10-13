# Pre-panel capacity-ceiling extension

The pinned Lean bound identifies K1798 as the last prefix width whose dense
packed-table request fits 2^29 slots (8 GiB of allocated slots) under
`states <= 83*K^2`; K1799 crosses this capacity threshold at the maximum
state count. This extension freezes one additional `K1798` dense arm **before
the 64-point panel is run**. It tests the largest base in that allocation
bucket, not an assumed monotonicity theorem for observed target probes.

Run `run_capacity_ceiling.sh <producer-binary> primary` after the four-arm
primary comparison, and `run_capacity_ceiling.sh <producer-binary> panel`
after same-point rho and the four-arm panel. The producer, seed 29, 12 rank
threads, one target thread, 900-second cap, public point input, and binary
hash are identical to the four-arm protocol. The panel run archives the full
point-defined base for independent four-point replay. Any cap or resource
failure remains an explicit censored row. Selection over all completed arms
uses the already frozen mean-probe rule on the same 64 public points; online
wall times remain exploratory on this shared host.
