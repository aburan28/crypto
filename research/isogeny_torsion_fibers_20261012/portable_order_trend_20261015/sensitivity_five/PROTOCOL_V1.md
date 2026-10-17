# Five-order missing-price sensitivity protocol

**Frozen before the derivation.** This extends the committed four-order
`W_break` calculation to the independently replayed selected order-1,201
cell, without rerunning a DLP or changing any public target, cap or arm.
Inputs are the five run/replay digests and the native five-order counted
trend in [`RUN_SPEC_V1.json`](RUN_SPEC_V1.json). The executable source
and release binary hashes are frozen there. The first four results are
historical controls; the fifth is the only new cell.

For every `r`, prefix `N=1/4/16/64`, and round `0/1/2`, recompute the
cold full-image and matched direct-rho group additions plus doublings
from the raw record. They include common work, setup, queries and final
verification. Confirm equality to the independently derived five-order
trend row. Sum all observed rho walk steps, including restarts. Return

`W_break = (G_full - G_rho) / rho_steps`

in additional group-addition equivalents per observed rho step. It is
the omitted rho price required to tie the counted full-image lower bound
**if all omitted full-image work costs zero**. The denominator includes
all observed steps and must be positive. Report 60 matched rows, the
median/minimum/maximum of the three N64 rounds per order, and preserve
the four earlier values exactly. Each run digest, independent replay,
84 unique children, 1,785 verified scalar outputs and 64-target public
prefix is a fail-closed input gate. Any mismatch exits before an output
is written. The output path must not exist before execution.

The reference is the same-target direct rho arm. This calculation is a
price **sensitivity threshold**, not a measured primitive price,
complete portable `S`, wall-speed claim, or fitted exponent. The
existing [complete-price protocol](../PRICE_PROTOCOL_V2.md) still
requires disjoint native counts, calibrated common conversion and
second-host verification. A threshold smaller than a future measured
rho overhead would require computing the actual complete totals; a
larger threshold would preserve the full-image disadvantage for that
row under the protocol's bound. The four-order V1 derivation stays
unchanged; this V1 is a new five-order versioned addition.
