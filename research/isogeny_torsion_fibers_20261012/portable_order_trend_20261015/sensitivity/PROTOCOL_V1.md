# Native price-sensitivity calculation, V1

**Frozen before calculation.** Inputs are the four already independently
checked child-run files: order 229 `../../different_j_order_batch_20261015/RUN_V1.jsonl`
(SHA-256 `d0e4616236aeb8978f5ac1549dde373be294bdf52503c954a6f964904f68faf8`),
and orders 421, 433 and 457 at `../BATCH_457_421_V1.run.jsonl`,
`../BATCH_409_433_V1.run.jsonl` and `../BATCH_421_457_V1.run.jsonl`
(SHA-256 `b1054cd665cae39d07a4f64adbb1132f4d5756b59e8b3e44e1b46ce4e32009fb`,
`b167728e7c0809803f387c6db8ddad4520a5fce603774a74c7d981947f1d4632`,
`d75247c817c17c0bb0966f06792fb4d58f74d24434d0e83dcba45396ec98caa0`).
The matching replay certificates in the parent study must be checked before
derivation. No child is dropped or rerun.

For every prefix `N=1,4,16,64` and round `0,1,2`, reconstruct all cold
group additions and doublings of the full-image and source negation-rho
arms, including common validation, setup, every target and final replay.
Sum `attempts[].steps` over rho targets and require it positive. Define

`W_break = (G_full - G_rho) / rho_steps` group-addition equivalents per
observed rho step. If `G_full <= G_rho`, record that fact instead of a
positive threshold. The formula sets *all* omitted full-image work to
zero, then asks what extra cost **only rho** would need per step to reach
the full-image group-count lower bound. Any positive omitted work in
full-image raises the required rho charge. It does not estimate the true
cost of a hash, table probe, root or row operation, and is not a complete
`S` or speed ratio. Keep the separate Fp M/I axes and memory in the
parent report, since they overlap group operations.

The native analyzer must verify run-file hashes, all 84 child keys per
cell, verified scalar outputs, matching target counts, the published
replay `run_sha256` values, and the `TREND_V1.json` group totals before
computing thresholds. Emit all 48 paired rows, per-cell N64 medians and
an SVG. A mismatch fails closed. The result report and rendered PDF must
name the worst and median thresholds, exact assumptions, and the
remaining need for pinned native-operation conversion and isolated
timing. This derivation is supplementary to the complete portable-price
gate in `../PRICE_PROTOCOL_V2.md`.
