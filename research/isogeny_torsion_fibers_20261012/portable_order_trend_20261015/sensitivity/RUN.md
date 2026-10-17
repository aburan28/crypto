# Native price-sensitivity derivation: run record

**Date and host:** 2026-10-15, Darwin arm64 shared host (L0). This is a
deterministic derivation of existing count records; its elapsed time is not
used as a performance measurement.

**Protocol:** [PROTOCOL_V1.md](PROTOCOL_V1.md), committed before this
calculation. The four input run SHA-256 values and the matching replay
certificates are pinned there. The executable source is
[`examples/different_j_price_sensitivity.rs`](../../../../examples/different_j_price_sensitivity.rs).

**Build and exact command:**

```sh
CARGO_TARGET_DIR=target/order-batch-check cargo build --release --example different_j_price_sensitivity
study=research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015
out="$study/sensitivity"
target/order-batch-check/release/examples/different_j_price_sensitivity "$study" "$out/SENSITIVITY_V1.json" "$out/SENSITIVITY_V1.svg" > "$out/SENSITIVITY_V1.stdout" 2> "$out/SENSITIVITY_V1.stderr"
```

The executable exited **0** ([exit record](SENSITIVITY_V1.exit)); stdout
and stderr are retained. The [JSON](SENSITIVITY_V1.json) contains all 48
paired rows and four N64 summaries. It checks the four complete run-file
digests, matching replay digests and certificates, the 84 child keys and
1,785 scalar outputs per cell, and all 48 corresponding group totals in
[`TREND_V1.json`](../TREND_V1.json). It refuses changed or missing files.

SHA-256: analyzer source
`356a7a52e6e069497c245585b7234b0963bcb4556e46505dd5a497a07bd67637`,
derived JSON
`a96ff6e33d77b5053723635bf5ec504afe9d31fe9e721231f34faa3e1ce0bc41`,
and SVG
`7c5da1fe0c4d8723f30ab3589d21fb367c0884de1f13539c4943e9fe78f2f79e`.

This run adds no new solver child and changes no frozen input. Its output
is an algebraic missing-price threshold, not a conversion measurement.
The [results report](RESULTS.md), [SVG](SENSITIVITY_V1.svg), [PNG](SENSITIVITY_V1.png),
and [PDF](RESULTS.pdf) carry the visual record. The complete portable
pricing and isolated-wall gates remain in
[`PRICE_PROTOCOL_V2.md`](../PRICE_PROTOCOL_V2.md).
