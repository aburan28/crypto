# Five-order missing-price derivation, 2026-10-16 01:53 UTC

The [protocol](PROTOCOL_V1.md) and [input/source/binary digest spec](RUN_SPEC_V1.json)
were committed as `52fcf45b6` before execution. On the shared Darwin
arm64 host, the frozen release executable ran:

```sh
target/release/examples/different_j_price_sensitivity_five \
  research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015 \
  research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/sensitivity_five/SENSITIVITY_FIVE_V1.json \
  research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/sensitivity_five/SENSITIVITY_FIVE_V1.svg
```

The process exited **0**; [stdout](RUN_V1.stdout) and
[stderr](RUN_V1.stderr) are retained, as is the [exit file](RUN_V1.exit).
The native analyzer checked five frozen run digests and their independent
scalar replay receipts, 420 unique child records, 8,925 scalar outputs,
and equality of every recomputed group count to the five-order trend.
It wrote [60 matched rows and five N64 summaries](SENSITIVITY_FIVE_V1.json),
SHA-256 `58d6d193cee0800f425b264893712c79a909abb0ce33a7cf757c4c989a18d250`,
and [SVG](SENSITIVITY_FIVE_V1.svg), SHA-256
`604c669ab96483e6b0f59813dbf2047564f6a61b2b0fa3514d209573e13c15c7`.
The previous 48 rows and four N64 summary objects were byte-equal after
canonical JSON serialization to the earlier V1 derivation. The new
order-1,201 N64 rows had `(G_full,G_rho,rho_steps)` of
`(3,398,567,9,030,2,437)`, `(3,396,221,8,676,2,150)`, and
`(3,400,936,9,297,2,678)` for rounds 0, 1 and 2. Their derived
thresholds are 1,390.865, 1,575.602 and 1,266.482 GAE per observed
rho step.

The [source-linked visual run report](../reports/run-20261016T0153Z-five-price/REPORT.md)
and [results report](RESULTS.md) state the assumption and limits. No
new solver or target was run here; this is a deterministic count
derivation from completed, independently replayed solves. Host CPU
timing is irrelevant to the threshold.
