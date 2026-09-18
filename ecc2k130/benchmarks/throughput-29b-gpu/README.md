# 29 B/s on this RTX PRO 6000

Written before the timings. Walk-rate engineering, not an ECDLP exponent claim.

The 20 B/s attempt on this SKU stopped at 17.298 B/s
([../throughput-20b-gpu](../throughput-20b-gpu/README.md)). This round asks
whether the leftover arithmetic in the product (hoisted 3-bit correction) and
the inverse chain (ONB multiply inside `inv131`), plus a compact-state
slot pipeline, can clear **29 B complete scalar updates/s**.

Boundaries, floor arithmetic, and the falsification target:
[THROUGHPUT-29B.md](../../THROUGHPUT-29B.md).

## Acceptance

Success is a verified median **> 29.0 B/s**. Same 5.3125 products/update.
300-report replay required before timing each binary. Inadmissible: two GPUs
as one, a different SKU, dropping identity, changing the product count.

```sh
bash benchmarks/throughput-29b-gpu/run.sh
```
