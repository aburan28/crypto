# Can one RTX PRO 6000 clear 18 B scalar updates/s?

Written before the timings. This is walk-rate engineering, not an ECDLP
exponent claim.

## Boundaries

- **Floor:** one-addition-per-step walk, ≈ 23–25 B/s
  ([THROUGHPUT-30B.md](../../THROUGHPUT-30B.md),
  [ITERATION-FUNCTION.md](../../ITERATION-FUNCTION.md)). 18 is under the
  floor, so the hardware does not forbid it.
- **Reference:** shipping product on the same SKU, 15.115792 B/s matched
  control from [top-clmad](../top-clmad/summary.json).
- **Unit:** billions of completed scalar updates per second.
- **Priced leftover:** table walk (`WALK_TABLE=1`) measured 16.56 B/s
  against 14.41 on an earlier pair; `TABLE_PIVOT_BYTES=1` is the one
  unmeasured ALU cut, priced at about −2.4% if the ALU still binds
  ([THROUGHPUT-20B.md](../../THROUGHPUT-20B.md)). Scaling the 16.56/14.41
  ratio onto today's 15.116 control predicts ≈ 17.4 B/s; plus the pivot,
  ≈ 17.8. That is the top of what this tree can price. It is still
  short of 18.
- **Class:** engineering. The table walk is a different iteration and
  cannot join the live campaign.

## Acceptance, written before the run

Success is a verified median **> 18.0 B/s** on one RTX PRO 6000, identity
matching the requested walk, every sample valid, same 5.3125 products per
update. Inadmissible: changing the GPU SKU, counting two GPUs as one,
dropping identity, changing the product count, or quoting a model as a
rate.

If every priced leftover stays ≤ 18.0, this thread stops. The remaining
distance is the product and the reduction, which are the floor.

## Result (2026-09-17)

Two successive Modal RTX PRO 6000 allocations, driver 580.95.05, automatic
workers (96,256 threads), 50,465,865,728 updates per sample, identity
matched (`packed table walk: 1`, pivot 0 then 1, shared 48,508 / 48,732
bytes).

| variant | median B/s | / 18 | / shipping 15.116 | / 23 B floor | correctness | class |
|---|---:|---:|---:|---:|---|---|
| shipping product | 15.116 | 0.840 | 1.000 | 0.657 | top clmad 0 | reference |
| table walk | 16.601 | 0.922 | 1.098 | 0.722 | walk 1, pivot 0 | engineering |
| table walk + byte pivot | 16.667 | 0.926 | 1.103 | 0.725 | walk 1, pivot 1 | engineering |

Neither arm clears 18. The byte pivot is +0.4% on the table walk, not the
priced +2.4%. The campaign default stays the shipping walk. This thread
stops: the leftover that could still be priced has been measured.

## Recipe

Automatic workers: that is the geometry that produced 16.56 B/s.

```sh
make bench-rtx-pro6000 RTX_PRO6000_WALK_TABLE=1 RTX_PRO6000_WORKERS=0
make bench-rtx-pro6000 RTX_PRO6000_WALK_TABLE=1 RTX_PRO6000_TABLE_PIVOT_BYTES=1 RTX_PRO6000_WORKERS=0
```
