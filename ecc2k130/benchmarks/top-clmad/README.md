# RTX PRO 6000 `PACKED_TOP_CLMAD` paired bench

Measure the shipping clmad product against `PACKED_TOP_CLMAD=1` on one
Modal RTX PRO 6000, using the documented Make targets in
[TOP-CLMAD.md](../../TOP-CLMAD.md). This is an engineering it/s
comparison of the same iteration function. It is not an ECDLP exponent
claim and does not belong on the index-calculus scoreboard.

## Boundaries, stated before the run

- **Floor:** one-addition-per-step walk, about 23–25 B scalar updates/s
  on this GPU ([THROUGHPUT-30B.md](../../THROUGHPUT-30B.md),
  [THROUGHPUT-20B.md](../../THROUGHPUT-20B.md)). The knob cannot cross
  that by rearranging the existing product.
- **Reference:** public Make/Modal audit of the shipping preset,
  14.637530 B/s complete-scalar median and 14.106673 B/s DP34 median
  ([SHARED-SIGMA.md](../../SHARED-SIGMA.md)). The matched control from
  this pair is the number the ratio column uses.

Unit: billions of completed scalar updates per second (`B/s`). Class:
engineering unless the ratio to the one-add floor moves.

## Acceptance, stated before the run

Promote the knob into the default RTX PRO 6000 preset only if the
candidate median is at least 1% above the matched control on the
complete-scalar bench, every control/candidate sample pair is faster,
identity gates match the requested build, and `[k]P = Q` is not in
scope here (this is a walk-rate bench). Otherwise leave the default at
0. A later DP34 audit can tighten that; it cannot loosen it.

Inadmissible: changing batch, workers, steps, launches, or the rest of
the preset; dropping identity checks; comparing against a different GPU
or a historical rate instead of the matched control; treating wall time
as the metric.

## Run

From `ecc2k130/`:

```sh
bash benchmarks/top-clmad/run.sh
```

The runner writes scratch logs under `/tmp/top-clmad-receipts` so Modal's
source snapshot does not see the tree change mid-build, then copies the
JSON receipts here. `summarize.py` cites those files and does not
recompute rates from a model.

## Result (2026-09-17)

| variant | median B/s | / control | / 23 B floor | correctness | class |
|---|---:|---:|---:|---|---|
| shipping product | 15.115792 | 1.000 | 0.657 | top clmad 0 | reference |
| `PACKED_TOP_CLMAD=1` | 12.859376 | 0.851 | 0.559 | top clmad 1 | engineering |

Acceptance failed. The default stays 0. Numbers from [summary.json](summary.json).
