# Order-1201 seven-arm cold batch: run record

The [committed run spec](RUN_SPEC_V2.json) froze `(p,r,h)=(1237,1201,1275)`,
the 64 public targets and checker-only scalars, source and binary hashes,
the seven arms, prefixes, three rounds, caps and the rotated order before
any measured child. The preceding [screen](SCREEN_RUN.md) and
[full-fiber audit](FIBER_RUN.md) passed independent Rust replay.

The first target-freeze attempt is preserved as
[`BATCH_1237_1201_V1.freeze.stderr`](BATCH_1237_1201_V1.freeze.stderr),
exit 101: `j_pair()` did not include the newly selected prime. It failed
before writing target files. The corrected
[V2 target freeze](BATCH_1237_1201_V2.freeze.json) exited zero; its
[public file](BATCH_1237_1201_V2.public.json) contains no answer scalars,
and the [checker file](BATCH_1237_1201_V2.checker.json) holds them separately.
The unused V1 batch script is retained for provenance; V2 is the executed
script.

Exact child command, from the checkout root:

```sh
zsh research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/order1000/run_order1201_v2.zsh
```

The script verified all hashes in `RUN_SPEC_V2.json` before launching.
The [84 child records](BATCH_1237_1201_V2.run.jsonl),
[meter records](BATCH_1237_1201_V2.meter.jsonl), and
[`raw-v2/`](raw-v2/) hold the exact stdout, stderr, meter and exit
for every arm/round/prefix. All **84/84** child exits were zero,
and the batch exit was zero. The
[independent replay](BATCH_1237_1201_V2.replay.json) exited zero,
recomputed **1,785/1,785** returned scalar identities and
**12,486** relation witnesses, checked every row and map identity,
and accepted all 84 child keys. The run-file SHA-256 is
`d314764e87988b2dc12675d199612f50abc538d1782761c562c2c74be5e34d89`;
the replay receipt SHA-256 is
`54e782e4df58cfb117ad87816e9a1a2dfd2593aa210d01526ca8623dc15ee41d`.

The native [five-order analyzer](../../../../examples/different_j_trend_five_summary.rs)
read all five run files and their matching replay certificates. Its
[V1 output](TREND_FIVE_V1.json) and figure remain as the first complete
calculation. V2 only moved one overlapping figure label; the raw
counts and conclusions did not change. The [V2 JSON](TREND_FIVE_V2.json),
[SVG](TREND_FIVE_V2.svg), [PNG](TREND_FIVE_V2.png),
[run PDF](BATCH_RUN.pdf), and [results PDF](RESULTS.pdf) are the current
visual record. The [results report](RESULTS.md) gives the one-unit table.

Host: Darwin arm64 shared L0. The child CPU and wall meters describe this
run only; they are not an isolated wall admission or an independent-host
performance replay. Complete portable `S` remains unpriced under
[`PRICE_PROTOCOL_V2.md`](../PRICE_PROTOCOL_V2.md).

![Five-order counted group-operation trend](TREND_FIVE_V2.svg)
