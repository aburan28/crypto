# P-256 Dickson quadratic-S3 round 3: result

Date run: 2026-10-04

The direct quadraticization hypothesis is **falsified**.  All four variable
layouts retained F4 solving degree 4 at depths 4 and 5, on both terminal-zero
and the round-2 winning nonzero fibres.  Lowering maximum input degree from 4
to 2 did not lower solving degree.

The quadratic systems did reduce field operations in several cells.  At depth
5, grouped variables cut the median from 36.7 million to 18.3 million field
multiplications on terminal zero, and from 31.6 million to 16.8 million on
terminal 369.  This is solver engineering, not a regularity result: matrices
were wider and the degree ratio stayed 1.00.

All 60 systems completed, all were certified, and every verdict agreed with
exhaustive signed-point addition.  There were no timeouts or pairs above the
degree bound.

| depth | terminal | arm/layout | targets | correct | complete | input degree | solving degree | columns median | field ops median | degree / quartic S3 | classification |
|---:|---:|:--|---:|:--:|:--:|---:|:--|---:|---:|---:|:--|
| 4 | 0 | quartic S3 | 3 | yes | yes | 4 | 4/4/4 | 309 | 1,183,930 | 1.00 | reference |
| 4 | 0 | grouped | 3 | yes | yes | 2 | 4/4/4 | 360 | 1,150,381 | 1.00 | relabelling |
| 4 | 0 | aux-first | 3 | yes | yes | 2 | 4/4/4 | 357 | 1,095,382 | 1.00 | engineering only |
| 4 | 0 | level-interleaved | 3 | yes | yes | 2 | 4/4/4 | 367 | 1,083,807 | 1.00 | engineering only |
| 4 | 0 | reverse-blocks | 3 | yes | yes | 2 | 4/4/4 | 390 | 900,802 | 1.00 | engineering only |
| 4 | 782 | quartic S3 | 3 | yes | yes | 4 | 4/4/4 | 309 | 1,105,052 | 1.00 | reference |
| 4 | 782 | grouped | 3 | yes | yes | 2 | 4/4/4 | 360 | 1,052,502 | 1.00 | engineering only |
| 4 | 782 | aux-first | 3 | yes | yes | 2 | 4/4/4 | 357 | 1,048,909 | 1.00 | engineering only |
| 4 | 782 | level-interleaved | 3 | yes | yes | 2 | 4/4/4 | 367 | 985,811 | 1.00 | engineering only |
| 4 | 782 | reverse-blocks | 3 | yes | yes | 2 | 4/4/4 | 390 | 836,223 | 1.00 | engineering only |
| 5 | 0 | quartic S3 | 3 | yes | yes | 4 | 4/4/4 | 647 | 36,725,681 | 1.00 | reference |
| 5 | 0 | grouped | 3 | yes | yes | 2 | 4/4/4 | 714 | 18,309,849 | 1.00 | engineering only |
| 5 | 0 | aux-first | 3 | yes | yes | 2 | 4/4/4 | 737 | 28,698,701 | 1.00 | engineering only |
| 5 | 0 | level-interleaved | 3 | yes | yes | 2 | 4/4/4 | 722 | 16,792,103 | 1.00 | engineering only |
| 5 | 0 | reverse-blocks | 3 | yes | yes | 2 | 4/4/4 | 804 | 14,237,319 | 1.00 | engineering only |
| 5 | 369 | quartic S3 | 3 | yes | yes | 4 | 4/4/4 | 642 | 31,566,344 | 1.00 | reference |
| 5 | 369 | grouped | 3 | yes | yes | 2 | 4/4/4 | 706 | 16,778,687 | 1.00 | engineering only |
| 5 | 369 | aux-first | 3 | yes | yes | 2 | 4/4/4 | 720 | 25,099,432 | 1.00 | engineering only |
| 5 | 369 | level-interleaved | 3 | yes | yes | 2 | 4/4/4 | 719 | 15,412,143 | 1.00 | engineering only |
| 5 | 369 | reverse-blocks | 3 | yes | yes | 2 | 4/4/4 | 802 | 13,397,495 | 1.00 | engineering only |

The frozen winner prioritizes matrix width after degree, so `aux-first` wins
both depth-4 cells, while `grouped` wins both depth-5 cells.  If field
multiplications were prioritized instead, `reverse-blocks` would win, but
changing that rule cannot create a degree reduction.

## Reproduction and evidence

```bash
cargo run --release --bin p256_dickson_quadratic_s3 -- \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_quadratic_s3_round3_20261004/degree-result.json
```

`degree-result.json` is 34,746 bytes with SHA-256
`380afc1d1e6c6521e818024b54cc1760107a21b5029c07a29b98de334137ff15`.

The next iteration should branch on the terminal one or two levels of the
Dickson chain and run the quadratic system on every resulting component.  That
is the first tested change that can remove quadratic constraints from each
component rather than merely rename them.  It must charge all branch systems,
report the maximum deciding degree across them, and keep total field
operations beside the degree.

No P-256 Gröbner degree, ECDLP solve, or end-to-end speedup is claimed.  The
best verified P-256 inventory remains round 2's `FB1h2f8621cda105`.
