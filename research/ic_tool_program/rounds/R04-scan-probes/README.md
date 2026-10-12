# R04: where a scanned summand's time goes — results

**A stage diagnostic: accounting, no speedup, no new baseline.** R04 was
declared in [`PROTOCOL.md`](PROTOCOL.md) (#1128) with amendments 1 and 2
(#1146, #1150), before any run. It priced each stage of the collection
scan inside the scan, with time-stamp-counter probes compiled in only by
the `scan-probes` feature. Both arms were built from `9a48b389` (v0′ plus
the probes), and ran on `M1`'s 22 rows, three rounds, ABAB, isolated.

## The answer

**The batched subtraction holds the most scan time at the top sizes.**
The protocol names the leading stage at "the top four sizes,
`2^44.3`–`2^47.2`". The suite has three sizes in that range, and
`analyse.py` takes the four largest, which adds `2^39.0`. Both readings
give the same answer:

| curve | `log₂ r` | leading stage | its share of the scan |
|:--|--:|:--|--:|
| `icv1-f2m41-tm2308219-7f48b14a` | 39.0 | subtract | 37% |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | admitted | 42% |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | subtract | 41% |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | subtract | 40% |

Three of four, or two of three, name the subtraction.

**Read with a caveat the protocol states.** It trusts a size's shares
where the probes' overhead interval lies inside that size's A/A band.
That held at no size (below), so every share here is reported with the
overhead beside it, as the protocol directs, and none is "trusted".

## Every stage, every size

Each cell is the stage's median share of the scan's cycles, with its
nanoseconds per scanned summand at the measured counter rate. The scan
is `subtract + key + filter + admitted`; `trial` is the rest of each
trial (the walk step and the bookkeeping), outside the scan. Six
processes a size, every one clean.

| curve | `log₂ r` | ns a summand | subtract | key | filter | admitted | admitted keys a summand | trial, ns a summand |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| `icv1-f2m19-tm797-9c54981b` | 18.0 | 507.8 | 14% (74.3) | 9% (46.5) | 2% (12.2) | 74% (375.4) | 0.331 | 78.7 |
| `icv1-f2m23-tm5197-1f85e9e1` | 22.0 | 198.3 | 30% (58.1) | 18% (36.0) | 4% (9.1) | 47% (93.9) | 0.168 | 39.7 |
| `icv1-f2m45-tm6236725-40939294` | 24.8 | 156.7 | 33% (51.6) | 14% (22.7) | 5% (7.2) | 48% (75.1) | 0.139 | 41.1 |
| `icv1-f2m37-tm534059-32aad96b` | 27.8 | 76.5 | 60% (46.0) | 17% (12.8) | 8% (6.5) | 14% (11.3) | 0.135 | 24.1 |
| `icv1-f2m43-tm998717-e2e742b0` | 32.1 | 58.0 | 48% (27.8) | 27% (15.8) | 10% (5.7) | 15% (8.5) | 0.153 | 13.0 |
| `icv1-f2m47-t22705043-f4e44623` | 36.6 | 39.2 | 39% (15.1) | 27% (10.8) | 13% (5.1) | 21% (8.2) | 0.131 | 4.3 |
| `icv1-f2m57-tm747311035-c1f545af` | 38.0 | 44.1 | 31% (13.5) | 28% (12.5) | 13% (5.7) | 28% (12.3) | 0.199 | 3.5 |
| `icv1-f2m41-tm2308219-7f48b14a` | 39.0 | 31.5 | 37% (11.6) | 30% (9.2) | 14% (4.3) | 20% (6.3) | 0.118 | 2.8 |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 54.9 | 24% (13.1) | 20% (10.8) | 14% (7.5) | 42% (23.3) | 0.198 | 1.3 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 78.6 | 41% (32.5) | 15% (12.3) | 10% (8.0) | 33% (26.3) | 0.206 | 1.5 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 82.8 | 40% (33.3) | 16% (13.2) | 11% (9.0) | 32% (27.0) | 0.170 | 1.3 |

What the table shows:
- **The two wide-tail sizes pay for the scalar subtraction.** At
  `2^44.5` and `2^47.2` the subtraction costs 32.5 and 33.3 ns a
  summand, against 13.1 at `2^44.3`, where the 8-lane kernel runs. That
  is the premium R02's kernel removes; R02b, declared in #1157, re-tests
  it.
- **The presence filter admits a key in five to eight** from `2^36.6`
  up: 0.118–0.206 admitted keys a summand. At `2^44.3` the admitted
  stage is the largest, at 23.3 ns a summand.
- **The key and the filter are steady:** 9–16 ns and 4–9 ns a summand
  across the top six sizes.
- **Below `2^28` the admitted stage dominates,** at 47–74% of the scan.
  The bases there are small, and the scan is short.

## The probes' overhead

Default over probes in cold time, so a ratio below 1 means the probes
cost time. Six pairs a size; R01's A/A bands are from ten.

| curve | default over probes | R01's A/A band | inside |
|:--|--:|--:|:--|
| `icv1-f2m19-tm797-9c54981b` | 0.976 [0.874, 1.090] | [0.930, 1.118] | no |
| `icv1-f2m23-tm5197-1f85e9e1` | 0.960 [0.882, 1.044] | [0.984, 1.059] | no |
| `icv1-f2m45-tm6236725-40939294` | 1.040 [0.962, 1.124] | [0.916, 1.060] | no |
| `icv1-f2m37-tm534059-32aad96b` | 0.974 [0.872, 1.088] | [0.903, 1.055] | no |
| `icv1-f2m43-tm998717-e2e742b0` | 0.961 [0.935, 0.987] | [0.936, 1.129] | no |
| `icv1-f2m47-t22705043-f4e44623` | 0.962 [0.891, 1.038] | [0.922, 1.127] | no |
| `icv1-f2m57-tm747311035-c1f545af` | 1.049 [0.881, 1.249] | [0.935, 1.095] | no |
| `icv1-f2m41-tm2308219-7f48b14a` | 0.975 [0.784, 1.214] | [1.005, 1.149] | no |
| `icv1-f2m53-tm56619371-dac20a85` | 1.034 [0.952, 1.124] | [0.951, 1.088] | no |
| `icv1-f2m59-tm943548413-98844ecc` | 1.025 [0.967, 1.087] | [0.934, 1.049] | no |
| `icv1-f2m61-t158598901-ab42b6c5` | 0.969 [0.918, 1.023] | [0.926, 1.049] | no |

- **No interval lies inside its band.** With six pairs, each overhead
  interval is as wide as the band it must fit inside, or wider. The
  protocol's three rounds could not have met its own test at most sizes.
  That is a design fault of the protocol, recorded here; the rule is
  applied as written.
- **Ten of eleven intervals contain 1.** At `2^32.1` the probes cost
  3.9% [1.3%, 6.5%]. The geometric means run 0.96–1.05.
- **So the shares are reported, not trusted.** The ranking at the top
  sizes rests on gaps of 7–18 points of share, much larger than a few
  per cent of overhead could move. That is a reading, not the protocol's
  test.

## The pin and the accounting

- **The pin held on all 22 rows.** The probe arm's outputs equal the
  default arm's: the counts, both arms' scalars, rho's counts and the
  verification flags. The probe arm reported `scan_probes`; the default
  arm reported nothing.
- **138 processes ran:** 132 planned, six of them contended and retried.
  None failed. The contended runs are kept in the archive and enter no
  figure.
- **The host matched R01's,** so the A/A bands are R01's.

## A check after the run: the filter's admitted keys

This check was not declared. It was made after the analysis, from the
same archive. [`filter_check.py`](filter_check.py) predicts the presence
filter's pass rate from how the table sizes it, and compares that with
the admitted fraction R04 measured.

The filter is one bit per stored pair under one hash. Its size is the
next power of two above four bits per stored pair, so 4–8 bits a pair.
A key that is not in the table passes with the fraction of bits set,
`1 − (1 − 2^{−bits})^{stored}`.

| curve | `log₂ r` | stored pairs | filter bits a pair | predicted pass rate | measured admitted a summand | measured over predicted | admitted stage, ns an admitted key |
|:--|--:|--:|--:|--:|--:|--:|--:|
| `icv1-f2m47-t22705043-f4e44623` | 36.6 | 150,024 | 6.99 | 0.1333 | 0.1306 | 0.980 | 62.5 |
| `icv1-f2m57-tm747311035-c1f545af` | 38.0 | 237,120 | 4.42 | 0.2024 | 0.1994 | 0.985 | 61.8 |
| `icv1-f2m41-tm2308219-7f48b14a` | 39.0 | 265,680 | 7.89 | 0.1190 | 0.1176 | 0.988 | 53.8 |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 3,707,880 | 4.53 | 0.1983 | 0.1976 | 0.996 | 117.9 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 3,881,728 | 4.32 | 0.2066 | 0.2056 | 0.995 | 127.7 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 6,265,920 | 5.36 | 0.1703 | 0.1697 | 0.996 | 158.8 |

- **From `2^36.6` up, the prediction holds to within 2%.** True hits
  are one summand in 7,200 to 221,000 there. So nearly every admitted
  key is a false positive of the filter, and the admitted stage, 20–42%
  of the scan at these sizes, is spent on them.
- **An admitted key costs 54–63 ns up to `2^39.0`, and 118–159 ns
  above.** The three largest tables hold 3.7–6.3 M pairs.
- **Below `2^36.6` the model does not hold.** The measured fraction is
  0.78–0.94 of the prediction at four sizes, and 2.15 times it at
  `2^18.0`, where a tenth of the summands are true hits. Those tables
  hold 1,368–11,696 pairs. `filter_check.json` has all eleven sizes.

## What follows

- **The decision, by the protocol: the subtraction.** The next scan
  round targets it. That round is R02b (#1157), already declared, which
  re-tests R02's 8-lane kernel at the two sizes where the subtraction
  runs scalar.
- **A second lever, seen here, not decided here: a sharper filter.** At
  the three largest sizes the admitted stage costs 23.3–27.0 ns a
  summand, nearly all of it on false positives. A filter that passes
  fewer absent keys would cut that, if its own probe cost little more.
  No such filter was built or measured here. It is a candidate for a
  later scan round, if plan §11 leaves the scan open after R02b.

## Files

- [`PROTOCOL.md`](PROTOCOL.md): the declaration and amendments 1–2.
- [`run.py`](run.py) and [`analyse.py`](analyse.py): the runner and the
  analysis, as merged before the run.
- [`analysis.json`](analysis.json): every figure above except the
  filter check's. It reproduces byte for byte from the archive.
- [`filter_check.py`](filter_check.py) and
  [`filter_check.json`](filter_check.json): the check after the run,
  which also reproduces byte for byte from the archive.
- `runs.tar.xz` and `runs.tar.xz.sha256`: the run tree, every attempt
  kept, contended ones included.
