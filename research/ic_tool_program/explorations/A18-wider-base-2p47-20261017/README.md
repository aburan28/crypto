# A18: a wider factor base at 2^47.2, with and without the huge-page filter

**What this is.** This is a timed exploration on v8. It ran on 2026-10-17 from 06:00 to 06:13 UTC on the
Emerald Rapids host R12 ran on. It is a record, not a round's measurement. It changes the recipe, which the plan's A12 classes as an algorithm change, so no round would take it without a declaration of its own.

**Why.** At v8 the 2^47.2 set-up is unbalanced: collection takes 82% of
cold time and the build 13%
([stage probes](../scan-stages-v8-20261017/README.md)). A simple model (build ∝ |F|²,
collection ∝ r/|F|) puts the cheapest base about 1.47× wider, at about 470 orbits rather than 320.
The model predicts about 12% off cold time. It is a model, not a measurement. Huge pages for the filter
([A17](../A17-filter-huge-pages-20261017/README.md)) matter more on a larger table, so the two
levers are crossed here.

**The outcome: the wider base buys little, because the build grows faster than the table.**
- **On v8:** 400 and 470 orbits make cold time 1.020 [1.002, 1.037] and 1.023 [0.994, 1.054] faster.
- **On A17:** no gain beyond A17's own.
- **Why:** from 320 to 470 orbits, collection falls 27% (1316 → 959 ms), but the
  build rises 2.28× (217 → 495 ms) for 2.17× the stored pairs (6.27M → 13.62M).
  Selection and the other phases grow too.
- **A17 replicates:** 1.065 [1.040, 1.090] at 2^47.2 on the suite's recipe, a second independent
  measurement after A17's 1.051.
- **The build, about 31–35 ns per stored pair at both sizes, is the next lever.** A cheaper build
  helps directly at 2^44.3 and 2^44.5, where it is 30% of cold time, and it is
  what would make a wider base pay at 2^47.2.

## What ran

- **The rows:** `M1`'s two rows at `icv1-f2m61-t158598901-ab42b6c5` (2^47.2).
- **The recipes** ([`params/`](params/)): the suite's file with only the factor base's
  `points` changed:
  - 39,040 points (320 signed-Frobenius orbits of 122). This file equals the suite's, apart from its name.
  - 48,800 points (400 orbits).
  - 57,340 points (470 orbits).
- **The binaries:** `v8` (`ic-r12-v7-23138d20`) and `a17` (`ic-a17-v8-fb4284d8`),
  hashed in [`binaries.sha256`](binaries.sha256).
- **The schedule:** six arms, four rounds, the order rotating, 48 processes, all completed and
  none contended.
  - Each process ran `ic price --single-target` with the suite's rho seed and
    `RAYON_NUM_THREADS=1`, isolated on CPU 2
    ([`explore.sh`](explore.sh), [`explore.log`](explore.log)).
  - The busy gate refused a run 21 times before it started, and each refused run was asked for
    again.
- **Logarithms:** every arm recovered each row's logarithm, the suite's
  ([`pairs-*.tsv`](.)).
- **The smoke records** ([`smoke/`](smoke/)): one process each at 320 and 470 orbits, before the timed runs.

## Results

Cold time, paired against `v8-c320` (above 1, the arm is faster). The intervals are 95% t intervals of the
geometric mean over 8 pairs ([`summarise.sh`](summarise.sh), `summary-*.tsv`):

| arm | cold time | 95% interval |
|:--|--:|:--|
| `v8-c400` | 1.020 | [1.002, 1.037] |
| `v8-c470` | 1.023 | [0.994, 1.054] |
| `a17-c320` | **1.065** | [1.040, 1.090] |
| `a17-c400` | 1.055 | [1.019, 1.093] |
| `a17-c470` | 1.058 | [1.042, 1.074] |

Each phase is the geometric mean over an arm's 8 processes of the median repetition, in ms
([`phases.sh`](phases.sh), [`phases.tsv`](phases.tsv)):

| arm | cold | collection | build | rest |
|:--|--:|--:|--:|--:|
| `v8-c320` | 1629.1 | 1315.7 | 216.8 | 96.6 |
| `v8-c400` | 1597.6 | 1123.2 | 365.0 | 109.4 |
| `v8-c470` | 1591.9 | 959.1 | 495.0 | 137.7 |
| `a17-c320` | 1530.1 | 1217.9 | 226.3 | 85.9 |
| `a17-c400` | 1543.8 | 1006.2 | 420.1 | 117.6 |
| `a17-c470` | 1540.1 | 866.6 | 546.0 | 127.5 |

## Graph and report

- **Graph:** [`report/cold-ratio.svg`](report/cold-ratio.svg).
- **Report:** [`report/report.pdf`](report/report.pdf).
- **Copies:** in [`benchmarks/20261017-ic-a18-wider-base/`](../../../../benchmarks/20261017-ic-a18-wider-base/).

## Scope

This covers one host and one core, one size, and two rows. Only the base's width changed;
the collection window and the trial budget are the suite's.
