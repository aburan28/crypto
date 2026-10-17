# A19: the huge-page filter with the folded build's 64 partitions

**What this is.** A timed exploration on v8. It ran on 2026-10-17 from 06:44 to 07:02 UTC on the
Emerald Rapids host R12 ran on. It is a record, not a round's measurement.

**Why.**
- **[A17](../A17-filter-huge-pages-20261017/README.md)'s huge-page filter** gains about 5% at 2^47.2 and nothing clear at 2^44.
- **The build is about 30% of cold time at 2^44.3 and 2^44.5.** Step timers in a diagnostic build
  ([`build-timing-diagnostic.patch`](build-timing-diagnostic.patch),
  [`partbits-sweep.tsv`](partbits-sweep.tsv)) split its cost per stored pair:
  - the addition and the key, about 8 ns;
  - the append into the partition streams, about 12–13 ns;
  - the place step (count, place, filter), 9.5–11.6 ns.
- **The append writes two streams for each of 256 partitions.** With 64 partitions it costs about 30% less, and the build about 8–10% less, in two processes per setting. Fewer partitions make the place step dearer.

  The table does not depend on the partition count: within a bucket, entries keep row order either way.

This candidate puts the two levers together, to see whether every target size clears a round's 1.03.

**The outcome: not at 2^44.5, so no round is declared.**

| curve | `log₂ r` | cold time, v8 / a19 | 95% interval |
|:--|--:|--:|:--|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | **1.067** | [1.010, 1.126] |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 1.014 | [0.978, 1.051] |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | **1.054** | [1.012, 1.097] |

The phases, as the geometric mean over 12 pairs of each pair's ratio of median repetitions
(`v8 / a19`; [`phases.sh`](phases.sh), [`phases.tsv`](phases.tsv)):

| curve | collection | build | selection |
|:--|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 1.045 | 1.096 | 1.091 |
| `icv1-f2m59-tm943548413-98844ecc` | 1.003 | 1.048 | 1.025 |
| `icv1-f2m61-t158598901-ab42b6c5` | 1.071 | 1.014 | 1.015 |

- **2^44.5's collection does not respond, as in A16 and A17.** It is the one size where
  three explorations in a row saw none of the filter's gain. The VM has no performance
  counters, so the cause is not measured.
- **A round accepts only if each target size's interval has its lower end above 1.03,**
  on the suite rows and on the holdouts (`icprog`'s `bound_reasons`). Here 2^44.5's lower end is 0.978,
  so a round on this candidate would very likely be rejected there.

## The candidate

`6c998fb0` is A17's commits plus one, on v8 (`23138d20`), frozen as
[`candidate.patch`](candidate.patch). Its one new line is the folded build's
`part_bits`, 6 instead of 8. Its tests were not run, because no round was declared.

## What ran

- **The rows:** `M1`'s two rows at the three target sizes.
- **The arms:** `v8` (`ic-r12-v7-23138d20`) and `a19` (`ic-a19-v8-6c998fb0`),
  hashed in [`binaries.sha256`](binaries.sha256).
- **The schedule:** six rounds with the order alternating, 72 processes, all completed and none
  contended ([`explore.sh`](explore.sh), [`explore.log`](explore.log)).
- **Refusals:** the busy gate refused a run 40 times before it started, and each refused run was asked for again.
- **Logarithms:** every pair recovered the same logarithm in both arms ([`pairs.tsv`](pairs.tsv)).

## Graph and report

- **The graph:** [`report/cold-ratio.svg`](report/cold-ratio.svg).
- **The report:** [`report/report.pdf`](report/report.pdf).
- **Copies:** in
  [`benchmarks/20261017-ic-a19-filter-pages-build-parts/`](../../../../benchmarks/20261017-ic-a19-filter-pages-build-parts/).

## Scope

One host and one core, the three target sizes, two rows a size.
