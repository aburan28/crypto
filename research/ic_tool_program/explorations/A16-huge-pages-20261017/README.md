# A16: the folded table's scan arrays collapsed onto 2 MiB pages

**What this is.** This is a timed exploration on v8. It ran on 2026-10-17 from 04:24
to 04:54 UTC on the Emerald Rapids host R12 ran on. It is a record, not a
round's measurement. None of its runs is pooled with any round's.

**Why.** The stage probes at v8
([record](../scan-stages-v8-20261017/README.md)) show the scan's filter
stage as the largest. [A15](../A15-filter-slots-20261017/README.md) shows
that the stage waits on memory.
- **The host's THP setting is `madvise`, and the tool never advised:** the
  4–8 MiB filter and the pair table sat on 4 KiB pages.
- **The smoke records** ([`smoke/`](smoke/)): single processes at
  `icv1-f2m61-t158598901-ab42b6c5` (2^47.2), isolated, before any timed
  run. These are observations, not measurements.
  - **Every allocation on huge pages** (glibc's `glibc.malloc.hugetlb=1` on
    v8): collection 1395 → 1177 ms, but the build 213 → 565 ms, with minor
    faults 88,040 → 6,009.
  - **A16's collapse of every scan array:** collection 1332 → 1209 ms, the
    build 207 → 264 ms, and cold time 1626 → 1558 ms.

**The outcome: the filter alone is the lever, and this form of it falls
short.**
- **Collapsing the filter alone:** cold time is **4.6% faster at 2^47.2**,
  with the interval wholly above 1. At 2^44.3 and 2^44.5 it is
  inconclusive.
- **Collapsing every array:** no gain anywhere. The build pays for copying
  and compacting the rests.
- **Why the filter-only gain is uneven:** `MADV_COLLAPSE` covers only the
  whole 2 MiB pages inside an allocation. The allocator's block starts just
  past a boundary, so only 3 of an 8 MiB filter's 4 pages collapse, and 1 of a
  4 MiB filter's 2.
- **No round is declared on this form.**
  [A17](../A17-filter-huge-pages-20261017/README.md) puts the whole filter on
  huge pages from its first touch.

## The candidate

`e7136eb8`, two commits on v8 (`23138d20`), frozen as
[`huge-pages.patch`](huge-pages.patch):
- **`on_huge_pages`** runs after the folded build. It advises and collapses
  (`MADV_HUGEPAGE`, then `MADV_COLLAPSE`) the whole 2 MiB pages inside:
  - the filter;
  - the bucket starts;
  - the stored keys.
- **`KIC_HUGE_PAGES`** chooses the arm: `0` collapses nothing, `filter`
  collapses the filter alone, and unset collapses every array.
- **The values are unchanged.**

The tests were not run on its tree, because no round was declared.

## What ran

- **The rows:** `M1`'s two rows at the three target sizes.
- **The arms:** `off`, `filter` and `all`, all on one binary
  (`ic-a16-v8-e7136eb8`, [`binaries.sha256`](binaries.sha256)).
- **The schedule:** six rounds with the order rotating, 108 processes, all
  completed. Each process ran `ic price --single-target` with the suite's rho
  seed and `RAYON_NUM_THREADS=1`, isolated on CPU 2
  ([`explore.sh`](explore.sh), [`explore.log`](explore.log)).
- **Refusals:** the busy gate refused a run 61 times before it started, and each
  refused run was asked for again. Some of these came while the session merged
  main into the PR branch (04:25–04:32 UTC). The heavy steps of that merge ran
  under the benchmark lock.
- **Contended runs:** one run (`off`, 2^44.3, round 1) was marked contended and is left out.
- **Logarithms:** every pair recovered the same logarithm in all three arms.

## Results

Cold time, `off` over each arm. Above 1, the arm is faster. The intervals are
95% t intervals of the geometric mean over the clean pairs
([`summarise.sh`](summarise.sh), `summary-filter.tsv`, `summary-all.tsv`).

| curve | `log₂ r` | `filter` | 95% | `all` | 95% |
|:--|--:|--:|:--|--:|:--|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 1.017 | [0.968, 1.068] | 0.970 | [0.910, 1.034] |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 0.981 | [0.936, 1.029] | 0.979 | [0.941, 1.020] |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | **1.046** | [1.026, 1.067] | 1.012 | [0.976, 1.049] |

The phases are the geometric mean of the ratio of median repetitions over all 12 pairs,
`off` over each arm ([`phases.sh`](phases.sh), `phases-*.tsv`):

| curve | collection, `filter` | build, `filter` | collection, `all` | build, `all` |
|:--|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 1.027 | 1.016 | 1.012 | 0.963 |
| `icv1-f2m59-tm943548413-98844ecc` | 1.018 | 0.927 | 1.024 | 0.899 |
| `icv1-f2m61-t158598901-ab42b6c5` | 1.062 | 0.995 | 1.061 | 0.774 |

## Graph and report

- **Graph:** [`report/cold-ratio.svg`](report/cold-ratio.svg).
- **Report:** [`report/report.pdf`](report/report.pdf).
- **Copies:** in
  [`benchmarks/20261017-ic-a16-huge-pages/`](../../../../benchmarks/20261017-ic-a16-huge-pages/).

## Scope

This covers one host and one core, the three target sizes, and two rows a size. Huge pages
depend on the host:
- **The THP setting:** here `madvise`, with defrag `madvise`.
- **Free memory for 2 MiB pages:** a host with `never` gets no change.
