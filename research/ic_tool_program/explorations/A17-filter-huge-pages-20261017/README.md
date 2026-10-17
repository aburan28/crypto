# A17: the presence filter on huge pages from its first touch

**What this is.** This is a timed exploration on v8. It ran on 2026-10-17 on the Emerald Rapids host R12 ran on:
- the timed runs, from 05:08 to 05:35 UTC;
- a stage diagnostic after them.

It is a record, not a round's measurement. None of its runs is pooled with any round's.

**Why.** [A16](../A16-huge-pages-20261017/README.md) found that collapsing only the
filter onto 2 MiB pages gains 4.6% at 2^47.2, but it collapsed only the whole
pages inside an allocation that starts off a 2 MiB boundary. This candidate
puts the whole filter on huge pages from its first touch.

**The outcome: a gain at 2^47.2 only, so no round is declared.**
- **2^47.2:** cold time is 1.051 [1.029, 1.075] faster than v8, and collection 7.6% faster.
- **2^44.3:** 1.047 [0.999, 1.098].
- **2^44.5:** 0.976 [0.944, 1.009].

The stage probes put the effect in the filter stage at 2^47.2 (6.43 → 5.32 ns), with no consistent change at 2^44.3 or 2^44.5. A round needs every target size above 1.03 on the suite rows and the holdouts, so it would very likely be rejected at 2^44.5.

## The candidate

`fb4284d8` is two commits on v8 (`23138d20`), frozen as
[`filter-huge-pages.patch`](filter-huge-pages.patch). It adds `FilterWords`, the
filter's storage, which dereferences to `[u64]`:
- **For a filter of at least 2 MiB:** a fresh anonymous mapping, trimmed to
  start on a 2 MiB boundary and advised `MADV_HUGEPAGE` before the words are
  copied in. The first touch then faults in whole huge pages.
  - **Why not the heap:** the first commit used an aligned heap allocation, and a diagnostic run
    showed it got one huge page, not two. The allocator handed back
    memory it had already touched, and advice does not convert existing pages.
- **`KIC_FILTER_HUGE=0`** keeps the copy on the heap: a same-binary control.
- **The values are unchanged.** Each table is built as in v8, and its filter is then
  copied once.

**Diagnostics** (unisolated, `AnonHugePages` sampled from
`/proc/<pid>/smaps_rollup`):
- the whole filter sat on huge pages, 4,096 kB at 2^44.3 and 2^44.5 and
  8,192 kB at 2^47.2;
- v8 had none.

The tests were not run on its tree, because no round was declared.

## What ran

- **Timed** ([`explore.sh`](explore.sh), [`explore.log`](explore.log)):
  - **The rows:** `M1`'s two rows at the three target sizes.
  - **The arms:**
    - `v8` (`ic-r12-v7-23138d20`);
    - `a17` (`ic-a17-v8-fb4284d8`);
    - `heap`, the same binary with `KIC_FILTER_HUGE=0`.
  - **The schedule:** six rounds with the order rotating, 108 processes, all completed.
  - **Refusals:** the busy gate refused a run 55 times before it started, and each refused run was asked for again.
  - **Contended runs:** four were marked contended and are left out:
    - `v8` at 2^44.5, round 1;
    - `v8` at 2^44.3, round 5;
    - `heap` at 2^47.2, round 4;
    - `heap` at 2^44.3, round 1.
  - **Logarithms:** every pair recovered the same logarithm in all three arms.
- **Stage probes** ([`probe.sh`](probe.sh), [`probe.log`](probe.log)):
  - **The binaries:** R04's probes on v8 (`786c06f7`) and on `a17` (`a610cdab`).
  - **The schedule:** interleaved, four rounds, 48 processes, all completed.
  - **Refusals:** 28.
  - **Contended runs:** two, both at 2^44.3, and they are left out.

## Results

Cold time, with 95% t intervals of the geometric mean over the clean pairs
([`summarise.sh`](summarise.sh), `summary-*.tsv`):

| curve | `log₂ r` | v8 / a17 | 95% | v8 / heap | heap / a17 |
|:--|--:|--:|:--|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 1.047 | [0.999, 1.098] | 1.014 | 1.023 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 0.976 | [0.944, 1.009] | 0.981 | 0.994 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | **1.051** | [1.029, 1.075] | 0.991 | **1.054** [1.026, 1.082] |

The phases, as the geometric mean of the 12 pairs' ratios of median repetitions
([`phases.sh`](phases.sh), `phases-*.tsv`):

| curve | collection, v8 / a17 | build, v8 / a17 | collection, heap / a17 |
|:--|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 1.067 | 1.064 | 1.041 |
| `icv1-f2m59-tm943548413-98844ecc` | 1.003 | 0.966 | 1.014 |
| `icv1-f2m61-t158598901-ab42b6c5` | 1.076 | 0.940 | 1.081 |

The stages are in ns per scanned summand, each the mean over the clean processes
([`stages.sh`](stages.sh), `stages-v8.tsv`, `stages-a17.tsv`):

| curve | filter, v8 | filter, a17 | scan, v8 | scan, a17 |
|:--|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 4.75 | 4.78 | 15.69 | 16.34 |
| `icv1-f2m59-tm943548413-98844ecc` | 4.92 | 4.56 | 16.48 | 16.08 |
| `icv1-f2m61-t158598901-ab42b6c5` | 6.43 | **5.32** | 17.37 | **16.50** |

## What it shows

- **The copy costs nothing measurable** (`v8 / heap`). The gain at 2^47.2 is
  the pages' (`heap / a17`, 1.054).
- **The gain follows the filter's size.**
  - At 2^47.2 the filter is 8 MiB, 2,048 base pages: as many as this core's
    second-level TLB holds, and the TLB also serves everything else. There, huge
    pages cut the filter stage by 1.1 ns.
  - At 2^44.3 and 2^44.5 the filter is 4 MiB, 1,024 base pages, and the
    probes show no consistent change.

  This explanation is a hypothesis. The VM exposes no performance counters, so TLB misses
  were not counted.
- **The cold-time readings at 2^44.3 (1.047) and 2^44.5 (0.976) disagree with each
  other** and with the probes there. They are read as noise around a small effect.
- **The lever should grow with the problem.** The filter doubles with the
  table, so beyond 2^47.2 a larger share of filter reads would miss the TLB.
  That is where the distance to rho is largest (11.51× at 2^47.2, rule v7).

## Graph and report

- **The graph:** [`report/cold-ratio.svg`](report/cold-ratio.svg).
- **The report:** [`report/report.pdf`](report/report.pdf), with its source [`report/report.typ`](report/report.typ).
- **Copies:** in [`benchmarks/20261017-ic-a17-filter-huge-pages/`](../../../../benchmarks/20261017-ic-a17-filter-huge-pages/).

## Scope

One host and one core, the three target sizes, two rows a size. Huge pages
depend on the host's THP setting (`madvise` here) and on free 2 MiB pages.
