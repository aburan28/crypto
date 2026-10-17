# A14: R12's run prefetch gated on the table's size

**What this is.** This is a timed exploration on v8, R12's accepted baseline. It ran on 2026-10-17 from
03:10 to 03:28 UTC on the Emerald Rapids host R12 ran on, before any round
was declared on the gate. It is a record, not a round's measurement. None of
its runs is pooled with any round's.

**Why.** R12 made cold time 1.055–1.063× faster at the three target sizes
([R12](../../rounds/R12-admitted-prefetch/README.md)). From 2^24.8 to 2^38.0
its point estimates read 0.956–0.987, and at 2^27.8 its interval lay
below 1. These differences were inside the declared A/A rule. R12's write-up put that down to the extra
pass where the table is small. The candidate here skips the pass below
2^20 stored pairs (16 MiB). The suite's tables hold 0.15M–0.27M pairs at
2^36.6–2^39.0 and 3.7M–6.3M at 2^44.3–2^47.2, so the target sizes run v8's
code unchanged.

**The outcome: no round is declared.** Below the gate, the gated arm is no
faster than v8. Cold time, v8 over the gate, reads 0.992 at 2^27.8 and 0.961
at 2^38.0, and both intervals include 1. So R12's small-size readings are
not the pass's cost, or they are too small to measure with 12 pairs a size. The gate
cannot change the target sizes, so a round on it could not meet the 1.03
acceptance line. v8 stays the baseline and the gate is not landed.

## The candidate

`ebfd2d2f`, one commit on v8 (`23138d20`), frozen here as
[`size-gate.patch`](size-gate.patch) (13 lines in
`koblitz_index_calculus.rs`):
- `RUN_PREFETCH_MIN_PAIRS = 1 << 20`.
- The second pass (`prefetch_run` over the admitted keys) runs only when
  the table holds at least that many pairs.
- The lookups, the matches and their order are unchanged.

The tests were not run on its tree: no round was declared, and the code is not
landed.

## What ran

- **The rows:** `M1`'s two rows at three sizes:
  - `icv1-f2m37-tm534059-32aad96b` (2^27.8), where R12's interval lay below 1;
  - `icv1-f2m57-tm747311035-c1f545af` (2^38.0), R12's lowest point estimate;
  - `icv1-f2m53-tm56619371-dac20a85` (2^44.3), above the gate.
- **The arms:** `v8` (`ic-r12-v7-23138d20`, the binary R12 measured) and
  `gate` (`ic-r13-v8-ebfd2d2f`). Both are hashed in
  [`binaries.sha256`](binaries.sha256).
- **The schedule:** six rounds with the order alternating. Every process ran
  `ic price --single-target` with the suite's rho seed and
  `RAYON_NUM_THREADS=1`, isolated on CPU 2: 72 processes, all completed
  ([`explore.sh`](explore.sh), [`explore.log`](explore.log)).
- **Refusals before a run started:** the machine-busy gate refused a run 201 times before it
  started (PSI avg10 above 5). The refused runs never started, and each was asked for again:
  - 185 refusals were skipped and picked up by the first passes;
  - after commit `35ad12a5`, the script waited 15 s and retried each refusal (16 retries).
  The refusals ran through the whole exploration (155 from 03:10 to 03:19, 45
  after 03:20). Most followed the previous run, while its own pressure
  decayed.
- **Other work during the runs:** from about 03:11 to 03:15 UTC the session
  merged main into the PR branch. During that time it ran the site's
  generators and checks, which do not take the benchmark lock, and
  `cargo check`, which did take it (`isolated_bench busy`). Runs that
  completed then carry the harness's contention check, like every other run.
  Nothing else excludes them.
- **Contended runs:** seven of the 72 were marked contended, all `v8`: two at
  2^27.8 and five at 2^38.0. They are left out of the intervals, and the
  all-pairs figures are given beside them.
- **Logarithms:** every pair recovered the same logarithm in both arms, all
  36 ([`pairs.tsv`](pairs.tsv)).

## Results

Cold time, `v8` over `gate` (above 1, the gate is faster), with 95% t
intervals of the geometric mean over the clean pairs
([`summarise.sh`](summarise.sh), [`summary.tsv`](summary.tsv)):

| curve | `log₂ r` | prefetch in `gate` | clean pairs | cold time | 95% interval | all pairs |
|:--|--:|:--|--:|--:|:--|--:|
| `icv1-f2m37-tm534059-32aad96b` | 27.8 | skipped | 10 | **0.992** | [0.964, 1.022] | 0.991 |
| `icv1-f2m57-tm747311035-c1f545af` | 38.0 | skipped | 7 | **0.961** | [0.854, 1.083] | 0.987 |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | runs | 12 | **1.011** | [0.978, 1.046] | 1.011 |

The phases, as the geometric mean of all 12 pairs' ratios of median repetitions,
`v8` over `gate` ([`phases.sh`](phases.sh), [`phases.tsv`](phases.tsv)):

| curve | collection | build | selection |
|:--|--:|--:|--:|
| `icv1-f2m37-tm534059-32aad96b` | 1.003 | 0.975 | 0.993 |
| `icv1-f2m57-tm747311035-c1f545af` | 1.029 | 0.998 | 1.001 |
| `icv1-f2m53-tm56619371-dac20a85` | 1.010 | 1.027 | 0.993 |

The median cold time over each arm's 12 processes, in ms:

| curve | `v8` | `gate` |
|:--|--:|--:|
| `icv1-f2m37-tm534059-32aad96b` | 3.09 | 3.09 |
| `icv1-f2m57-tm747311035-c1f545af` | 38.97 | 38.65 |
| `icv1-f2m53-tm56619371-dac20a85` | 406.37 | 398.84 |

- **At 2^38.0, skipping the pass makes collection 2.9% faster** (all pairs).
  Collection is about 17 ms of the run's 39 ms. That would be about 1% of
  cold time, which this exploration cannot resolve: its clean interval is ±11%.
- **At 2^27.8 collection is unchanged** (1.003), and the whole run takes 3.1 ms.
- **R12's 0.971 [0.950, 0.992] at 2^27.8 is not recovered here** (0.992
  [0.964, 1.022]).

## Graph and report

[`report/cold-ratio.svg`](report/cold-ratio.svg) shows each pair's ratio and
the clean intervals. [`report/report.pdf`](report/report.pdf) is the report,
and [`report/report.typ`](report/report.typ) is its source. Copies are in
[`benchmarks/20261017-ic-a14-prefetch-size-gate/`](../../../../benchmarks/20261017-ic-a14-prefetch-size-gate/).

## Scope

This covers one host and one core, three sizes and two rows a size. It shows that
skipping the pass recovers no measurable cold time at these sizes. It does
not show that R12's small-size readings were noise. R12's own A/A rule
already placed them inside its band.

## Files

- `runs.tar.xz` holds every run's `price.json` and isolation record.
- `SHA256SUMS` hashes the files here.
