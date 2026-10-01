# R04: where a scanned summand's time goes, priced inside the scan

**Declared 2026-10-01, before any R04 code.** This is Track A in
`research/notes/index-calculus/IC_TOOL_PROGRAM.md` (§8, A2 and A11). It
is a stage diagnostic and claims no speedup. Nothing below changes after
the first run, except by a dated amendment appended at the end.

## The question

**The collection scan dominates the set-up at the top sizes, and half of
a scanned summand's cost has no native price.**
- **The scan's share.** It is 65–89% of the set-up from `2^39` to
  `2^47.2` (plan §8).
- **The cost per summand.** At `n = 53` a scanned summand costs about
  40 ns (R01: 1.49 units).
- **What is priced.** The A11 exploration priced two stages outside the
  scan, warm:
  - the AVX-512 key: 5.6–10.2 ns a key;
  - the 8-lane batched subtraction: 8.3–8.4 ns a summand.

  Together they are about half.
- **What is not.** The rest has no native price:
  - the presence-filter probe;
  - the admitted keys' bucket probes and pair recovery;
  - the loop's bookkeeping.
- **Why the usual tools cannot answer.** Valgrind 3.22 hides AVX-512, so
  callgrind sees only the portable paths. The container has no hardware
  counters, so `perf` cannot split the time either.

This round prices each stage inside the scan, on the suite's own rows.
The next scan round's lever is chosen from its answer.

## The instrument

- **A cargo feature, `scan-probes`, off by default.** With it, the scan
  (`PairSumTable::witnesses_fast_scan`) and the collection loop
  (`collect_walked`) read the time-stamp counter at stage boundaries.
  They add the cycles to per-stage totals:

  | stage | what it holds |
  |:--|:--|
  | `subtract` | the batched subtraction (`add_many_lazy`, or `add_many` on an unfolded table) |
  | `key` | `keys_of`, block by block |
  | `filter` | the presence-filter loop, its prefetches included |
  | `admitted` | the admitted rests: `finish_lazy`, `pairs_for_key` and the sink |
  | `trial` | the rest of each trial in `collect_walked`: the walk step and recording the attempt, as the trial's whole time less the four above |

- **Counts beside the cycles:** trials, summands scanned, blocks,
  admitted keys and pairs returned.
- **The report.** With the feature on, `ic price` writes the totals as
  `scan_probes` in its report. The counter's rate is measured against the
  monotonic clock in the same process.
- **The default build.** Without the feature, the probes are not in the
  source the compiler sees, so the default build's code is unchanged.
- **The cost of a probe.** A counter read costs tens of cycles, and the
  probes read about seven per trial. A trial scans hundreds of summands
  at about 40 ns each, so the probes are expected to cost well under 1%.
  Measurement 3 checks that rather than assuming it.

## Arms

- **default**: the newest accepted baseline's commit plus R04's change,
  built without the feature.
- **probes**: the same commit built with `--features scan-probes`.

Both are built from one commit, recorded with their SHA-256 in the
round's manifest. R04 runs after R02, R03, B0 and B1, in the order of
the benchmark lock's queue.

## Rows

- **The rows.** `M1`'s 22 rows: two targets at each of the 11 sizes.
- **The runs.** Three rounds, ABAB, isolated, with the programme's
  runner: 132 processes.
- **The command.** `ic price --single-target`, as the suite runs it.

## Measurements

1. **The pin.** The probe arm's outputs must equal the default arm's on
   every row:
   - the counts;
   - both arms' scalars;
   - rho's counts;
   - the verification flags.
2. **The shares.** Per size and per stage:
   - the cycles per scanned summand;
   - the stage's share of the scan's cycles.

   Each is the median over the size's rows and rounds, and is also given
   in nanoseconds at the measured counter rate.
3. **The overhead.** The paired cold-time ratio, probes over default,
   per size, with its 95% interval.

## What it decides

- **The answer.** It names the stage that holds the most scan time at
  the top four sizes, `2^44.3`–`2^47.2`. The next scan round's lever
  targets that stage.
- **Whether to trust it.** The shares are trusted at a size where the
  overhead's interval lies inside that size's A/A band. Elsewhere they
  are reported with the overhead beside them.
- **Class: accounting.** It is a stage diagnostic: no speedup is
  claimed, and no baseline changes.

## Inadmissible

- Reading a share as a speedup, or a stage's price as the method's.
- Dropping a size.
- Counting a contended run.

## Cost

132 timed processes, about 40 minutes.
