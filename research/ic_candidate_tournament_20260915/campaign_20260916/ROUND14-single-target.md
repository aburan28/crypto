# Round 0014 pre-registration: the worker as one optimisation unit

Written before any measured stage of this round was run.

## Where the job's time sat after round 0013

Round 0013 left the IC arm at 0.7697 of rho's instructions and 0.9098 of
its native time. To find what the native clock is spending, this round's
development work instrumented a copy of the round-0013 winner with a
monotonic clock and `getrusage` at every phase boundary, and measured both
arms on three confirmation fixtures per cell, twelve repetitions each:

- **Most of the process wall is not the job.** On `n13a0` the phases sum
  to 324 µs of a 810 µs process wall; the rest is the spawn, the image
  page-in and the reap that the evaluator's protocol charges to every
  trial. Both arms pay it.
- **The phase with the worst throughput is shared**: `curve_and_targets`
  runs at 1.10 ns per instruction on `n13a0` against 0.17–0.29 for the
  job's own phases, because it contains the first `blake3` hash (whose
  runtime CPU-feature detection costs 30–37 µs on this machine, measured
  separately) and 19–23 of the job's minor page faults.
- **The job faults about 50 pages** in and touches a 2.75 MB text segment.

Two changes were tried in development and **dropped**, which is why they
are not in this round:

- *No runtime CPU probing* (`blake3`'s `pure` feature and a compile-time
  `pclmulqdq`): 3.4% fewer instructions but 10 µs *more* on the clock,
  because portable `blake3` costs more than the `CPUID` trap it removes.
- *Pre-faulting with `madvise(MADV_POPULATE_READ|WRITE)`*: the phases do
  get much faster (324 µs → 252 µs on `n13a0`, and the faults vanish from
  them) but the two system calls cost about 380 µs, so the process wall
  rose from 810 µs to 1119 µs. The faults are real and the cure is dearer
  than the disease.

## What changes, and for whom

**Baseline source** `lto` (`--source-root`; [round14-lto.patch](round14-lto.patch)
and [round14-arena-tests.patch](round14-arena-tests.patch) on the round-0013
source), shared by every arm including rho:

- The release profile builds the worker as **one optimisation unit**: fat
  LTO, a single codegen unit, and dead sections dropped at link time. No
  source of the library or the worker changes, so every report byte is
  unchanged. Text falls from 2,749,393 to 2,067,481 bytes (25%), data from
  167,592 to 111,080.
- Tests only: the three arena tests took the process-global bump pointer
  for granted while the parallel test harness allocated from it too
  (reported by a review bot on PR #403). They now hold one lock and retry
  the sequences whose claim is about the most recent block. The measured
  executable contains none of this.

**This round is expected to move the headline ratio slightly against the
IC arm**, and it is registered that way before the run. Cross-crate
inlining makes *both* arms faster, and in development it took 1.6–4.8% of
their instructions — but slightly more of rho's (0.952–0.967 of its
round-0013 count) than of the IC arm's (0.964–0.984), because rho runs a
longer stream of the same single-word arithmetic. The absolute job gets
faster for both: natively, the IC arm 0.754–1.069 ms against round 0013's
0.814–1.133, and rho 0.871–1.198 against 0.898–1.257. `beats_rho_strict`
is expected to hold on every cell; a *worse* winner/rho ratio than round
0013 with both arms faster is the expected and fully reported outcome, and
the strict-win record stays with round 0013 unless this round's ratios
beat it.

**Challengers**, both IC-only, both byte-identical to the incumbent on all
304 frozen round-0013 fixtures:

- `lto_canon` ([round14-canon.patch](round14-canon.patch)): the orbit name
  of an abscissa found from the longest circular run of zeros rather than
  by trying all `n − 1` rotations, with the rotation loop kept as the test
  reference. Round 0013 measured it at 0.9758 of the incumbent; under one
  optimisation unit development puts it at 0.950–0.997 by cell, so it sits
  **on** the 0.98 no-regression line and is expected to fail it on the
  small cells.
- `lto_ic` ([round14-fastio.patch](round14-fastio.patch) on `lto_canon`):
  adds the report serialisation change (decimal digits two at a time from
  a table, appended without a UTF-8 validation pass). Development: 0.906–
  0.945 of the incumbent's instructions. It has now failed the native gate
  in rounds 0011, 0012 and 0013 on a clock that could not resolve 5–8%;
  this round asks the same question of a faster baseline, and a third
  *retained* verdict is a fully reported outcome.

## Objective and gates

`--objective rho`, as rounds 0007–0013. Incumbent: the round-0013 winner
configuration (`batch_trials: 1`) on the baseline above. A challenger is
promoted only if it passes the no-regression gate against the incumbent
(instruction ratio at most 0.98 with the paired upper limit below one,
native upper limit below one, no cell more than 10% worse) and the strict
rho gate on confirmation and replay. The decision records the winner's
`winner_over_rho`, `beats_rho_strict` and `rho_parity`.

## Parent, seed, budget

Fresh seed 2026091614; target count 1; pilot profile; 2,400 paired-job
budget; one pinned CPU; 8 GiB cap; 60-second watchdog; Valgrind 3.22.0
`Ir`; native progress recorded; the round-0010 evaluator (spawn without a
fork, caps before job delivery).

## Boundary, floor, class, honesty

Unit and boundary unchanged. Base support, `m = 3`, no direct relations,
every verification obligation, the worker's phase dumps and report fields
are unchanged; the final check is still `[d]G = Q` in the general
arithmetic on every recovered scalar. Class: engineering, and of the
build rather than the algorithm: nothing the compiler was asked to do
changes what is computed. Absolute times and instruction counts are
comparable with rounds 0010–0013 (same evaluator, same host, same
compiler) and not with earlier rounds. No arithmetic-complexity,
family-wide or cryptographic-size claim follows. Fresh fixtures from the
new seed; every failure retained; panels stay separate.
