# B6: fuzzing and differential checks

**Declared 2026-10-01, before any B6 measurement.** This is Track B's
B6 in `research/notes/index-calculus/IC_TOOL_PROGRAM.md` §9:

> A seeded generator of parameter files, valid and corrupted, runs in CI
> with a time budget. Answers are checked against `oracle.py` for
> `n ≤ 61`. Larger fields are checked against a slow reference
> implementation.

It is done when there is "no panic and no wrong answer within the
budget". Nothing below changes after the first measurement, except by a
dated amendment at the end.

**Disclosure.** A development version of the generator ran before this
declaration, against B2's development build on 2026-10-01 (seeds 1–8,
about 2,400 documents). It found no panic, no undocumented exit and no
wrong answer. That was a development check, not B6's measurement.

## What B6 does

1. **The generator**, [`../../fuzz/fuzz_v2.py`](../../fuzz/fuzz_v2.py),
   frozen with this declaration. It mixes three kinds of document:
   - random documents across every field kind, curve form, target form
     and method, with corrupted keys and values;
   - the conformance cases' parameter files, mutated lightly: only the
     target or the method changes, so most still run;
   - the same files mutated anywhere.
2. **The checks**, on each document, under `ic check` or `ic price` with
   a wall budget of 1–3 s:
   1. The exit status is 0–4: never a panic (101) or a signal.
   2. A JSON report exists whenever the exit status is 0, 2, 3 or 4.
      Exit 1 is a run its budget stopped, or a failed run.
   3. No "panicked" on stderr.
   4. **The answer is checked independently.** Every `complete` report's
      certificates are replayed as `[d]G = Q` in the generator's own
      Python arithmetic, which shares nothing with the Rust. A known
      logarithm must equal the reported scalar.
   5. **The resolved document gives the same answer.** For a quarter of
      the complete reports, the report's `resolved` document, run again,
      must give the same scalar (design §6).

   The plan names `oracle.py` for the answers. The replay in check 4
   serves the same purpose at every width: it computes `[d]G` directly,
   needs no size limit, and shares no arithmetic with the Rust.
3. **CI.** A workflow on pull requests that touch the `ic` tool. It
   builds `ic` and runs the conformance suite at the accepted steps. It
   then runs the generator for 600 documents, with a seed taken from the
   run number. The time budget is ten minutes.
4. **The campaign**, which is B6's measurement: seeds 1–20, 1,000
   documents each, 20,000 in all. It runs on the build that carries
   every accepted step through B2. It is untimed, and runs under the
   benchmark lock so that it never shares the machine with a timed
   process.

## Class

**Robustness.** B6 changes no pipeline. A finding is a bug in the step
that introduced it. It is fixed in that step's code before that step is
measured, or in a B6 fix that carries the pin and the timing check of
the rounds before it.

## When B6 runs

After B2 is accepted. CI needs the steps through B2 on `main`. Until
then, the conformance suite's cases for B0–B3 would fail there, since
their code is not merged.

## Acceptance

B6 is **accepted** when the campaign finds none of the following:
- a panic;
- an undocumented exit;
- a missing report;
- a wrong answer;
- a resolved document that disagrees.

CI must also pass. Each finding is recorded with its seed and document.
B6 is not accepted until every finding is fixed and the campaign run
again.

## Inadmissible

- Changing the generator after a campaign run to avoid a finding.
  `SHA256SUMS` pins it.
- Excluding a document or a seed from the campaign.
- Counting a budget stop (exit 1 without a report) as a wrong answer, or
  a wrong answer as a budget stop.

## Cost

- The campaign: 20,000 processes. Most take well under 0.1 s and none
  more than 3 s, so it takes under an hour.
- CI: about five minutes a pull request.

## Amendment 1 (2026-10-01, before the campaign): the generator's method mutation

- **The defect.** A development run of the generator, on B4's build with
  seeds 101–110 rather than the campaign's, crashed the generator itself
  at seed 106.
  - A heavy mutation's "junk" pass had made a seed document's
    `method.index_calculus` a non-object.
  - A later "method" pass then wrote into it, and Python raised a
    `TypeError`.
  - The tool was not at fault, and the run found nothing in it: 6,300
    documents on the other nine seeds, with no panic, no undocumented
    exit, no missing report and no wrong answer.
- **The fix.** The "method" mutation now replaces a non-object
  `method.index_calculus` with an object before writing into it.
  - It draws nothing more from the generator's random stream, so every
    document the old generator made is made again.
  - On seeds 101 and 104, 300 documents each, the old and the new
    generator gave identical outcome histograms.
  - Seed 106 now runs its 1,000 documents, and none is bad.
  - No seed document has a non-object `method`, so that is the only path
    to the crash.
- **`SHA256SUMS`** pins the amended generator (`f150f910…`).
- **Unchanged:** the campaign's seeds and counts, the checks, the
  acceptance and the inadmissible list. The campaign has not run.
