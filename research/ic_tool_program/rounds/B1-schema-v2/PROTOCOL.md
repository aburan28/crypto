# B1: parameter schema v2, and the Koblitz pipelines on imported instances

**Declared 2026-10-01, before any B1 code.** This is Track B's second
round in `research/notes/index-calculus/IC_TOOL_PROGRAM.md` §9. It is
built to the design in
[`../../design/schema-v2.md`](../../design/schema-v2.md), which fixes:
- the schema;
- the checks and their codes;
- the routing;
- the conformance cases.

Nothing below changes after the first measurement, except by a dated
amendment at the end.

## What B1 does

B1 is the design's §10 row:
- **Schema v2's parser** (design §3), strict at every level.
- **The checks for binary instances:**
  - §4.1 and §4.2's binary rows;
  - §4.3 and §4.4;
  - §4.5's first two methods: the Koblitz recurrence, and counting over
    a subfield.
- **The disclosures for binary fields** (§4.6).
- **An importer.** It builds the pipeline's `KoblitzCurve` from a
  validated document: the modulus, `a`, `b`, `r`, `h` and the
  generator, with `λ` computed. `kic`, `rho-koblitz` and `trivial`
  then run on it.
- **Routing among those three pipelines** (§5.1–§5.2), with the gate
  codes.
- **`ic check`**, which reads v2 and inspection v1 files, and
  `ic check --translate`, which writes a v1 workflow file's v2
  translation (§8).
- **`ic price` reading v2**, single-target by construction.

v1 files keep their parser and their code path. B1 adds a path beside
it and changes nothing on v1's.

## Class

**Robustness.** B1 changes what the tool accepts, not what it computes
on any input it already accepted. It claims no speedup.

## Arms

- **The base**: the newest accepted baseline when B1 runs, as a commit,
  recorded in B1's manifest. B0 comes before B1, because B1's refusals
  build on B0's.
- **B1**: that commit plus B1's change, built with
  `IC_BUILD_COMMIT=$(git rev-parse HEAD)`.

## Measurements

1. **Tests.** `cargo test --release --bin ic`, and the library tests of
   every module B1 touches.
2. **Conformance.** `conformance/v2/run.py --through B1` on B1, which is
   C001–C031. The same runner also runs on the base, and its failures
   are recorded as what the base could not do.
3. **The pin, untimed.** B1 runs all 90 suite rows from their v1 files.
   Every output must equal v0's from R01's profile pass: the counts,
   both arms' scalars, rho's counts and the verification flags.
4. **The translation, untimed, on all 90 rows.**
   - B1's `ic check --translate` must give the same JSON value as a
     translation written in Python by the design's §8 rules. That
     translation is the one `conformance/v2/make_cases.py` uses for
     C009, applied to each row.
   - `ic price` on each translation must give outputs identical to that
     row's v1 output in step 3.

   That is 180 processes.
5. **No slowdown on v1 files, timed.** The base against B1 on `M1`'s 22
   rows from their v1 files: five rounds ABAB, isolated, 220 processes,
   with the programme's runner. The figure is the paired cold-time
   ratio per size.
6. **No slowdown on the v2 path, timed.** B1 on `M1`'s 22 rows, v1 file
   against v2 translation: five rounds ABAB, 220 processes. Real use
   goes through the importer, so its cost must match the v1 path's.

## Acceptance

B1 is **accepted** when all of the following hold:
- every test passes;
- C001–C031 pass on B1;
- the pin holds on every row;
- both parts of the translation check hold on every row;
- no size regresses beyond its A/A band, in either timed comparison. A
  regression is a size whose interval's upper end lies below the lower
  end of that size's A/A interval: R01's, or the current round's own
  when the host has changed.

It is **stopped** if any v1 output differs. It is **rejected** if a test
or a case fails, the translation check fails, or a size regresses. A
rejected B1 is fixed and declared again, by amendment.

## Inadmissible

- Loosening a case after a run. The only exception is the design's
  `until` rule, applied by a later step.
- Changing a file under `conformance/v2/` after this declaration.
  `SHA256SUMS` pins the B1 files.
- Counting a contended run.

## The registry

B1's results PR teaches `scripts/build_curve_registry.py` to read v2
documents. It then registers the curves B1's cases name and the registry
lacks, as AGENTS.md §11 requires:
- `y² + xy = x³ + 1` under `z^31 + z^6 + 1` (C011);
- the degree-32 curve (C030);
- the gate curve's EC1 representation with its frozen generator
  (C027, C028).

## Cost

- The pin and the translation: 270 untimed processes.
- The timing checks: 440 timed processes, about 90 minutes.
- The tests and the conformance suite: a few minutes.
