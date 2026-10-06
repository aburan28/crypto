# R02b: the wide-tail kernel, re-tested on fresh holdouts with power

**Declared 2026-10-01, after R02's rejection (#1152) and before
any R02b build or run.** This is Track A in
`research/notes/index-calculus/IC_TOOL_PROGRAM.md` (§8, A2). Nothing
below changes after the first run, except by a dated amendment appended
at the end.

## What R02 found

R02 gave the AVX-512 batched addition the two suite fields whose
polynomial has `n + deg t = 66`. Its results are in
[`../R02-wide-tail-kernel/README.md`](../R02-wide-tail-kernel/README.md).
- **The mechanism held.** With identical outputs on all 90 rows,
  collection ran 1.29–1.31× faster at the two sizes.
- **Cold time:**

  | curve | suite, 40 pairs | holdouts, 10 pairs |
  |:--|--:|--:|
  | `icv1-f2m59-tm943548413-98844ecc` | 1.182 [1.151, 1.213] | 1.161 [1.089, 1.239] |
  | `icv1-f2m61-t158598901-ab42b6c5` | 1.267 [1.242, 1.292] | 1.262 [1.216, 1.310] |

- **The rule** asked for a lower end above 1.10 on each set separately.
  The holdouts at `icv1-f2m59-tm943548413-98844ecc` missed it by 0.011,
  so R02 was rejected. That decision stands, and R02's data are not
  reused here.
- **The callgrind control also failed.** Under valgrind both arms run the
  portable paths, which R02's amendment 3 required to be
  instruction-identical within `10⁻³`. Every arithmetic kernel matched
  to the instruction, and the candidate ran 0.085–0.102% fewer
  instructions in total, all in the compiler's code for two scan
  functions whose source the patch does not touch
  (`witnesses_fast_scan`, `pairs_for_key`).

## Why a re-test, and what limits it

**The holdout set was too small for its rule.** It had two rows per
size, ten pairs, against the suite's forty. At
`icv1-f2m59-tm943548413-98844ecc` its interval's half-width was 6.7%,
against the suite's 2.6%. The two sets' point estimates agree, at 1.161
and 1.182.

R02b asks the same question on new rows, with the holdout set's power
fixed before anything runs:
- **The same claim, on new data.** It does not re-analyse R02's runs.
- **The last test of this kernel.** If R02b is rejected, the kernel is
  retired, and no third test is declared.
- **Plan §11 applies as written.** R02b is the second consecutive
  declared round on the collection scan. If it fails, the scan is set
  aside, and Track A moves to the next phase in the profile.

## Hypothesis

The same as R02's. At `icv1-f2m59-tm943548413-98844ecc` and
`icv1-f2m61-t158598901-ab42b6c5` the scan's batched addition runs the
scalar path because the 8-lane kernel refuses `n + deg t = 66`. The
candidate carries `H·t`'s bits above `z^63` into the second fold, so the
kernel serves both fields, with the same field elements lane for lane.

## Class

Engineering, if it gains: the counts are unchanged by construction.

## Arms

- **The base**: the newest accepted baseline when R02b runs, as a
  commit recorded in the manifest. That is R03's candidate if R03 is
  accepted, and v0′ (`c1a2e5f8`) otherwise.
- **The candidate**: the base plus R02's
  [`candidate.patch`](../R02-wide-tail-kernel/candidate.patch),
  unchanged. It touches only `src/cryptanalysis/koblitz_fast.rs`. R03
  touches only `factorise_u64` in `koblitz_index_calculus.rs`, so the
  two apply without interaction.

## Before any timed step

1. **Tests**: `cargo test --release --lib -- koblitz_fast` on the
   candidate's tree, which includes the patch's wide-tail and
   scalar-equivalence tests.
2. **The pin, untimed**: the candidate on all 90 suite rows. Every output
   must equal v0's from R01's profile pass: the counts, both arms'
   scalars, rho's counts and the verification flags. The names must
   agree, as in R02's pin.

## The callgrind control

R02's control, run again after the timed steps (`run.py callgrind`):
`ic workflow` under callgrind, both arms, on `M1`'s first target at
`icv1-f2m61-t158598901-ab42b6c5`, `icv1-f2m59-tm943548413-98844ecc` and
`icv1-f2m41-tm2308219-7f48b14a`. Valgrind hides AVX-512, so both arms
take the portable paths.

**It holds when, on every row:**
- every arithmetic kernel is identical to the instruction: every
  function of `koblitz_fast`, `Gf2`'s field operations, and the bulk
  key (`PairSumTable::keys_of`);
- the whole run's instructions are within `5·10⁻³` of the base's;
- both arms recover the same logarithm.

Every other function whose count differs is named in the analysis, with
both counts.

**This restates R02's condition, and says so.** R02 asked for the total
within `10⁻³`. It is written here knowing R02's figure, 1.00102 at
`icv1-f2m41-tm2308219-7f48b14a`, and before any R02b run. The reason:
the control exists to show that the candidate's work is the kernel's.
Identity of the kernels function by function shows that directly. A
total that moves by a tenth of a per cent through codegen in unchanged
callers does not bear on a 1.18× effect, and R02's nine narrow-kernel
sizes, which run those callers with the same codegen, read 0.973–1.027
in time. The total stays bounded, at five times R02's shift. R02's
verdict on its own control stands.

## Rows

- **The suite rows (40 pairs at each target size):**
  - `M1`'s 22 rows: two targets at each of the eleven sizes;
  - at the two target sizes, the other six suite rows each (`M2`–`M4`).

  That is 34 rows. Five rounds, ABAB, isolated: 340 processes.
- **The fresh holdouts (40 pairs at each target size):** at the two
  target sizes only, eight rows each, by suite v1's own construction
  (`make_suite.py`), as R02's holdouts were built:
  - recipe seeds 206, 207, 208 and 209;
  - two targets per seed, `T103` to `T110`;
  - `public_hash_seed` 23103 to 23110, and rho seeds `0x230000 + 103`
    to `+ 110`.

  Neither the seeds nor the targets have been run by any round. Five
  rounds, ABAB, isolated: 160 processes.
- **The extension.** At a target size, a set whose interval half-width
  (`hi / geomean − 1`) exceeds 3% after five rounds gets rounds 6–10, and
  its figure pools all ten. This depends only on the interval's width,
  never on its position.

## Power

R02's per-pair spread at `icv1-f2m59-tm943548413-98844ecc` gives a 95%
half-width of 2.6% to 3.4% at forty pairs. If the true ratio is 1.16,
the expected lower end is near 1.12–1.13. If it is 1.18, as R02's suite
read, the lower end is near 1.14–1.15.

## Prediction

- **`icv1-f2m59-tm943548413-98844ecc`: 1.18×** [1.15, 1.21], on both
  sets.
- **`icv1-f2m61-t158598901-ab42b6c5`: 1.27×** [1.24, 1.29], on both
  sets.
- **The nine other sizes: 1.00×.** Their kernel is the narrow one, whose
  instructions the patch leaves unchanged.

## Success and stop

**Accepted** if all of the following hold:
1. The tests pass, and every pinned output is identical.
2. At both target sizes, the paired cold-time ratio's 95% interval lies
   above 1.10, on the suite rows and on the fresh holdouts separately.
3. No size regresses beyond its A/A band. A regression is a size whose
   interval's upper end lies below the lower end of R01's A/A interval,
   or the round's own if the host has changed.
4. The callgrind control holds, as stated above.

**Rejected** otherwise. The kernel is then retired.

**Stopped** if a test fails, an output differs, or a verification fails.

**If accepted**, the candidate becomes the next baseline, with its
binary hash, host, and `S` from this round's runs. It also gets the
rule's comparison that R02's protocol required: §23's protocol at its
six sizes, with 64 targets, on the candidate (384 processes).

## Inadmissible

- Reusing R02's runs, pooling them with R02b's, or choosing between the
  two rounds' figures.
- Changing the field polynomial, the suite, the recipes or the targets,
  or drawing new holdout seeds after a run.
- Pooling contended runs.
- Quoting collection or build as the speedup: the speedup is cold time.
- Claiming a gain at a size whose kernel did not change.

## Cost

- The tests and the pin: a few minutes, and 90 untimed processes.
- The callgrind control: six untimed runs under valgrind, about an hour
  on one core.
- The suite rows: 340 timed processes. The holdouts: 160, and up to 160
  more if they are extended. About two and a half hours together.
- If accepted, the rule's comparison: 384 processes, about 40 minutes.

## Order

R02b runs after R03's decision and before B0. B0's base is then the
newest accepted baseline after R02b. B0 touches neither `koblitz_fast.rs`
nor `factorise_u64`, so it merges with both without interaction.

## Amendment 1 (2026-10-01, before any R02b build or run): R05 runs first

- **The order.** R02b runs after R05's decision, not after R03's. R05,
  declared in the same pull request
  ([`../R05-presence-filter/PROTOCOL.md`](../R05-presence-filter/PROTOCOL.md)),
  tests a sharper presence filter on the same scan.
- **Why.** R04's check after its run found the filter lever: the admitted
  stage, 32–42% of the scan at the three largest sizes, is spent almost
  entirely on false positives. An exploration measured the lever at
  1.27–1.37× in cold time, against R02's 1.16–1.27× for this kernel.
  Plan §11 sets the scan aside after two consecutive failed rounds on
  it, so whichever of the two runs second would not run if the first
  failed. The larger lever goes first.
- **If R05 is rejected,** the scan has failed twice running (R02, R05). It
  is set aside, and R02b does not run. This declaration stays on record,
  unrun, and the kernel stays as R02 left it.
- **If R05 is accepted,** R02b runs as declared, with R05's candidate as
  its base. A rejection then retires the kernel, as declared, but is not
  a second consecutive failure, so the scan stays open.
- **The base and the patch do not interact.** R05's patch touches the
  filter's code in `koblitz_index_calculus.rs`, `gpu/ecc2k/` and
  `examples/load_fold_table.rs`. R02's touches only `koblitz_fast.rs`.
- **Unchanged:** the hypothesis, the rows and holdouts, the callgrind
  control, the prediction, the rule and the cost. On R05's base the
  subtraction is a larger share of cold time, which if anything raises
  the ratio R02b measures. The prediction stays as declared.

## Amendment 2 (2026-10-02, before any R02b build or run): native tooling, v2, and the host

**The base and the candidate.**
- **The base is v2:** R05's candidate, `edcb0bec`, whose binary has
  SHA-256 `76a2a2fd…` (`baselines.json`). Amendment 1 named it.
- **The candidate** is `edcb0bec` plus R02's
  [`candidate.patch`](../R02-wide-tail-kernel/candidate.patch), unchanged
  (SHA-256 `2d3fb1b5…`). It applies to v2 without conflict.
- **A tree on v1 was committed and abandoned.** It is `c88cf863`: R03's
  candidate plus the same patch, made just before amendment 1 moved
  R02b after R05. It was never built or run.

**The tooling is native.** AGENTS.md's rule against Python tooling
(#1180) merged after this declaration.
- **Each step runs on `icprog run r02b <step>`,** with `isolated_bench`
  (plan §10a, N2): `manifest`, `pin`, `aa` (below), `compare`,
  `holdout`, `extend` and `callgrind`.
- **The analysis is `icprog analyse r02b`.**
- **Each piece reproduces the declared script's frozen outputs, byte
  for byte:**
  - the pin reproduces R05's `pin.json` from R05's run tree;
  - the callgrind phase split reproduces R02's six profiles and R01's
    two;
  - the runner resumed R05's holdouts;
  - the analysis is R05's, whose statistics and run-tree reading are
    the ones that reproduce R03's frozen analysis.
- **The callgrind control is new code.** It has no frozen output to
  match, so it was checked on R02's profiles against R02's README
  instead. The base ran 1.00102, 1.00085 and 1.00098 times the
  candidate's instructions, and every kernel was identical.
- **`run.py` and `analyse.py` stay** as the record of what was declared.

**The host has changed, so the A/A is the round's own.**
- R01's manifest names the kernel build `6.18.44-fc-v50`, and this
  container runs `-fc-v51`. The declared manifest compares that field,
  so it reports the host as changed.
- "Success and stop" 3 then measures a regression against the round's
  own A/A. Plan §5 asks for one every round in any case.
- **The A/A:**
  - the base against a byte-identical copy of itself;
  - on `M1`'s 22 rows, five rounds, ABAB, isolated;
  - 220 processes, as R01's;
  - run after the pin and before the comparison.
- **The declared `analyse.py` read R01's band whatever the host.** The
  native analysis reads the round's own band when `aa-source.json`
  reports the host as changed. A size without a band is then a reason
  to reject.
- If the manifest finds the host matching R01's after all, the bands are
  R01's as declared, and the A/A runs as a check only.

**The holdouts' files are checked by `SHA256SUMS`.**
- Their sums were written from suite v1's construction when the files
  were committed (#1157).
- The declared runner rebuilt the files with `make_suite.py` and
  compared. That script stays unported (plan §10a).

**The analysis adds** these, beside the declared figures:
- the extension rule's five-round test, checked against
  `extended.json`;
- each process's isolation tool;
- the own A/A's record.

**Unchanged:** the hypothesis, the rows, the callgrind control's rule,
the prediction and the success rule. The cost grows by the A/A's 220
processes, about an hour.

## Amendment 3 (2026-10-05, after its pin and part of its A/A): suspended, because the base moved to main

**What ran.** No comparison row ran, so no candidate timing exists. Both
stopped run trees are archived in [`stopped/`](stopped/), checked by
`SHA256SUMS`.
- **Attempt 1, on `6.18.44-fc-v51`, 2026-10-02 from 02:34 UTC:**
  - the manifest, and the pin, which held;
  - one A/A pair, before a container rebuild ended the run.
- **Attempt 2, on `6.18.44-fc-v70`, 2026-10-05 from 17:35 UTC:**
  - a fresh manifest, and the pin, which held again;
  - 52 of the A/A's 220 processes, one of them a retry of a run the
    isolation tool marked contended;
  - one container restart on the same host build, at 17:47, recorded in
    `host-resumed.json` and `resumes.log`;
  - stopped deliberately at 19:27 UTC.

**Why it stopped.**
- **Main's head has the Gf2 two-fold reduction** (#1242, merged
  2026-10-02 as `809a7318`). It speeds the scalar subtraction this
  kernel replaces.
- **A diagnostic** ran one process per arm on `M1`'s first target, under
  the benchmark lock but not isolated
  ([R07's PROTOCOL.md](../R07-main-head/PROTOCOL.md), "The diagnostic").
  It found:
  - collection fell from 4276 to 2615 ms at
    `icv1-f2m61-t158598901-ab42b6c5`, and from 790 to 509 ms at
    `icv1-f2m59-tm943548413-98844ecc`;
  - on main a scanned summand cost 35–38 ns at all three of the largest
    sizes, the narrow 8-lane kernel's cost at
    `icv1-f2m53-tm56619371-dac20a85`;
  - the outputs were identical.
- **So the lever R02b tests has largely moved under it.** R02b's base
  rule names the newest accepted baseline, and R07 now measures main's
  head as that baseline.

**What happens next.**
1. **R07 runs first.** R02b waits for R07's decision.
2. **Then a go/no-go exploration on the newest baseline:**
   - the base itself, against the base plus R02's
     [`candidate.patch`](../R02-wide-tail-kernel/candidate.patch);
   - if the patch no longer applies, it is ported, and the port is
     disclosed with its diff;
   - the tests are this protocol's;
   - `M1`'s two rows at each target size, three rounds, ABAB,
     isolated.
3. **The rule.** If the cold-time ratio's geometric mean is below 1.05
   at both target sizes:
   - R02b is withdrawn without running;
   - the kernel is retired, as a rejection would have retired it;
   - the withdrawal counts as neither an acceptance nor a failed round
     on the scan (plan §11).

   Otherwise R02b runs as declared on that base, with a fresh manifest,
   pin and A/A, and the exploration is disclosed.
4. **The stopped trees** are never pooled with any later run.

**Unchanged:** the hypothesis, the rows and holdouts, the callgrind
control's rule, the prediction and the success rule.
