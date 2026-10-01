# R02: the AVX-512 batched addition for fields with a wide reduction tail

**Declared 2026-10-01, before any candidate code and before R01's
results were read.** This comes from ledger §23's reports and a reading
of the code; R01's profile checks it rather than sources it. The plan
is `research/notes/index-calculus/IC_TOOL_PROGRAM.md` (Track A, A2).
Nothing below changes after the first run except by a dated amendment
appended at the end.

## Hypothesis

**The scan is slower at the two largest sizes because the 8-lane kernel
is off there.** At `icv1-f2m59-tm943548413-98844ecc` and `icv1-f2m61-t158598901-ab42b6c5` the collection
scan's batched addition falls back to the scalar path. The SIMD kernel
`Simd512` in `src/cryptanalysis/koblitz_fast.rs` refuses those fields:
- It reduces a product `L + z^n·H` through the tail `t` of the field's
  polynomial (`z^n ≡ t`).
- It needs `H·t` to fit in one word, which means `n + deg t ≤ 65`.
- Both fields have `n + deg t = 66`:
  - `GF(2^59)` is `z^59 + z^7 + z^4 + z^2 + 1`;
  - `GF(2^61)` is `z^61 + z^5 + z^2 + z + 1`.
- Every other suite size qualifies, at 44–61.

**The cost per summand follows that split exactly.** In §23's reports,
at `0bf67f16`, one thread, isolated, a scanned summand costs:

| curve | `n + deg t` | kernel | ns per summand | units per summand | collection's share of the set-up |
|:--|--:|:--|--:|--:|--:|
| `icv1-f2m47-t22705043-f4e44623` | 52 | on | 38.7 | 1.44 | 42% |
| `icv1-f2m57-tm747311035-c1f545af` | 61 | on | 44.8 | 1.48 | 31% |
| `icv1-f2m41-tm2308219-7f48b14a` | 44 | on | 33.4 | 1.37 | 67% |
| `icv1-f2m53-tm56619371-dac20a85` | 59 | on | 48.7 | 1.59 | 65% |
| `icv1-f2m59-tm943548413-98844ecc` | 66 | **off** | 71.7 | 2.38 | 71% |
| `icv1-f2m61-t158598901-ab42b6c5` | 66 | **off** | 68.9 | 2.30 | 89% |

**The candidate.** The kernel carries the bits of `H·t` that pass `z^63`
into its second fold. Those bits are `h >> (64 − t_i)` for each tail
term `t_i ≥ 1`. The second fold's input becomes
`(f >> n) | (o << (64 − n))`. That input still has degree at most
`deg t − 2`, and the fold's output lands below `z^n` while
`2·deg t ≤ n + 1`, which is the kernel's existing second condition.

So the field elements are the same as the scalar path's, lane for lane.
The condition `n + deg t ≤ 65` is dropped.

- **Fields that already qualify keep their exact instruction stream.**
  The carry is a const generic, monomorphised away when off.
- **The kernel serves four callers:** the scan, the folded pair-table
  build, the full-scan variants and the descent's walk. All four take it
  at the two sizes.

## Class

Engineering, if it gains: the counts are unchanged by construction, so
the ratio to the floor is flat.

## A control first, on v0, before any candidate is timed

v0 at `icv1-f2m53-tm56619371-dac20a85`, `M1`'s two rows, five rounds ABAB. One arm sets
`KIC_SCAN_SIMD=0`, the existing same-binary switch to the scalar path;
the other is the default. That is 20 isolated processes.

- **Confirmed** if the scalar arm's collection time per scanned summand
  is at least 1.3× the kernel arm's (geometric mean over the 10 pairs).
- **Otherwise the mechanism is wrong.** R02 stops before the candidate
  is timed and records the control.

## Pinned outputs

- **Every row, untimed.** The candidate runs on all 88 S rows and both
  smoke rows. Every row's outputs must equal v0's from R01's profile
  pass: counts, both arms' scalars, rho's counts, and the verification
  flags.
- **Tests:**
  - the kernel against the scalar field on every degree 2–62 the
    relaxed condition accepts;
  - the batched addition against the scalar batched addition on the
    `icv1-f2m59-tm943548413-98844ecc` and `icv1-f2m61-t158598901-ab42b6c5` curves themselves;
  - the two suite fields marked wide.

## Timed comparison

- **Rows.** v0 against the candidate on all 88 S rows, five rounds,
  with the arms' order alternating (IC_TOOL_PROGRAM.md §5). That is 880
  isolated processes.
- **Holdouts.** Two holdout rows per size, on seed 205 with targets
  `T101` and `T102`: `public_hash_seed` 23101 and 23102, rho seeds
  `0x230000 + 101` and `+ 102`. Five rounds give another 220 processes.
- **The A/A.** R02 runs in R01's container session. It uses R01's A/A
  if the host manifest's CPU, kernel, memory and THP settings are
  unchanged, and runs its own on `M1`'s rows otherwise.
- **Callgrind, untimed.** Both arms run at `icv1-f2m61-t158598901-ab42b6c5` and
  `icv1-f2m41-tm2308219-7f48b14a` `M1-T01`, `--repeats 1`, as in R01. The figures are
  the instruction counts and the collection's share.
- **If R02 is accepted,** the new baseline gets the rule's comparison:
  §23's protocol at its six sizes with 64 targets, on the candidate.

## Prediction

- **`icv1-f2m61-t158598901-ab42b6c5`: 1.37×** [1.18, 1.53], the paired cold-time ratio,
  v0 over the candidate. Collection, 89% of the set-up, falls from 2.30
  units a summand to `icv1-f2m53-tm56619371-dac20a85`'s 1.59. The range runs from 1.9 to
  1.4 units. The build's and the descent's gains are left out, so the
  prediction is conservative.
- **`icv1-f2m59-tm943548413-98844ecc`: 1.31×** [1.15, 1.45], by the same arithmetic at
  71%.
- **The other nine sizes: 1.00×.** Their kernels are unchanged.

## Success and stop

**Accepted** if all of the following hold:
1. The control confirms the mechanism.
2. Every pinned output is identical.
3. At `icv1-f2m59-tm943548413-98844ecc` and `icv1-f2m61-t158598901-ab42b6c5`, the paired cold-time ratio's
   95% interval lies above 1.10, on the suite rows and on the holdouts
   separately.
4. No size regresses beyond its A/A band. A regression is a candidate
   interval whose upper end lies below the lower end of that size's A/A
   interval.

**Rejected** otherwise. The code is kept on record and the numbers go
in the ledger.

**Stopped early** if the control fails, if any output differs, or if a
verification fails.

## Inadmissible

- Changing the field polynomial: that changes every coordinate, so it
  is a different instance.
- Changing the suite, the recipes or the targets.
- Pooling contended runs.
- Quoting the stage figures (collection, build) as the speedup: the
  speedup is cold time.
- Claiming a gain at a size whose kernel did not change.

## Cost

- The control: 20 processes.
- The pin: 90 untimed.
- The comparison: 1,100 processes, about 2.5 hours.
- Callgrind: four runs, about 40 minutes.
- If accepted, the rule's comparison: 384 processes, about 40 minutes.

## Amendment 1 (2026-10-01, after R01's timed steps, before any R02 run)

**What R01 said.** R01's protocol, under "What follows", said: "a size
whose A/A interval is wider than ±5% gets ten rounds instead of five."
R01's A/A intervals are ±4–11% at ten of eleven sizes.

**Why the rule is miscalibrated.** R01's A/A interval rests on ten
pairs: two rows × five rounds. R02's comparison has 40 pairs per size:
eight rows × five rounds. Its interval is therefore expected to be about
half the A/A's width (`√(40/10) = 2`), so ±2–5%. The intent of R01's
rule is precision near ±5%, and five rounds of R02's design meet it.
Ten rounds would double the cost, about four more hours, and add
nothing to that intent.

**The amendment.** R02 keeps five rounds. Any size whose comparison
interval has a half-width above 5% then gets five more rounds, and its
figure pools all ten.

Nothing else changes: the rows, the order, the control, the pin, the
success condition and the stop rules. The figure is fixed now, before any
R02 process runs. One untimed sanity check ran before this amendment:
the candidate on `icv1-f2m61-t158598901-ab42b6c5` and `icv1-f2m59-tm943548413-98844ecc` `M1-T01`, outputs only,
under `taskset` on a machine that was compiling. Both outputs equal v0's.
Its collection phase read about 1.32× faster than v0's R01 profile row
at both sizes. That figure was seen before this amendment was written.
It is not evidence, it is not used, and the amendment's argument does
not depend on it.

## Amendment 2 (2026-10-01, before any R02 run): what callgrind measures

**What R01 found.** R01's callgrind of `ic price --repeats 1` showed that
most of the process is the pricer's own calibration, not the index
calculus. At `icv1-f2m41-tm2308219-7f48b14a`:
- `UnitBench`'s timing loop is 1.56 × 10⁹ instructions;
- the pass up to its first online boundary is 0.12 × 10⁹.

A cross-check on that process would therefore mostly count calibration
code. Some of it runs through the very kernel R02 changes.

**The amendment.** R02's instruction cross-check runs callgrind on the
pipeline alone.
- **The command:** `ic workflow --params <row> --dir <fresh dir>`, with
  `RAYON_NUM_THREADS=1` and no rho baseline in the file.
- **The arms and rows:** both arms, at `icv1-f2m61-t158598901-ab42b6c5`, `icv1-f2m59-tm943548413-98844ecc`
  and `icv1-f2m41-tm2308219-7f48b14a` `M1-T01`.
- **The options:** the same callgrind options as R01.
- **The figures:**
  - the total instruction ratio, v0 over the candidate, per row;
  - the per-function instructions of the batched addition (scalar and
    vector paths), the scan and the build.

A workflow run must also recover the same logarithm in both arms. This
replaces the `ic price` callgrind runs in "Timed comparison". The timed
comparison itself is unchanged.

## Amendment 3 (2026-10-01, before any R02 run): the baseline arm, two controls, the names and the run tree

**What moved.** After this protocol was written, `main` merged AGENTS.md
§11, which names every curve by its ICV1 slug, together with the
`src/` change that carries it out. From v0 (`46ae2014`) to this round's
base (`c1a2e5f8dda5613bf797c4c38afc8b0b433ea348`), `src/` changed in eight files:
- a new naming module, `curve_id.rs`, with its alias table and
  `mod.rs`;
- the labels: a report's `curve` field is now the slug;
- name matching in `ic rho`'s replays and in `ic experiment`, now through
  the registry;
- the unit calibration's pinned-ratio lookup, likewise;
- the text of one error message;
- a new strong-rho module, which the pricer does not call.

No pipeline step changed. The calibration runs before the timed
intervals start, and the labels are written after they end.

**The baseline arm is v0′.** v0′ is `main` at `c1a2e5f8`, the commit
the candidate is built on. The candidate is v0′ plus the kernel change
in `src/cryptanalysis/koblitz_fast.rs` and its tests, so the two arms
differ in that one file. Wherever this protocol names v0 as an arm (the
control, the comparison, the holdouts, callgrind), read v0′.

**v0′ is pinned like the candidate.** Before any timed step, v0′ runs
the pin on all 90 rows, untimed. Its outputs must equal v0's from R01's
profile pass. If any differs, R02 stops, because the baseline itself has
changed. The pin also checks the names: every report's `curve` must
equal the registry's slug for its row, and v0's old label must resolve
to that slug.

**v0 against v0′, timed: accounting.** v0 (R01's binary, `c776cc04…`)
and v0′ run on `M1`'s 22 rows, five rounds ABAB: 220 processes. The
figure is the paired cold-time ratio per size. It says whether `main`'s
renaming moved the cost. It gates nothing, because R02's speedup is v0′
over the candidate. A ratio outside the A/A band is recorded as
`main`'s own change, not R02's.

**The control switches off two kernels, not one.** `KIC_SCAN_SIMD=0`
turns off the AVX-512 batched addition (`Simd512::new`) and also the
AVX-512 canonical key (`canon_simd_enabled`). At
`icv1-f2m53-tm56619371-dac20a85` both are on by default, so the control
measures the two together. Its declared rule changes as follows:
- **A ratio below 1.3 still stops R02.** If the two kernels together are
  worth less than 1.3× in collection per summand, the addition alone is
  too.
- **A ratio at or above 1.3 no longer confirms the mechanism.** It shows
  only that the two kernels together pass the bar.
- **The mechanism is tested by the comparison itself.** At the two
  target sizes the candidate changes only the addition, and the AVX-512
  key runs in both arms.

**Callgrind becomes a control.** Valgrind 3.22 hides AVX-512: under it,
CPUID reports `avx512f`, `vpclmulqdq` and `gfni` as absent (R01's
correction). Every callgrind run therefore takes the portable paths, in
both arms, and cannot see the kernel. Amendment 2's step runs as
declared, with a new reading:
- **The two arms' portable paths must be instruction-identical:** a
  total ratio of 1.000 per row and the same logarithm.
- **A difference fails the step.** It would mean the candidate changed
  more than its kernel.
- **It yields no speed figure.** The speed evidence is the timed
  comparison and the holdouts.

**Names.** This protocol now writes every curve by its ICV1 slug, as
§11 requires. The rewrite is the commit before this amendment, made by
`scripts/check_curve_names.py --fix`. It changes no figure, row,
target, seed or rule.

**The run tree.** A run directory is named `<slug>/<recipe>-<target>`,
for example `icv1-f2m61-t158598901-ab42b6c5/M1-T01`, not by the suite's
row id. The holdouts' parameter files are named the same way, while the
`name` field inside them stays as suite v1's construction writes it. The
suite's own frozen files keep their names. The tree is committed as
`runs.tar.xz` with its SHA-256, as R01's is.

**R01's outputs** for the pin come from R01's archive (`runs.tar.xz`,
SHA-256 `ae20e0bd…`). It is checked against its hash and extracted in
place, where `.gitignore` keeps it out of git.

**Amendment 1's extension, made exact.** A size's interval half-width
is `hi / geomean − 1` of its paired cold-time ratio. A size above 5%
gets rounds 6–10, and its figure pools all ten. The extension applies to
the suite rows, not to the holdouts.

**Order and cost.**
1. manifest;
2. the control;
3. the pin, both arms;
4. R02's own A/A, only if the host differs from R01's;
5. the comparison, then its extension;
6. the holdouts;
7. v0 against v0′;
8. callgrind.

This adds 90 untimed processes and 220 timed ones, about 40 minutes.

**The kernel change is held back from the declaration's pull request.**
The repository's approval bot merges pull requests once CI passes,
drafts included, so a declaration carrying candidate code would land
that code unmeasured. The candidate was written after this protocol's
first commit, as local commits `76efcd1e` and `7c356fac` (formatting).
It lands with R02's results.

**What ran before this amendment.** Both arms were built:
- v0′, SHA-256 `a89e1be0…`;
- the candidate, SHA-256 `19240609…`.

The candidate's `koblitz_fast` tests pass, 28 of 28, including the
wide-tail and scalar-equivalence tests. The runner's pin logic was then
exercised untimed on two small rows: the smallest suite row's `M1-T01`
and the first smoke row. On both, both arms' outputs equal v0's, and
the names agree. Nothing was timed.
