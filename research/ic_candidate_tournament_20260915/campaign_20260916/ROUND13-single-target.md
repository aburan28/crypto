# Round 0013 pre-registration: one inversion per scalar product, and what it does to rho

Written before any measured stage of this round was run.

## Where the job's cost sat after round 0012

The archived round-0012 receipts put the IC arm's median complete job at
1.67M instructions and 1.155 ms against rho's 3.94M and 1.448 ms. A
line-level profile of the winner (a release build with line tables, phase
costs identical to the frozen executable's) shows that the two shared
phases — the curve construction with the target lift, and the general
final check `[d]G = Q` — spend most of their instructions in field
**inversions**: every scalar product in the library, the general
`scalar_mul` of the final check and the single-word `FastCurve::mul` of
the generator search (`[h]P` and the order test `[r]P` on every
candidate), is an affine double-and-add with one inversion per step. On
`n23a0` that is 0.53M of the IC arm's 2.98M instructions in the final check
(55% of it in the Itoh–Tsujii inversion) and 0.48M in the construction
(63% of it in `Gf2::inv`).

The same profile of the rho arm shows something the earlier rounds did not
look for: rho's own solve, 6.39M of its 7.47M instructions on `n23a0`,
begins every restart by drawing its jump table and its walk starts as
`[a]G + [b]Q` — two of those affine scalar products per jump and per walk —
and that setup, not the walk, is most of what rho pays on this panel.

## What changes, and for whom

**Baseline source** `ld` (`--source-root`; [round13-ld.patch](round13-ld.patch)
on the round-0012 source), shared by every arm including rho: the general
`scalar_mul` and the single-word `FastCurve::mul` and `mul_u64` run in
López–Dahab projective coordinates (Hankerson–Menezes–Vanstone algorithms
3.24 and 3.25, mixed addition, the degenerate `Q = ±P` cases through the
affine formulas), one field inversion for the whole product. Every result
is the same affine point: the affine double-and-add is kept as
`scalar_mul_affine` / `mul_affine` and tests compare the two on every
catalogue curve with scalars from 0 to 40, the order and its neighbours and
pseudo-random scalars of nine widths, and on nine Koblitz cells with random
points, `O`, the subgroup order, the group order and the cofactor. Nothing
else in the library moves; the worker and the cargo config are those of
round 0012.

Rho is the shipped solver, untouched; it runs on the same single-word curve
and its setup uses the same `mul_u64`, so it gets the same improvement.
**This round is expected to move the ratios toward one**, because rho's
setup was a larger share of rho's job than the shared phases were of the
IC job. Development evidence before freezing (not a claim; Callgrind on
one confirmation case per cell, native on two fixtures per cell with
twelve repetitions under the evaluator's spawn):

| cell | IC instructions (round 0012 → `ld`) | rho instructions | IC/rho native (round 0012 → `ld`) |
|---|---:|---:|---:|
| n13a0 | 0.874M → 0.759M | 1.736M → 1.047M | 0.82 → 0.89 |
| n17a1 | 1.196M → 0.986M | 3.052M → 1.538M | 0.75 → 0.85 |
| n19a0 | 1.810M → 1.399M | 5.048M → 2.814M | 0.71 → 0.80 |
| n19a1 | 2.060M → 1.631M | 4.383M → 2.035M | 0.76 → 0.83 |
| n23a0 | 3.103M → 2.540M | 6.936M → 2.996M | 0.80 → 1.06 |

So the incumbent is expected to keep the strict win in instructions on
every cell (winner/rho 0.50–0.85 by the development profiles) and **to
lose it on the native clock on `n23a0`**, where the IC job, with fewer
instructions than rho's, takes longer to run them (a cache simulation
shows no memory stall in either arm; the difference is per-instruction
throughput). A retained incumbent with `beats_rho_strict: false` is a
fully reported outcome of this round, and the strict-win record stays
with round 0012, which holds under its own baseline; the RESULTS note
will say that the margins of rounds 0007–0012 included rho's affine setup.

**Challengers**, both IC-only, both byte-identical to the incumbent on
every frozen fixture:

- `ld_canon` ([round13-canon.patch](round13-canon.patch)): the orbit name
  of an abscissa (the least rotation of its normal-basis coordinate word)
  found from the longest circular run of zeros — the runs by doubling,
  their maximal length by a binary search over those masks, and only the
  one or two rotations that put such a run at the top compared — instead
  of trying all `n − 1` rotations; the loop is kept as the test reference
  and compared on every width the pipeline runs, random and sparse words,
  the all-ones word and periodic words. Development: 0.757M–2.403M
  instructions on the five cells, 0.95–0.99 of the baseline.
- `ld_ic` ([round13-fastio.patch](round13-fastio.patch) on `ld_canon`):
  adds the report serialisation change of rounds 0011 and 0012 (decimal
  numbers two digits at a time from a table, appended without a UTF-8
  validation pass). Those rounds measured that change alone at 0.965 and
  0.954 of the incumbent's instructions but could not resolve it natively;
  here it rides with `ld_canon`. Development: 0.720M–2.293M instructions,
  0.90–0.95 of the baseline; native, within the noise of the development
  clock (the twelve-repetition medians differ from the baseline's by less
  than their standard deviations on every cell).

The whole library test suite was run on the `ld_ic` tree (the superset of
the three): 2,428 pass, 97 ignored and 5 fail, the same five general
pair-table tests of `koblitz_index_calculus` that fail identically on the
unmodified round-0010 snapshot (recorded in the round-0011
pre-registration); the three new reference tests (projective against
affine products in the general and single-word curves, orbit name against
the rotation loop) pass. The worker's own four tests pass on `ld_ic`.

A rotation-invariant prefilter on the pair-table lookups was also tried in
development and dropped: it rejects 90% of the probes on `n23a0` but the
orbit naming it skips is a small part of each probe's cost, and marking
the invariant of every stored key made the table build dearer than the
scans got cheaper (`n23a0` 2.540M → 2.542M).

## Objective and gates

`--objective rho`, as rounds 0007–0012. Incumbent: the round-0012 winner
configuration (`batch_trials: 1`) on the baseline above. A challenger is
promoted only if it passes the no-regression gate against the incumbent
(instruction ratio at most 0.98 with the paired upper limit below one,
native upper limit below one, no cell more than 10% worse) and the strict
rho gate on confirmation and replay. `ld_canon` is expected to fail the
0.98 line on the small cells and `ld_ic` to pass it on instructions; on
the native clock neither is expected to resolve, and neither is expected
to restore the strict win on `n23a0`. Either way the decision records the
winner's `winner_over_rho`, `beats_rho_strict` and `rho_parity`.

## Parent, seed, budget

Fresh seed 2026091613; target count 1; pilot profile; 2,400 paired-job
budget; one pinned CPU; 8 GiB cap; 60-second watchdog; Valgrind 3.22.0
`Ir`; native progress recorded; the round-0010 evaluator (spawn without a
fork, caps before job delivery).

## Boundary, floor, class, honesty

Unit and boundary unchanged. Base support, `m = 3`, no direct relations,
every verification obligation, the worker's phase dumps and report fields
are unchanged; the final check is still `[d]G = Q` in the general
arithmetic on every recovered scalar, in the library's general curve
module, with the projective formulas tested against the affine ones it
replaces. Class: engineering. Absolute times and instruction counts are
comparable with rounds 0010–0012 (same evaluator, same host, same
compiler) and not with earlier rounds. No arithmetic-complexity,
family-wide or cryptographic-size claim follows. Fresh fixtures from the
new seed; every failure retained; panels stay separate.
