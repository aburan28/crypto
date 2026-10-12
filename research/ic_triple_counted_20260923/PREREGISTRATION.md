# Pre-registration: sizing the triple collector by counting four-sums

Written and committed **before the counted sizing exists as code and before any
base other than two orbits has been run**, per `AGENTS.md` §4.  Nothing below is
a result.  Follows `research/ic_triple_table_20260923` (PR #668), whose arm this
changes in one function.

## Why

The triple arm sizes its base with W₃, the pair model's units carried over.  Its
rate `λ₃·cov₃` a rest, `λ₃ = C(|F|+2, 3)/r`, is modelled over the `|F|` rests of
a target, which counts about four times the witnesses that exist.  A target is
found **exactly** when it is a sum of four base points with at least one summand
in a built row, so the chance a fresh target hits is

```text
p(K, t) = 1 − exp(−[C(|F|+3, 4) − C(|F|−2nt+3, 4)] / r),      |F| = 2nK
```

and the expected work, in the solver's units, is

```text
U(K, t) = t·|F|(|F|+1)/2  +  (K+1)·|F| / p(K, t)
```

table entries plus rests scanned: `K` relations for rank and one descent, a
target costing `|F|` rests.  `model.py` writes this down.

**Checked, not fitted.**  Against the triple arm's committed trial counts
(`confirm.json`, 32 fixtures a cell), the count predicts mean trials of 5.9,
37.1 and 396.1 at `n23a1`, `n37a0` and `n43a1`; the arm measured 6.1, 36.3 and
387.6.  Nothing in `U` was adjusted to get there.

**Where it is optimistic, and in which direction.**  The count treats distinct
multisets as distinct points.  Exact enumeration on a real `n23a1` base
(`table_completeness.py`, the checker's own curve arithmetic,
`table-completeness-n23a1.json`) finds 1,988 distinct three-sum orbits where
the multiset count gives 2,538, and the arm's own table stores exactly 1,988.
(A 1,989th is reachable only through the row `C_0 + (−C_0)`, whose "triples" are
bare base points; the arm skips it by design.)  So **the table is complete** and
the shortfall is sums that coincide, which at `n37a0` and `n43a1`
(the arm's reported `sum_orbits` against the multiset count) is 12% and 10%.
The count therefore overstates coverage somewhat at two orbits.  Coincidences
should be rarer at three, where more sums cross orbits; if so the error favours
a lower ratio below.  That is expected, not verified.

## What changes

`U` chooses the base (`model.py`, output in `model-output.txt`):

| cell | `r` | W₃ chooses | **the count chooses** | `U` ratio | whole job, counted/triple |
|---|--:|--:|--:|--:|--:|
| `n23a1` | 4,196,903 | K=2, t=1 | K=2, t=1 | 1 | identity |
| `n37a0` | 230,603,167 | K=2, t=1 | K=2, t=1 | 1 | identity |
| `n43a1` | 4,644,189,029 | K=2, t=1 | **K=3, t=1** | 0.791 | **0.795–0.890** |
| `n59a0` | 10,063,074,221 | K=2, t=1 | K=2, t=1 | 1 | identity |
| `n61a1` | 11,514,943,771 | K=2, t=1 | K=2, t=1 | 1 | identity |

**One cell changes.**  At `n43a1` two orbits give a target a 0.76% chance of
being a four-sum, so the arm draws about 400 targets and spends four-fifths of
its work scanning them.  Three orbits cost more table (33,411 entries against
14,878) and cut the targets to about 125.  At the holdouts the Frobenius orbits
are larger (118 and 122 points), so two orbits already reach and the count
keeps them.

The whole-job range applies `U` at the triple arm's own measured 1,132
instructions a unit, plus 2.03M of fixed phases that do not depend on the base
(`phase-profile-probe.json`).  The low end is `K+1 = 4` hits.  The high end
allows one extra relation for rank, which at three orbits should be needed
about one time in twenty.

**The sizing search itself also gets cheaper.**  W₃'s search evaluates every
`t ≤ K ≤ 64`, twice over.  The profile puts `triple_rows` alone at 120,528
instructions.  The counted search stops once one row's table exceeds the best
total so far, which is a few dozen evaluations.  That difference is fixed, so
it shows at the identity cells as a ratio just under one, largest at `n23a1`,
whose whole job is only 4.0M.

## The arm

`counted-sizing.patch`, applied on top of `triple-table.patch`:

- `Collector::TripleCounted`: the triple collector exactly, with `(K, t)` from
  `U` in place of W₃.  Same sampler, same table, same scan, same witnesses.
  The search is a fixed function of `(n, r)`, decided before any target is drawn.
- The worker routes `solver: "triple_counted"` at `summands: 4` to it, and only
  that pairing.
- `Collector::Pair` and `Collector::Triple` are unchanged, so every committed
  measurement of either still describes its code.

## Measurement

**Arms:**

- `triple`, the committed arm;
- `counted`, the new one;
- `rho`, the worker's matched rho (`automorphism_order = 2n`), run from the
  `triple` build. Neither patch touches the rho path.

**Cells:** `n23a1`, `n37a0`, `n43a1` and the round-0023 holdouts `n59a0` and
`n61a1`. Round 0023 chose the holdouts because rho finishes there at
`max_trials` 65,536. Round 0023's configuration otherwise.

**Unit:** callgrind `Collected:`.

**Design:** paired. Every fixture runs on all three arms with the same target
and algorithm seeds, and every report goes through `oracle.py` (summands 4 for
the IC arms, rho mode for rho).

- **Probe** (does not count): 4 fixtures a cell from `random.Random(20260925)`.
- **Confirmation:** 32 fixtures a cell from `random.Random(20260926)`, drawn as
  `(target_seed, algorithm_seed)` pairs of `getrandbits(64)`, cells in the
  order above. Both streams are disjoint from the triple study's 20260923 and
  20260924.

**Statistic:** per cell, the geometric mean of paired ratios, with a 95%
percentile bootstrap interval (10,000 resamples, seed 0) and the median beside
it. This is `check.py`'s `summarise()` from the triple study, unchanged.

## Decision rules, fixed now

**Gate, absolute.** Every report from every arm verifies, and all three arms
recover the same logarithm on every fixture. Any failure is a defect, not a
measurement. A run that does not finish, in any arm, removes that fixture's
cell from every comparison below, and the cell is reported as having none. An
unfinished run is never evidence.

**1. Identity at the four unchanged cells** (checks the patch).

- At `n23a1`, `n37a0`, `n59a0` and `n61a1`, `counted` and `triple` report the
  same factor-base orbits, trials, columns, accepted relations and logarithms
  on every fixture.
- Their instruction ratio lies in `[1 − 150,000/Ī − 0.005, 1.005]`, where `Ī`
  is the triple arm's mean at that cell. That puts the band at about
  `[0.957, 1.005]` at `n23a1`, `[0.986, 1.005]` at `n37a0`, and tighter at the
  holdouts.
- Failing either is a defect in the patch, and it voids item 2 until explained.

**2. The one change: `counted`/`triple` at `n43a1`.** `counted` chooses K=3,
t=1, and the predicted ratio is 0.795–0.890.

- **Supported:** geometric mean ≤ 0.90, with the interval's upper end below one.
- **Refuted:** geometric mean ≥ 0.95. That would mean the count's advantage at
  three orbits does not survive the per-entry cost of a larger table, and the
  note will say by how much.
- **Otherwise:** partial.

**3. End to end against rho — registered this time.** The last study measured
this only as an unregistered secondary.

| cell | predicted `counted`/rho | registered |
|---|--:|---|
| `n23a1` | ≈ 1.03 (the triple arm's, unchanged) | not predicted below one |
| `n37a0` | ≈ 0.81 (unchanged arm, new fixtures) | below one |
| `n43a1` | 0.74–0.83 | below one |
| `n59a0` | ≈ 0.77 | below one |
| `n61a1` | ≈ 0.76 | below one |

- **Per cell:** *below rho* if the interval's upper end is below one; *above
  rho* if its lower end is above one; otherwise *not resolved*.
- **Supported:** below rho at all four of `n37a0`, `n43a1`, `n59a0` and
  `n61a1`.
- **Refuted:** above rho at any of them.
- **Otherwise:** partial, cell by cell.

The holdout predictions carry `n43a1`'s measured triple/rho (0.933) by the
ratio of `U` and by rho's expected steps, `√(r/n)`. That assumes instructions
per unit and per rho step grow alike with `n`, since both use the same field
arithmetic and the same canonical form. It is the least-grounded number here,
and it is registered as a direction, with the point value for reference.

## Scope, fixed now

- **Instructions only.** The tournament also gates on native time, and none of
  its stages is run here.
- **One machine, one configuration, five cells.**
- **The arm has four summands.** The `OPERATIONS.md` summand-matching question
  raised on #668 is unchanged, and it is the tournament lane's call.
- **Class: engineering.** A sizing rule with a derivation. Item 3, if
  supported, is a measurement: this pipeline used fewer instructions than
  matched rho at the cells named, and nothing beyond them. It does not change
  the verdict at 131 bits: this collector family's IC/rho still grows as
  `r^{1/10}`.
- **Not in this arm:** a pair path for `n23a1`. That is the switch named on
  #668, and it needs its own design, because the checker fixes the summand
  count per job.

## Addendum, before the confirmation run (after a 4-fixture probe)

Added after `counted-sizing.patch` was built and tested (12 of 12,
`unit-tests.txt`) and a 4-fixture-a-cell probe from `random.Random(20260925)`
had been seen (`probe-does-not-count.json`), and before any confirmation seed
was run.  **The probe does not count.**

**What it showed.**

- Every report verified, with the same logarithm on all three arms.
- The identity held on every fixture at `n23a1`, `n37a0`, `n59a0` and `n61a1`:
  the same base, trials, columns and relations, and ratios of 0.968, 0.990,
  0.999 and 0.997, each inside its band.
- At `n43a1`, four ratios of 0.74, 1.17, 0.74 and 1.60, with trials of 67,
  110, 49 and 147 against the triple arm's 339, 152, 289 and 161.  That is a
  geometric mean of 1.01 [0.75, 1.37], which says nothing at four.

**Two things fixed now.  Both are quantified from `model.py` alone
(`model_gm.py`, `model-gm-output.txt`), not from the probe.**

1. **The band was computed on means, but the rule scores a geometric mean.**
   The two arms draw different bases, so a fixture's two trial counts are
   independent, and the triple arm's count is the heavier-tailed one.
   Simulating trials as the count says they are, the ratio of means
   reproduces the registered 0.795 and the geometric mean of ratios is
   **0.859**.  That point sits inside the registered band, so **the band, the
   point range and the rules stand unchanged**.
2. **At 32 fixtures the rule cannot reliably pass even if the prediction is
   exactly right.** The spread of the log ratio is 0.51, so the expected
   interval around 0.859 is [0.72, 1.02]. Its upper end would sit above one
   about half the time. **`n43a1` is extended to 128 fixtures**, expected
   interval [0.79, 0.94].
   - The registered 32 come from `random.Random(20260926)` exactly as
     registered, in the registered order. The other 96 come from a separate
     stream, `random.Random(20260927)`, appended after all five cells
     (`check.py --extend n43a1:96:20260927`).
   - Items 2 and 3 at `n43a1` are scored on all 128. The registered 32 alone
     are reported beside them, and nothing is decided on them.
   - Every other cell stays at 32.

**Noted, not changed.** The probe's `counted`/rho at `n61a1` (0.40) is far
below the registered ≈ 0.76. The triple arm reads the same there (0.40),
since the two arms are identical at that cell, so this is about rho's cost at
`n61a1`, not about the sizing. The holdout point values were registered as
the least-grounded numbers here, and item 3 is scored on direction only.
