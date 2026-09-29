# Results: counted sizing, and what "below rho" was measured against

Two pre-registered measurements on the same 256 fixtures, instructions only
(callgrind `Collected:`, valgrind 3.22.0).  All three arms' reports and the
control's verified in `oracle.py`, with the same logarithm on every fixture,
and no run went unfinished.

| | class | outcome |
|---|---|---|
| Counted sizing (`PREREGISTRATION.md`) | engineering | **supported**: 0.773× the triple arm at `n43a1`, an identity at the other four cells |
| End to end against the worker's rho (same, item 3) | measurement | **supported as registered**: below that rho at `n37a0`, `n43a1`, `n59a0` and `n61a1` |
| Rho with the IC arm's own canonical form (`PREREGISTRATION-rho-control.md`) | **accounting** | **the IC arm is above it at all five cells**, 1.19× to 2.88×.  By the rule's letter this is *partial*: two cells are unscored, and both read the same way |

**The headline is the third row.**  The worker's rho names each Frobenius class
by walking the orbit with field squarings, about 9,100 of its 9,930
instructions a step at `n43a1`.  The IC arm has named the same classes with a
normal-basis rotation since round 0007.  Give rho that one piece of the IC
arm's own code, and every cell where this pipeline "beat rho" becomes a cell
where it loses.

## 1. Counted sizing

`counted-sizing.patch` sizes the triple collector by counting the four-sums a
target can be, in place of W₃ (`PREREGISTRATION.md`).  The count predicted the
triple arm's committed mean trials without fitting (5.9, 37.1 and 396.1
against 6.1, 36.3 and 387.6).  It chooses three orbits at `n43a1` and two
everywhere else.

`confirm.json`: 32 paired fixtures a cell from `random.Random(20260926)`, and
at `n43a1` 96 more from `random.Random(20260927)`.  The addendum fixed that
extension before any confirmation seed ran.  Geometric mean of the paired
ratios, 95% bootstrap interval.

| cell | counted/triple | registered | |
|---|--:|--:|---|
| `n23a1` | 0.968 [0.968, 0.969] | identity, [0.959, 1.005] | **met**: every fixture inside (0.966–0.973) |
| `n37a0` | 0.991 [0.991, 0.991] | identity, [0.986, 1.005] | **met** (0.989–0.993) |
| `n43a1`, 128 fixtures | **0.773** [0.708, 0.843] | ≤ 0.90 with the upper end below one | **supported** |
| `n43a1`, the registered 32 alone | 0.742 [0.623, 0.881] | reported, not decided on | |
| `n59a0` | 0.999 [0.999, 0.999] | identity, [0.995, 1.005] | **met** (0.9985–0.9995) |
| `n61a1` | 0.997 [0.997, 0.997] | identity, [0.994, 1.005] | **met** (0.9965–0.9969) |

At the four identity cells the two arms reported the same base, trials,
columns, relations and logarithm on every fixture.  What remains of the
ratio is the cheaper sizing search, as registered.

At `n43a1`:

- **The whole-job ratio is below its predicted range.**  It measured 0.773,
  against 0.795–0.890, or 0.859 under the geometric mean.
- **Trials fell from 337.8 to 111.5** a fixture, against the count's 396.1 and
  124.9.
- **Three of 128 fixtures needed a relation beyond K**; about six were
  expected.

The count overstates coverage at two orbits, where sums coincide, and less at
three (`PREREGISTRATION.md`). The ratio falling below its prediction is the
direction that bias implies, but this run does not isolate that as the cause.

## 2. End to end against the worker's rho, as registered

| cell | counted/rho | registered | |
|---|--:|---|---|
| `n23a1` | 1.210 [1.121, 1.303] | not predicted below one; point ≈ 1.03 | above rho; **the point missed** |
| `n37a0` | 0.795 [0.679, 0.935] | below one | below |
| `n43a1`, 128 fixtures | 0.824 [0.752, 0.904] | below one, 0.74–0.83 | below |
| `n59a0` | 0.823 [0.763, 0.889] | below one, ≈ 0.77 | below |
| `n61a1` | 0.714 [0.580, 0.874] | below one, ≈ 0.76 | below |

**Supported, as registered.**  The same fixtures give triple/rho of 1.249 at
`n23a1`, where the last study read 1.029 on other fixtures. So fixture
variation at that cell is wider than either interval shows.

**§3 is why this table does not mean what it appears to.**

## 3. The control: rho with the IC arm's canonical form

### What was found, and when

The profile was taken while §2's confirmation was running; the control was
pre-registered after it finished and before it existed as code
(`PREREGISTRATION-rho-control.md`, commit `e11b59a`). Where rho's
instructions go at `n43a1` (`rho-lines-n43a1.txt`: line tables, callgrind self
cost by source line over one exposed fixture's `rho_solve`):

- 54.7% in field squaring and reduction;
- 15.0% in the reduction-table index;
- about 21% in the orbit loop's compare and range;
- under 1% in the point addition.

That is `FastRhoWalk::canonicalize`, which runs `x` through all `n` Frobenius
images by squaring and then squares `y` to the winning power.

### The control

`rho-normal-basis.patch` replaces that loop, and nothing else:

- It calls the IC arm's own `NormalBasis::canon`, now `pub(crate)` and built
  from the rho field's squaring.
- It recovers `φ^k(x)` and `φ^k(y)` with an inverse basis change.
- It keeps rho's sign rule and coefficient scaling as they were.

**Tests** (`rho-control-unit-tests.txt`):

- The IC module passes 12 of 12, including a new test of the inverse basis
  change and of the rotation direction.
- The library's rho tests pass 39 of 40. They include log recovery, step
  counts against the birthday bound, fruitless-cycle escape and closed charge
  ledgers.
- **The one failure is `single_word_rho_walk_matches_the_reference_step_for_step`, by construction.**
  That test pins the fast walk's representative to the reference walk's, and
  the control names a different member of the same class. It found and
  verified the right logarithm in a different number of steps. The control is
  a patch file for measurement, and the library is unchanged.

### Result

`rho-control.json`: the control once on every §1 fixture.

| cell | fixtures | **counted/control** | predicted | control/rho | walk additions, control/rho |
|---|--:|--:|--:|--:|--:|
| `n23a1` | 32 | **1.501** [1.433, 1.574] | ≈ 1.85 (1.2–3) | 0.806 | 1.251 [0.947, 1.649] |
| `n37a0` | 32 | **1.723** [1.569, 1.892] | ≈ 3.1 (2–5) | 0.461 | 0.989 [0.822, 1.208] |
| `n43a1` | 128 | **2.882** [2.680, 3.095] | ≈ 5.2 (3–8) | 0.286 | 0.979 [0.852, 1.122] |
| `n59a0` | 32 | **1.192** [1.148, 1.242] | ≈ 1.25 (1.05–1.6) | 0.690 | 1.056 [0.810, 1.392] |
| `n61a1` | 32 | **2.736** [2.306, 3.258] | ≈ 3.9 (2.5–6) | 0.261 | 1.178 [0.947, 1.464] |

**Scored exactly as registered:**

- **Item 1, walk lengths in [0.85, 1.15]:** met at `n37a0`, `n43a1` and
  `n59a0`. **Missed at `n23a1` (1.251) and `n61a1` (1.178)**, so the rule
  leaves item 2 unscored at those two cells.
- **Item 2:** above one at all three scored cells. By the rule's letter that
  is **partial**, since "supported" needed all five.

What the two misses are and are not:

- **Sampling spread, by every check available.** All five bootstrap
  intervals on the walk ratio contain one (`walk_lengths.py`, added after the
  run). At both missed cells the control walked closer to the birthday-bound
  expectation than the worker's rho did: 402 and 12,994 against expectations
  of 379 and 12,176, where rho walked 322 and 11,033. A [0.85, 1.15] band was
  about ±1.2 standard errors at 32 fixtures, too tight for that sample. That
  was my mistake in the registration.
- **The caveat cuts in the IC arm's favour, not against it.** Where the
  control walked longer, it spent more instructions, which lowers
  counted/control. The two unscored cells read 1.50 and 2.74 despite that.

**The magnitude was over-predicted, at `n37a0` and `n43a1` most.**

- The control costs about 3,000 instructions a walk addition at `n43a1`
  (`rho-control-phases-probe.json`), against the old walk's 9,930 and an IC
  probe's roughly 1,100. The prediction assumed about 1,400.
- So **the control is still not a lean rho.** It costs about twice what the
  prediction assumed a step, and this study did not profile where the
  difference goes. A rho as lean per step as the IC arm's probe would widen
  every ratio in the table.

### What this changes

- **Every IC-below-rho result compared against the worker's rho is a
  comparison with a baseline carrying a canonicalisation the IC arm had
  engineered away.** That covers §2 of this note, the secondary rho table of
  `research/ic_triple_table_20260923`, and the IC tournament's
  `beats_rho_strict` rounds 0007–0020, each read at the cells and seeds it was
  scored on. Against a rho given the IC arm's own orbit naming, the triple
  pipeline with counted sizing loses at every cell measured here, by 1.19× to
  2.88×.
- **Class: accounting** (`AGENTS.md` §3). No algorithm got worse. What the
  ratios were counted against changed.
- **It says nothing against index calculus in general**, and nothing new
  about ECC2K-130, where the verdict already stood against IC.
- **What it does say** is that this lane has not yet shown an IC pipeline
  beating a matched rho at any cell. The cheapest way to keep that honest is
  for the tournament's rho to carry the same orbit naming as its IC arm.

## Related: #650

Another lane's open draft, #650 ("Matched rho references"), builds matched rho
references for the scoreboard's `ic bench` and `ic boundary` figures, and
re-prices those figures against them.  Its Koblitz classes are named the same
way as this control's, by least normal-basis rotation.  The two were found
independently, and they agree on the direction.

#650 does not touch the tournament worker's `koblitz_signed_frobenius_rho`,
which is the rho every IC/rho ratio in this note and in the tournament was
taken against.  This control is the corresponding check for that code path.

## Scope

- **Instructions only, one container.** Native time was not measured.
- **Five cells, the fixtures named above, one configuration.**
- **The rho control changes one function** and is not offered as the best
  rho. Other asymmetries may remain, in either direction.
- **Nothing here is a tournament round.** Whether the tournament's matched rho
  should change is that lane's decision.

## Reproducing

```sh
research/ic_triple_counted_20260923/build_arms.sh /abs/work   # triple, counted and the rho control
python3 research/ic_triple_counted_20260923/check.py TRIPLE COUNTED \
    --seed 20260926 --fixtures 32 --extend n43a1:96:20260927 --out confirm.json
python3 research/ic_triple_counted_20260923/rho_control.py TRIPLE CONTROL --out rho-control.json
python3 research/ic_triple_counted_20260923/walk_lengths.py
python3 research/ic_triple_counted_20260923/model.py            # the count and its predictions
python3 research/ic_triple_counted_20260923/model_gm.py         # the addendum's statistic and sample size
```
