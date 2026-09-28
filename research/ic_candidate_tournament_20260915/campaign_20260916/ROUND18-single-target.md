# Round 0018 pre-registration: the scan block is oversized and the cofactor is counted twice

Round 0017 restored `beats_rho_strict` on all eight cells by shrinking the
certificate. It left the n23a1 collection untouched and said so, and the
campaign's next-step record named that collection as round 0018's target.
This round went there, found the hypothesis waiting for it was wrong, and
took what the measurement pointed at instead.

## What was killed first, and what it cost

Round 0017 measured its `orbits_rows` arm — the pair-table row rule
re-weighted at the build/scan ratio round 0015 measured, 680/380 — and found
it badly wrong: it moved n23a1 from three rows to two, predicting 0.973x, and
two rows cost **1.272x** three. The natural reading is that the shipped rule
is wrong in the same direction and that the real optimum is elsewhere.

It is not. `campaign_20260916/round18_row_sweep.py`, with
`round18-rows-probe.patch` making the row count settable and the probe first
checked byte-identical to the frozen `orbits` worker, measures every row count
at every cell, three fixtures each:

| cell | K | \|F\| | shipped | measured best | best / shipped |
|---|---:|---:|---:|---:|---:|
| n13a0 | 7 | 182 | 1 | 1 | 1.0000 |
| n17a1 | 8 | 272 | 1 | 1 | 1.0000 |
| n19a0 | 8 | 304 | 1 | 1 | 1.0000 |
| n19a1 | 8 | 304 | 1 | 1 | 1.0000 |
| n23a0 | 8 | 368 | 2 | 3 | 0.9899 |
| **n23a1** | 8 | 368 | **3** | **3** | **1.0000** |
| n29a1 | 8 | 464 | 1 | 1 | 1.0000 |
| n31a0 | 8 | 496 | 1 | 1 | 1.0000 |

The shipped rule sits at the measured optimum at seven of eight cells, n23a1
among them. The one miss is n23a0, where three rows beat two by 1.0%, inside
the +-4-5% per-cell spread the A/A control shows. The sweep also reproduces
round 0017's number from the other side: at n23a1 two rows cost 1.307x three
here against 1.272x there. **The row rule is not the lever**, and the round
that would have been built on it was not built. Cost of finding out: one
probe build and about fifteen minutes.

## What the sweep exposed instead

`decompose` scans the base in blocks. Each block pays one batch inversion and
is discarded on an early exit, so the rule sizes a block at about one expected
witness: `chunk = (1/hit + 1).clamp(8, 64).min(|F|)` with
`hit = lambda * cov(t)`. That ceiling of 64 binds at exactly two cells of the
panel, the two with the largest subgroups — the rule asks for 71 at n23a0 and
**102 at n23a1** — and nowhere else.

The obvious reading is that the cap is too low. Measurement says the opposite
(`round18_block_sweep.py`): at n23a1 the value the rule asks for costs
**1.040x** what 16 costs, and the cap has been quietly protecting the code
from its own rule. One expected witness per block is the wrong target — the
inversion a block amortises saturates after a dozen or so rests, while the
work discarded on an early exit keeps growing with the block.

A single fixed block is nevertheless a wash. `round18_block_grid.py` sets the
block outright — which the ceiling sweep cannot do below the rule's own value,
so it could only ever measure 8 at the three cells where the rule already asks
for the floor — over a grid at every cell, six fixtures each. The best constant
on the six development cells is 12, at 0.9910, and it makes n13a0 (1.0166) and
n29a1 (1.0088) worse. What does work is leaving the rule's shape alone and
capping it lower: at **16** every cell whose value is already at or below 16 is
untouched and every capped cell improves.

| cell | rule asks | block 16 | block 24 | block 32 |
|---|---:|---:|---:|---:|
| n13a0 | 8 | unchanged | unchanged | unchanged |
| n17a1 | 8 | unchanged | unchanged | unchanged |
| n19a0 | 13 | unchanged | unchanged | unchanged |
| n19a1 | 25 | 0.9833 | 0.9990 | 1.0100 |
| n23a0 | 71 (capped 64) | 0.9627 | 0.9654 | 0.9691 |
| n23a1 | 102 (capped 64) | 0.9808 | 0.9788 | 0.9798 |
| n29a1 | 8 | unchanged | unchanged | unchanged |
| n31a0 | 50 | 0.9723 | 0.9762 | 0.9827 |

The ceiling was chosen on the **development cells only**
(`campaign_20260916/round18_ceiling_choice.py`, which reads the grid's saved
raw rows and scores each candidate ceiling by the block it would give each
cell, `min(rule value, ceiling)`):

| ceiling | 8 | 12 | **16** | 24 | 32 | 48 | 64 |
|---|---:|---:|---:|---:|---:|---:|---:|
| gm, development cells | 0.9911 | 0.9869 | **0.9859** | 0.9866 | 0.9885 | 0.9941 | 1.0000 |

16 is the best of them, and **it is barely the best**: 12, 16 and 24 sit
within 0.1% of each other, so this is a choice among three near-equivalent
values rather than a sharp optimum, and the round should not be read as having
located one. What is not marginal is the direction — anything at or below 32
beats the shipped 64 — and the property that decided it: under a ceiling of 16
no cell regresses, because a ceiling can only ever lower a block and every cell
already at or below 16 is untouched. The holdout columns above are what the
tournament will test; they took no part in the choice.

## The second lever: the cofactor is applied twice

`TinyIc::new` samples a curve point, projects it into the order-`r` subgroup
with `[h]`, and closes the abscissa under Frobenius. Every orbit
representative `rep` is therefore **already in the subgroup**. The column is
then formed as `[h]rep` — a second cofactor multiplication of a point that is
already where it needs to be. It scales every column by the constant `h`, and
the row's right-hand side (`h*a`) and the descent (`invmod(h*b)`,
`sum - h*a`) each carry a matching `h` to compensate.

Dropping all four leaves exactly the same logarithms and saves `K` scalar
multiplications by the cofactor per job. That is nothing where `h` is 2 or 4
and most of the base phase where it is not: `h` is 12,646 at n29a1 and 1,492
at n31a0, which are the panel's two large-cofactor cells.

Round 0017 deferred this deliberately, because the checker hard-codes `[h]P`
as a base point's column and `h*a` as its row. That is the amendment below.

## The checker amendment

`campaign_20260916/round18-oracle-convention.patch`, applied to `oracle.py`,
is additive in the same sense round 0017's was. A report may declare
`"column_convention": "representative"`, and then a base point is located by
**itself** and a relation's row reads `sum coeff*log == a`. A report without
the key is read exactly as before — `[h]P` and `h*a` — and an unknown value is
refused. The checker recomputes every point and every row with its own group
arithmetic under either convention. Under `representative` every base point
must carry a column, where the old convention excused a point whose `[h]`
image is the identity; the new reading is therefore strictly stronger, not
weaker.

The amendment was admitted only after it was attacked. At n29a1, where
`h = 12,646` and so the two conventions differ enormously:

| attack | checker |
|---|---|
| a `representative` report with the label stripped | refused (uncovered base column) |
| a legacy report relabelled `representative` | refused (uncovered base column) |
| a `representative` report relabelled `cofactor` | refused (uncovered base column) |
| an unknown convention string | refused (unknown column convention) |
| a column logarithm moved by one | refused (bad column log) |
| a relation scalar moved by one | refused (incorrect point relation) |
| a recovered logarithm moved by one | refused (incorrect scalar) |
| an orbit representative's coordinates swapped | refused (point does not lift to curve) |
| the two reports as emitted | both accepted |

and after a regression check that it changes nothing for what already exists:
280 frozen round-0016 and round-0017 profiles sampled at random still verify
unchanged.

## Objective, gates, arms, seed, budget

Objective `rho`, unchanged. Promotion over the incumbent requires the
instruction ratio at or below 0.98 with a 95% upper limit below 1, a native
95% upper limit below 1, and every cell within 1.10, on confirmation and on
replay. `beats_rho_strict` requires the winner over rho below 1 in both
metrics in every cell, on both stages. No gate, margin, bootstrap, repetition
count, confidence level or panel changes in this round. The one protocol
change is the checker's second accepted column convention, stated above.

Parent: round 0017 (`runs/round-0017`). Four arms plus the implicit rho
reference:

- `incumbent` — the round-0017 winner `orbits`, exactly as sealed.
- `block` — `orbits` + `round18-block.patch`: the clamp ceiling 64 -> 16.
- `column` — `orbits` + `round18-column.patch`: the column is the
  representative; declares `column_convention: representative`.
- `both` — both patches.

Panel: the round-0016/0017 eight cells, `13a0, 17a1, 19a0, 23a0, 23a1, 31a0`
with `19a1, 29a1` held out to confirmation and replay. Seed 2026091818, pilot
profile, one target per job, CPU 3, 60 s timeout, `--max-processes 3000`.
Trials, counted from `stage_arms` with four IC arms and rho: 36 aa, 90 smoke
(6 cells x 3 repetitions x 5 arms), 270 development (18 cases x 3 x 5), 216
selection (the two best challengers, the incumbent and rho), 864 confirmation
and 864 replay (incumbent, the provisional challenger and rho) — **2,340
trials**, inside the 3,000 budget.

## Boundary, floor, class, honesty

Unchanged from round 0017. The unit is Valgrind 3.22 `Ir` over the complete
worker process; native timings are paired cold-process wall under the same
blocking-reap protocol and carry their own scope. The class is accounting, not
a complexity claim, and nothing here is a statement about any curve outside
the eight cells named above.

## Development measurement

Every arm verified against the oracle on all 96 round-0017 confirmation
fixtures, each returning the incumbent's own solutions and base hash, with
`block` declaring the legacy convention and `column` and `both` declaring
`representative`. Ir by callgrind, wall the median of five interleaved cold
runs, four fixtures a cell (`campaign_20260916/round18_measure.py`):

| cell | block/inc | column/inc | both/inc | inc/rho |
|---|---:|---:|---:|---:|
| n13a0 | 1.0005 | 0.9755 | 0.9755 | 0.7154 |
| n17a1 | 1.0004 | 0.9928 | 0.9928 | 0.6050 |
| n19a0 | 1.0003 | 0.9897 | 0.9897 | 0.6856 |
| n19a1 (holdout) | 0.9826 | 0.9896 | 0.9720 | 0.6585 |
| n23a0 | 0.9575 | 0.9927 | 0.9501 | 0.6649 |
| n23a1 | 0.9820 | 0.9949 | 0.9767 | 0.7561 |
| n29a1 (holdout) | 1.0002 | **0.9227** | **0.9227** | 0.8588 |
| n31a0 | 0.9707 | 0.9499 | **0.9205** | 0.5786 |
| **gm, development cells** | 0.9851 | 0.9825 | **0.9672** | |
| **gm, all eight** | 0.9867 | 0.9757 | **0.9621** | |

**The native side of that script does not resolve this round, and the round
says so before running rather than after.** `block` is a no-op by
construction at n13a0, n17a1, n19a0 and n29a1 — the rule's own value there is
already at or below 16, so the executable cannot behave differently — and the
instruction column confirms it to four decimals (1.0005, 1.0004, 1.0003,
1.0002). The same script's native column at those four cells reads 1.0427,
1.0190, 1.0283 and 1.0438. That is the native instrument's floor at five
repetitions, not an effect, and it is the size of everything this round is
trying to measure. Native is left to the tournament's own protocol.

That accident is also this round's best control. **The four cells where
`block` is a no-op are an A/A control carried inside a live arm**, and
whatever `block`/`incumbent` reads there at confirmation — in either metric —
is the floor against which every other number in this round should be read.

## Predictions, in the order they will be checked

1. **Every arm answers exactly what the incumbent answers.** On every
   confirmation and replay fixture, all three arms return the incumbent's
   solutions and factor-base hash. Any difference is a bug and ends the round.
2. **`block` is a no-op at n13a0, n17a1, n19a0 and n29a1.** Its instruction
   ratio to the incumbent is 1.000 to within 0.001 at each, and its relation
   count and trial count are identical. If this fails, the patch is not what
   it claims to be.
3. **`both` beats the incumbent on instructions panel-wide**, point estimate
   between 0.95 and 0.975, with the largest gains at n31a0 and n29a1 (both
   near 0.92) and the smallest at n17a1 (above 0.98).
4. **The `column` arm's gain is ordered by the cofactor.** Ranking the eight
   per-cell `column`/`incumbent` instruction ratios, n29a1 (h = 12,646) and
   n31a0 (h = 1,492) are the two largest gains, and the six cells with
   h in {2, 4} are the six smallest. A cell out of order falsifies the
   mechanism as stated even if the arm still wins.
5. **The holdout cells behave as the development cells chose.** `block` at
   n19a1 lands near 0.983 and is a no-op at n29a1; neither holdout was allowed
   to choose the ceiling.
6. **Native is the binding question, and it may refuse the round.** The
   instruction gain is 3.8% panel-wide and the native gate needs a 95% upper
   limit below 1. Prediction: `both` passes the instruction gate
   (<= 0.98, upper < 1) on both stages, and its native point estimate against
   the incumbent is below 1. Whether the native *upper limit* clears 1 is
   genuinely open, and if it does not, the arm is not promoted and this
   pre-registration will have said so first.
7. **`beats_rho_strict` still holds on all eight cells**, for whichever arm is
   carried to confirmation, in both metrics and on both stages.

Prediction 6 is the one most likely to fail, and failing it is not a failure
of the round: a refused promotion with a measured 3.8% instruction gain and a
native instrument whose own floor is +-4% is a statement about the instrument,
which is worth recording either way.
