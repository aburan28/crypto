# Round 0016 pre-registration: how many curves the gate needs, and the two cells nobody has run

Round 0015 produced the largest IC-only instruction win the campaign has
measured — `scan_io` at 0.9059 of the incumbent on confirmation, 0.6679 of rho —
and the native gate refused it, because the paired native interval was
[0.9785, 1.0073] about a point estimate of 0.9925. That was the fifth round in a
row to fail the same way. Before spending another round on the code, this round
asks what the gate can resolve at all, and answers it with arithmetic that is
written down here, before anything is measured.

## What the round-0015 receipts already say

Every number in this section comes from the frozen round-0015 confirmation
receipts (`runs/round-0015/runs/confirmation/**/receipt.json`, 540 trials) and
from re-running the evaluator's own nested bootstrap — `comparison()` in the
frozen `runs/round-0015/evaluator/tournament.py` — over them. The reproduction
of the published interval is exact: [0.97852, 1.00727] against the recorded
[0.9785179757065511, 1.0072720504423012].

That bootstrap has two levels. It resamples the C curve cells with replacement,
and then resamples the fixtures inside each chosen cell. Only the inner level
shrinks when a round buys more fixtures. Holding round-0015's data fixed and
varying only the number of fixtures per cell gives the gate's resolution budget
for `scan_io` against the incumbent on native wall:

| fixtures per cell | half-width | upper limit |
|---:|---:|---:|
| 3 | 0.02570 | 1.01927 |
| 6 | 0.02010 | 1.01338 |
| 12 (round 0015) | 0.01438 | 1.00727 |
| 30 | 0.01090 | 1.00373 |
| 100 (`--profile standard`) | 0.00828 | 1.00069 |
| 400 | 0.00759 | 0.99987 |
| infinite | 0.00719 | 0.99952 |

The standard profile is 8.3x the jobs of the pilot and it does not promote this
challenger. Nothing does: with five cells the upper limit converges to 0.99952,
so the whole remaining margin at infinite replication is four parts in ten
thousand. **Replication is not the lever.** That disposes of the plan this round
was going to be.

Neither is the measurement floor, which is the other thing five rounds have
blamed. A worker instrumented with cross-process `CLOCK_REALTIME` marks
(`campaign_20260916/round16-floor-probe.patch`, applied to the frozen round-0015
incumbent source, driven by `round16_floor_probe.py`) splits a measured wall
into two spans that really are serial -- the time inside the worker's `main`,
and everything else -- and then splits the first by the worker's own clock. It
does not split the second, because the child execs and can reach `main` while
the evaluator is still inside Python: those spans overlap and do not partition
the wall. An earlier version of this table presented them as if they did and
produced a negative segment, which is how the overlap was found.

| cell | arm | wall | in `main` | outside | algorithm | algorithm share |
|---|---|---:|---:|---:|---:|---:|
| n13a0 | IC | 1057.3 | 363.8 | 693.4 | 305.4 | 0.289 |
| n13a0 | rho | 899.0 | 386.1 | 512.8 | 345.7 | 0.385 |
| n17a1 | IC | 882.4 | 330.7 | 551.7 | 295.9 | 0.335 |
| n17a1 | rho | 960.9 | 473.1 | 487.8 | 435.1 | 0.453 |
| n19a0 | IC | 970.7 | 397.6 | 573.2 | 361.1 | 0.372 |
| n19a0 | rho | 1071.9 | 564.9 | 507.0 | 524.2 | 0.489 |
| n19a1 | IC | 953.8 | 410.6 | 543.2 | 374.6 | 0.393 |
| n19a1 | rho | 1041.7 | 552.4 | 489.2 | 507.4 | 0.487 |
| n23a0 | IC | 1179.0 | 627.6 | 551.3 | 587.3 | 0.498 |
| n23a0 | rho | 1160.9 | 678.3 | 482.6 | 639.6 | 0.551 |

All times are microseconds, medians over 12 fixtures and three repetitions.
Between 29% and 55% of what the tournament calls a native measurement is the
algorithm; the rest is process creation, exec and page-in, stdin, the report
write, exit and reap, plus the evaluator's own Python.

Two things follow, and the second was a surprise.

First, **the remainder is not arm-independent.** The IC arm carries 60 to 180 us
more of it than rho in every cell of the probe, and the frozen round-0015
receipts show the same sign and ordering without any instrumentation at all:
subtracting each receipt's in-worker `elapsed_seconds` from its recorded
`process_wall_seconds` gives 684.4 / 676.4 / 716.2 / 708.9 / 718.9 us for the
incumbent against 656.9 / 631.9 / 651.4 / 651.4 / 662.1 us for rho, cell by
cell. The IC arm emits a factor base, a relation list and a column-log table
where rho emits one scalar, and it dirties far more arena memory retiring them.
That is a real cost of producing the larger certificate and the process wall is
right to charge it -- but it means a native ratio is not an algorithm ratio, and
no round may read it as one. This round records the effect rather than removing
it.

Second, removing floor would not help the challenger anyway. Subtracting a
constant from round-0015's walls and re-running the same bootstrap:

| floor removed | `scan_io`/incumbent point | upper | incumbent/rho point | upper |
|---:|---:|---:|---:|---:|
| 0 us | 0.9925 | 1.0073 | 0.9246 | 0.9503 |
| 108 us | 0.9915 | 1.0080 | 0.9168 | 0.9455 |
| 250 us | 0.9899 | 1.0093 | 0.9036 | 0.9373 |
| 404 us | 0.9872 | 1.0113 | 0.8835 | 0.9255 |

Removing floor moves the challenger's point estimate the right way and its
**upper limit the wrong way**, because shrinking both denominators inflates the
per-case log-ratio spread faster than it separates the arms. A floor fix would
flatter the headline incumbent-versus-rho number and make the challenger gate
harder. This round therefore proposes no floor fix, and records the table so no
later round proposes one believing it will help.

What is left is the outer level: the number of curve cells. The outer variance
is exactly sigma_b^2 / C, so every half-width above scales by sqrt(5/C).

## The panel is smaller than the curve family allows

The panel has been five cells since round 0006 and was never chosen; it was
inherited. A census of every `(degree, curve_a)` pair the worker accepts
(odd degrees 5..31, both `a`), asking the frozen round-0015 worker for a fixture
and keeping the pairs that yield a usable subgroup, finds **17 usable cells**:

| cell | subgroup order | | cell | subgroup order |
|---|---:|---|---|---:|
| n5a0 | 11 | | n15a1 | 211 |
| n5a1 | 11 | | n17a1 | 65,587 |
| n7a0 | 29 | | n19a0 | 130,873 |
| n7a1 | 71 | | n19a1 | 262,543 |
| n9a0 | 127 | | n23a0 | 2,095,853 |
| n9a1 | 37 | | **n23a1** | **4,196,903** |
| n11a1 | 991 | | n29a1 | 42,457 |
| n13a0 | 2,003 | | **n31a0** | **1,439,393** |
| n15a0 | 751 | | | |

Two of those have never been run: **n23a1**, whose subgroup is twice the largest
cell in the current panel, and **n31a0**, the highest usable degree. `n21`,
`n25` and `n27` have no usable subgroup at either `a`, so degree 31 is not the
start of a ladder — it is the top.

This round takes the panel from five cells to eight by adding `n23a1` and
`n31a0` to development and selection, and `n29a1` alongside the existing `n19a1`
holdout in confirmation and replay. The legacy five remain, so every round-0015
comparison is recomputable from this round's receipts as a subset.

## What this round claims, and how it can fail

The headline claim of this campaign is that the IC winner beats rho — 0.7372 of
rho's instructions and 0.9246 of its native wall on round-0015 confirmation,
with every cell below the parity margin. That claim is scoped to five cells
topping out at a subgroup of 2,095,853. The single weakest cell in it is the
largest one: at n23a0 the winner's native ratio against rho was 0.9642 against a
panel mean of 0.9246, and its instruction ratio 0.9012 against 0.7372. The
margin closes as the curve grows.

Extending the panel to a subgroup of 4,196,903 and to degree 31 is therefore not
a formality. It is the first real test of whether the campaign's central result
survives at the top of the family it is stated over.

Pre-registered predictions, in the order they will be checked:

1. **`scan_io` still fails the promotion gate.** At eight cells the half-width
   scales by sqrt(5/8) = 0.79, giving roughly 0.0114 about a point near 0.9925
   and an upper limit near 1.004. Promotion needs an upper limit below 1. This
   prediction is against the challenger the campaign has spent two rounds on,
   and it is recorded so that a promotion would be a genuine surprise rather
   than a rationalised one.
2. **The non-algorithm remainder stays arm-dependent, in the same direction.**
   In every cell of the widened panel, the IC winner's `process_wall_seconds`
   minus its in-worker `elapsed_seconds` will exceed rho's, as it did in all
   five round-0015 cells. If the sign reverses anywhere, the explanation offered
   above -- that the IC arm pays for emitting and retiring a larger certificate
   -- is wrong and the remainder is something else that needs finding.
3. **The winner still beats rho panel-wide**, but by less than round-0015's
   0.9246 native, because the two added large cells sit where the margin is
   thinnest. A panel-wide native ratio above 0.95 would be consistent; above
   1.00 would falsify the headline.
4. **The margin narrows with subgroup size.** The winner's native ratio against
   rho at n23a1 and n31a0 will exceed its value at n17a1 (0.9148).

Falsification, stated before the run: if the winner's native or instruction
ratio against rho exceeds the 1.10 parity margin in **any** cell, `rho_parity`
fails and the campaign's rho claim is rescoped in the decision record to the
degrees where it holds, rather than restated at panel level. If prediction 2
fails, this round reports no native verdict at all.

## Objective and gates

Objective `rho`, unchanged. Promotion over the incumbent still requires the
instruction ratio at or below 0.98 with a 95% upper limit below 1, a native 95%
upper limit below 1, and every cell within 1.10. `beats_rho_strict` still
requires the winner over rho below 1 in both metrics, every cell, on
confirmation and replay. No gate, margin, bootstrap, repetition count or
confidence level changes in this round. The only protocol change is that the
panel is now named on the command line instead of being hard-coded, and the
defaults reproduce the rounds 0006-0015 panel exactly.

## Parent, seed, budget

Parent: round 0015 (`runs/round-0015`). Both arms are the byte-identical frozen
trees that round sealed — incumbent from `runs/round-0015/source`, `scan_io`
from `runs/round-0015/source_candidates/scan_io/source` — so no source changes
in this round at all and every difference measured here belongs to the panel.
Seed 2026091616, pilot profile, one target per job, CPU 3, 60 s timeout.

## Boundary, floor, class, honesty

Unchanged from round 0015. The unit is Valgrind 3.22 `Ir` over the complete
worker process; native timings are paired cold-process wall under the same
blocking-reap protocol and carry their own scope. The class is accounting, not
a complexity claim: nothing here is a statement about any curve outside the
eight cells named above, and adding degree 31 to the panel widens the tested
scale, not the claim's family. Cross-round ratios continue to carry fixture
variation at least as large as the differences between recent rounds, as
`WINNER-single-target.json` records.
