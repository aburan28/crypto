# Experiment proposals after the MMT / Lemma 6 correction

2026-10-14. Five proposals were drafted from the corrected magic-number audit
(`docs/audits/GHS_MAGIC_NUMBER_MMT_CORRECTION_20261013.md`), ranked by
information per unit cost, and then run. This file is the plan and the index:
each result lives in its own PR with a frozen protocol, controls, a report, a
graph and a PDF. Nothing here is a new measurement.

| # | question | status | where |
|:--|:--|:--|:--|
| 1 | Can a short 3-isogeny walk leave a blocked twist tower? | run | PR #1843 |
| 2 | Is the cost model `10 g²` behind `formula_advantage` calibrated? | run (as a counted-operations measurement) | PR #1846 (merged) |
| 3 | Walk-distance spectrum in binary isogeny classes | already answered; not rerun | `research/weak_isogeny_classes_20261013/` |
| 4 | Smith's Prop 2 / Prop 4 conditions, natively | run | PR #1848 |
| 5 | Which odd-characteristic descents land on genus 3? | run | PR #1845 |

## 0. Pilot: where condition (8) bites on the catalog

MMT Lemma 6 condition (8) blocks a tower `(N, l, n = N/l)` only when
`Tr_{K/F_2}(a) = 1`, `n` is even and `Tr_{K/k}(b) = 0`. On the 91 binary catalog
curves (179 GHS towers) no tower is blocked as recorded; 113 have `Tr(a) = 1`.
A quadratic twist keeps `b` and flips `Tr(a)`, and 13 towers on 5 curves are
blocked on the twist, including `m = 7` towers on `c2pnb208w1` and
`c2pnb368w1`. The native tally is `tally` mode of
`examples/twist_walk_blocked_towers.rs` (PR #1843), output
`research/twist_walk_blocked_towers_20261014/pilot_towers.csv`. An earlier
Python version of the same count is superseded by it.

Correction recorded here: an early version of this pilot assumed every catalog
curve has `a ∈ {0, 1}`. That is false (27 towers have `Tr(a) = 1` with even
`N`, normal form `a ∈ {0, γ}`). The twist logic was fixed to `Tr(a_twist) =
1 − Tr(a)`; the counts above are the corrected ones.

## 1. Walk past the blocked twist towers (PR #1843)

*Hypothesis as drafted:* if the trace were an independent fair coin on each
neighbour, a walk of `d` steps escapes with probability `1 − 2^−d`.
*That null was wrong, and the protocol says so.* In characteristic 2 only
degree-3 isogenies walk (`Φ_2 mod 2` is Frobenius), and the blocked status at
`l = 1` is constant on a trace class, so the fair-coin null does not apply.

*Result (exact computation, 13 towers):* 3 escape at depth 1 (`N = 24, l = 2`;
`c2pnb368w1` at `l = 2` and `l = 46`, `m' = 12, 8, 8`), 7 have no rational
3-isogeny, 3 stay blocked within depth 4. *Class ceiling (exact, `N = 8…22`):*
for `l = 1` no trace class mixes blocked and unblocked curves, so no isogeny of
any degree helps there; for `l ≥ 2`, 60% to 99.9% of blocked curves lie in a
mixed class. Random blocked starts escape at 0% (`l = 1`), 36–40% (`l = 2`) and
about 75% (`l ≥ 13`). The catalog is not tested against those rates: five
independent graphs cannot reach the corrected level.

## 2. Calibrate the cost model behind `formula_advantage` (PR #1846, merged)

*As drafted* this was an end-to-end comparison on toy towers. The repository's
`ghs_full_attack` runs only `m = 1` end to end, so the study counts operations
of textbook Cantor instead, and compares them with the model's `10 g²`.

*Result (exact computation):* textbook Cantor is 2.04 bits above `10 g²` at
genus 2 and 3.32 bits above it at genus 32; the fitted exponent is 2.35, not 2.
Of 848 `formula_advantage` verdicts, 14 lie within 1 bit of the 8-bit margin
and at most 24 would drop below it if the whole route cost scaled with the
group operation (an upper bound). *Not shown:* that the model is wrong; textbook
Cantor is not the cheapest algorithm, and `g = 2, 3` have explicit formulas.

## 3. Walk-distance spectrum in binary isogeny classes (not rerun)

`research/weak_isogeny_classes_20261013/` already enumerates trace classes
exhaustively and measures class and 3-walk coverage against direct coverage.
The proposal restated that study, so it was not rerun.

## 4. Smith's Prop 2 / Prop 4 conditions, natively

*Result:* all four of Smith's §8 rates reproduce natively on 10^6 random degree-8
models at `p = 2^31 − 1` (seed 20261014): tractable subgroup 0.50187
(prediction 0.50159), Prop. 2 0.49999, Prop. 4 given Prop. 2 0.50102, and
`full_success` 0.18638 (Theorem 2: 0.18574). Each 98.75% Clopper–Pearson interval
lies inside prediction ± 0.01. The positive control (Smith's Example 1 over
`F_37`) reproduces exactly, and an independent point count of the construction
agreed with Prop. 4 in 18,381,760 cases and disagreed in none.

*Class level (exhaustive, curve plus twist):* with `D` the direct `full_success`
share, `C` the share whose class contains such a curve and `U` the share whose
class contains a tractable member, the frozen statistic `r = (C − D)/(U − D)` is
0.811 at `p = 7` and 0.938 at `p = 11`. The rule fixed in advance calls
`r ≥ 0.75` "near the ceiling", so the walk extension is worth costing. This
corrects the hypothesis as drafted here, which expected 50% against 18.6%: the
class-level rates are `C` = 0.430 and 0.475, not 50%, because the two conditions
are correlated inside a class.

*Post hoc, descriptive:* `C` is below `U`, so the reachable ceiling for a walk is
`C`. *Not shown:* no walk or explicit isogeny between class members is built, so
no walk length or cost is measured; only `p = 7, 11`, hyperelliptic models with
`deg f = 8`, nothing about characteristic 2. Smith's Prop. 2 "iff" has a
small-`p` exception (27,645 subgroups at `p = 7`, 535,435 at `p = 11` pass the
square test and have no degree-3 map); the frozen rule counts them as failures.
The problem record `problems/smith-kernel-isogeny-class-walk.md` stays open with
the walk as its next step.

## 5. Which odd-characteristic descents land on genus 3 (PR #1845)

*Result:* yes. For `n = 3` there is a descent onto a genus-3 hyperelliptic curve
for every odd `q` (`y² = e·g(x)(x−α)(x−α^q)`; Momose–Chao, Joux–Vitse §4.1),
and the native tool confirms Joux–Vitse's model on 81 instances. For
`F_q(j) = F_{q^n}` with `n ≥ 4` every `x` gives genus at least 5, and the
extension degrees 6, 10 and 14 used by the catalog's `cover_descent` routes give
minimum genus 9, 49 and 41. The problem record
`problems/ghs-to-smith-pipeline-characteristic-mismatch.md` is **not** closed as
an obstruction.

**Consequence for existing labels (acted on 2026-10-15).** 846 of the 848
`formula_advantage` verdicts were extension-field `cover_descent` routes priced at genus `k`, which no GHS
cover has: the GHS minimum genus at `k = 6, 10, 14` is 9, 49 and 41 when `j` generates the
extension. The rows are now priced at that floor and relabelled; 599 of the 846 are still
`formula_advantage` there, 41 are `formula_marginal` and 206 are `boundary`. See
[`../cover_descent_genus_floor_20261015/`](../cover_descent_genus_floor_20261015/REPORT.pdf).

## Order and what is left

Proposal 3 was dropped as already answered. Follow-ups worth a protocol of their
own: build and price the isogeny walk between Smith class members (proposal 4's
next step); relabel or re-derive the 846 `cover_descent` verdicts; extend the class
ceiling beyond `N = 22` or prove that `Tr_{K/F_2}(b)` is constant on trace
classes (conjectured, not proved); price Cantor against explicit genus-2/3
formulas before using `10 g²` for genus above 3.
