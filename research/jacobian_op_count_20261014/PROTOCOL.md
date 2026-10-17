# Protocol: field-multiplication count of a Jacobian group operation against the cost model

Version 1, 2026-10-14. Frozen before the full run; its sha256 is recorded in
`manifest.json`. Skill: `.claude/skills/controlled-experiment/SKILL.md`.
Instrument: `examples/jacobian_op_count.rs` (native Rust). Conductor is not
installed in this environment, so no task was attached.

## Question and claims

`docs/curves/jacobian-transfers.json` prices a Jacobian group operation at
`jacobian_group_op_mults_per_genus_squared · g² = 10 g²` field multiplications and
requires an advantage of `advantage_margin_bits = 8` bits over rho. Of its 848
`formula_advantage` curve verdicts, 846 are extension-field `cover_descent`
routes (genus 6, 10 or 14 for most; `log2_genus` up to 7) and 2 are binary GHS
routes. How many multiplications does a textbook Cantor group operation use as
`g` grows, and could a gap to `10 g²` move a verdict?

| claim | class (§1 of the skill) | licensed by |
|:--|:--|:--|
| "textbook Cantor addition/doubling over `F_p` uses `M(g)` multiplications, `I(g)` inversions on average, for `g = 2 … 32`" | exact computation, averaged over random divisor pairs | the counted rows, each checked against the repository's group law |
| "`log2 (M(g) + 8 I(g)) = a + e log2 g` with `e` in [lo, hi]" | exact computation on a sample, bootstrap interval | the fit over genera |
| "the number of advantage verdicts within `δ` bits of the margin" | exact computation on the committed JSON | `docs/curves/jacobian-transfers.json` |
| "the cost model is wrong / right" | **not made** | see Decision rule |

## Instrument

- Textbook Cantor over `u64` residues (`p < 2^31`): two extended-Euclid gcds,
  composition, reduction until `deg u <= g`. One counter each for field
  multiplication, addition/subtraction and inversion, and one for reduction steps.
  An inversion is reported separately and weighted `INV_WEIGHT = 8`
  multiplications in the summed cost (a convention, written here before the run).
- Inputs: random smooth monic curves `y² = f(x)`, `deg f = 2g + 1`, over `F_p`
  for `p ∈ {1009, 65521}` (`p ≈ 2^10` matches the extension routes' base field);
  divisors are sums of `g` random points, resampled until `deg u = g`
  (generic position). Addition uses two independent divisors, doubling one.
- 200 trials per `(g, p, op)`; `g ∈ {2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24, 32}`;
  seed `20261014`.

## Controls

- **Positive control (correctness of the instrument).** Every counted result
  `(u, v)` is compared with `MumfordDivisorP::add` of the repository on the same
  inputs. A single mismatch fails the run (exit code 1). Result of the dry run:
  0 mismatches in 50 rows; the full run must also give 0.
- **Counts are not a measurement of time.** No timing enters the result.
- **Null / second population.** The two primes give two independent samples of
  curves and divisors; the counts should not depend on `p` beyond the inversion
  term. Both are reported.

## Decision rule

Reported per genus: mean counts, `log2` of mean total over `10 g²`, and the fitted
exponent with a 95% percentile bootstrap interval (2000 resamples of trials within
each genus). Verdict wording, fixed before the run:

- "within one bit of the model" if `|log2 ratio| <= 1` at every `g` in the range;
- "model below textbook Cantor by `r` bits at genus `g`" if `log2 ratio > 1`;
- "model above textbook Cantor" if `log2 ratio < -1`;
- "exponent 2 is inside / outside the fitted interval".

Sensitivity (exact count on the JSON): for each of the 848 `formula_advantage`
verdicts, `slack = ratio_log2 − 8`. A verdict is counted as able to drop below the
margin if `slack <= Δ(g)`, where `Δ(g) = log2(fitted measured cost) − log2(10 g²)`
at the route's genus and `Δ(g) > 0`; this assumes the whole route cost scales
with the group operation, so it is an upper bound on fragility. The counts of
verdicts with `slack <= 1, 2, 3` bits and the number of `formula_marginal`
verdicts are reported beside it.

No test against a null is made: the comparison is between a measured count and a
model constant, both exact, and a significance level would have no meaning.

## What the dry run disclosed

On scratch output (5 trials, `p = 1009`, `g ≤ 8`, not committed): the counted
additions use about `125 / 294 / 505 / 1253 / 2487` multiplications at `g = 2 / 3 /
4 / 6 / 8` against `10 g² = 40 / 90 / 160 / 360 / 640`, a ratio of about 3 to 4.
The textbook algorithm is not the cheapest known (explicit formulas exist for
`g = 2, 3`), so a ratio above 1 says the model is below *this algorithm*, not
below the best one. This was seen before freezing; it did not change the design.

## Pre-flight runtime estimate

Dry run on scratch output: 3 trials at `g = 32` take 0.6 s, i.e. about 0.2 s per
trial dominated by the repository's divisor generation. Estimate for the full run:
about 3 to 4 minutes in total. Record the actual time in `manifest.json`.

## What does not count

Changing `INV_WEIGHT`, the genera, the primes, the trial count or the sensitivity
rule after seeing results. Any such change is a labelled new protocol version, and
both are reported.
