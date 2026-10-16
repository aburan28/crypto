# Protocol: fixing the D₄ + double-large-prime regression, k = 4 IC on E(F_{p⁴})

IC improvement round 2 (`../RULES.md`), item 04. Frozen 2026-10-15 before
any data was collected. Prior work: `../../ic_improve_20261014/01-d4-2lp/`
(item 01: tuned D₄+2LP / D₄ whole 1.025 [1.016, 1.037], REGRESSES).

## Question

Item 01 showed that a tuned D₄+2LP regresses against D₄ because the relation
phase decides the cost and 2LP does not cut the relations needed. Which of the
three proposed fixes removes the regression, and does any of them improve on
the D₄ incumbent?

## Derivations to be checked before measurement (stated in the README)

1. **Large primes inside the decomposition base B.** On a shared trial
   stream, every full relation the large-prime graph (or any elimination of
   large-prime columns) produces is a linear combination of the collected
   rows. So the logarithm is determined by those rows no later than by the LP
   combinations. A stop at the *first linear dependency* of the collected base
   rows therefore needs at most as many trials as any LP scheme over F ⊂ B
   that sees the same rows. Shrinking F cannot lower the relation phase below
   that. Prediction: the gain item 01's Gaudry arm hinted at (≈ 0.92 ×
   columns relations) comes from stopping at the dependency rather than at
   columns + 1. It is available without large primes.
2. **Large primes outside B** (one summand from a larger set Λ ⊃ B). Such a
   relation is useful only after a collision in Λ. Derive the trial count and
   compare it with D₄'s.
3. **Pre-filtering trials with ≥ 3 large primes.** Say whether the
   large-prime count can be bounded before the C₄ solve, and why.
4. **Adaptive arm.** A cost model fitted on training curves chooses α per
   size, with a 5 % margin before leaving plain D₄.

## Arms (one shared trial stream per curve, same targets, same LA seeds)

| arm | role | definition |
|:--|:--|:--|
| `d4` | **incumbent** | item 01's `d4`, unchanged. First LA attempt at ncols + 1 full relations, retry every +5 %, square-core filter, then Wiedemann. |
| `d4-ed` | **candidate A** (paths 1–2) | `d4`'s base, columns, pick rule and solver. The only change is the stop: square-core filtering is attempted after *every* new relation, and Wiedemann runs as soon as the filtered core is square (rows = columns incl. d). A failed solve retries at the next relation. No large primes, so no rejected trials. This is the large-prime graph taken to its limit (every column is a "large prime" and a hypergraph dependency is the cycle). |
| `d4-adapt` | **candidate B** (path 3) | Choose α ∈ {plain, 0.95, 0.90, 0.85, 0.80, 0.75} at the curve's log₂ #E (public before solving) from the frozen cost model. Leave plain only if the modelled whole ratio is < 0.95. Otherwise this is `d4`, the incumbent schedule. |
| `d4-ed-all` | exploratory | `d4-ed`, but it feeds every distinct row of a decomposed trial, picked row first. |
| `d4-2lp` | reference | item 01's frozen D₄+2LP (α = 0.95, cap ∞, uniform, ω = 1). Expected to reproduce the regression. |
| `d4-2lp-null` | **null A/A** | the 2LP path with s = ncols. Must equal `d4` in trials, full relations, LA counts and output on every curve. |
| `d4-late110` | **positive control** | `d4` with the first LA attempt at ⌈1.10·(ncols + 1)⌉. Known regression of about +10 % relations. |

Diagnostic (uncharged): the exact first-dependency relation index FD of
`d4`'s row sequence, by incremental elimination mod r in hindsight. This is
the minimum any stop rule could reach on that row sequence.

## Training and freeze

- **Training curves:** seeds 4, 5 at p ∈ {101, 181, 269, 397, 521, 769, 1033}
  (14 curves, item 01's training split). Not used for testing.
- **Cost model for `d4-adapt`.** For each α ∈ {0.75, …, 0.95}, use
  configuration `a<α>-c0-uniform-w1.0`, plus `d4`. Fit OLS of log₂(relation
  GAE) and of log₂(LA GAE) on log₂ #E over the 14 training curves. The α
  rows come from item 01's committed `raw/train.jsonl`: same curves, same
  stream, same code, kept for frugality. This run's own training pass
  re-runs `d4` and α = 0.95 on the same curves and must reproduce their
  trial counts exactly. Otherwise the model is refitted on this run's data
  only.
- `d4-ed` has no tuned parameter. On training it is only checked for
  correctness and run once before the freeze.
- Frozen file: `raw/frozen.json`. No change after test data exists.

## Test (held out)

Seeds 1, 2, 3 at p ∈ {101, 181, 269, 397, 521, 769, 1033, 1301}: 24 curves,
item 01's test split. Stretch: p = 1601 (seeds 1–3) if it fits 3 GB and
60 min. For that the oracle stores 48-bit bucketed fingerprints instead of
64-bit ones. This is output-identical, because the pair table is exact and
every decomposition is re-added in the group. `d4` trial counts must still
equal item 01's on all 24 shared curves.

Curves: `quartic_folds::generate_translation_instance(p, seed, 16)`, with one
planted target per curve. The trial stream is `trial_ab(seed, t)`, as in
item 01.

## Metric (item 01's charge, unchanged)

S_charged = (factor base + trial scalar mults + trials·C₄/97 + re-adds +
LA mod-r mults·16/97 + final check)/√r, in GAE, with C₄ = 2.54·10¹⁰.

Phases:

- relation = scalar + trials·C₄/97 + re-adds;
- LA = (Wiedemann + merge) mod-r mults·16/97;
- whole = the sum.

The integer bookkeeping of square-core filtering is uncharged in every arm,
as in item 01. It is counted (index operations) and reported, with a
sensitivity check that charges it at 1 F_p mult per index operation.

## Statistics and decision rule (item 01's)

- **Paired ratio.** Candidate/incumbent per curve. Geometric mean per size
  and pooled, with a 95 % bootstrap over curves (2,000 resamples).
- **Exponent.** OLS slope of log₂ ratio on log₂ #E over all test curves
  (≥ 4 sizes), with a bootstrap CI. Per-arm, per-phase exponents use OLS of
  log₂ cost on log₂ #E.
- **Constant level.** ADVANCES if the CI < 1, REGRESSES if the CI > 1,
  MATCHES if it covers 1. A ratio identical to 1 on every curve is MATCHES.
- **Exponent level.** ADVANCES if the CI < 0, REGRESSES if the CI > 0,
  MATCHES otherwise.
- Each is applied to whole, relation and LA. TRADES when phases or levels
  disagree in sign.
- Verdicts are also given per size (3 curves each). With 3 curves, a per-size
  bootstrap CI is coarse, and this is stated.
- **Headline.** `d4-ed` whole-cost verdict at the constant and exponent
  levels. `d4-adapt` must be at least MATCHES at every measured size.

## Controls and gates

- **Correctness.** Every arm, every curve: [d]G = Q and d = the planted value.
  Every decomposition is re-added in the group. Any failure fails the arm.
- **Null.** `d4-2lp-null` must be identical to `d4`.
- **Positive control.** `d4-late110` must be called REGRESSES at the constant
  level, with the point in [1.04, 1.16].
- **Reproduction.** `d4` trial counts must equal item 01's `raw/test.jsonl`
  (24 curves) and `raw/train.jsonl` (14 curves).
- **Dominance check of derivation 1.** On every curve, compare the trials of
  `d4-ed-all` and `d4-ed` with `d4-2lp`'s, and on training with every item-01
  α configuration. Report any curve where an LP arm used fewer trials.

## ecbench

Add `base = yline-d4-ed` to `ic.k4_gaudry`: D₄ Y-line with the eager stop.
This is an additive edit to `fp4.rs` (a new field `eager_la` in `K4Params`)
and to `methods.rs` (a new base value). Existing specs' parameters and hashes
are unchanged. Run p = 101, 181, 269, 397 (seed 1, 3 targets) with arms
`yline-d4` and `yline-d4-ed`, audit it with replays, and report
`ecbench compare`.

## Budget

`nice`, single thread, `prlimit --as=3000000000`, `timeout 3600` per job, at
most two jobs at once.

## Inadmissible

Changing C₄, the charges, the split, the arms, the margin or the frozen model
after test data; dropping a curve or a failed arm; quoting a model handover
as measured.
