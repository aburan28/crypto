# Protocol: D₄ symmetrisation stacked with double large primes, k = 4 IC on E(F_{p⁴})

IC improvement round, item 01 (`../RULES.md`). Frozen 2026-10-14 before any
data was collected. Prior work: `../../ic_next10_20261014/01-fp4-ladder/`
(its PROTOCOL.md, README.md and `raw/ladder.jsonl`).

## Question

Does adding double large primes (2LP) to the D₄ Y-line arm lower the whole
charged cost of a k = 4 index-calculus logarithm, and at which level
(exponent, constant)? Does a *tuned* 2LP still lose at the relation phase,
as the untuned x-line 2LP did in item 1 (relation exponent 0.238 → 0.333)?

Rho is a reference column only; no verdict depends on it.

## Incumbent and candidate

| arm | definition |
|:--|:--|
| `d4` (**incumbent**) | item 1's D₄ arm unchanged: Y-line base {P : Y(P) ∈ F_p}, one column per {±P, ±(P+T)}, one decomposition per decomposed trial drawn uniformly, first LA attempt at (columns + 1) full relations, retry every +5 %, filtering to a square core then Wiedemann |
| `d4-2lp` (**candidate**, tuned and frozen) | same base, columns, trial stream and charge. Small base F' = the first s = ⌈ncols^α⌉ columns; the remaining columns are large primes; a decomposition is admissible with ≤ 2 large-prime columns; partials go to `gaudry_cubic::LargePrimeEliminator` with merge cap `cap` (0 = uncapped); pick rule `uniform` (uniform over admissible decompositions) or `minlp` (uniform over the admissible decompositions with the fewest large primes); first LA attempt at ⌈ω·(s + 1)⌉ full relations, retry every +5 % |
| `d4-2lp-gaudry` (untuned stacked arm, reported) | the candidate at item 1's DLP rule: α = 0.75 (Gaudry's \|F\|^{1−1/k}), cap 0, uniform, ω = 1 |

### The four tuned 2LP parameters and their grid (60 configurations)

- large-prime bound α ∈ {0.75, 0.80, 0.85, 0.90, 0.95} (s = ⌈ncols^α⌉);
- cycle-finding cut-off, merge cap ∈ {0 (uncapped), 2, 4};
- relation selection ∈ {uniform, minlp};
- relations before the first LA attempt, ω ∈ {1.0, 1.1}.

## Split (curves are disjoint)

- **Training:** seeds 4, 5 at p ∈ {101, 181, 269, 397, 521, 769, 1033}: 14 curves.
  All 60 configurations plus the incumbent ride the same trial stream per
  curve.
- **Test (held out):** seeds 1, 2, 3 at p ∈ {101, 181, 269, 397, 521, 769, 1033}
  (item 1's 21 curves, so the incumbent can be checked against the committed
  data), plus the stretch size p = 1301 (seeds 1–3) if it fits the 3 GB cap
  and 60-minute job limit. Arms: `d4`, `d4-2lp` (frozen), `d4-2lp-gaudry`,
  and the controls below.

Curves: `quartic_folds::generate_translation_instance(p, seed, 16)`; one planted
target per curve; every arm sees the same (a, b, pick) trial stream
(`trial_ab(seed, t)`), the same base points and the same target. Each arm has
its own LA random stream seeded identically.

### Tuning rule (fixed now)

Freeze the configuration minimising the geometric mean over training curves
with p ≥ 521 (6 curves) of the paired whole-cost ratio S_charged(config) /
S_charged(d4). On training curves the stream stops once the incumbent is
solved and the trial count reaches 8× the incumbent's; an unsolved
configuration is scored at its cost at the stop (a lower bound) and is
flagged. No change to the frozen configuration after test data exists.

## Metric and charge (item 1's, unchanged)

S_charged = (factor base + trial scalar multiplications + trials × C₄(D₄)/97
+ decomposition re-adds + LA mod-r multiplications × 16/97 + final check) / √r,
GAE. C₄(D₄) = 2.54·10¹⁰ F_p multiplications. The exact enumeration oracle
supplies decompositions and is uncharged (as in item 1; the D₄ solver matched
it on 6 of 6 targets there). Phases:

- relation phase = trial scalar + trials × C₄/97 + re-adds;
- linear algebra (LA) = (Wiedemann + filtering + large-prime merge mod-r
  multiplications) × 16/97;
- whole = S_charged · √r.

## Statistics

- Per arm and phase: OLS of log₂(cost) on log₂ #E over the test curves, 95 %
  bootstrap over curves (2,000 resamples).
- Paired: per test curve, ratio candidate / incumbent for whole, relation and
  LA cost. Per size: geometric mean over the curves at that size with a 95 %
  bootstrap over those curves. Pooled: geometric mean over all test curves,
  bootstrap over curves. Trend: OLS slope of log₂(ratio) on log₂ #E with
  bootstrap CI (this is the exponent difference).
- Handover: the log₂ #E where the fitted log₂(whole ratio) line crosses 0
  (extrapolation, labelled as such); also where the per-phase power-law
  models (relation fit + LA fit per arm) cross.

## Decision rule (per level, per phase)

- **Constant level** (pooled paired ratio R, 95 % CI): ADVANCES if CI < 1;
  REGRESSES if CI > 1; MATCHES if it covers 1.
- **Exponent level** (slope of log₂ ratio, 95 % CI, ≥ 4 sizes): ADVANCES if
  CI < 0; REGRESSES if CI > 0; MATCHES if it covers 0.
- Applied to whole, relation and LA separately. **TRADES** when one phase or
  level advances and another regresses. The headline verdict is the whole-cost
  verdict at each level.

## Controls

- **Null A/A** `d4-2lp-null`: the 2LP code path with s = ncols (no large
  primes possible), cap 0, uniform, ω = 1. Must equal `d4` in trials, full
  relations and LA counts on every curve, train and test.
- **Positive control** `d4-half`: the incumbent that discards every
  odd-indexed trial after paying for it. Expected whole ratio ≈ 2. The
  instrument must call it REGRESSES at the constant level with the CI
  excluding 1 and the point in [1.6, 2.5].
- **Reproduction:** `d4`'s trial counts on test seeds 1–3 must equal item 1's
  committed `raw/ladder.jsonl` (same stream; the y-base RNG is consumed after
  the x-line lifts, as item 1 did).
- **Correctness gate:** every logarithm is checked [d]G = Q and against the
  planted d; every decomposition is re-added in the group. An arm with a
  wrong output fails.

## ecbench

If the arm can be registered with an additive edit (a new `base` value in
`ic.k4_gaudry` that allows large primes on the Y-line), run a 4-size session
(p = 101, 181, 269, 397, seed 1, 3 targets) with `d4` and the frozen
`d4-2lp`, audit it with replays, and report `ecbench compare`.

## Budget

Every job single-threaded, `nice`, `prlimit --as=3000000000`, `timeout 3600`.

## Inadmissible

Changing C₄, the 97/16 charges, the grid, the split, the tuning objective or
the frozen configuration after test data; dropping a curve or a failed arm;
quoting the handover as measured.
