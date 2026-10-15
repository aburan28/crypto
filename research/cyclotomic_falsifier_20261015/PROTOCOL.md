# PROTOCOL — the cyclotomic γ-weight base's falsifier (encyclopedia §8, D4; lead A2)

**Frozen:** 2026-10-15, committed before any measured cell of this study.
**Lane:** `cyclotomic-falsifier` (branch `claude/leads-cyclotomic-falsifier`).
**Conductor:** the control plane in this container refuses connections
(`conductor check` → `dial tcp 127.0.0.1:8080: connection refused`); recorded,
work proceeds without a task card.

## 1. The hypothesis and its own falsifier, quoted

`docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md`, item A2 (L187–221):

> **Falsifier.** At `n=23`: build the `γ`-weight-`w` base for `w = 2, 3`,
> report `|F|`, the subgroup-usable count, and natural 2-sum and 3-sum yields
> against a Hamming-weight ONB base of matched size; then compare F4 and
> WDSat solving time and solving degree on the two systems for 16 targets.
> Fail if the cyclotomic system's solving degree is not lower, or its yield
> per column is lower.

The encyclopedia (`docs/ic/factor-base-encyclopedia/README.md` §5 D4) records
D4 as "an open engineering lead; the falsifier is not established". The prior
census (`research/ic_techniques_20261012/a2_cyclotomic_and_delta/`) measured
coverage on 32 targets with witness counts capped at 256, against one ONB
base per cell and no null population; its solver cells at `n = 23` timed
out.

This study establishes the **yield clause** exactly, with a null population
and a positive control, and measures the **end-to-end cost** in ecbench. It
does not re-run the solver clause (F4/WDSat solving degree): the brief says
not to depend on the timed-out solver cells, and the `n = 23` systems did not
finish in the prior budget. The `n = 11` solver record is cited as prior
evidence only (§6).

## 2. Instrument (yield): exact sumsets in discrete-log space

Example `examples/cyclotomic_falsifier_yield.rs` (Rust; one thread).

- Curves: `K_a / F_{2^n}` from ecbench's `CurveSpec::Koblitz {a, n}`
  (`ic_boundary::koblitz_instance`), the curves `ic.pipeline` runs on.
- Every base is a set of points of the order-`r` subgroup `G`, closed under
  negation and Frobenius (`P ↦ λP`), i.e. a union of `K` projected signed
  Frobenius orbits ("columns"), each of `2n` points.
- A full table of `[k]G` gives the exact log of every base point, so a base
  becomes `B ⊂ Z/r`. Self-checks (abort on failure): `log φ(P) = λ·log P`;
  for plug-in bases, `log(point) = coef · log(column representative)` for
  every dumped point.
- Metrics, each exact over **all `r` targets** (no sampling):
  - `cov2 = |B + B|` — targets with at least one ordered 2-decomposition
    (relation probability for `m = 2` is `cov2 / r`);
  - `cov3 = |B + B + B|` — the same for `m = 3`;
  - `energy = #{(a,b,c,d) ∈ B⁴ : a + b = c + d}` — additive energy, the
    collision count behind `cov2`.
- Identity (derived, no control needed): the mean number of ordered
  `m`-decompositions per target is `|B|^m / r` for **every** base of `|B|`
  points, so at matched columns the yield per column can differ only through
  sum collisions — exactly what the three metrics measure.

### Bases (arms), all at the same column count `K`

| arm | construction | role |
|:--|:--|:--|
| `cyclotomic_w{w}` | plug-in `cyclotomic:w=W,columns=K` via `ecbench::methods::dump_factor_base` | target |
| `compact_orbit_scan` | plug-in `compact-orbit-scan:columns=K,raw_x_cap=2^n−1` (same dump path) | comparator (incumbent D2) |
| `onb_weight` | least ONB Hamming weight `v` reaching `K` projected columns; first `K` columns in the set's order; same projection/dedup rule as the plug-ins | A2's own comparator |
| **null** | `N = 1999` bases: abscissae `x` uniform in `F_{2^n}^*` (splitmix64, seed `20261015 ⊕ draw·0x2545f4914f6cdd1d`), each lifted, projected by `[h]`, deduplicated by signed Frobenius orbit, until `K` columns (nested prefixes of one draw serve every `K`) | negative control: random Frobenius-stable base of the same size and fold |
| `positive_lambda_progression` | columns with representatives `log = 1, 2, 3, …` (deduplicated by orbit) | positive control (planted additive structure) |

The null is the right one: it differs from `C_w` only in how the abscissae are
chosen (uniform vs γ-weight), with identical lifting, cofactor projection,
Frobenius/negation fold and column count.

**Positive control, expected detection (derived before data).** In the
λ-progression base, pairs inside the same Frobenius power `λ^i·{±1..±K}` sum
into an interval of `4K+1` values, so energy gains about
`n·(2/3)(2K)³` over a random base's `≈ 2|B|² + |B|⁴/r` (at `n = 23, a = 1,
K = 15`: about +40%), and `cov2` loses about a `1/n` share of its pairs.
Expected: energy above every null draw and `cov2` below the null on the
`n = 23` cells.

### Cells (pre-registered by rule, not by data)

For each rung and each `w ∈ {1, 2, 3}` whose cyclotomic base exists
(`K_full ≥ 1` projected columns):

- transition sizes `K_m = round((m!·r·ln 2)^{1/m} / 2n)` — the column count at
  which a random base covers half the targets with `m`-sums (the only regime
  where yield can differ; saturated cells have coverage 1 for every base);
- cell A: `K = min(K_2, K_full)`, metrics `cov2` and `energy`;
- cell B: `K = min(K_3, K_full)`, metric `cov3` (merged with cell A when the
  two `K` coincide).

Rungs: **primary** `n = 23`, `a = 0` and `a = 1`; **secondary (cheap)**
`n = 11`, `a = 0` and `a = 1` where the curve has a usable subgroup. A rung
where `C_w` has no projected column records "no base" (an existence result).

**`n = 83`: no honest cheap cell.** The plug-in and `CyclotomicField` are
single-word (`n < 64`); `ecbench fb` refuses wide curves; exact enumeration
needs `r ≈ 2^81` table entries. Recorded as out of scope, not run.

### Decision rule (yield clause)

- For each test (base × metric × cell): empirical two-sided add-one rank
  p-value against the `N = 1999` null draws,
  `p = min(1, 2·min(#{null ≤ v}+1, #{null ≥ v}+1)/(N+1))`.
- Family: the `n = 23` cyclotomic tests — at most 5 bases (`a = 0`: w = 2, 3;
  `a = 1`: w = 1, 2, 3) × 3 metrics = **15**; Bonferroni, `α = 0.05`,
  critical `p < 0.05/15 = 0.00333`. Power check: smallest attainable
  `p = 2/2000 = 0.001 < 0.00333`; `N + 1 = 2000 ≥ 4·15/0.05 = 1200` ✓. The
  same critical value is applied (descriptively) to comparators, controls and
  `n = 11`.
- **Lower yield** = `cov2` or `cov3` significantly below the null median
  (energy significantly above is reported as the collision mechanism but does
  not decide by itself).
- **Positive control gate:** the λ-progression base must be significant on
  `energy` (above) in both `n = 23` rungs' cell A. If not, the instrument is
  blind and the yield verdict is **inconclusive**.
- **Yield-clause verdict:**
  - **FALSIFIED** — any `n = 23` cyclotomic `cov2`/`cov3` test is
    significantly below the null;
  - **SURVIVES (not lower)** — none is, with the positive control passing;
    sub-labelled *yield-neutral* if no test is significant in either
    direction, *yield-advantaged* if some `cov` test is significantly above;
  - **INCONCLUSIVE** — positive control fails, or no `n = 23` cell runs.
- Comparators (`compact_orbit_scan`, `onb_weight`) are one deterministic draw
  each; their rank in the null and their difference from `C_w` are reported
  descriptively.

## 3. Instrument (cost): ecbench end to end

Specs `research/cyclotomic_falsifier_20261015/specs/*.json`, written by the
rule below once the census gives `K_full` (the census is a count of the base,
not a measurement of yield or cost):

- reference `rho-frob` = `rho.signed_frobenius_full`; IC arms
  `ic.pipeline_full` with `oracle = mitm-frobenius-counted:m=2`,
  `targets = stride-renew`, `max_trials = 1000000000` — the configuration of
  the incumbent study `research/ecbench_frobenius_orbit_vs_rho_20261014`;
- `a = 1, n = 23`: `cyclotomic:w=1` vs `compact-orbit-scan:columns=1`
  (the prior smoke's base, now matched), and `cyclotomic:w=2` (all columns,
  `K_full`) vs `compact-orbit-scan:columns=K_full`;
- `a = 0, n = 23`: `cyclotomic:w=2` (all columns) vs
  `compact-orbit-scan:columns=K_full`, and `cyclotomic:w=3,columns=3` vs
  `compact-orbit-scan:columns=3` (the incumbent's configuration there);
- `raw_x_cap = 100000000` for compact; 16 public targets per curve
  (`target_kind: public`), 2 rounds, 1 warmup, `isolation_required: L0`,
  timeout 3600 s; operation counts `S` only (wall time is not a result here).
- Audit: `ecbench verify --replay-all`; every replay must be identical, or the
  session is no result.
- Comparison: `ecbench compare --a <compact arm> --b <cyclotomic arm>` and each
  against `rho-frob`.
- **Cost verdict per pair:** *cyclotomic costlier* if the paired-ratio 95%
  interval (cyclotomic / compact) lies above 1; *cheaper* if below 1;
  *unresolved* otherwise. Not a speed claim; L0.

## 4. Overall verdict for D4 (what the encyclopedia entry may say)

- **falsified** — yield clause FALSIFIED;
- **survives** — yield clause SURVIVES; then the lead's status depends on the
  untested solver clause and the cost verdicts, which are reported beside it;
- **inconclusive** — yield clause INCONCLUSIVE.

## 5. What does not count

Changing `N`, the seeds, the metric, the cells' `K` rule, the null
construction or the critical value after seeing results. Any change is a new,
dated amendment below, and both versions are reported.

## 6. Prior evidence cited, not re-run

`a2_summary_n11.md` (2026-10-12): at `n = 11, a = 1, w = 1, m = 2` the
cyclotomic membership system's F4 max degree was 4 against 3 for the matched
ONB system, F4 mean 754.6 ms vs 7.3 ms, and the Macaulay probe resolved 0/16
vs 3/16 targets at degree 4. The `n = 23` solver cells timed out. This study
neither re-runs nor reinterprets them.

## 7. Pre-flight

Disclosed before the measured run: a dry run with `N = 19` on scratch output,
used only to check that the code runs, the self-checks pass and to time one
null draw. Its timings and anything it reveals are recorded in the amendment
below before the measured run.

## Amendments

### A1 (2026-10-15, after the disclosed dry run, before any measured cell)

Dry run: `cyclotomic_falsifier_yield --a 1 --n 23 --nulls 19` on scratch output
(not committed). It showed, and this amendment records without changing any
rule, `N`, seed, metric, cell or critical value:

1. **Census (counts of the bases, not yields).** `K_full` (projected columns):
   `a = 1, n = 23`: w = 1 → 1, w = 2 → 30, w = 3 → 338; `a = 0, n = 23`:
   w = 1 → none ("no gamma-weight-1 abscissa carries a subgroup point"),
   w = 2 → 17, w = 3 → 367. `K_0 / F_{2^11}` has no usable prime-order
   subgroup in the repository's constructor, so the `n = 11` rung is `a = 1`
   only. Transition sizes at `a = 1, n = 23`: `K_2 = 52`, `K_3 = 6`.
2. **The `K = 1` cells are degenerate (identity).** Every one-column base is
   `l·{±λ^j}` for a unit `l`, and `|B+B|`, `|B+B+B|` and the energy are
   invariant under multiplying `B` by a unit, so every 1-column base — null
   draws, `C_1`, compact, ONB and the positive control alike — has the same
   three numbers (seen: 921, 10534, 7314 at `a = 1`). The `w = 1` tests stay in
   the family of 15 (conservative) and return `p = 1` by construction; the
   `w = 1` base is therefore yield-identical to every other one-column base,
   a derived fact rather than a test.
3. **Speed.** About 0.9 s per null draw with `cov3` computed at every `K`;
   the null loop now computes `cov3` only at the `K` of a `cov3` cell (it is
   never read elsewhere; comparators still compute it at every `K`). No rule
   changes.
4. **Values seen in the dry run** (N = 19, not evidence): cyclotomic effects
   within ±0.4% of the null median on every metric; the positive control at
   +279% energy and −40% `cov2` at `K = 30`.
5. **Cost-leg specs** written by the §3 rule with these `K_full`:
   `specs/cost-a1-n23.json`, `specs/cost-a0-n23.json` (`make_specs.sh`).

### A2 (2026-10-15, cost leg, before any verified cost record)

The first cost session (`sessions/void-a1-n23-full-unpriced/`, spec
`specs/void/cost-a1-n23-full.json`) errored on every record before it was
stopped: `src/cryptanalysis/ecbench/full_prices.json has no row for
icv1-f2m23-tm5197-1f85e9e1` (the `*_full` methods need a measured price row,
and `K_1 / F_{2^23}` has none). Kept as evidence, no values. Change: both cost
specs use the operations-only methods `rho.signed_frobenius` and `ic.pipeline`
(same oracle, targets, factor bases, seeds, sizes), one instrument on both
rungs; `S` then counts group operations without pricing native field work, as
in the 2026-10-13 smoke. A one-target, one-round smoke of the new `a = 1` spec
(scratch, not committed) verified all five arms; its values are not evidence.
Decision rule unchanged.
