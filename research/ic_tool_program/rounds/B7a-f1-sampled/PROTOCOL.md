# B7a: the F1 sampled level at one-word sizes

**Declared 2026-10-01, before any B7a code.** This is Track B's B7a, as
designed in [`../../design/f1-sampled.md`](../../design/f1-sampled.md)
with its amendment 1. B7a's cases are C071–C077, frozen with this
declaration in
[`../../conformance/v2-b7a/`](../../conformance/v2-b7a/cases.json). The
constants F1 carries are frozen in [`carried.json`](carried.json).
Nothing below changes after the first measurement, except by a dated
amendment at the end.

## What B7a does

### 1. `fidelity: F1` runs

- **`kic`** is extrapolated on any instance it admits at F0: Koblitz
  curves, and curves over `GF(2^k)` with `1 < k ≤ 8` once B2b is in.
  This works alone or paired.
- **Every rho pipeline** is extrapolated too. Its online cost is its
  step cost, measured on the instance, times the expected steps
  `√(πr/2A)`. So F1 runs rho at any size its pipeline admits: the m = 83
  gate's field on `rho-koblitz`, or secp256k1 on `rho-bignum`.
- **`ic-binary-s4` and `ic-prime-s3`** have no F1 model. At F1 they are
  refused with a new code, `no-f1-model`, class `unsupported`, exit 3.
- **The size gates** are F0's. F1 runs the real field's kernels, so a
  pipeline F0 gates is gated at F1 too. (F2 alone skips the five
  word-width gates.)

### 2. What F1 runs, for `kic`

This is the design's §3 and §4, made exact.

- **In full:** the curve, the factor base, its projection and the pair
  table.
- **Collection:**
  - The yield is counted. With `K` the table's stored pairs and `e` the
    orbit length, `μ = r/(2e·K)`.
  - Each stored pair is a key up to colliding sums, which are negligible
    while `K` is far below `r/2e`. The carried `κ` was computed the same
    way.
  - The cost of a summand is sampled: 2^12 summands of warm-up, then
    2^16 timed. They are scanned by collection units, aimed as F0 aims
    them, from unit 0.
  - The sample's relations are reported beside the count's prediction,
    as the check `κ_sample`.
- **The total summands:**
  `S = max(ρ·κ·c·μ, units_planned · unit_trials · window)`.
  The second term is the unit F0 always runs, which dominates at the
  smallest sizes.
- **The descent:**
  - Its probes come from `μ` by §20's `probes(μ, granule)`, with
    granule 64 for `m = 2` and `min(1024, |F|)` for `m = 3`.
  - The cost of a probe is sampled: 8 probes of warm-up, then 64 timed,
    on probes of the real target.
- **Linear algebra:** the block Wiedemann solver runs on a random sparse
  system. It has `c` columns, `⌈ρ·c⌉` rows of weight 3 over the signed
  orbits, and residues modulo `r`. Its first 8 iterations are timed, and
  a solve is taken as `2c/B + 8` iterations, `B` the block size.
- **The unit:** measured before and after, as F0 measures it.

### 3. What F1 carries

The constants are in [`carried.json`](carried.json), each with its
source and its interval.

- **`ρ`, relations per column at the end of collection.** It is the
  median of v0's `ρ` at the six largest suite sizes (R01's A/A runs):
  1.0135, with the range `[1.000, 1.075]` as its interval. One constant
  serves every size, so no size is extrapolated from its own F0.
- **`κ`, the count's correction.** It is 1, with the interval
  `[0.93, 1.07]`: the spread of the three v0 sizes with more than 200
  relations, rounded outward.

### 4. The intervals

Both come by Monte Carlo, from 10,000 draws over the constants'
intervals (uniform within each) and the timed samples (resampled).
- **The expectation's interval** draws only those.
- **The predictive interval, for one run,** also draws:
  - the collection's summands, as Gamma with shape `ρ·c` and mean `S`;
  - one target's descent probes, as exponential with mean `probes(μ, g)`.

### 5. The report

- **The run's keys:**
  - `status: extrapolated` and `fidelity: F1`;
  - `label`: "an extrapolation from samples (F1); not a measurement; it
    does not discharge the m = 83 gate";
  - `speedup_eligible: false`.
- **`phases.ic`**: for each of `select`, `build`, `collect`,
  `linear_algebra` and `descent`, these keys:
  - `how`: `in_full`, `counted_and_sampled`, `sampled` or `synthetic`;
  - the sample's counts and times;
  - the constants, with their intervals;
  - the formula.
- **`phases.rho`**: the step sample, the expected steps and the formula.
- **`extrapolated.ic` and `extrapolated.rho`**:
  - `cold` and `online`, each in `units` and `seconds`, with `median`,
    `expectation_lo`, `expectation_hi`, `predictive_lo` and
    `predictive_hi`;
  - `complete: true` once every phase has been priced.
- **`carried`**: the constants of §3, with their sources.

### 6. The level

- **`fidelity: auto` still means F0 or a refusal.** A run asked for as
  a measurement never becomes an extrapolation on its own.
- **The over-budget refusal** now also gives F1's own estimated cost.
  F1 is asked for by name.
- **`fidelity: F1`** is refused `over-budget` when F1's own estimated
  cost exceeds the budget. That cost covers the curve, selection and
  build by B2's model, and the samples at their declared sizes.

### 7. The runner and the steps

B7a and B7b are added to `../../conformance/run.py`'s `STEPS`. B7a runs after B2b; the order of `STEPS` is only how reports list the steps.

## The cases

C071–C077 are in
[`../../conformance/v2-b7a/cases.json`](../../conformance/v2-b7a/cases.json).
The generator `make_cases.py` copies earlier steps' frozen documents and
changes their method, and `SHA256SUMS` pins the files. An F1 figure
varies from run to run, so the cases check the report's shape and its
codes, not its numbers.

| case | input | expected |
|:--|:--|:--|
| C071 | B1's C010 (curve A) at `fidelity: F1`, paired | exit 0, `extrapolated`; `phases.ic.collect.how: counted_and_sampled`; `extrapolated.ic.complete` and `extrapolated.rho.complete` true |
| C072 | C010 at F1, `solve: rho` | exit 0, `extrapolated`; `rho.pipeline: rho-koblitz` |
| C073 | C010 at F1, `solve: index_calculus` | exit 0, `extrapolated`; `operation: ic_single_target` |
| C074 | the m = 83 gate file at F1, `solve: rho` | exit 0, `extrapolated`, on `rho-koblitz` at `n = 83` |
| C075 | B2's C047 (secp256k1) at F1, rho alone | exit 0, `extrapolated`, on `rho-bignum` |
| C076 | B2's C038 (a generic curve, `ic-binary-s4`) at F1 | exit 3, `no-f1-model` |
| C077 | B2b's C059 (a curve over `GF(4)`) at F1 | exit 0, `extrapolated`; `kic` and `rho-negation` |

## Class

**Robustness and fidelity.** F1 adds a level. It changes nothing at F0,
so the pin and the timing check apply as in every Track B step, and it
claims no speedup.

## Arms

- **The base**: the newest accepted baseline when B7a runs. B2b comes
  before B7a.
- **B7a**: that commit plus B7a's change, built with
  `IC_BUILD_COMMIT=$(git rev-parse HEAD)`.

## Measurements

1. **Tests.** `cargo test --release --bin ic`. They cover:
   - the count against an enumerated table at small `n`;
   - the Monte Carlo intervals against closed forms where those exist;
   - a partial table's key count against `g²` times the full one's.
2. **Conformance**, `--steps <accepted>,B7a`.
3. **The pin and the translation**, as in B1 (`bround.py pin` and
   `translate`). F0 is unchanged.
4. **No slowdown at F0, timed**: the base against B7a on `M1`'s 22
   rows, five rounds ABAB, isolated (`bround.py timing`). The B7a arm's
   F0 runs are also the F0 reference of measurement 5.
5. **F1 against F0.** F1 on `M1`'s 22 rows with B7a's build, three
   rounds, isolated: 66 processes. Per size, each figure set against
   measurement 4's B7a F0 cold times on the same rows:
   - F1's expected cold cost over F0's median: the ratio;
   - whether F0's median lies in F1's predictive interval;
   - F1's wall time over F0's.
6. **The count's check and the carried constants.** Per size:
   `κ_sample` with its Poisson interval, F0's `ρ` and `κ` from
   measurement 4's reports, and a fit of `κ` against `n`.
7. **Partial tables.** At the six largest sizes, F1 forced to partial
   tables at `g = 1/2` and `g = 1/4`, one round of `M1`'s rows each: 24
   processes. The figures:
   - `b`, `u` and the descent's probe cost, against the full table's;
   - the key count, against `g²` times the full one's.

## Acceptance

B7a is **accepted** when all of the following hold:
- every test passes;
- every case the runner selects passes;
- the pin holds on every row;
- the translation check holds on every row;
- no size regresses beyond its A/A band at F0;
- **the falsification target** (the design's amendment 1):
  - in measurement 5, F0's median lies in F1's 95% predictive interval
    at 9 or more of the 11 sizes;
  - at `n = 53`, `59` and `61`, F1's expected cold cost lies within
    `[0.85, 1.18]` of F0's.

Measurements 6 and 7 are reported, not gated. They are B7b's evidence.

If the target is missed, B7a is **revised once**:
- the sample sizes and nothing else change, by dated amendment, before
  the revision runs;
- measurement 5 then runs again on the revised build.

If the target is missed again, B7a is **abandoned**, and the record says
by how much. Any other failure **rejects** B7a, as for every step.

## Inadmissible

- Changing `carried.json` after any B7a run.
- Choosing a size's constants from that size's F0 figures.
- Widening an interval after seeing a miss.
- Dropping a size.
- Reporting an F1 figure as a measurement or a speedup.
- Counting a contended run.

## Cost

- The tests and the conformance suite: about ten minutes.
- The pin and the translation: 270 untimed processes.
- The timing check: 220 processes, about 90 minutes.
- F1 against F0: 66 processes. F1 runs the full table, so each costs
  20–40% of its F0 row.
- Partial tables: 24 processes.
