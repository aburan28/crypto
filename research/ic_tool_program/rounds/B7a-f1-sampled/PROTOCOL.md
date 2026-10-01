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

## Amendment 1 (2026-10-01, before any B7a run)

Written while implementing F1, before any measurement this protocol
declares. F1 was run only on documents no measurement uses: the
conformance documents (curve A, the suite's smoke size
`icv1-f2m31-tm90707-c95f16f5`; C059's curve over `GF(4)`; the m = 83
gate's field, rho alone; secp256k1, rho alone). No F1 run touched any of
measurement 5's eleven sizes. Two findings at the smoke size shaped
items 1 and 9, and each rests on an argument that does not depend on
them.

Everything not named here stands: the carried `ρ` and `κ`, the
measurements, the acceptance rule, the falsification target, the one
revision and the abandonment.

1. **The count uses the table's distinct keys.** `K` is
   `PairSumTable::distinct_keys`, not the stored entries.
   - A folded table enumerates each sum of two points of one signed
     orbit from both offsets `g` and `g⁻¹` of its row, so it stores
     those sums twice, and it stores the identity `P + (−P)` once per
     orbit. That is about `1/(c + 1)` of its entries: 11% at `c = 8`,
     0.3% at `c = 320`. Two pairs whose sums share an orbit also share a
     key.
   - So `μ = r/(2e·K)` on a folded table, and `r/K` on an unfolded one,
     whose keys are points.
   - `κ` was computed from stored entries at the three largest sizes,
     where the two counts differ by 0.3–0.4%, inside its interval. It is
     unchanged.
   - A test enumerates every pair sum of two small bases and checks the
     count against them.
2. **Summands per relation are `probes(κμ, w)`**, not `κμ`. `w` is the
   window one collection trial scans, and a trial yields one relation
   at most. When `w ≪ κμ` this is `κμ + w/2`; at the smallest sizes it
   is up to twice `κμ`. So the total summands are
   `S = max(ρ·c·probes(κμ, w), units_planned · unit_trials · w)`.
3. **The collection sample uses F0's unit layout.**
   - The trials run from 0, each unit's part in one call. Each unit is
     aimed as F0 aims it, and the coverage is updated at its end.
   - The first call follows the collector's construction, as F0's first
     unit does. It is timed apart, and what it costs beyond the timed
     rate is charged once.
   - The collector and the coverage are charged twice when F0 would
     extend past its planned units, since F0 builds them again then.
4. **The verification is priced per relation.**
   - Relations are pushed one at a time: the collection sample's own
     relations when it found 16 or more, else synthetic ones. A
     synthetic relation fails after the same group arithmetic, so only
     its row's entry into the system goes unpriced, and the report says
     which kind was used.
   - The relations priced are `ρ·c` when the count decides how long
     collection runs. Otherwise they are what F0's planned units find.
   - The first push is timed apart and its premium charged once.
   - Pinned to one CPU, as measurement 5 runs, a batch and single pushes
     cost the same per relation.
5. **The linear algebra is solved in full on a synthetic system.** This
   replaces "the first 8 iterations timed".
   - The system has `c` columns and `⌈ρ·c⌉` rows of weight 3. Each of
     its first `c` rows holds a column of a random permutation, so every
     column is covered and the generic rank is full. Its coefficients
     are uniform and nonzero modulo `r`, and its right-hand sides come
     from a planted solution.
   - F0's own solver solves it in full, with F0's options: the filter,
     the fold, block Wiedemann, reconstruction and the row check.
   - The table's certification, one scalar multiplication a column, also
     runs in full, on random logarithms.
   - Why: at the suite's sizes v0's filter leaves 0–63 core columns, so
     the solve takes microseconds and the certification is most of F0's
     `la` clock. Eight iterations would have priced neither the filter,
     the approximant basis nor the certification.
6. **The log solver and the descent's solver are built in full.** The
   descent's is built on a placeholder table, every column's logarithm 0,
   which indexes as the real table does.
7. **The descent.** Its online cost is
   `start + probes(κμ, g)·d + one witness`, plus the first start-up's
   premium.
   - `m = 2`: `d` comes from differencing walks of 1,024 and 17,408
     probes on the real target, through the placeholder table, five
     times. A probe that decomposes then fails its recovery check. The
     timed probes' witnesses are counted on one observed walk, and their
     cost (a decomposition and a recovery check, timed) is subtracted.
     This replaces "8 probes of warm-up, then 64 timed", which is one
     round of the walk and too short to time.
   - `m = 3`: `d` is per summand, from 64 full scans of the base (after
     8), on states `[a]G + [b]Q`.
   - The start-up and the witness: 64 timed after 8.
8. **Rho.**
   - Its reusable set-up runs in full for each walk.
   - Its step cost and per-target start-up are fitted by least squares
     to seven capped walks on the target, under the seeds `seed` to
     `seed + 6`. The caps are 4,096 and 16,384 steps in turn, 65,536 in
     all, after a warm-up walk of 4,096.
   - The expected steps are `√(πr/2A)`, `A` the walk's own; one run's
     steps are Rayleigh-distributed.
   - The unit is the one-word batched addition. Rho alone on
     `rho-koblitz` uses the two-word one its walk runs on. A prime field
     has no unit in the tool, so its figures are in seconds only.
9. **A carried context `χ`.**
   - F1 times each part once, as F0's first repetition does. F0 reports
     the median of its repetitions.
   - In v0's A/A runs, the first repetition's set-up is 0.883 to 1.115
     times the median over the eleven sizes. That makes `χ`
     `[0.88, 1.12]`, with 1 at its centre.
   - `χ` multiplies `kic`'s cold and online totals in both intervals.
   - It is added to `carried.json` with its source, one figure for each
     size. `ρ` and `κ` are unchanged.
10. **The predictive interval** also draws the relations of F0's planned
    units, binomially, when those units decide the collection.
11. **Partial tables (measurement 7)** are forced by
    `ic price --f1-partial g`.
    - The table is built over a random fraction `g` of the base's signed
      orbits, rounded and at least two, from a fixed seed. The samples
      run on that base.
    - The full table's entries come from the build's row structure:
      `(c·|F| + |F|)/2` folded, `|F|(|F| + 1)/2` unfolded. Its keys are
      the partial table's divided by the fraction squared.
    - The build is priced per stored pair, times the full table's
      entries.
12. **`conformance/v2-b7a/params/README.md`** is added. Some cases name
    their documents as `{here}/../../v2-b2/params/…`, a path through
    `params/`. A directory with nothing in it cannot be committed, so
    the path did not resolve. No case and no expectation changes.
13. **The sample sizes** are `carried.json`'s amended `samples` block.
14. **The runner and the analysis** are committed with this amendment,
    before any run: [`run.py`](run.py) and [`analyse.py`](analyse.py).
    - Measurements 2–4 are `harness/bround.py`'s steps. `run.py` adds
      `f1` and `partial`: each M1 row's v2 translation at `fidelity: F1`,
      run through the programme's isolation and retries.
    - `analyse.py` fixes the statistics.
      - F0's cold cost per size is the median over measurement 4's B7a
        processes, each in its own units.
      - F1's per size is the median over its six processes of each
        figure.
      - "Inside" means F0's median lies in F1's predictive interval.
      - The ratio is F1's median over F0's.
    - A size missing either figure makes the decision incomplete, not a
      rejection.
