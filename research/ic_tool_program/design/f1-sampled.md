# F1, the sampled level: design (B7)

**Written 2026-10-01, before any B7 code.** This designs Track B's B7
(`research/notes/index-calculus/IC_TOOL_PROGRAM.md` §9) for `kic`, the
programme's index calculus. B7 is "done when `n = 83` and `131` are
reported as labelled extrapolations, with their samples". It is split
into two steps:

- **B7a** builds F1 where F0 also runs, at one-word sizes (`n ≤ 62`). It
  measures F1's error there, size by size, against F0.
- **B7b** runs F1 at `n = 83` and `131`, once B3b and B4 give `kic`'s
  kernels two and three words. B7a is its evidence.

F1 never discharges the `m = 83` gate (AGENTS.md §8a; plan §3).

## 1. The level

Plan §3 defines F1 as "the real field and curve, with each phase's
kernel run on a sample: relation yield per trial, build cost per row,
linear algebra per non-zero, descent per trial. The totals are
extrapolated." An F1 report:
- claims an extrapolation, labelled as one;
- gives its samples, its formulas, the constants they measured, and an
  interval for every total;
- names nothing a speedup.

## 2. What F0 costs, by phase

§20's phase model (`research/ic_exponent_20260926/predict.py`,
`ic_phases`) prices one target in units, the batched addition. With
`|F|` base points, `n` the orbit length (`e = n/k` on a curve over
`GF(2^k)`), `c = |F|/2n` columns and `m` descent summands:

| phase | cost |
|:--|:--|
| select | `|F| · s` |
| build | `(|F|²/4n + |F|) · b` |
| collect | `ρ · c · μ · u` |
| descent | `probes(μ, granule) · d` |
| linear algebra | not in §20's model: under 1% of the set-up at the suite's sizes (plan §8, A4) |

- **The constants.**
  - `s`: units per selected point.
  - `b`: units per stored pair.
  - `u`: units per scanned summand.
  - `d`: units per descent probe.
  - `μ`: summands scanned per relation. §20 writes it as `μ = Y · r/|F|²`,
    where `Y` is a yield constant (2.02 at `n = 41`).
  - `ρ`: relations per column at the end of aimed collection (1.026 at
    `n = 41`).
- **The probes.** `probes(μ, g) = g / (1 − e^{−g/μ})`: the probes to the
  first witness, paid `g` at a time.

§20 took every constant from one older run. F1 measures each cost on
the instance it reports, and it counts the yield. That is the difference
between F1 and F2.

**The yield is a count, not a rate to sample.**
- A scanned summand `R − P_k` is a relation when its canonical key is
  one of the table's keys. A key stands for a signed Frobenius orbit of
  `2n` points. So with `K` distinct keys, a summand hits with
  probability `2n·K/r`, and `μ = r/(2n·K)`.
- The folded table stores about `|F|²/4n` keys, so `μ ≈ 2r/|F|²`. That
  is `Y = 2`.
- §20 measured `Y = 2.0214` at `n = 41`: 883,200 summands for 197
  relations at `|F| = 15,744`. The count predicts 4,436 summands a
  relation, against 4,483 measured. The ratio
  `κ = μ_measured/μ_counted` is 1.011.

## 3. What F1 runs

The curve and the factor base are built in full, as F0 builds them. They
are the instance, and every later sample needs them. Then:

| phase | at one word (B7a) | measured |
|:--|:--|:--|
| select | in full, timed | `s` |
| build | in full when the table fits its byte budget, timed. Otherwise a **partial table** (below). | `b` |
| collect | **the yield by counting**: the table's distinct keys `K`, so `μ = r/(2n·K)`; and **a short sample of the scan**, aimed as F0 aims it, `S_s` summands, to price one | `μ`, `u` |
| descent | the probes to a witness from `μ`, by §20's `probes(μ, granule)`; and **a short sample of probes** on the real target, to price one | `d` |
| linear algebra | the solver on a **synthetic system** with the instance's column count and row weight, a fixed number of iterations | cost per iteration, so per non-zero |
| ρ | not sampled: it is a property of finishing collection. It is carried from F0 at the nearest suite size, and named. | — |

- **The sample sizes.**
  - `S_s = 2^16` summands, after a warm-up of `2^12`: enough to time a
    summand to about ±2%.
  - At the six largest suite sizes, a collection scans far more:
    `ρ·c·μ` is `2^18.8` to `2^26.1` summands.
  - At the five smallest, the whole collection scans `2^5.5` to
    `2^16.3`: no more than about the sample, so F1 is no cheaper there.
    It is validated there, not needed.
  - The descent's sample is 64 probes after a warm-up of 8.
  - Each sample also has a wall budget, a declared fraction of the
    run's. A sample that hits it reports a wider interval and says so.
- **The count's own check.** The collection sample's relations are
  counted too, and compared with the count's prediction `S_s/μ`. Over
  B7a's rows they give `κ` with a Poisson interval. A `κ` far from 1
  would mean the count misses something, such as duplicate pair sums,
  and B7a reports it.
- **The units.** Every timed sample is converted at the unit measured in
  the same process, as F0's are (`price::UnitBench`).
- **The extrapolation.** It is §2's table, with F1's constants and
  intervals. Each interval comes from its sample: Poisson for counts,
  the order statistics of repeated timings for costs. The total's
  interval follows by Monte Carlo over the constants' intervals.

## 4. Partial tables

A full pair table holds about `|F|²/4n` pairs. At `n = 83`, `|F|` is
about `r^{1/3} ≈ 2^27`, so the table would hold about `2^45` pairs: far
beyond memory. F1 therefore defines a **partial table**:
- the pairs `{i, j}` whose orbits both lie in a random subset of the
  base's orbits, a fraction `g` of them;
- so a fraction about `g²` of the full table's pairs.

What a partial table can measure:
- **The build constant `b`.** A stored pair costs the same however many
  are stored, provided the table and its scratch exceed the last-level
  cache, as the full one does. F1 chooses `g` so that they do.
- **The scan constant `u`.** The same proviso applies: the presence
  filter and the buckets must miss the cache as the full table's do.
  Huge pages, if a later round adds them, change `u`, and the report
  names the build it measured.

What it does not need to measure:
- **The yield.** The count gives it from the full table's key count,
  which is `|F|²/4n` up to colliding pair sums. Collisions are
  negligible while `|F|²/4n` is far below `r/2n`: at `n = 83`, `2^45.6`
  against `2^73.6`.
- No sample could measure the yield there anyway. Through a partial
  table, relations arrive `g²` times as often, so one relation would
  need about `μ/g² ≈ 2^45` scanned summands.
- What B7b carries is `κ`, the count's correction. B7a's measurement 3
  measures `κ` at every suite size with its interval, and fits it
  against `n`. If it drifts, B7b carries the fitted drift, with its
  exponent named, and widens the interval to cover the fit's residuals.

At one word, B7a also forces partial tables at `g = 1/2` and `1/4`. It
compares `b` and `u` with the full table's, and the partial table's key
count with `g²` times the full one's, at the same sizes. So partial
tables' own error is measured before B7b depends on them.

## 5. Linear algebra

LA is under 1% of the set-up at the suite's sizes. It grows faster than
the other phases (plan §8, A4; AGENTS.md §5), so F1 prices it.
- **The solver.** F1 runs the block Wiedemann solver
  (`koblitz_sparse_la.rs`) on a random sparse system. The system has the
  instance's column count, its relation rows' weight (`m = 3` summands,
  folded to signed orbits), and residues modulo the instance's `r`.
- **The cost.** A fixed number of iterations is timed, giving the cost
  per iteration. A full solve takes about `2c/B` iterations of block size
  `B`.
- **The extrapolation.** The cost per iteration is about the non-zeros
  times a constant, and F1 checks that law at two column counts.
- **Past one word.** At `n = 83`, `r ≈ 2^81` needs multi-word residues.
  B3b gives the solver those, and F1 then times them. It does not carry
  the one-word cost per non-zero across.

## 6. The report

`ic price` with `fidelity: F1`:
- **`status: extrapolated`**, `fidelity: F1`.
- **`phases`**: for each phase, what ran (`in_full`, `sampled`,
  `partial_table`, `synthetic`, or `carried`), the sample's counts and
  times, the measured constants with intervals, and the formula.
- **`extrapolated`**: the cold cost and the online cost in units and in
  seconds on this host, each a median and a 95% interval.
- **`carried`**: every constant not measured on this instance, with its
  source.
- **`label`**: "an extrapolation from samples (F1); not a measurement;
  it does not discharge the m = 83 gate".

There is no rho arm. Rho's own F1 is its step cost times the expected
steps `√(πr/2A)`. The report gives that too, with the step cost measured
on the instance, so the two arms are extrapolated alike.

## 7. B7a's measurement and its falsification target

On the suite's eleven sizes, `M1`'s 22 rows, one target each:
1. **F1 against F0.** The ratio of F1's extrapolated cold cost to the
   newest baseline's measured F0 cold cost (the same rows, the same
   build), per size.
   - **The target**, declared now: every size's ratio in `[0.80, 1.25]`,
     and F0's value inside F1's 95% interval at 9 of 11 sizes or more.
   - **The model's error, for comparison.** B2's F2 model, with no
     samples, measured 0.34–1.98 against v0 in development. F1 is
     worth having only if it is much tighter.
2. **F1's own cost.** Its wall time as a fraction of F0's, per size. At
   one word the full table bounds it from below, since build is 9–28% of
   the set-up.
3. **The count's correction.** `κ` at every size, with its Poisson
   interval, and a fit of `κ` against `n`. This is B7b's evidence for
   carrying it.
4. **The carried `ρ`.** F0's `ρ` at every size, against the value F1
   carried for it.
5. **Partial tables.** As §4 says: `b`, `u` and the key count at
   `g = 1/2` and `1/4`, against the full table's.

**B7a is abandoned** if measurement 1 misses its target at more than two
sizes after one revision of the sample sizes. That revision is
declared, by dated amendment, before it runs.

**Inadmissible:**
- choosing the carried `ρ` after seeing F0's figures;
- dropping a size;
- widening an interval to cover a miss;
- reporting an F1 figure as a measurement or a speedup.

## 8. What B7b adds

- **The kernels.** B3b and B4 give `kic`'s select, build, scan, descent
  and LA kernels two and three words.
- **The instances.** F1 then runs at `n = 83`, the m = 83 gate's field
  `x^83 + x^45 + x^2 + x + 1` (AGENTS.md §8a), and at `n = 131`,
  ECC2K-130's field. Both use partial tables, the counted yield, and the
  carried `κ` and `ρ`.
- **The report.** The two extrapolations go on the scoreboard, marked as
  extrapolations. Their exponents and carried constants are named,
  beside rho's F1 on the same instance.

## 9. Order

1. B7a's protocol, with its cases, declared before code.
2. B7a's code: the sampled phases, partial tables, the synthetic LA, and
   the report. The F0 path is unchanged, and the pin holds.
3. B7a's measurement (§7).
4. B3b (two-word `kic` kernels), then B7b.

## Amendment 1 (2026-10-01, before B7a's protocol and before any F1 run)

v0's own runs (R01's A/A) show that §7's target ignored F0's sampling
noise. B7a's protocol uses the corrected target below.

- **F0 is one draw.**
  - A collection that needs `R` relations scans a Gamma-distributed
    number of summands, with relative spread `1/√R`.
  - At the five smallest suite sizes, `R` is 8 to 16, so F0's collection
    cost scatters by ±25–35% on its own. Both `M1` rows of a size share
    one set-up, so each size has a single draw.
  - The `m = 2` descent's probes are exponential: one target is one
    draw, with a spread equal to its mean.
- **The count's correction `κ` in those runs** (summands per relation
  measured, over counted) is 0.88–1.33 across the eleven sizes.
  - Most of that is the same Poisson noise.
  - The three sizes with 256–336 relations give 0.96–1.04.
- **So F1 reports two intervals.**
  - The interval of its expectation, from its constants' uncertainty.
  - A predictive interval for one run. It adds the Gamma spread of the
    collection's relations and the exponential spread of one target's
    descent.
- **The target is restated, before any F1 run:**
  1. F0's measured cold cost lies inside F1's 95% predictive interval at
     9 or more of the 11 sizes.
  2. At the three sizes where F0's collection holds 200 relations or
     more (`n = 53`, `59`, `61`), F1's expected cold cost lies within
     `[0.85, 1.18]` of F0's.

  The rule for abandoning B7a is unchanged.
