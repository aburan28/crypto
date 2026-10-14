# Greedy coverage-driven factor bases: relation yield at equal unknowns

Run `2026-10-14`, seeds 1–20, exhaustive prime-order subgroup targets, pair sums (`m = 2`).
Report: [`REPORT.pdf`](REPORT.pdf) (source [`REPORT.typ`](REPORT.typ)).

## Answer

At equal unknowns (signed Frobenius orbits, the columns of the expected-trial formula)
and equal table size (points), a greedy choice of orbits covers more prime-order
subgroup targets than a uniform random choice of the same orbits, when coverage is
only partial. At the best column count the expected trials fall by:

| curve | greedy `E` at `U*` | random `E` (mean of 20) | ratio |
|:--|--:|--:|--:|
| `K₀/2¹³` (`r = 2003`) | 4.46 (`U* = 3`) | 5.29 | 1.18 |
| `K₀/2¹⁵` (`r = 751`) | 3.00 (`U* = 2`) | 3.25 | 1.08 |
| `K₁/2¹⁷` (`r = 65587`) | 16.35 (`U* = 12`) | 18.18 | 1.11 |

`E = (U + 1) T / covered`, from the census. The gain vanishes once coverage saturates.
It is a constant-factor effect at small and medium degree. The construction enumerates
the subgroup, so it does not yet build cryptographic-size bases.

## Files

| path | what it is |
|:--|:--|
| `REPORT.pdf`, `REPORT.typ` | the results report, with the four figures |
| `data/rows.json` | one row per base (2776 rows): arm, seed, columns, points, covered, coverage, expected trials, build and census times, scalar-model count and whether it matches the census |
| `data/summary.json` | 56 column counts, with means over seeds and 95% bootstrap intervals on the coverage ratio |
| `data/run.log` | the run's stderr |
| `graphs/plot_yield.py` | the figure source; reads only the two JSON files above |
| `graphs/fig1_coverage_vs_columns.{svg,png}` | coverage against column count, three curves |
| `graphs/fig2_coverage_ratio.{svg,png}` | greedy over random coverage ratio, with bootstrap band |
| `graphs/fig3_expected_trials.{svg,png}` | expected trials against column count |
| `graphs/fig4_construction_time.{svg,png}` | construction wall time on `K₁/2¹⁷` |

Benchmark copies (report PDF, SVG and 150 dpi PNG): [`benchmarks/20261014-greedy-factor-base-yield`](../../benchmarks/20261014-greedy-factor-base-yield/).

Code: `src/cryptanalysis/koblitz_greedy_base.rs` (builder, control, tests),
`FactorBaseSpec::GreedyCoverage` in `src/cryptanalysis/koblitz_factor_base_search.rs`,
`examples/koblitz_greedy_factor_base_yield.rs` (the run).

## Reproduce

```sh
# run (about 10 minutes, release build, 4 cores)
cargo run --release --example koblitz_greedy_factor_base_yield -- \
    research/greedy_factor_base_yield_20261014/data 9,13,15,17 20

# figures and report, from the recorded values
python3 research/greedy_factor_base_yield_20261014/graphs/plot_yield.py
typst compile research/greedy_factor_base_yield_20261014/REPORT.typ
```

The example writes `rows.json` and `summary.json` into the `data/` directory. Those
two files are the record, and re-running overwrites them. The coverage, covered-count
and expected-trial fields should reproduce exactly for the same seeds; the timing
fields will differ between hosts and runs.

## Limits

- Degrees 9, 13, 15, 17 only. `K₁/2⁹` and `K₁/2¹⁵` are saturated by one orbit and do not discriminate.
- Pair sums (`m = 2`) only. Per-trial oracle cost is not measured; bases are matched in points, so the pair-table cost is matched by construction.
- The greedy is a heuristic, deterministic in full-scan form, with no optimality claim.
- The construction enumerates the prime-order subgroup and is capped at `r ≤ 2²²`. Scaling it needs sampled candidates and a hashed set of covered sums; that is follow-up work.
- The construction time is not comparable with the shipped sampler, which does not enumerate the subgroup.
