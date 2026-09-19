> This archive records the measured source variants. This PR publishes their
> evidence and UI; it does not replace the runtime Gaudry or GPU implementations.

# Summation-polynomial setup and specialization ablations

This round pursues three paths in `gaudry_cubic::SymmetrisedS4` independently:
borrowing cached symbolic expansions, fixed-slot coefficient accumulation, and
per-monomial Horner evaluation in the residual abscissa. It also evaluates the
compatible combined implementation. Results and the retention decision are in
[RESULTS.md](RESULTS.md).

## Frozen boundary and comparison

[contract.json](contract.json) is written before measurement. All variants must
produce the same specialized component polynomials, complete decompositions,
relation streams and recovered scalars. The summation-polynomial degree and
relation yield are unchanged. Classification is engineering, with a separately
identified accounting correction. No normalized-cost or generic-group advance
is implied by a faster specialization routine.

Six implementations are frozen under [sources/](sources):

| Variant | Change relative to the previous allocation-optimized Gaudry code |
|---|---|
| legacy | Previous implementation, including its historic specialization charge |
| reference | Same arithmetic, with actual specialization multiplication accounting |
| setup | Reference plus borrowing the cached symmetric expansion instead of cloning it |
| fixed | Reference plus fixed-slot accumulation, emitting each nonzero output coefficient once |
| horner | Reference plus once-per-curve grouped coefficients and Horner evaluation |
| combined | Borrowed setup plus grouped Horner evaluation |

`fixed` and `horner` are alternative specialization strategies, not cumulative
per-target passes. The combined strategy incorporates fixed-slot grouping at
curve setup and Horner evaluation per target. It does not run both specializers.

The fixed accumulator has 125 slots indexed by `(a,b,c) -> 25a+5b+c`, with only
the 35 monomials of total degree at most four emitted. It replaces repeated
hash-map updates with field additions in the array, preserving products exactly.
Horner groups up to five target-power coefficients for each occupied monomial
once per curve. It initializes from the highest nonzero coefficient, then performs
one extension-field multiplication per remaining degree. All missing intermediate
coefficients are represented as zero; cancellation and zero/base-field targets
retain their exact semantics. It needs no per-target power table.

## Separate accounting correction

The old solver charges `15 * number_of_terms` for specialization. The actual
power-table evaluator performs four target-power multiplications plus one per
term. `Fp3::mul` counts eleven base-field multiplications each, so the charge is
`11 * (number_of_terms + 4)`. Neither formula was an optimization of the other.
The reference measures `f.muls()` immediately before and after specialization,
including power preparation; it does not reset or disturb the caller's counter.
This remains accurate when Horner changes the products performed.

All optimization comparisons use the corrected reference. Legacy is retained to
replay the prior exact results and quantify the accounting adjustment separately.
Setup and fixed must preserve the reference's field-operation counts exactly.
For Horner and combined, every solver-counter change must equal the independently
measured per-specialization reduction times the number of solver calls. All other
solver counters, complete solution sets and relation statistics must match.

## Validation and benchmark scope

The original WDSat binary-ANF protocol is inapplicable to this prime-cubic-field
S4 solver. Under the equivalent-suite exception in `AGENTS.md` §8, this round
replays all four previous 60-residual cases and all three previous cold full-DLP
cases, plus new curve seeds. The [parent accounting contract](../index_calculus_baseline_20260914/ec_index_calculus_contract.json)
still applies; it does not replace or modify the historical WDSat evidence.

The microbenchmark measures thirty complete symbolic setups, and sixteen passes
over 256 specialization inputs per curve. Targets include zero, one, base-field
and arbitrary extension-field abscissas. Inputs can repeat, especially in the
small field; these are evaluation counts, not counts of independent targets.
Canonical sorted polynomial coefficients are hashed with BLAKE3 for exact
cross-variant comparisons. The independent termwise unit test also exhausts
all 343 abscissas over F_(7³), plus zero/base targets and samples at p=67,271.

Stage runs retain full decompositions for sixty targets per curve and check them
against the independent MITM oracle. MITM-only two-term solutions are excluded
from the S4 triple contract, following the existing test; returned solutions are
also verified directly by curve arithmetic. Cross-variant comparisons retain
all outputs without filtering. Cold full-DLP runs build their own instance and
factor base, collect relations, solve the matrix, recover and verify the scalar,
and run matched rho. Every residual, including empty answers, is MITM-audited.

All six variants use identical isolated dependency versions and harnesses. Cases
run in randomized six-way blocks, three repetitions per case, pinned to one CPU.
Microbenchmark timers separate setup from specialization. Stage/full-DLP process
times include cold setup, verification, output, process startup and MITM/rho
validation work. No warm cache is shared between processes. Raw commands, errors,
timeouts, result hashes and per-run JSON are retained. Prior completed legacy
results are replayed against the earlier allocation benchmark, ignoring only
explicit timing fields.

Paired timing ratios use log ratios and bootstrap whole input clusters, keeping
all repetitions together. These intervals describe the small frozen corpus on a
shared host. A local setup/evaluation gain does not establish a full-DLP gain.
Full common-operation totals, S, cost/rho and cost/floor remain null because cold
setup and runtime bookkeeping are not completely calibrated in one operation
unit. Partial counted Fp products and existing solver totals remain explicit.

## Reproduce

Make an isolated copy of the previous Gaudry allocation candidate, with the
saved Cargo manifest and lock. Copy the three saved files from `harnesses/` into its `examples/` directory.
For each variant in `contract.json`, copy `sources/<variant>.rs` onto
`src/cryptanalysis/gaudry_cubic.rs`, then build:

```sh
cargo build --release --locked --example summation_polynomial_bench \
  --example gaudry_allocation_stage --example gaudry_cubic_bench
```

Save those three executables under `<build-root>/bin/<variant>/` before building
the next variant. Keep `<build-root>/Cargo.lock`, Cargo.toml and the harness
sources. Then run:

```sh
python3 research/summation_polynomials_20260914/run.py \
  --build-root /path/to/isolated-build --output /tmp/s4-fresh-run
cargo test --release --locked --lib cryptanalysis::gaudry_cubic::tests -- --test-threads=1
```

The output directory must not exist. The runner preserves failures and rejects
missing completions, polynomial mismatches or unexplained counter differences.

## Performance visualizations

[Interactive dashboard](../../docs/performance-gains.html) with metric switches,
paired confidence intervals and individual confirmation runs.
[PNG](../../docs/performance-gains/summary.png),
[SVG](../../docs/performance-gains/summary.svg), and
[PDF](../../docs/performance-gains/summary.pdf) exports are available.

Regenerate from frozen evidence (requires Matplotlib and NumPy):

```sh
python3 research/summation_polynomials_20260914/visualize.py
```

The generator reads saved benchmark results; it does not rerun measurements.

The dashboard also includes exact tables without JavaScript and a CSV export.
“Copy this view” preserves the selected metrics, controls, and section in the URL.
Local file previews provide a copyable file URL instead of claiming a public link.
Set `CRYPTO_VISUALIZATION_DIR` to optionally copy the generated page and assets
to another directory; the default generator writes only to repository `docs/`.

Validate the publishing routes and dashboard behavior:

```sh
python3 -m unittest discover -s scripts/site -p 'test_*.py'
node scripts/site/test_performance_ui.js
```
