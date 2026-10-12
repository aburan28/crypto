# F4 cached-layout resource-limit audit

This is a correctness and accounting correction to cached matrix construction.
It makes no performance, relation-yield, or discrete-logarithm claim.

## Contract and correction

A cached layout and a fresh layout must obey the same current matrix limits.
Previously, fresh construction checked `F4_F2_MAX_COLS` in `macaulay_columns`,
but cache hits bypassed that function. Warming a layout under a larger cap and
then lowering the cap could therefore admit a matrix that a fresh build refused.
This affects inherited fused packing, flat fused packing, and materialized
packing through the shared layout cache.

All three lookup sites now use `cached_f4_layout`, which checks the current
column cap before returning a layout. An oversized cached layout is a cache
miss, not an immediate refusal of the input: a new system with smaller support
can still fit, rebuild, and replace the old layout. The existing row-limit
checks remain in the builders. Rejected matrices retain the existing `None`
resource-limit outcome; this is not a mathematical infeasibility certificate.

The regression uses two Boolean variables and at most four monomials. Four
fresh child processes exercise inherited/flat storage with fused/materialized
packing. Each child checks a cold build, an exact-cap hit, a lowered column cap,
a smaller-support fallback, a lowered row cap, and restored caps. Environment
changes are confined to a process running only that test.

## Parallel test isolation

The old cache-hit test compared process-wide hit/miss totals even though the
cache itself is thread-local. Other tests can increment those totals or reset
them independently. The test now checks its own thread's retained `Rc` layout
identity and exact matrix output. It clears only its own cache and makes no
assertion about global counters.

## Validation and claim boundary

Run the nonignored algebra tests, including the isolated cap cases, with:

```sh
cargo test --release -q --lib cryptanalysis::koblitz_groebner::rref_tests:: -- --test-threads=8
cargo test --release -q --lib cryptanalysis::inherited_f4::tests::
```

Local validation: the first command passed 10 tests with 4 timing benchmarks
ignored; the second passed all 5 inherited-basis tests. The cap regression
passed in all four fresh-process packing configurations.

These are small algebraic correctness controls, not a timing experiment. The
fixed-cap historical timing records are left intact; this audit does not prove
or retract their timing ratios. A new performance claim would still require
fresh matched measurements. No full-pipeline cost or ratio to rho is inferred.
