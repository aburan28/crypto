# Bounded Boolean symbolic-product schedule study

This standalone experiment implements the proposed **complete generator-support
signature** contract on small Boolean-polynomial systems. It imports no curve,
decomposition, relation-collection, or solver modules. Its worker is compiled
directly with `rustc` and uses only the Rust standard library.

## Mathematical boundary

In `F2[x0,...,x(n-1)]/(xi^2-xi)`, every polynomial has a unique squarefree
monomial expansion with coefficients in `{0,1}`. Its normalized support is
therefore the whole polynomial. An ordered list of complete generator supports,
the variable count, matrix degree, and multiplier mask specify an identical
ordered matrix construction problem. An exact signature hit cannot represent
a changed canonical system.

The candidate precomputes surviving product rows as output column indices.
This saves repeated monomial unions, sorting and cancellation. A packed-matrix
cache is an essential competing baseline: on precisely the same key, it can
copy the already known matrix. Comparing only against repeated direct
construction would miss that alternative.

## Fixed experiment

`protocol.json` fixes all dimensions, discovery and holdout seeds, batch sizes,
run order, caps and acceptance criteria before timing. The four variants are
direct construction, verified layout-only reuse, exact product schedules and
an exact packed-matrix cache. Families repeat the same input or alternate it
with a changed constant, a removed term, or a generator with lower degree.
All cache variants retain one initial object; a mismatch falls back and does
not silently populate a larger cache.

Every cold batch pays preprocessing, lookups, fallbacks, fresh output storage,
and complete validation. The independent reference uses exhaustive submask
enumeration and set parity rather than the measured sort/cancel implementation.
Tests exhaust all polynomial supports through three variables, then compare
product evaluation, stale-key fallback and resource-cap behavior.

This is an equivalent **matrix-construction** control, not the repository's
full IC or WDSat regression: it builds no SAT formula and solves no system.
There is no common-unit end-to-end IC or rho measurement. Runtime ratios are
construction diagnostics only. The study cannot establish a cryptanalytic gain.

## Reproduce

```sh
python3 research/boolean_product_schedule_20260922/run.py --out /tmp/boolean-schedule-new-run
python3 /tmp/boolean-schedule-new-run/analyze.py /tmp/boolean-schedule-new-run
```

Use a new output directory each time. The harness freezes its source and
protocol, records compiler identity and actual executable hashes, launches
each fixed cell as a fresh metered process, retains every sample, verifies
stdout hashes and the full expected grid, and writes a SHA-256 manifest.
The report keeps cold and retained costs visible even if the candidate loses.

`run_01` is the original experiment. `run_02` repeats the unchanged grid with
only the schedule's retained column capacity compacted, declared in the
protocol amendment before execution. Both runs reject performance promotion;
the compacted schedule trades a little space for slower replay at 12 variables.
See [the decision](CONCLUSION.md) for both tables and the mathematical result.

Run `python3 research/boolean_product_schedule_20260922/test_evidence.py` to
replay the retained summaries and test fail-closed evidence checks.
