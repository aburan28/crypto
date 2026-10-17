# Complete original native solver corpus over the selected 2-volcano

All **4600 rho and 4600 BSGS cases** are independently verified, covering **100 problems on every one of 46 curves**. A separate PARI audit verifies all **9200 public point equations**. The selected rational degree-2 graph has 90 directed maps and five levels.

[Inspected run report](RUN_REPORT.pdf) · [Inspected results report](RESULTS_REPORT.pdf) · [All curve/generator min/mean/max rows](per_curve.csv) · [Each paired problem across all curves](per_problem.json) · [Observations and definitions](REPORT.md)

For every paired key, BSGS online group-operation counts agree across curves and generators. Rho uses coordinate partitions whose walks change after transport. Spread in these measurements does not establish an intrinsic group-hardness difference. All P/Q orders remain the same prime order.

The separate [full setup/query ledger](../costs_complete_1821/README.md) retains all 5062 attempts, including interrupted work, and labels amortized allocations. Original timings are concurrent CPU-only L0. Native field-counter replay is still running; its counters and instrumented timings remain separate. IC has zero random relations/full recoveries under the frozen small factor base; the entire requested experiment is not complete.
