# Stage 179 protocol: six-column BlockTables on the frozen target

## Hypothesis

The restored current F4 engine performs `147,794,583,858` table-assisted word
XORs on the frozen target. The repository already carries an opt-in
`f4-wide-tables` implementation that tabulates six pivot columns at a time;
five columns are the measured general default. Although six columns regressed
on smaller generic cases, the single `n=59, ell=9, m=3` Phase B workload has
large matrices and may amortize the larger tables.

This stage changes no source algorithm and does not revive either rejected
reducer index. It compares two binaries built from one exact source commit:

- control: default five-column `BlockTables`;
- candidate: `--features f4-wide-tables`, six-column tables.

## Correctness and accounting

- Run the nine Boolean-F4 tests for both feature configurations.
- Both binaries must return exhaustive UNSAT on the same authenticated target,
  visit all 512 masks, complete all 242 systems, and reproduce equation
  fingerprint
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- Basis, pair, matrix, divisor, row-equivalent XOR, and extraction counters
  must agree exactly. Actually performed XORs and table memory may differ.
- Run both builds under `scripts/process_meter.py`; charge build wall,
  core-seconds, and RSS separately.
- Run three interleaved target pairs in fixed order:
  `five, six, six, five, five, six`.
- Query controls are twelve Rayon workers, X1 batch 512, 300-second internal
  budget, 360-second watchdog, and one thread for every external BLAS/OpenMP
  runtime.

## Decision rule

Select six columns for this Phase B adapter only if all correctness gates pass
and the median paired six/five ratios are both below `0.97` for wall and total
core-seconds. Report performed-XOR and RSS ratios; do not hide a memory cost.
The repository-wide default remains a separate decision because the existing
general suite selected five columns.

Charge every build, test, candidate, control, timeout, and failure. Single-core
time remains null here. A target-specific solver-stage gain is engineering and
cannot establish relation yield, a full unknown-scalar run, a rho crossover,
independent reproduction, novelty, or SOTA.
