# Rejected measurement

Do not use this run to claim performance. It is retained to preserve the failed
measurement attempt.

- The shared Cargo target directory reused the baseline executable for both
  stage/DLP labels. Their identical SHA-256 hashes invalidate those comparisons.
- The kernel test did contain both implementations, but its linear xorshift
  stream produced synthetic dense matrices with rank at most 64. Their names
  therefore overstated the intended workload. Macaulay inputs were unaffected.
- Run 002 forces clean package builds and rejects identical binary hashes. Its
  SplitMix64 generator produces the intended dense ranks; the comparator rejects
  dense ranks more than four below the smaller matrix dimension.

The production optimization and acceptance threshold were unchanged. All raw
outputs and the original contract remain here. The corrected comparator marks
this run rejected in `summary.json`.
