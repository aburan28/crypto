# 2026-10-16: branch-free selections in production `Fp64` (`fp_backend` and `prime_fast`)

Copies for browsing. The source of truth is
[`research/fp64_branchfree_20261016/`](../../research/fp64_branchfree_20261016/),
which holds the protocol, the raw run records, the plotting source and the
Typst sources.

- `report.pdf`: results report, with claims, findings, the cost noted and limits.
- `run_report.pdf`: descriptive report of runs A1 B1 A2 B2 (started 20261016T2319Z).
- `fig_kernel`: `hwcost` rho kernel time per field, baseline against candidate binaries.
- `fig_ratios`: candidate/baseline ratio for every arm, with its drift margin.
- `fig_pfrho`: every sample of the production `prime_fast::rho_parallel` (44-bit curve, one thread).
- `fig_disasm`: conditional jumps and `cmov` in the kernel step and in `prime_fast::rho_worker`.

Every graph is here as SVG and as a 150 dpi PNG, both from the same plotting source.
