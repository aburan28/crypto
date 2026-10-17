# 2026-10-16: `hwcost` round 2 (native ISA, `fp64` branches, GPU profile, IC in time and dollars)

Copies for browsing. The source of truth is
[`research/hwcost_round2_20261016/`](../../research/hwcost_round2_20261016/),
which holds the protocol, the raw run records, the plotting source and the
Typst sources.

- `report.pdf`: results report, with claims, interpretation, limits and next steps.
- `run_report.pdf`: descriptive report of every run (started 20261016T1654Z).
- `fig_isa_kernel`: rho kernel time per field, baseline against `-C target-cpu=native`, runs A1 B1 A2 B2.
- `fig_isa_ratios`: native/baseline ratio for every measured arm, with its drift margin.
- `fig_fp64_branches`: `fp64` against `fp64bf`, the masks LLVM compiled back to branches (F1, F3).
- `fig_fp64_cm`: `fp64` against `fp64cm` (`select_unpredictable`), the a priori model, and the disassembly counts (F1′, F2′, F3).
- `fig_gpu`: predicted against measured update rate for 46 ECC2K-130 GPU builds (G1–G3).
- `fig_ic_koblitz`, `fig_ic_prime`: IC against matched rho in wall time, implied §5a weight, and both projections. The projections are extrapolated, not measured.

Every graph is here as SVG and as a 150 dpi PNG, both from the same plotting source.
