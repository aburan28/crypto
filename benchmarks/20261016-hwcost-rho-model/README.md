# 2026-10-16 — `hwcost` calibrated cost model, validated on the Pollard rho iteration

Copies for browsing; the source of truth is
[`research/hwcost_model_20261016/`](../../research/hwcost_model_20261016/)
(protocol, raw run records, plotting source, Typst sources).

- `report.pdf`: results report (claims, interpretation, rho application).
- `run_report.pdf`: descriptive run report of `runs/final` (started 20261016T0459Z).
- `fig_validation`, `fig_errors`: measured vs predicted ns per rho iteration and the error per field and batch size.
- `fig_roofline`: measured read bandwidth vs the bandwidth the kernel drew.
- `fig_opcost`: field-op cycles against limb count, with the size fits.
- `fig_rho`: predicted wall time and cost of parallel rho against group size (extrapolated; illustrative price).

Each graph as SVG and as a 150 dpi PNG from the same plotting source.
