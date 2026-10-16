#set page(paper: "a4", margin: 2cm)
#set text(font: "New Computer Modern", size: 10pt)

= Inverse-fiber primitive run record

*2026-10-16 · Native Rust · Darwin L0 and GitHub-hosted Linux L1*

The prespecified experiment is `PROTOCOL_V1.md`. The public input is
`order1000/BATCH_1237_1201_V2.public.json`, SHA-256
`24f2c56a9ba56f91c954f38b7752edf24fee5ed027b49b1ac0747771945099c4`.
The historical solver is `inverse_velu/FROZEN_RUNNER_V2.rs`, SHA-256
`7587fc43ffbfe6042c3915f90d4c7dc164337fce61d0599776f2b607150660c3`.
The instrumented source byte copy is `FROZEN_RUNNER_V1.rs.txt`, SHA-256
`0731ddc40f26397fdd30bed5098613a3bff91ec3c74c3863c7079d0313b8a250`.

Run command from the repository root:

`bash research/isogeny_torsion_fibers_20261012/primitive_price_20261016/run_panel_v1.sh research/isogeny_torsion_fibers_20261012/primitive_price_20261016/RUN_V1`

The shell launcher ran five arms at prefixes 1, 4, 16 and 64 over three
rotated rounds. `RUN_V1/raw/` retains each child stdout, stderr, exit,
meter stdout and meter JSON. `RUN_V1/build-sha256.txt` binds the binary,
sources, inputs and script. `RUN_V1/ARTIFACT_FILES.sha256` verified after
the run. All 60 exits were zero; `RUN_V1/replay.json` reports 1,275
scalar and 6,126 relation-witness checks passed. The native reduction
command was:

`isogeny_fiber_primitive_price audit RUN_V1 inverse_velu/RUN_V2 PRIMITIVE_GRAPH_V2.svg`

Its V2 output compared all 60 deterministic old/new records. The first
audit attempt failed because its reducer assumed rho had a
`query_group_ops` field. The actual field is `rho_group_ops`; the
failed attempt and corrected V2 output are both preserved beside this
report. The reducer error did not change the measured children.

Independent Linux artifacts came from PR #1939 runs 38961930762
(60-child inverse panel), 38961930725 (108-child nine-arm panel), and
38961930860 (hard-capped order-2,437 full-image probe). All copied
artifact manifests verified. The first two replay checkers passed.
The full-image child exited 134 after an allocation failure under the
512 MiB address-space cap, leaving a resource-censored row.

#image("PRIMITIVE_GRAPH_V2.svg", width: 82%)

This graph records disjoint top-level Fp² multiplication calls at N1.
Square-root internals, polynomial/table control flow, allocations,
I/O and independent-host calibration remain unpriced. The result PDF
and `RESULTS_V1.md` give the per-arm counts and exact limits.
