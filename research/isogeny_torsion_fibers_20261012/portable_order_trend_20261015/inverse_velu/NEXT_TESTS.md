# Next discriminating tests after direct inverse-Vélu construction

**Current boundary.** On the frozen `p=1237,r=1201` public cell,
inverse-Vélu construction selects the same eight columns as the full
source census and all 255 paired relation rows match. At one target,
`S_group=149.645` inverse versus 5.107 matched signed rho; at 64
targets, 151.595 versus 4.071. These are counted-group lower bounds,
not complete portable `S`. The [V2 result](RESULTS.md) and
[open cost problem](../../../../problems/nonrational-isogeny-fiber-dlp-cost-gap.md)
state the remaining gap. Each experiment below requires a versioned
frozen spec before execution and belongs in a focused follow-on PR.

**2026-10-16 execution status.** The [four-gate frozen
protocol](FOUR_GATE_PROTOCOL_V1.md) is committed. Gate 1 has a
[60-child disjoint Fp² call panel](../../primitive_price_20261016/RESULTS_V1.md):
1,275 scalars and 6,126 relation witnesses replayed, with 60/60
deterministic matches to the older V2 rows. Square-root internals,
polynomial and table work, allocation, I/O, calibration and an
independent counter checker still prevent a complete portable `S`.
Gate 2 has an
[exhaustively checked direct χ resultant](../../chi_polynomial_20261015/different_j_resultant_20261016/RESULTS_V1.md)
and a [nine-arm held-out solve
panel](../../chi_polynomial_20261015/different_j_resultant_20261016/HELDOUT_RESULTS_V4.md)
after exhaustive source and quotient S3 pair checks. All 108 scalars
and 996 relation witnesses replayed, with identical six-arm relation
support. This completes Gate 2's matched order-229 row controls;
the complete price remains required. Gate 3's independent Linux V2
and nine-arm V4 workflows returned passing L1 host artifacts with
frozen source checks; L2 access remains open. Gate 4 has
[four checked structural cells and 96 complete four-arm
children](trend_four_cell_v1/TREND_RESULTS_V1.md), while portable
pricing remains open. The [Linux p2473 N1 full-image
probe](trend_four_cell_v1/FULL_IMAGE_LINUX_SPEC_V1.json) ran under its
512 MiB address-space cap and exited 134 after an allocation failure;
the resource-censored child has no scalar-replay row. A memory-bounded
full-image constructor is needed for that cell. These are partial
gate statuses, not promotions.

1. **Complete the price of the working inverse oracle.** Extend the
   verified inside/outside Fp² call counts with disjoint counters for
   `sqrt_mod` branch and exponentiation work, polynomial gcd and split,
   table build/lookups, modular row operations, allocation, rejected
   candidates and final replay. Calibrate each native category into a
   common group-addition equivalent with intervals on Darwin and an
   independent Linux host. Recompute whole-child cold `N=1` and reuse
   `N=64` `S` for all five arms from raw counters. Accept this gate
   only when an independent Rust checker recomputes every phase and
   total and the uncertainty interval is reported against same-target
   rho and BSGS. A result above controls closes only this
   implementation's cost question at the selected order.
2. **Completed at order 229: test relation yield with a genuinely different oracle.** Derive
   the nonfactorable `chi`-coordinate decomposition resultant and a
   root procedure that returns all factor-base pairs, multiplicities,
   exceptional fibers and failure codes. First exhaustively compare
   its returned set on the smallest frozen different-`j` cell with an
   independent direct pair enumerator, including zero-decomposition
   cases. Then freeze four unseen public targets and run full relation
   collection, row solve and independently checked scalar recovery
   against direct `S_3`, quotient-invariant, inverse-Vélu and full-image
   controls, plus same-target signed rho, BSGS and rho A/A, with equal
   factor base, random stream, trial cap and memory where the relation
   arms share those settings. Charge every
   polynomial call and failed root. The discriminant is *useful rank
   or verified relations per complete priced work*, not only roots
   found or image-base setup cost. If support is identical, optimize
   cost rather than claim a new relation effect; if support differs,
   independently replay each differing witness before a larger panel.
   The [V4 result](../../chi_polynomial_20261015/different_j_resultant_20261016/HELDOUT_RESULTS_V4.md)
   met the exhaustive and held-out solve gates: all six relation arms
   had identical trials, relations and useful ranks in paired cells.
   A claimed portable cost still requires Gate 1.
3. **Admit wall timing.** Run the frozen five-arm panel interleaved on
   a reserved physical core with host, binary hash, affinity, thermal
   and RSS receipts and same-target rho A/A. Apply the repository's
   A/A noise admission rule before quoting a wall ratio. Current
   screened AWS On-Demand and Spot vCPU quotas are zero in five
   regions; an authorized isolated runner or quota change is the
   concrete resource exit. The inverse V2 and nine-arm V4 panels now
   have independent GitHub-hosted Linux L1 replays, linked from the
   [primitive report](../../primitive_price_20261016/RESULTS_V1.md);
   those virtual cores do not satisfy the L2 physical-core gate.
4. **Measure a credible order trend.** Freeze at least four more
   geometrically spaced prime-to-characteristic subgroup orders in
   this exact construction family, one unseen target per cell and a
   secondary 64-target reuse panel, with the same memory and cap
   contract. Screen curve, kernel, `j`, full-fiber and subgroup
   identities before selecting targets. Run direct inverse, full
   census where feasible, BSGS, signed rho and rho A/A with independent
   replay. Preserve censored time/memory rows. Four additional inverse
   orders now have 96 complete children and independent replay, but
   their selected sizes and heterogeneous field parameters do not
   support an exponent fit. The Linux full-image order-2,437 cell is
   censored at the frozen cap, so a memory-bounded constructor is the
   next implementation requirement for its paired cost.
5. **Test the factor-base scaling mechanism.** Preregister a supported
   grid of base sizes, beginning with 8, 16, 32 and 64 where the
   subgroup has enough signed orbits. Use the same public target and
   selected projected columns for direct-pair and χ decomposition
   controls at each size, with destination-direct and inverse-Vélu
   construction priced separately. Recompute pair support by independent
   enumeration, collect enough relations for useful rank, solve each
   target, and replay every scalar and witness. Charge polynomial calls,
   failed roots, base setup, matrix work and target descent in one cold
   `S` alongside same-target rho, BSGS and rho A/A. The test asks whether
   a growing base changes *priced useful relations per unit work*;
   relation count alone cannot establish an ECDLP crossover. Freeze the
   exact orders, targets, seeds, caps and selection rule in a new protocol
   before running this grid.

The larger-picture target is a **cold unknown-scalar ECDLP result**:
an explicit map and subgroup-preserving transfer, a useful
decomposition or destination solve, independent scalar replay, and a
complete `S=work/sqrt(r)` below the best same-target control over a
credible size trend, followed by admitted isolated wall evidence.
Multi-target reuse is reported separately. Pairing transfer and
hidden-witness access are separate routes with their own
construction, public-equivalent-route and finite-field DLP costs;
these inverse-cubic measurements do not price either route.
