# Stage 185 protocol: selected dense-pair default replay

## Purpose

Stage 184 accepted dense exact pair selection for current five-column
`BlockTables`. Commit `8014149a2` makes dense selection the default while
retaining `F4_F2_DENSE_PAIR_SELECT=0` as the quadratic control. This stage
validates that selected policy from an exact clean build before publication.

## Frozen validation

1. Clean release-build `koblitz_pdp_backend` from selection commit
   `8014149a2` under the process meter.
2. Run all Boolean-F4 tests with the environment unset (selected dense default)
   and with `F4_F2_DENSE_PAIR_SELECT=0` (quadratic control).
3. Run all backend tests in both modes.
4. Run the frozen `n=59, ell=9, m=3` target once with no pair-selector
   environment variable, twelve Rayon workers, X1 batch 512, and current
   BlockTables. It must report dense selector counters, zero full-M4RI
   matrices, exact current logical/performed XORs, exact exhaustive UNSAT, and
   the frozen equation fingerprint.
5. Run direct MITM from the same binary once for an updated same-host stage
   ratio. This is a decomposition reference, not a full attack or rho result.

Charge the build, four validation commands, selected replay, and direct arm.
Report wall, core-seconds, RSS, conflicts, and single-core validity. Any routing,
correctness, hash, timeout, or test failure rejects the selected default and
requires reverting the selection commit.

Passing this stage validates only the repository implementation default on one
public target. It does not satisfy natural relation yield, full unknown-scalar
recovery, the external solver panel, matched rho, independent reproduction,
novelty review, or SOTA.
