# IC live-trace integration is blocked by existing main compilation failures

- **Status:** `open`
- **Kind:** implementation gap in the integration baseline
- **Affected work:** [IC live relation output](../docs/ic/README.md#reports-and-resource-accounting)
- **Next action owner:** unassigned; recent quartic/ecbench integration maintainers
- **Last updated:** 2026-10-17

## Problem

The live relation/verbosity change cannot complete validation on main
`c6d64436` because the unchanged library already fails to compile outside
the trace paths. The original requirement remains live, selectable relation
output and factor-base plug-in documentation, with the ordinary solver
results, query streams, counters and JSON reports preserved.

## Evidence and uncertainty

- **Observed:** `cargo check --locked --lib --bin ic` on that integration
  baseline reports 25 errors. `ecbench/workload.rs` defines `Instance::Fp4`
  twice and refers to an absent `CurveSpec::PrimePublicFixture`;
  `ecbench::fbm` and `OnbWeightBase` imports are missing;
  `gaudry_quartic.rs` uses an undefined `sparsest` and
  `echelon_by_lead_with`. Remaining type/match failures follow those paths.
  These files were not changed by the live-trace implementation.
- **Control:** the earlier verbose-output PR #1880 head
  `14bf033e2ee85d05b465f6b99e685243d6b8a6fb`, without the new patch,
  passes `cargo check --locked --lib --bin ic` locally on Rust 1.99.1.
  Applying this PR's patch to that revision permits separate trace
  validation; that is not integrated-main validation.
- **Derived or proved:** no mathematical claim is involved.
- **Suspected explanation:** incomplete integration of recent independent
  changes; the exact responsible merge has not been established here.
- **Missing capability:** a compiling integrated main revision.
- **What this evidence does not establish:** passing checks on the earlier
  control cannot establish that the final current-main combination builds.

## What would unblock it

- [ ] Restore the intended quartic types, imports, dispatcher cases and
  elimination entry points on main without dropping supported methods.
- [ ] Apply the live-trace change and pass
  `cargo test --release --lib --bin ic --test ic_framework --test ic_progress`
  on that exact integrated revision, plus applicable lint/review gates.

These conditions establish integration and regression validation, not a
performance improvement or a security result. Trace I/O has no speed claim.

## Next discriminating step

- **Action:** restore the missing integration wiring in the affected modules
  and rerun the exact checks above.
- **Inputs and controls:** main `c6d64436`; unchanged and patched #1880 head;
  deterministic synthetic fixtures in the existing IC test suites.
- **Execution limits and stop condition:** local compilation and regression
  tests only; no cloud campaign or mathematical search.
- **Outcomes:** a compiling base permits integrated validation; persistent
  unrelated errors keep this record open; trace-specific failures must be
  fixed before the feature is considered validated.
- **Recheck trigger:** a repair to the affected main modules.

## Resolution and history

Unresolved. Creating or merging the feature PR does not satisfy the current-main
validation condition. Its PR description records the exact checks completed
on the separate control revision and any remaining failures.

2026-10-17: with the feature patch on the #1880 control revision, the release
library suite passed 4,297 tests (117 ignored). The IC binary unit tests
passed 10/10, CLI integration 25/25, and pipeline progress 5/5. The initial
verbosity-filter fixture covered every target and therefore had no misses;
it was replaced by a saved dimension-4 base containing both hits and misses,
and the CLI tests passed on rerun. The manual preview emits 9 successful
relations at `-v`/`-vv` and all 64 attempts at `-vvv`, with identical recovered
scalar and counted work. These are correctness and UI checks, not a speed
measurement, and do not resolve the integrated-main build gate.

Strict Clippy for the library, IC binary and both integration targets also
passes on CI's pinned Rust 1.98.0. Rust 1.99's new
`mismatched_bit_width_type` lint instead stops in unchanged
`src/symmetric/modes/pmac.rs:112`; no lint was suppressed for this feature.
