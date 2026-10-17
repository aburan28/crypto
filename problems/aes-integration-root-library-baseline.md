# AES integration: root library does not compile at the base

- **Status:** open
- **Kind:** existing implementation/build failure; separate regression-check gap
- **Affected work:** [AES integration](../research/aes_distinguisher_20261017/INTEGRATION.md), required root-library regression check, and root-linked follow-up drivers
- **Next action owner:** unassigned owner of the existing ecbench/Gaudry changes; AES integration session records evidence and keeps its independent checks scoped
- **Last updated:** 2026-10-17

## Problem

The required `cargo test --release --locked --lib` cannot reach test execution
because the unchanged `crypto@c6d64436a33359320ad9d5b42c4acebee2968b9a` root
library has compilation errors in its existing ecbench and Gaudry modules.
The new standalone AES crate does not depend on this root library.

## Evidence and uncertainty

- **Observed:** the [root check log](../research/aes_distinguisher_20261017/checks/crypto-lib-attempt-2.log)
  records duplicate `CurveSpec::Fp4`, missing `ecbench::fbm`, missing
  `OnbWeightBase`/`sparsest`, and inconsistent `PrimePublicFixture`/curve APIs.
  Rust 1.99.1 stable, release, locked dependencies, two build jobs and two test
  threads were used. No AES research trial is involved.
- **Observed:** the first attempt stopped at a missing sparse-checkout path
  dependency. That existing dependency and compile-time fixture files were
  materialized before the substantive second attempt; that access issue is
  distinct from the recorded Rust source errors.
- **Observed:** [base input identity](../research/aes_distinguisher_20261017/checks/base-library-identity.log)
  records zero diff for the root manifest, lockfile, `src`, `isogeny_algos` and
  `volcano-cli`, plus matching base/worktree Git blob identities for the named
  failing files. The integration HEAD was still the unchanged base. New AES
  files live in an independent crate with its own workspace and lockfile.
- **Derived:** these source errors are present in the unchanged base inputs;
  the AES mirror does not edit or link those modules. This does not certify
  that the rest of the root suite would pass after they are corrected.
- **Missing capability:** a compiling root-library baseline with the intended
  ecbench/Gaudry interfaces restored and its ordinary regression evidence.
- **What this evidence does not establish:** full root-library regression
  success, scientific correctness of the existing modules, a novel AES attack,
  or authority to skip a required CI/review gate.

## What would unblock it

- [ ] Correct the base ecbench/Gaudry interface and declaration errors in their
  owning change, without changing or suppressing tests to force a pass.
- [ ] On the repaired baseline, run `cargo test --release --locked --lib` to
  completion, preserving failures/ignored tests and linking its exact head.
- [ ] Recheck the AES PR against that baseline and all applicable current-head
  CI/review requirements before merging.

These conditions establish the missing regression check only. They do not
implement the blocked AES follow-ups or discharge scientific admission.

## Next discriminating step

- **Action:** repair the named base declarations/imports and replay the root
  library check; keep AES source parity and independent native tests intact.
- **Inputs and controls:** unchanged source identities and the complete error
  log above; compare the repaired diff with the intended ecbench interfaces.
- **Limits and stopping:** ordinary local native build/unit-test envelope; no
  cloud experiment or research spending cap is introduced. A new compiler
  error or test failure remains evidence rather than a passed check.
- **Outcomes:** a complete pass closes the regression gap; further failures
  identify the next concrete base defect. Infrastructure failure is separate.
- **Recheck trigger:** a new baseline commit fixing the recorded interfaces.

## Resolution and history

2026-10-17: native AES formatting, Clippy, fourteen tests and eleven-file source
parity pass. All nine blocked backlog proposals pass structural validation.
The required root check remains failed on unchanged base inputs. Unresolved.
