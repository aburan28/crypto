# Goal: function-first decomposition beyond the Semaev controls

Build a function-coefficient solver enforcing H divides L_V that reaches larger
fields at lower total cost per verified relation than the strongest matched
Semaev baseline. A breakthrough is an outcome to establish, not an assumption.

## Acceptance gate

At least 20% lower total cost per verified relation at three increasing sizes,
with identical subspace factor bases, targets, and three-point decomposition
size, and zero correctness failures. Include setup, unsuccessful targets,
branching, extraction, verification, and duplicate/rejected candidates.
Until calibrated operation accounting exists, timing and completion are
preliminary evidence. Full ECDLP costs and rho ratios are a separate gate.

## Milestones

1. Modular-squaring support circuit: implemented and toy-validated in solver_02.
2. Reliable eleven-bit solving: open. Test fresh targets against exhaustive
   ground truth, measuring first relation and complete enumeration separately.
3. Scaling and baseline gate: open. Include chained S3 and symmetrized S4 on
   the identical factor base; aggregate all target costs, including failures.

## Current evidence and next obstacle

solver_02 and solver_03 show that cubic and quadratic support circuits can be
correct while harder for SAT than explicit-root norm equations. The nine-bit
norm results motivate continued experiments; they do not establish root-free
superiority. Eleven-bit norm completion remains incomplete.

solver_04 tests the proved linear-image characterization of quadratic support
conditioned on its linear coefficient. Its coefficient branching and setup
must fit within one total budget. Prior samples are development data; fresh
samples must be identified as such and not silently selected for decomposability.

The repository's existing symmetrized-S4 implementation uses the transformed
coordinate u=1/(x+1) and a different factor base. Directly comparing that
implementation against the x-subspace experiments would violate this goal.
A matched adapter must enforce x in V, preserve exceptional coordinate charts,
and verify exact sums rather than admitting a sum shifted by two-torsion.
The strongest-baseline gate remains open until that adapter is validated.

Frozen runs remain immutable under research/nagao_relations. Each successor
records its contract before execution, source hashes, raw outcomes, scoped
verdict, and scoreboard update. A timeout remains unknown.

## Measured update: solver_04 through solver_06

The hybrid direct quadratic solver cleared all eight eleven-bit development
slots in both first-relation and complete-enumeration modes under three
seconds, with exact oracle equality. A matched elementary-S4 circuit and
transformed-S4 adapter are now implemented and validated. The adapter gap
above is closed for the tested proper x-subspaces; strongest implementation
and operation-accounting claims remain open. Elementary S4 wins at smaller
sizes, so the three-size 20% acceptance gate is not met.

See goal_round_20260913.md for the 288-trial ledger and remaining gates.
The goal is active, not achieved. This file records the research objective;
it does not configure an unattended background job.

## Scaling update: solver_07 through solver_09

Image-space support clears all tested 23- and 29-bit target/mode slots within
five seconds, with a dimension-six base and exact oracle agreement. The
matched S3/S4 controls and field-operation reductions are recorded in
scaling_23_29.md. This establishes a finite-panel field-degree improvement.
The overall goal remains open: larger factor bases, fresh post-development
seeds, and calibrated costs against the Semaev controls are still required.

## Subfield and larger-base checkpoint — 2026-09-14

The [subfield extension](subfield_dimension_results.md) supports non-F2 F4
coefficients and even ambient degrees 18 and 30. Dimensions now increase
through 6, 7 and 8, retaining the exact targets across dimensions. All eight
d8 supplemental enumerations complete in 11.28–20.73 seconds. Five-second
matched coverage is 35/48 for the hybrid and 0/48 for each of chained S3 and
symmetric S4. No correctness failures were found in exhaustive tiny-field
checks or the larger independent oracle/certificate replay.

The five-second full-enumeration gate fails at d8 for both fields. The
quadratic dependence on base size remains. The goal is **still unmet**:
these are stage diagnostics, with no calibrated SAT/field conversion,
required broad regression, complete ECDLP pipeline, or measured rho ratio.
Growing support changes the counting boundary and is not an advance against
a fixed boundary. See the frozen contract for the explicit scope and holdouts.
