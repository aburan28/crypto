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
