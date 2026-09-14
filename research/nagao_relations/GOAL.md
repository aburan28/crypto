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

## Structured-base and stronger-baseline checkpoint — 2026-09-14

The [structured scaling campaign](structured_01/RESULTS.md) extends prefix
bases through d10 and tests genuine F4-linear, Frobenius4-stable d8/d10
bases at 18 and 30 bits. Its 240 matched cold trials compare the original
hybrid, an early support filter, a direct S3 pair-invariant table, chained
S3 SAT and symmetric S4 SAT. Frozen d8 targets and fresh holdouts are kept
separate. Ten additional eight-target batches charge pair-table setup once.

The stronger Semaev baseline changes the conclusion: in three-second cold
enumeration trials the direct S3 table completes 14/24, while both hybrid
variants complete 0/24. First-relation resolutions are 14/24 for the table
and 13/24 for each hybrid; those counts include proved empty targets. At
18 bits specifically, the hybrids resolve more first-relation trials than
the table, so the enumeration finding is not a universal first-hit result.
Both SAT controls remain at 0/24 in each mode.

The parity filter reduces rejected-function work but does not improve this
panel's completion counts. The hybrid still visits quadratically many
coefficient branches. Ordinary setup caching cannot deliver a 20% saving:
previous setup was below 0.25% of total time. Frobenius-stable support cannot
be quotiented for a fixed arbitrary target without also transforming it.

All ten S3 batches complete, including d10 at 30 bits, but 19 of their 20
uniform 30-bit targets have no admissible relation. Supported targets must
not be used to estimate random-target yield. The pair table retains
quadratic setup and is itself a summation-polynomial algorithm. Its success
does not meet the function-first goal. The 20% calibrated total-cost gate,
broader regression and full-DLP gates remain **unmet**.

## Bound checkpoint and next candidate — 2026-09-14

The [exact bound audit](bounds_01/RESULTS.md) now separates signed-triple mass,
target support, and the current solvers' mandatory branch work. It verifies
265 exhaustive tiny target/space cases, all 104 prior larger inputs, and 60
fresh uniform inputs. Chart exclusions are counted exactly. A group-trace
quotient into E(F64) produces no additional global ceiling improvement and
no empty target class on any of the ten larger bases.

For complete enumeration of the exact same eight-target batches, the hybrid's
branch-only field-multiplication floor exceeds the measured total S3-table
multiplications by 15.73–53.03 times. At 30 bits and prefix d10 the numbers
are at least 58,605,624 versus 1,391,380. This is a component bound on
complete enumeration; it is not a runtime ratio or a first-hit result.

The [next proposal](bounds_01/NEXT.md) is to reject entire coefficient blocks
before conditioned-root solving. A general-coefficient normalized pullback
identity is proved and checked in 2,048 residual tests plus 756 exhaustive
GF64 target checks. Block pruning itself remains unimplemented and unproved.
With the current per-branch multiplication floor and zero pruning overhead,
the 30-bit prefix d10 batch would require at least 98.1007% branch removal
to reach 20% fewer multiplications than the measured S3 table. A new
per-branch formula requires re-deriving this threshold. Full-cost comparison
and the original three-size acceptance gate remain open.

## Certified block experiment — 2026-09-14

The [coefficient-block experiment](blocks_01/RESULTS.md) now implements a sound
binary-linear relaxation of the bilinear support equation at fixed b. Every
rejected block carries an independently checked dual separator. Exhaustive
GF64 validation and 192 matched cold trials find no correctness failures.
This is a coefficient/support hybrid with an auxiliary w witness, not a pure
root-free solver.

The new unpruned bilinear circuit completes 12/16 cold enumerations, versus
0/16 for the old filtered hybrid and 14/16 for the direct S3 table. Pruning
completes 6/16: it helps n30,d8 but loses all six n18,d8 completions. Both SAT
controls remain at zero. First-mode coverage is 15/16 unpruned, 14/16 pruned,
11/16 old hybrid and 14/16 direct S3; these counts include proved empty targets.

The predeclared pruning screen passes on six complete n30,d8 pairs, with
22.6% fewer field API calls and fewer selected binary operations. A fresh
eight-target batch confirms a 25.0% API-call reduction against the unpruned
circuit, but the pruned solver still uses 10,980,028 field API calls against
580,690 for S3. All four arithmetic variants complete the same batch with
four verified relations; the four uniform targets are empty. These API sums
are uncalibrated component diagnostics, not total-operation speedups.

The [new bound](blocks_01/BOUNDS_AND_NEXT.md) applies to circuit construction:
530,432 mandatory multiplications in the pruned batch, versus 207,401 total
S3 multiplications. Perfect free pruning cannot remove that setup floor.
Symmetry alone is insufficient; the next hypothesis is reusing coefficient
maps and elimination work across b values, charging all table/XOR costs.
The original three-size cost goal, broad regression and full-DLP gates
remain **unmet**.

## Coefficient reuse — 2026-09-14

The [complete reuse audit](reuse_01/RESULTS.md) verifies the quadratic
interpolation identity, 288 matched cold cells and the frozen eight-target
n30 d8 batch. All available relation and rejection certificates pass.
An unresolved original runner stall is preserved and the affected S3 cell
is excluded from conservative completion counts and cost ratios; lost
post-maintenance records were rerun and separately identified.

Batch inversion plus reuse lowers multiplications from 655,536 to 87,552,
but field API calls only fall from 10,980,028 to 10,495,212. The strongest
S3 baseline remains at 580,690 API calls. All four variants complete the
same eight targets and four relations. These are uncalibrated components,
not an end-to-end operation gain.

The unchanged rejection traversal alone costs 9,794,908 field additions,
16.8677 times the whole S3 API count. Coefficient-only tuning cannot cross
that architecture-specific floor. The next frozen candidate caches nested,
prefix-independent column spans while preserving the exact relaxation;
its proof and falsification conditions are in [rank_01](rank_01/README.md).
Coefficient preparation must also get cheaper before a 20% API saving
against S3 is possible. The original calibrated three-size goal remains
unmet.
