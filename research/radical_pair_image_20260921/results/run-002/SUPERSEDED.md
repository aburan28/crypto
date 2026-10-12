# Superseded run 002

All 276 solver processes completed and its timing/F4 evidence remains valid.
Bugbot correctly observed that the binary certificate generator asserted
presentation equivalence after checking the direct chain and group law, but did
not independently enumerate the symmetric roots, recover quadratic roots, or
evaluate both binary polynomial presentations.

Run 003 adds those certificate checks and is the accepted run.  Curves, factor
bases, targets, solver inputs/options, repetitions, timeout and execution-order
rotation are unchanged.
