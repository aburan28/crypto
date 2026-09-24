# Diagnostic refinement after the first run

The original results.json is preserved. Its saturation degree is six in
every system. This metric has an ambient-dimension floor: any nonzero ideal
in the six-variable Boolean function algebra contains a point indicator of
degree six, while rows admitted before D=6 have degree at most D. Hence full
ideal saturation cannot occur before six. Its observed constancy is not
evidence that all systems have equally difficult Groebner computations.

Before the next run, add two separate measurements to the same fixed matrices:

1. Earliest D where the constant polynomial 1 belongs to the row span,
   recording null for satisfiable instances.
2. Earliest D where the dimension of all linear-or-constant consequences
   equals its final value. Compute intersection dimension as rank(M) minus
   rank(M projected onto columns of degree greater than one).

Freeze all curves, bases, monomial orders, and rows exactly as in PROTOCOL.md.
Save this refined run as refined.json. Require the original fields of every
record and every prior counter to match results.json. Verify contradiction
membership agrees with the exhaustive zero count. The new projection-rank
calculations are diagnostic overhead and excluded from the original row-XOR
counter; no total-cost or speedup claim is permitted. These are certificate
and linear-consequence thresholds for a specified Boolean filtration, not
F4/F5 solving degrees or degree of regularity.
