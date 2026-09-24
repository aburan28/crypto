# Prospective half-word syndrome experiment

The fixed four-minor prepass is too costly in the measured discovery cells and
leaves substantially more candidates than exact affine consistency. The evidence
does not justify integrating that implementation into a complete solver. Alternate
minor selection, shared cofactors and truth-domain determinant arithmetic remain
unmeasured; none is assigned a speedup by this note.

A direct representation experiment is available in the strongest retained kernel:
it stores equation syndromes in 32-bit lanes. The byte filters use 8-bit lanes,
and the 16-point schedule's name refers to its number of points, not to 16-bit
equation lanes. No native 16-bit syndrome representation is implemented there.

Projecting onto the first 16 equation coordinates is an exact F2-linear map.
The projected zero test is necessary; when more than 16 equations are present,
every hit must be checked against all original equations. No uniformity assumption
or projected-zero probability substitutes for that check.

A bounded discovery experiment should:

1. Implement scalar and native 16-bit representations with identical assignment
   order. Reuse the exact Gray identities and fixed compiled schedule, not any
   coefficient-dependent output from another system.
2. Verify every initial block and scheduled update against direct polynomial
   evaluation, including the 16/17-equation boundary, the highest retained bit,
   false projected hits, multiple hits, constants and complete UNSAT scans.
3. Verify every projected hit on the complete original system. A failed check
   continues the scan. Count all hits, rejected hits and successful checks. Caps
   remain UNKNOWN, including when a block contains an unchecked projected hit.
4. Charge encoding, narrowing, schedule construction, scanning, original checks,
   validation and workspace release. Compare against both the unmodified 32-bit
   kernel and the existing byte filters in the same binary.
5. Use discovery fixtures first. Only a frozen promising candidate merits the
   full retained corpus and fresh holdouts under the unchanged complete-solve
   dramatic-gain criterion. No source tuning follows holdout timing.

The representation is an unimplemented engineering hypothesis on generated Boolean
systems only. Its full cost and success rate are unknown; it carries no full-IC,
rho or asymptotic claim.
