# Non-rational torsion fibers: first bounded check

Date: 2026-10-12. Work item: T-2070. This protocol precedes the run.

## Question and hypothesis

The existing `RESEARCH_ISOGENY_FIBERS_ECDLP.md` says in its summary and
§2.7 that the degree `d` of a separable isogeny is always a divisor of the
cofactor `h = #E(F_q)/r` when `r` does not divide `d`. The statement needed
for relation projection is narrower: the **rational kernel** has order
dividing `h`. A Frobenius-stable order-3 kernel with no nonidentity rational
point should give `d = 3`, `h = 1` on a prime-order curve.

## Frozen input and reference

Enumerate every prime `p <= 199` with `p > 3` and `p ≡ 1 (mod 3)`, and
every nonsquare `b` in `1..p`. Use `E_b/F_p: y² = x³ + b` and the
Frobenius-stable kernel `M = {O, (0, sqrt(b)), (0, -sqrt(b))}`.
For each curve count `#E_b(F_p)` exactly by two independent native methods:
direct affine point enumeration and quadratic-character summation. Check
the Hasse interval and record rows for which the order is prime.
The kernel has order 3 because the tangent at each `x=0` point is horizontal;
Frobenius exchanges the nonzero kernel points because `b` is nonsquare.
Thus `M(F_p) = {O}`. The quotient isogeny exists over `F_p`; no isogeny
evaluation or DLP operation is needed for this counting control.

## Success, failure, and stop conditions

Success for this correction is at least one prime-order row, with both
counts equal, where `3` does not divide `h = 1`. A zero count leaves the
general group-theoretic correction intact but does not supply this witness.
Any counting disagreement or failed Hasse check is a validation failure:
retain raw output and stop promotion. Scan exactly the stated finite set;
do not extrapolate the frequency. Record the source revision and exact
command, compiler, output, and exit status.

This is a structural check. It supplies no solver, timing, destination
advantage, or end-to-end `S`. The comparison reference for any later
speed claim is matched single-target rho on the same subgroup, with setup,
transport, relations, linear algebra, target recovery and verification
charged. A follow-on performance experiment requires a separate frozen
protocol and native execution path.
