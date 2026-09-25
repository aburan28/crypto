# Protocol: oriented binary Vélu transport

Registered before execution, 2026-09-24; follows merged PR #703.

## Hypothesis and scope

For an odd cyclic subgroup of an ordinary binary curve
y²+xy=x³+a x²+b whose nonzero kernel abscissae are rational, paired
Vélu sums yield an oriented full-point map over the base field. A generator
may be supplied on the curve or a quadratic twist with the same b.

For half-kernel abscissae u and c=u/(x+u), the normalized formulas are:

    X = x + sum(c+c²)
    Y = y + sum(u*(x²+y)/(x+u)² + c³+c)
    b' = b+t+t², t=sum(u), a'=a.

Infinity and valid kernel inputs map to infinity. Off-curve or malformed
inputs and invalid cyclic kernels are rejected. Rational abscissae are a
declared limitation; arbitrary nonrational kernel polynomials are not supported.

## Frozen inputs and reference

- Existing two-seed degree-263 representatives from
  research/ecc2k130_direction_review_20260924/twist_torsion_results.json.
- Independent reference arithmetic and direct extension-field Vélu sums:
  redteam_velu_replay.py from PR #703.
- Exhaustive odd-degree binary toy fields selected before execution:
  F_2^3 (x³+x+1), F_2^5 (x⁵+x²+1), F_2^7 (x⁷+x+1).
- Source and twist families a=0,1 with b=1; enumerate each rational subgroup
  of odd order found, choosing the largest odd cyclic subgroup from each
  toy curve, and verify every curve point and all point pairs where affordable.
- Exact-target public P,Q; planted scalars 0,1,-1,2,3,17,263.
- Independent deterministic validation seed 2026092403.

## Success and stop conditions

Pass only with zero disagreements:
- constructor verifies subgroup order, distinct half-kernel abscissae and
  compatible curve/field parameters;
- infinity, kernel, negatives, the 2-torsion point and malformed inputs behave
  explicitly;
- every tested output lies on the stated codomain, and full point additivity
  agrees on exhaustive toy pairs;
- explicit dual/reverse composition on toy fields agrees with multiplication
  by the degree up to a documented normalization isomorphism when constructed;
- all eight saved exact-target representatives match independent direct Vélu
  full coordinates, and retain public P,Q subgroup order;
- replacing a Y coordinate by its other lift triggers the sign-control failure.

Stop a run at 180 seconds and preserve failures. Record source/input hashes,
exact test counts and failures, command, host, elapsed time and diagnostic
native counters if available.

## Accounting and interpretation

This is a correctness and reusable-interface prerequisite, not a performance
iteration. No solver, relation yield, matrix work or DLP solve is measured.
Full ECDLP cost and speedup remain null. Existing frozen solver performance
gates do not apply to this new arithmetic interface; a downstream comparison
must supply its own matched end-to-end protocol.

## Reproduction

    python3 research/ecc2k130_oriented_transport_20260924/validate.py --out /tmp/oriented-velu-validation.json

The final PR includes code, replayable tests, compact validation receipt and
the outcome. Older saved preflight evidence is not overwritten.

## Additional generic-interface controls (registered before their execution)

The initial fixed corpus passed. Before final validation, also cover an odd
composite kernel on a non-unit b: enumerate b=2,...,31 on a=0 over F_2^5,
select the first b whose rational group has an odd composite cyclic subgroup,
and verify that subgroup's quotient against direct rational full Vélu sums
on every point and against all source point-pair additions. This selection is
solely a coverage rule, not a performance or curve-quality search. Add explicit
incompatible-field, incompatible-b, wrong-order and composition-model rejection
checks. The frozen exact-target corpus remains unchanged.
