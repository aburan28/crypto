# Experiment Contract: ECC2K-130 rational twist torsion preflight

Written before execution, 2026-09-24. This is a representation feasibility
experiment, not an index-calculus speed experiment.

## Hypothesis

For E: y²+xy=x³+1 over q=2^131, q-Frobenius is -1 on E[263].
Consequently E_tw: y²+xy=x³+x²+1 has all 263-torsion rational over F_q.
Projecting random twist points by H=#E_tw(F_q)/263² gives a directly
constructible two-dimensional 263-torsion space, whose 264 lines split
under 2-Frobenius into orbit lengths 1,1,131,131. Kernel abscissae give
262 distinct descending codomain j-invariants and two horizontal loops.

## Null hypothesis and falsification

The integer recurrence, group arithmetic, independent basis, Frobenius
action, kernel or codomain checks disagree. Any such disagreement fails
the run; it does not become evidence of a new curve structure.

## Parameters and boundary

- Exact public target q=2^131, reduction z^131+z^13+z²+z+1, a=0, b=1.
- Twist a=1, b=1; deterministic seeds 20260924 and 20260925.
- Degree 263; at most 180 seconds per run.
- Existing FastGF2m and Koblitz Python arithmetic imported read-only.
- Independent extended-Euclidean inversion checked against existing
  exponentiation inversion before use.
- Reference: prior proposed route constructs/factors the degree-34584
  division polynomial. This preflight measures no comparative speed.
- There is no PDP factor base, relation collection, linear algebra, or
  recovered challenge logarithm. Those costs and speedup stay null.

## Controls and checks

1. Confirm Frobenius recurrence agrees with Certicom's subgroup order.
2. Verify v_263(#E_tw)=2, H coprime to 263, and points satisfy a=1 curve.
3. Verify two nonzero order-263 points are independent by exhaustive
   enumeration of the first cyclic subgroup.
4. Compute 2-Frobenius matrix on that basis by bounded lookup; verify
   its quadratic characteristic relation and permutation on all 264 lines.
5. Enumerate each line's 131 distinct abscissae, check kernel closure,
   and check abscissa squaring maps exactly to the predicted next line.
6. If quotient parameters are computed, use the existing characteristic-2
   Vélu derivation b'=1+t+t² with t the half-kernel abscissa sum. Verify
   Frobenius covariance of b', count distinct j, and test the explicit
   x map at several planted scalar multiples of the public generator.
7. Negative controls: generator belongs to a=0, not a=1; q-Frobenius
   characteristic on twist is X²-X+2 (not X²+X+2); kernel orbit labels
   are not claimed to be distinct codomains without the quotient checks.

## Metrics and success criterion

Preserve seed, timestamp, dependency hashes, source hash, exact integer
constants, point coordinates, Frobenius matrix, orbit lengths, kernel
hashes, optional codomain parameters and verification results, operation
counters, and elapsed time. Success requires all checks on both seeds.
This earns only exact-target structural evidence; no speed claim follows.

## Reproduction

```sh
python3 research/ecc2k130_direction_review_20260924/twist_torsion_preflight.py --out /tmp/twist_torsion_replay.json
```
