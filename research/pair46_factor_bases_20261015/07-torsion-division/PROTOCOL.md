# 07-torsion-division: preregistered protocol

Pair-46 factor-base campaign, item 07. Frozen before any data is generated; not
edited afterwards (sha256 recorded in `PROTOCOL.sha256`).

Source idea: IDEA-20260807-0a3272 (2-torsion root count r as a class-varying,
gauge-invariant quantity) and EV-ICINV-6b6264 / EV-ICINV-8afb64 (stopped at the
SR3 backend gate: no Groebner/Singular/Macaulay2 backend; later run found the only
varying family was the 2-torsion factorisation type). This item asks the
factor-base version: can torsion points or division polynomials define a factor
base whose membership degree is far below its size, or whose relations are
better than an unstructured base of the same size?

## Backend gate: route-around (minimum engineering)

No Sage/Singular/M2/Magma on this host (checked: `which sage Singular M2 magma gp`
empty). We do not need a Groebner backend: the relation ideal for m = 2 is
I_R = <S_3(x1, x2, x(R)), g(x1), g(x2)> with g univariate. Its F_p-points are
found exactly by the univariate route (F_p-roots alpha of g, then the quadratic
S_3(alpha, x2, x(R)) in x2, then g(x2) = 0). The membership polynomial g is
built explicitly from division polynomials in a standalone Rust binary (AGENTS.md:
no Python for computation). Groebner-basis invariants (regularity, Betti tables)
are NOT measured; degrees are reported exactly.

## Claim class

- Exact statements (proved in RESULTS.md, checked computationally):
  E1. If #E(F_p) = N is prime, E(F_p)[m] = {O} for N not dividing m; torsion
      subgroups give no factor base on E.
  E2. The F_p-roots of psi_m are x-coordinates of E(F_p)[m] and of the quadratic
      twist E^t(F_p)[m]; for pair 46 all such roots lie on the twist (none on E),
      so V ∩ x(E(F_p)) is empty for every m < N.
  E3. A sum of l-torsion points is l-torsion: a torsion base on a curve with
      N = l*q covers only O in the order-q subgroup.
  E4. For gcd(m, N) = 1 the preimage base [m]^{-1}S is a group-automorphism
      transport of S: its relation set is exactly [m]^{-1} of S's, coverage is
      identical, and its membership polynomial prod_s (phi_m - s psi_m^2) has
      degree m^2 |S| with |S| F_p-roots on E (ratio size/degree = 1/m^2).
- Empirical (toy): coverage of each arm, size/degree ratio of psi_m level sets.

## Metric (same instrument on every arm)

- M1 coverage C(V) = |{ ±P1 ± P2 : P1, P2 in V }| / N, exact enumeration over
  discrete logs on the whole group (m = 2 decompositions, P1 = P2 allowed).
- M2 rho(V) = |V| / deg g_V (membership-polynomial degree): explicit list
  prod (x - x_i) has rho = 1.
- M3 (torsion arm) subgroup coverage |T ∩ <G_q>| / q.
- M4 (ideal check) number of F_p-solutions of I_R by the univariate route equals
  the log-table relation count for R, on 20 targets per curve (consistency).

## Arms (per toy prime-order curve; |V| = B = ceil(sqrt(N)/2))

- SX: small-x base (first B on-curve x values from x = 0 upward).
- RND: B random on-curve x values (seeded) -- NEGATIVE/NULL control.
- KL: known-log progression {kG : k = 1..B} -- POSITIVE control (structure that
  collapses coverage to about 4B/N; shows the instrument sees structure).
- PRE_m, m in {2, 3, 5, 7}: x([m]^{-1} P) for P in SX, built from division
  polynomials by evaluating phi_m/psi_m^2 at every x in F_p (cross-checked
  against the log table).
- LVL_m, m in {3, 5, 7}: level sets {x on E : psi_m(x) = c}, 50 seeded c per curve.
- TOR_l, l in {3, 5, 7}: on separate toy curves with N = l*q (q prime), base =
  F_p-roots of psi_l on E; POSITIVE control for the root-count detector.

## Sizes and seeds

- Toy prime sizes: 12, 14, 16 bits; 5 prime-order curves per size (seeds 0..4);
  TOR curves: 3 per (size, l). PRE_m explicit polynomial g built at 12 bits only.
- Pair 46: crater and floor curves; F_p-root count of psi_m for m in
  {2, 3, 5, 7, 11, 13} via deg gcd(psi_m, x^p - x); twist exponent check with 20
  seeded twist points: [5v^2]Q = O (order) and whether [5v]Q = O (crater
  predicted Z/5v x Z/v: always; floor predicted cyclic: generically not).
- Bootstrap: 2000 resamples over curves, seed 7, 95% percentile intervals.

## Predictions (stated before data)

- PRE_m coverage minus SX coverage = 0 exactly on every curve (E4).
- RND vs SX: no material difference (|diff| within bootstrap CI of 0 or small).
- KL coverage << SX coverage (positive control).
- LVL_m mean size about 0.5 to 1 (random-polynomial root count, half on E),
  rho << 1.
- TOR_l: #roots = (#E(F_p)[l] - 1)/2 exactly; subgroup coverage = 1/q.
- Pair 46: psi_m F_p-roots = 0 for m in {2,3,7,11,13}; psi_5 has exactly 2 roots,
  both twist x-coordinates.

## Stop rule

Run all arms once (deterministic seeds). Hard cap 40 min wall on 1 core; if
exceeded, cut the 16-bit size first, then m = 13 on pair 46, and report the cut.
No reruns to change outcomes.
