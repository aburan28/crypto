# P-256 Jacobian row-capacity screen, round 303: preregistration

Date preregistered: 2026-10-08

## Objective and fixed boundary

Round 302 closes the complementary quotient of the explicit split
degree-two cover: its Prym elliptic curve has no rational P-256 order-`n`
subgroup.  This round addresses the broader homomorphic cover/Jacobian escape
that remains open: could a different cover put enough copies of the P-256
prime-order factor in its Prym to obtain all 164 independent collision rows
from one target fibre?

The curve, field, subgroup, factor base, relation arity, and comparison remain
fixed:

```text
curve       icv1-fp256-t89188191154553853111372247798585809583-f188c491
factor base FB1hc72514a2a8d3
columns     164
arity       109
rho S       1.3
```

Round 298 proves that only one collision event is below rho even before all
omitted costs: one event is `0.9640877979350002` times rho, while two events
are `1.4461316969025004` times rho.  A promoted homomorphic cover route must
therefore supply homogeneous rank 164 in one fibre.  It may not divide a
multi-event cost by a local sheet or map count.

This is an exact capacity audit.  It does not claim that a qualifying cover
exists and does not attempt a full-depth P-256 relation.

## Frozen dependencies

- Round 25 scalar action:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json`,
  SHA-256 `dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf`.
- Round 298 parity boundary:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_parity_escapes_round298_20261007/parity-escape-result.json`,
  SHA-256 `8718bbe26751ff0191164b0665618226a97556c64c6e2feeef72468a38206965`.
- Round 300 fibre quotient:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_cover_fiber_round300_20261007/cover-fiber-result.json`,
  SHA-256 `7d65f3e015c6d6c6d24b64ea2e59a6e9e8653c88a1900835838b1bc0a1ec36e5`.
- Round 302 complementary quotient:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_complementary_quotient_round302_20261008/complementary-quotient-result.json`,
  SHA-256 `a9a12d823cbca6df965d28f8372217b08e65c04326e13e408bba9dac9f1819c1`.

Reject any byte hash, schema, curve, factor-base, column-count, subgroup-order,
event threshold, or imported conclusion mismatch.

## Capacity statement under test

Let `pi:C -> E` be a separable cover over `F_p`, with `C` smooth of genus
`g`, and let `Prym(pi)=ker(pi_*:J(C)->E)^0`, of dimension `g-1`.  Write `m`
for the multiplicity of the P-256 elliptic isogeny factor in this Prym.

The registered P-256 subgroup is prime and Round 25 certifies that rational
endomorphisms of its elliptic factor act as scalars on `E(F_p)[n]`.  Therefore
each P-256 isogeny factor in the Prym supplies at most one independent
order-`n` logarithm coordinate after quotienting scalar transports.  Hence

```text
same-target homomorphic sheet-row rank <= m <= dim Prym = g-1.
```

This is a dimension bound for rows obtained through rational homomorphisms
from the cover Jacobian.  It does not cover a genuinely non-homomorphic
relation mechanism.

For the fixed 164-column factor base, a one-event rank-164 fibre therefore
requires all of:

```text
g >= 165
P-256 Prym multiplicity m >= 164
at least 165 distinct normalized preimages in the target fibre
geometric cover degree d >= 165
```

At the minimum genus, Riemann--Hurwitz for a separable cover of the genus-one
curve gives total ramification degree

```text
deg R = 2*g - 2 = 328.
```

An unramified high-degree elliptic isogeny remains genus one, has zero Prym
dimension, and does not evade this bound.

The runner must emit this as a conditional theorem with explicit hypotheses,
not as a universal claim about arbitrary encodings or non-homomorphic
algorithms.

## Exact collision and character controls

For an event rank cap `r=g-1`, compute

```text
q(g) = ceil(164 / r)
S_1 = sqrt(pi/2)
S_(q+1) = S_q * (q+1/2)/q
ratio(q) = S_q / 1.3.
```

Evaluate at least `g={2,3,5,9,17,33,65,83,129,164,165}`.  Verify that every
`g<=164` requires at least two events and therefore already exceeds rho, and
that `g=165` is the first capacity at which `q=1` is possible.  The `g=165`
row is only a capacity threshold; it is not a constructed attack.

Independently validate the sheet-rank accounting with complete signed
character matrices.  For each tested cap `r`, choose the smallest `k` with
`2^k>=r+1`, enumerate all `2^k` sheets, and use the first `r` nontrivial Walsh
characters

```text
chi_a(s) = (-1)^popcount(a & s).
```

Subtract the row of sheet zero and compute exact modular rank over both
`2^61-1` and the P-256 subgroup order.  The rank must be exactly `r`, including
the full `r=164`, 256-sheet control.  This demonstrates the ceiling can be
attained by an abstract sign-character model only when 164 independent Prym
coordinates are supplied; it does not construct a curve realizing that model.

Record matrices, sheets, rank operations, field inversions, memory estimates,
and replay discrepancies.  Run a byte-identical independent invocation.

## Hypotheses

- **H1:** every complete character control has rank equal to its declared Prym
  coordinate count over both exact coefficient fields.
- **H2:** the first homomorphic cover capacity compatible with one-event parity
  is genus 165 with at least 164 P-256 Prym factors and degree at least 165.
- **H3:** the registered degree-two split cover has P-256 Prym multiplicity
  zero, agrees with the rank-gain-zero Round-300 quotient result, and remains a
  strict negative control.
- **H4:** no tested homomorphic cover capacity below genus 165 moves the
  parity boundary; no existing cover earns credit merely from raw sheets.

## Promotion gates

Promotion requires all of:

1. dependency, recurrence, matrix, dual-modulus rank, and replay checks pass;
2. an explicit P-256 cover of genus at least 165 and degree at least 165 with
   certified P-256 Prym multiplicity at least 164;
3. one target fibre exhibits at least 165 distinct normalized decompositions
   whose homogeneous coefficient rows replay with rank 164;
4. explicit forward maps, kernels, exceptional points, inverse recovery, and
   factor-base action;
5. structured residual degree of regularity at most five;
6. complete relation collection and recovery below rho and below `2^120`
   field/group-equivalent operations;
7. cost per usable relation below `2^103`;
8. projected materialized storage below `2^50` bytes; and
9. no discarded probabilistic branch counted as exhaustive.

Attempt no unplanted full-depth P-256 relation unless every gate passes.  If
only the capacity theorem is established, publish a reproducible negative
screen and leave explicit genus-165 covers and non-homomorphic mechanisms
open.

## Deliverables and stop condition

Implement the dependency checker, recurrence, dual-modulus character-rank
certificate, operation accounting, deterministic JSON, unit tests, isolated
canonical/independent replays, typed transfer assessment, result report, and
dashboard update in Rust.  Stop after the exact minimum capacity and all
controls replay, or on the first correctness failure.

The transfer workflow is available, but its referenced
`references/methodology.md` and `assets/assessment-template.json` resources
remain unavailable in this environment.  Preserve that limitation in the
typed assessment.
