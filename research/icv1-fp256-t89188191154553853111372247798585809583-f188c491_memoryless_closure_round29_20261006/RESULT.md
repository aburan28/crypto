# P-256 exact-memoryless factor-base closure, round 29: result

Date run: 2026-10-06

No proper signed fixed-S17 factor base in the registered exact-memoryless class
reaches rho parity.  The one-operation transition that coalesces is addition of
a representation-independent group element, which is ordinary partitioned
Pollard rho and does not preserve a proper signed factor base.  The transitions
which preserve the factor-base representation are either non-coalescing local
updates, scalar actions at least 228.548 times over rho's operation budget, or
coordinate correspondences that retain the independent-log quotient at least
392.155 times over rho.

This is a negative closure result for the registered class, not an
impossibility theorem for every conceivable future ECDLP algorithm.

## Why the one-operation row collapses to rho

Let `F=-F` be a nonempty signed factor base in the prime-order P-256 group.
If changing one represented column from `P` to `u(P)` is to factor through the
represented group sum, then every representation in one transition cell must
have the same delta

```text
u(P) - P = C.
```

For nonzero `C`, closure under the update makes `F` invariant under translation
by `C`.  Since the group has prime order, that translation has one orbit of
length

```text
115792089210356248762697446949407573529996955224135760342422259061068512044369.
```

It therefore forces the complete group, including the identity, rather than a
proper factor base.  For `C=0`, the update does not mix columns.  Replacing a
bounded tuple changes nothing: if the tuple's total delta is independent of
the representation, the represented state executes `S -> S+C`.  Choosing `C`
from a partition of `S` is Pollard rho with redundant tuple state.

The implementation exhaustively checked this implication on every nonempty
signed subset for prime groups 3, 5, 7, 11, 13, 17, 19, 23, 29 and 31:

| complete reference count | value |
|:--|--:|
| signed nonempty subsets | 52,068 |
| nonzero-translation cases | 1,501,168 |
| affine cases | 45,069,808 |
| invariant proper sets with nonzero translation | 0 |
| false positives / false negatives | 0 / 0 |

Ordered reference digest:
`b1ad06f3258f76c33dad7b089c40315eb1628429784143d4531f8ee222f65944`.

## Algebraic transport closes to scalar actions

An elliptic-curve morphism which sends the identity to the identity is a group
homomorphism.  On the prime-order P-256 subgroup it acts as a scalar.  Adding a
fixed point gives an affine map `f(x)=lambda*x+b`.  For negation `nu`, direct
composition gives

```text
f o nu o f^-1 : x -> -x + 2b
(f o nu o f^-1) o nu : x -> x + 2b.
```

Thus a signed base invariant under an affine action with nonzero `b` is also
invariant under a nonzero translation and must be the complete group.  A
proper signed invariant base has `b=0` and is a union of scalar orbits.

The binary independently replayed 4,096 hash-selected affine coefficient
controls, including the inverse, conjugated-negation and generated-translation
identities, with zero failures.  Control digest:
`985b3e7ccb864a204a064f2a8b6864a0c973ba793daa81ec0aed7129e6a75006`.

The imported complete P-256 scalar census contains 10,240 divisors, 79
eligible orders, 1,712,366 exact-order elements, 856,183 folded actions and
12,712 width-band combinations.  Across every action with at least 17 folded
columns:

| objective | folded width | operation floor | corrected ratio to rho |
|:--|--:|--:|--:|
| minimum binary chain | 26,483 | 340 before target correction | at least 328.839 for the corresponding width-band construction |
| minimum doubling floor | 9,301 | 236 + 1 target correction | **228.548** |

The doubling row is deliberately optimistic.  It ignores additions above the
bit-length floor and still misses parity by more than two orders of magnitude.

## Terminal boundaries

| terminal class | preserves fixed S17 | coalesces | exact quotient | non-generic | ratio / rho | decision |
|:--|:--:|:--:|:--:|:--:|--:|:--|
| partitioned group translation | no | yes | yes | no | 1.000 | Pollard rho reference |
| known-log `K=1` collision oracle | yes | no | yes | no | 0.964 | `2^133.326` materialized bytes; memoryless form is rho |
| one-addition scalar-orbit neighbour | yes | no | yes | yes | 0.964 oracle-only | distinct deltas prevent coalescence |
| proper signed scalar action | yes | yes | yes | yes | **228.548** | cost gate fails |
| best-width independent-log quotient | yes | no | no | yes | **392.155** | quotient, collection and storage fail |

The isogenous-model screen does not open another row.  All 4,097 models retain
the universal degree-four leading support of short-Weierstrass `S3`; the
degree-11 pullback is not a degree reduction, and no structured residual
degree of regularity at most five was measured.

The best imported fixed-row collection remains `2^137.363459` operations,
`2^120.288827` per usable relation and `2^142.502917` bytes.  The collection,
per-row, storage, degree, complete-cost and non-generic gates fail.  All
imported complete false-positive, false-negative, transport and replay checks
pass.

## Decision and scope

Reject exact-memoryless bounded-update and algebraic-transport selectors for a
proper signed fixed-S17 P-256 factor base.  Within this class, parity is
attained only by dropping the factor-base constraint and using Pollard rho.
Continuing toward a genuine index-calculus parity result now requires leaving
at least one registered assumption: use an unbounded/non-memoryless
representation solver, abandon fixed signed S17, or discover a relation
algorithm outside exact algebraic group transport.  Screening more geometric
factor bases without such a mechanism cannot fix the end-to-end obstruction.

No full-depth unplanted relation was attempted because the promotion gates did
not pass.

## Reproduction

```bash
cargo test --bin p256_memoryless_factor_base_closure
cargo run --release --bin p256_memoryless_factor_base_closure -- \
  --round24 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_negation_cross_colour_round24_20261006/negation-result.json \
  --round25 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json \
  --round26 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_selector_coalescence_round26_20261006/coalescence-result.json \
  --round27 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_invariant_round27_20261006/affine-invariant-result.json \
  --round28 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_isogenous_dickson_round28_20261006/isogenous-dickson-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_memoryless_closure_round29_20261006/memoryless-closure-result.json
```

Canonical result artifact: `memoryless-closure-result.json`, 14,905 bytes,
SHA-256
`9be6e5fc38644c34ad746dec9a6542da7e35ef2d8e11eb39a83de7af0f54f443`.
