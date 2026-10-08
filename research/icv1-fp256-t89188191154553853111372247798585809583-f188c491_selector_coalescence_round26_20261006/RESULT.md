# P-256 scalar-orbit selector coalescence, round 26: result

Date run: 2026-10-06

The scalar-orbit base can reach the rho collision constant only by ceasing to
be an S17 selector.  A local neighbour update costs one group addition and
preserves the 17-term representation, but the two colours have disjoint legal
deltas at every valid distinct-column collision, so their walks cannot
coalesce.  A global orbit action coalesces optimistically, but its cheapest
nontrivial scalar costs 344 binary group operations per sample.

This closes the remaining low-cost memoryless-selector route for
`FB1he6b6b6e25de6`.  No full-depth unplanted P-256 relation was attempted.

## Exact local obstruction

For the full signed orbit `A_j=[m^j]H`, with order `r=279184=2B`, the forward
one-atom move has delta

```text
D_j = A_(j+1)-A_j = [(m-1)m^j]H.
```

The complete census enumerated all 279,184 deltas:

| check | result |
|:--|--:|
| signed deltas | 279,184 |
| distinct deltas | 279,184 |
| duplicate deltas | 0 |
| antipodal pairs `D_(j+B)=-D_j` | 139,592 / 139,592 |
| delta blob SHA-256 | `9bef47854a134951a3970dc45b31af2153ffea5b37dca9562577310fde01f9a4` |

At a cross-colour collision, a left move adds `D_j`; moving a right base atom
forward changes `R=T-B` by `-D_k`.  Injectivity and antipodality give

```text
D_j = -D_k  iff  j = k+B,
```

which means the two atoms are the same folded column.  The S17 contract
requires all 17 columns distinct, so the legal left and right next-step sets
are disjoint exactly when a valid relation is found.

The deterministic panel constructed 4,096 hash-selected planted S17
collisions.  All 69,632 signed atoms used distinct folded columns.  Every
scalar relation and every P-256 group relation replayed, with zero false
positives, false negatives, or group failures.  The measured number of legal
delta intersections was **zero**.

Thus a one-addition local transition depends on hidden support.  Equal group
points do not take the same next step, so Pollard/Floyd cycle finding cannot
observe the collision.  Materialising the collision table restores
correctness but restores the round-25 `2^133.326`-byte projection.

## Global transition screen

All 279,184 subgroup elements were enumerated, then identity/negation pairs
were folded to 139,591 nontrivial actions.  The cheapest action is the original
orbit generator:

```text
lambda = 379483656589291393959605088584129393463128953089898693989044430844151895048
bits = 248
Hamming weight = 96
binary double/add operations = 344
unavoidable doubling lower bound = 247
folded permutation order = 139592
```

All 4,096 deterministic signed-S17 rotation controls replayed both as scalar
identities and P-256 group equations.  There were zero failures.  This screen
is optimistic for a target-shifted right colour: rotating `R=T-B` also needs a
target-affine correction, which is not charged below.

## Parity table

Round 24's known-anchor oracle constant is `S=1.253637`; matched rho is
`S=1.3`.  Parity therefore permits at most 1.036982 group-equivalent
operations per accepted sample.

| transition | operations/sample | coalesces | preserves fixed S17 | / rho | storage | class |
|:--|--:|:--:|:--:|--:|--:|:--|
| local neighbour | 1 | **no** | yes | 0.964336 oracle-only | `2^133.326` bytes when materialised | hidden-state selector |
| global rotation, binary cost | 344 | optimistic | yes | **331.732** | small | reject |
| global rotation, doubling-only floor | 247 | optimistic | yes | **238.191** | small | reject |
| unrestricted group translation | 1 | yes | **no** | 0.964336 | small | Pollard rho |

The last row is the apparent parity result.  It updates the group point and
tracked coefficients with a precomputed translation, but it no longer
maintains eight and nine distinct factor-base atoms.  It is ordinary Pollard
rho, not index calculus.

## Decision

No scalar-orbit selector passes.  The exact dichotomy is:

1. preserve fixed signed S17 locally, pay one addition, retain hidden support,
   and lose coalescence; or
2. make the transition a function of the group point, recover coalescence, and
   either pay at least 247 doublings for the orbit action or drop S17 and become
   Pollard rho.

Correctness gates pass.  Coalescence, operation, storage, degree, complete-cost
and non-generic gates fail.  The dominant obstruction is no longer an
unmeasured implementation constant: it is the incompatibility between
fixed-cardinality representation state and memoryless group collision
detection.

## Reproduction

```bash
cargo test --bin p256_selector_coalescence
cargo run --release --bin p256_selector_coalescence -- \
  --round25 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_selector_coalescence_round26_20261006/coalescence-result.json
```

Canonical artifact: `coalescence-result.json`, 9,269 bytes, SHA-256
`463c100dd951e2a7a7a419118f84bde90fc419358d8279a430a5a0e6191db029`.
