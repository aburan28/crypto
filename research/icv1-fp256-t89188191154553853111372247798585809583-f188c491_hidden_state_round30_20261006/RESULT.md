# P-256 non-memoryless hidden-state selector, round 30: result

Date run: 2026-10-06

Retaining representation state does not recover rho parity for the one-addition
scalar-orbit selector.  The 279,184 distinct oriented update deltas collapse to
exactly 139,592 sign-folded compatibility classes.  An exact lifted collision
therefore pays a `sqrt(139592)=373.620` hidden-state overhead.  Even after
crediting the optimistic round-26 oracle constant, the lower boundary is
**360.296 times rho**, or about `2^136.872` operations.

Compressing the hidden key cannot improve that boundary.  It produces more
same-bucket collision candidates, but the exact-replay survival probability
falls by precisely the factor which cancels the apparent collision gain.

No full-depth unplanted relation was attempted.

## Complete delta reconstruction

The binary reconstructed every local update coefficient

```text
D_j = (m-1)m^j mod n
```

from the frozen round-25 scalar and P-256 order.  The result exactly reproduces
round 26:

| check | result |
|:--|--:|
| signed atoms | 279,184 |
| distinct oriented deltas | 279,184 |
| duplicate deltas | 0 |
| modular multiplications | 558,368 |
| full scalar cycle closes | yes |
| half cycle is negation | yes |
| ordered delta SHA-256 | `9bef47854a134951a3970dc45b31af2153ffea5b37dca9562577310fde01f9a4` |

For exact compatibility modulo global sign, the hidden key is

```text
K_j = min(D_j,n-D_j).
```

This quotient has 139,592 classes, each containing exactly two opposite
deltas.  All 139,592 compatible pairs replay in scalar group arithmetic; none
are equal-oriented, all are opposite-oriented, and false positives, false
negatives and replay failures are zero.  The exact class inventory has digest
`ce60be2de1c1af0baaf0f79f774068298c8ff2214f00d695d8fb5869118c1191`.

The class entropy is 17.090857 bits.  A perfect rank needs an 18-bit fixed key
or 314,082 packed bytes for the complete base; storing every canonical
coefficient uses 4,466,944 bytes.  These small metadata tables are not the
failed storage gate—the materialized collision frontier remains
`2^133.326` bytes before the hidden-state penalty.

## Why compression cannot help

Let bucket `i` contain `w_i` exact delta classes and set

```text
S2 = sum_i w_i^2.
```

A bucket collision survives exact delta replay with probability
`B/S2`, where `B=139592`.  Bucket collisions arrive faster by the effective
state-space factor `B/sqrt(S2)`, but `S2/B` such collisions are required per
usable one.  Their product is `sqrt(S2)`.  Since `S2 >= B`, equality is reached
only when every exact class remains distinct.  Even a perfectly balanced
oracle partition therefore cannot beat the full exact quotient.

The deterministic SHA-256 prefix ladder confirms the accounting on the fixed
complete class set:

| prefix bits | occupied buckets | max width | false merged pairs | replay survival | measured lower bound / rho |
|--:|--:|--:|--:|--:|--:|
| 0 | 1 | 139,592 | 9,742,893,436 | 0.000716% | 134,613.658 |
| 8 | 256 | 607 | 38,066,211 | 0.183019% | 8,421.922 |
| 16 | 57,737 | 11 | 148,214 | 32.015045% | 636.769 |
| 17 | 85,907 | 8 | 74,231 | 48.460358% | 517.566 |
| 18 | 108,312 | 7 | 37,062 | 65.316588% | 445.807 |
| 24 | 139,007 | 3 | 588 | 99.164583% | 361.810 |
| 32 | 139,588 | 2 | 4 | 99.994269% | 360.306 |
| 34 | 139,592 | 1 | 0 | 100% | **360.296** |
| 48 | 139,592 | 1 | 0 | 100% | **360.296** |

The first collision-free hash prefix is 34 bits.  A perfect rank reaches the
same boundary with 18 fixed bits; the wider hash label is only a storage
encoding artifact.  No hash branch is discarded, and every same-bucket
candidate is charged through exact replay.

## End-to-end effect

Round 26's one-addition row alone had an oracle ratio of 0.964336 and failed
coalescence.  Lifting it with the exact hidden state changes the optimistic
time boundary to

```text
0.964336 * sqrt(139592) = 360.295518 rho
```

and adds 8.545428 bits to the oracle collision work:

```text
2^(128.326143 + 8.545428) = 2^136.871572 operations.
```

The exact one-class log transport and fixed signed S17 representation survive,
but the usable-collision time gate fails.  The imported structured degree of
regularity remains unknown, relation collection is not below `2^120`, cost per
usable relation is not below `2^103`, projected collector storage is not below
`2^50`, and the end-to-end construction is not non-generic.  Thus every
promotion gate except correctness, completeness, representation preservation
and log transport fails.

## Decision

Reject compressed non-memoryless state for `FB1he6b6b6e25de6`.  Exact state
resolves round 26's coalescence failure, but its entropy recreates essentially
the same hundreds-of-times-rho quotient penalty.  Shorter keys do not remove
that entropy; they move it into false candidate collisions and replay retries.

The remaining route must change the delta equivalence itself—not merely its
encoding.  A useful next candidate would need a factor base whose local moves
share only a constant number of exact group deltas while still spanning enough
distinct columns for S17 relations, or must abandon local fixed-S17 moves.

## Reproduction

```bash
cargo test --bin p256_hidden_state_selector
cargo run --release --bin p256_hidden_state_selector -- \
  --round25 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json \
  --round26 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_selector_coalescence_round26_20261006/coalescence-result.json \
  --round29 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_memoryless_closure_round29_20261006/memoryless-closure-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_hidden_state_round30_20261006/hidden-state-result.json
```

Canonical result artifact: `hidden-state-result.json`, 48,674 bytes, SHA-256
`8927791b3b149c60a0d4cfb43b77aca21ede9ff21c9c9451a534247ec1f56ad1`.
