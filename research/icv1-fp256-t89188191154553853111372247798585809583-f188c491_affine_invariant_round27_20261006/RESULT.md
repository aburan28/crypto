# P-256 affine-invariant factor-base screen, round 27: result

Date run: 2026-10-06

No proper affine-invariant P-256 factor base reaches rho parity.  The complete
screen reduces additive progressions and affine recurrences to a sharp
dichotomy: a nonzero translation forces the full prime-order group, while a
proper signed base is necessarily a union of scalar orbits.  Exhausting every
such action that fits between 131,458 and 139,592 columns leaves a rigorous
best-case floor of 237 group operations per target-corrected sample, 228.548
times the parity budget of 1.036982 operations.

No full-depth unplanted relation was attempted.

## Affine classification

Let the signed coefficient set be `S=-S` in `F_n`, with zero excluded, and
let `f(x)=a*x+b` permute `S`.  Conjugating negation `nu(x)=-x` by `f` gives

```text
f o nu o f^-1: x -> -x+2b.
```

Composing that reflection with negation gives `x -> x+2b`.  If `b` is
nonzero, this translation has orbit length `n` because the P-256 group order
is prime.  A nonempty invariant set is therefore the full group, not a
131,458--139,592-column factor base.  Hence every proper signed base with a
global affine action has `b=0` and is a union of scalar orbits on
`F_n^*/{+-1}`.

The explicit arithmetic progression `{[1+j]G : 0 <= j < 131458}` confirms
the boundary.  All 131,457 internal `+G` edges replay, but the last edge wraps
with nonzero defect `[131458]G`; it reaches neither the first point nor its
negative.  Correcting a wrap therefore depends on representation state.  The
proper prefix is not a cyclic memoryless walk, and exposing its coefficients
would in any case make it a generic representation encoding.

## Complete action census

The binary independently replays the Round 25 factorisation of `n-1`, its
10,240 divisors, and primitive root 7.  It retains all 79 nontrivial orders
whose negation-folded width is at most 139,592, then visits every exact-order
element of every retained subgroup.

| census item | exact count |
|:--|--:|
| exact-order elements visited | 1,712,366 |
| unique actions after folding `a` and `-a` | 856,183 |
| folded duplicates removed | 856,183 |
| width-band `(action width, orbit count)` combinations | 12,712 |
| nondominated frontier rows | 24 |

The ordered census SHA-256 is
`d302ee4eebc9ae48e54109b552bdf54cacf43333cc47e79c050f5156ad436e34`.
The enumeration is exhaustive, not a small-scalar sample.

The three relevant extrema are:

| objective | width x orbits | columns | scalar action | binary ops | doubling floor | target-corrected floor |
|:--|--:|--:|:--|--:|--:|--:|
| minimum binary cost | 26,483 x 5 | 132,415 | `154410347155025119293275570919859591050360234274381103831027912432220504873` | **340** | 246 | 247 |
| minimum rigorous floor | 9,301 x 15 | 139,515 | `152918900535311551617472835739352639811851243692460361273286367213104126` | 349 | **236** | **237** |
| minimum unknown anchors | 139,592 x 1 | 139,592 | `379483656589291393959605088584129393463128953089898693989044430844151895048` | 342 | 247 | 248 |

Round 27 uses the preregistered standard left-to-right binary count
`bits + Hamming weight - 2`; the target correction is one further group
addition on the right colour.  Even granting the minimum doubling-only floor
over every candidate, `237 / 1.036982 = 228.548` times rho parity.  The
minimum binary-cost action is 328.839 times the parity budget after its target
correction.

## Native factor base and exact replay

The minimum-binary-cost candidate was fully materialised from five
independent SHA-256-labelled anchors in the repository's native FB1 identity
and wide-dump format:

```text
factor base       FB1hd37fd8aa9f44
columns           132,415
signed points     264,830
orbit width       26,483
unknown anchors   5
points SHA-256    2e0f455497e49ba430fff5c6279cc30c834f3b09ba2a7e3bb76f603cfc211849
dump bytes        66,614,641
dump SHA-256      d81022821b00140b97cf0137b88fe2e8ae36a319a8affa9b446c8374ea38094e
```

All five first-choice anchors were disjoint.  Every point is nonidentity and
on curve, all 132,415 columns are unique up to sign, all 132,415 action edges
replay, and all 4,096 independent hash-selected controls replay.  False
positives, false negatives, transport failures, and control failures are all
zero.  The native wide dump round-trips to the same FB1 identity.

The charged construction and replay counts are 25,158,850 scalar additions,
65,148,180 scalar doublings, 389,120 control additions, 1,007,616 control
doublings, 262,915 progression additions, 662,070 normalisation
multiplications, and five batch inversions.

## End-to-end consequence

Unlike Round 26's cheap local update, a global scalar action does coalesce and
permutes the complete signed base.  Its arithmetic cost is fatal before any
summation-polynomial solve: the best rigorous floor over the complete width
band is 228.548 times the rho per-sample allowance.

The native minimum-cost base also retains five independent unknown anchor
logarithms.  Publishing them makes every factor-base log known and turns the
remaining target representation into a generic claw/rho search; withholding
them requires relation information to solve those global classes.  Neither
case is a non-generic index-calculus advance.

Coordinate membership remains the degree-132,415 support polynomial.  No
structured residual degree of regularity at most five was measured.  Thus the
operation, degree, relation-collection, per-row, storage, complete-cost, and
non-generic gates fail; only the exact correctness gates pass.

## Decision

Reject affine-invariant factor bases for rho parity on this P-256 curve.
Additive progressions cannot wrap at a proper width, nonzero affine
translations force the full group, and the exhaustive scalar-union remainder
misses the operation gate by more than two orders of magnitude even under its
best doubling-only floor.  Local degree claims cannot change that complete
action obstruction.

## Reproduction

```bash
cargo test --bin p256_affine_invariant_factor_base
cargo run --release --bin p256_affine_invariant_factor_base -- \
  --round25 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json \
  --round26 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_selector_coalescence_round26_20261006/coalescence-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_invariant_round27_20261006/affine-invariant-result.json \
  --dump target/p256-affine-invariant-round27.factor-base.json
```

Canonical result artifact: `affine-invariant-result.json`, 47,507 bytes,
SHA-256
`256271ebb6f1c3c736a01da4358f01a6822c4315bb2081adcba1dd05b1d847a5`.
