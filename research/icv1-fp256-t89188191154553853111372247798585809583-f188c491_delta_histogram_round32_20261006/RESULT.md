# P-256 folded-delta histogram bound, round 32: result

Date run: 2026-10-06

Adding exact local-delta classes does not improve the uniformly sampled
adjacent-edge selector.  A complete sweep over all 131,457 possible maximum
folded-delta multiplicities has no relaxed row at or below rho.  The global
minimum is the same two-class relaxation found in round 31:

```text
M = 124524, R = B-M = 6934
complete lower bound = 1.0164599772037055 rho.
```

This is still an optimistic lower bound.  It grants the support upper bound
full P-256 coverage and omits construction, verification, relation collection,
degree, storage and linear algebra.  It is not a measured attack cost or a
usable relation yield.

## Universal histogram calculation

For a cyclic known-log base, fold adjacent coefficient deltas modulo sign and
let `M` be the largest folded class.  Deleting the other `R=B-M` edges leaves
at most `R` paths whose steps are `+a` or `-a`.  Every coefficient on one path
is its path anchor plus an integer multiple of `a`.  Consequently

```text
|S17| <= U(R) = min(n, (2R)^17 * (34B+1)).
```

For uniformly sampled cyclic edges with class counts `N_j`, exact compatibility
has collision probability `C=sum(N_j/B)^2`.  Among every histogram capped by
`M`, convexity gives

```text
q = floor(B/M), s = B mod M,
C <= C_max(M) = (q*M^2+s^2)/B^2.
```

Round 32 simultaneously grants `C_max`, grants `U(R)`, and charges at least
`n/U(R)` target retries.  Those choices all favor the candidate.

| relaxed histogram landmark | M | R | cap classes | hidden lower bound / rho | support upper bits | complete lower bound / rho |
|:--|--:|--:|--:|--:|--:|--:|
| hidden-only parity edge | 126,680 | 4,778 | 2 | 0.999997 | 246.868955 | 560.682922 |
| **global minimum** | **124,524** | **6,934** | **2** | **1.016460** | **256.000000** | **1.016460** |
| near-balanced two-class cap | 65,730 | 65,728 | 2 | 1.363778 | 256.000000 | 1.363778 |
| near-balanced three-class cap | 43,820 | 87,638 | 3 | 1.670280 | 256.000000 | 1.670280 |
| near-balanced four-class cap | 32,865 | 98,593 | 4 | 1.928673 | 256.000000 | 1.928673 |
| near-balanced five-class cap | 26,292 | 105,166 | 5 | 2.156322 | 256.000000 | 2.156322 |

The sweep records all 724 changes in `floor(B/M)`.  Its ordered-row digest is
`82a9ef98f020d0afbbc37f065321d3405b964d01f380d84692a2b605f28bbabf`.
The best hidden-only parity row again misses the group by 9.131045 support bits.

## Exact references

The convex histogram cap was checked against every 7,336 integer partition for
totals 2 through 24.  There are 299 equality cases and zero violations.  The
partition digest is
`7c91819a39e5f84981f985a508c4a6b26715e924d945c42026688a1f14551117`.

The complete cyclic references enumerate every orientation and every ordering
of the full sign-class bases for `(p,B,m)=(7,3,2)`, `(11,5,3)` and `(13,6,3)`:

| exact reference total | count |
|:--|--:|
| oriented cycles | 49,968 |
| signed distinct-column sums | 7,680,576 |
| replayed edges | 295,824 |
| zero-exception cycles | 0 |
| edge failures | 0 |
| collision-cap violations | 0 |
| support-bound violations | 0 |
| false positives / false negatives | 0 / 0 |

The ordered cycle digest is
`ff86ba9359e84e0f9044268ab7951baeb8389986a31d609055d74a0bd88ab657`.
Round 31's P-256 overlap reproduces exactly.  Its valid `R=6935` primitive
cycle contributes another 131,458 exact coefficient-group edge replays with
zero failures and coefficient digest
`980917981827d813e60484abb0655e8bd527b0beb540ea146d2974ff17303a53`.

## Requirement status

| requirement | status | evidence or gap |
|:--|:--|:--|
| all 17 columns remain variable | verified | support theorem covers the complete signed distinct-column S17 set |
| arbitrary folded-delta count | verified within scope | every `M=1..B-1`, with the optimal histogram at each cap |
| exact deterministic replay | verified | complete partition/cycle references and imported native P-256 edge receipt |
| time at or below rho | failed | global optimistic lower bound is 1.016459977 rho |
| proved P-256 relation probability | not established | only a support upper bound is used |
| structured residual degree at most 5 | not established | imported degree remains unknown |
| collection, per-row and storage gates | not established | no end-to-end candidate survives the time gate |
| non-generic complete algorithm | not established | the screened object is a local-selector relaxation |

No factor base is promoted and no unplanted full-depth relation is attempted.

## Decision and exact scope

Reject ordered known-log factor bases whose selector samples cyclic edges
uniformly and requires equality of one adjacent folded-delta class.  Adding
delta classes cannot beat the two-class relaxed optimum.

This is not an impossibility result for arbitrary factor bases or ECDLP.  In
particular, it does not cover biased edge schedules, nonlocal multi-column
transitions, or selectors that avoid adjacent-edge compatibility.  A next
iteration must analyze one of those mechanisms explicitly rather than add more
classes to the same uniform local walk.

## Reproduction

```bash
cargo test --bin p256_delta_histogram_bound
cargo run --release --bin p256_delta_histogram_bound -- \
  --round31 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_low_delta_round31_20261006/low-delta-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_delta_histogram_round32_20261006/delta-histogram-result.json
```

Canonical result artifact: `delta-histogram-result.json`, 453,895 bytes,
SHA-256
`3672088ee692ac66f309d8a412cff7c1169519c12f2b78b86137d068f07ed736`.
