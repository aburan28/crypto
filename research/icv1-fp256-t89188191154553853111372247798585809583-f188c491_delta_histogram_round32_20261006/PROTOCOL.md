# P-256 folded-delta histogram bound, round 32: protocol

Date registered: 2026-10-06

## Requested scope and evidence plan

| continuing requirement | round-32 evidence |
|:--|:--|
| keep screening factor bases toward rho parity | classify every folded adjacent-delta histogram of an ordered known-log 131,458-column base |
| retain all 17 variable columns | bound the complete signed distinct-column S17 support; do not fix atoms or discard branches |
| exact replay and deterministic artifacts | import and cross-check round 31's native P-256 cycles; exhaust and replay complete small-prime cyclic bases; emit deterministic JSON and hashes |
| honest end-to-end comparison | charge hidden compatibility and target coverage against the frozen rho boundary; retain all degree, collection, per-row, storage and non-generic gates |
| no local-only promotion | attempt no unplanted full-depth relation unless every registered promotion gate passes |

## Question

Can three or more exact local-delta classes escape round 31's two-delta
coverage-versus-hidden-state obstruction?

Consider any cyclic ordering of nonzero known-log coefficients
`c_0,...,c_(B-1)` in the prime-order P-256 subgroup, unique up to sign, with
`B=131458`.  Fold every adjacent delta by global sign:

```text
d_i = c_(i+1)-c_i mod n
k_i = min(d_i,n-d_i).
```

Let the folded class counts be `N_1,...,N_t`, let `M=max N_j`, and set
`R=B-M`.  Round 32 sweeps every integer `M=1..B-1`; it does not choose a
delta count or histogram after seeing the result.

## Registered support theorem

Choose a representative `a` of a largest folded class.  Delete the `R` cycle
edges outside `{+a,-a}`.  The remaining graph has at most `R` path runs.  In
each run every coefficient has the form

```text
s_r + z*a mod n,  with integer |z| <= B.
```

Every signed 17-sum is therefore specified by at most 17 signed run labels and
one integer offset in `[-17B,17B]`.  Its support obeys

```text
U(R) = min(n, (2R)^17 * (34B+1)).
```

This overcounts orders, repeated run labels, duplicate columns, sign choices
and coincident offsets.  It is an optimistic upper bound on target coverage,
not an achieved relation count.

`R=0` cannot give a valid distinct cyclic base.  A closed walk whose every
edge is `+a` or `-a` has integer displacement zero because `B<n`; any
nearest-neighbour integer walk returning to its starting value revisits a
vertex before the cyclic endpoint.  The corresponding coefficients are not
distinct.

## Registered histogram optimum

For uniformly visited local edges, exact folded-delta compatibility has
collision probability

```text
C = sum_j (N_j/B)^2.
```

For a fixed maximum multiplicity `M`, convexity maximizes `C` by filling as
many classes of size `M` as possible and putting all remainder in one class.
With

```text
q = floor(B/M),  s = B-qM,
C_max(M) = (q*M^2+s^2)/B^2,
```

every histogram satisfies `C <= C_max(M)`.  Credit the candidate the optimistic
local ratio

```text
round26_local_oracle_ratio / sqrt(C_max(M)).
```

Charge at least `n/U(B-M)` independent target trials when `U<n`.  The product
is a lower bound on complete work relative to rho.  Sweep every `M=1..B-1`
and report the global minimizer, all change points in `floor(B/M)`, the
parity-side coverage deficit, and whether any row reaches rho.  A discarded
histogram, invalid cycle or unproved support value may not be credited as an
algorithm.

For `M>B/2`, the maximizing histogram is exactly the two-class distribution
screened in round 31.  Round 32 must reproduce the overlapping rows but is not
allowed to assume that the global optimum lies in this range.

## Complete references

1. Enumerate every integer partition of `B=2..24`.  For each partition compare
   its exact collision numerator with the registered multiplicity-cap formula;
   record equality cases, violations and a row digest.
2. Exhaust every cyclic ordering and every coefficient orientation of the full
   sign-class bases for prime-order scalar groups `(p,B,m)=(7,3,2)`,
   `(11,5,3)` and `(13,6,3)`.  For every cycle:
   - verify nonidentity and uniqueness up to sign;
   - replay every edge;
   - compute the exact folded-delta histogram;
   - enumerate every signed distinct-column `m`-sum;
   - compare exact support with `U(B-M)` and exact collision with `C_max(M)`.
3. Import round 31's P-256 dependency by hash.  Cross-check its curve,
   relation model, local oracle constant, universal overlapping minimum,
   primitive native landmark and exact edge-replay receipt.

Record cycle counts, edge replays, signed sums, exact support maxima, histogram
widths, collision probabilities, false positives, false negatives, theorem
violations and deterministic digests.  Small-prime checks validate the
implementation and theorem accounting; they are not P-256 exponent fits.

## Gates and stop condition

Promotion requires all of:

- zero partition, support, replay, false-positive and false-negative failures;
- fixed signed S17 and an exact known-log transport class;
- proved P-256 usable-relation probability;
- complete projected time at or below rho;
- structured residual degree of regularity at most 5;
- relation collection below `2^120` operations;
- cost per usable relation below `2^103`;
- peak materialized storage below `2^50` bytes;
- a non-generic end-to-end algorithm.

If the universal complete lower bound exceeds rho, reject every ordered
known-log factor base whose selector requires exact compatibility of one
adjacent folded-delta class.  This conclusion does not cover nonlocal
multi-column transitions, selectors that avoid edge compatibility, unknown-log
bases, or arbitrary ECDLP algorithms.  If any row reaches rho, preserve it for
a separately preregistered native construction rather than claiming success
from the histogram relaxation alone.
