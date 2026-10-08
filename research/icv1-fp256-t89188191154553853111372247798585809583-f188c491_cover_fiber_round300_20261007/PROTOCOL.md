# P-256 cover-fiber quotient screen, round 300: preregistration

Date preregistered: 2026-10-07

## Objective

Test the remaining suggestion that a lift to a cover or Jacobian can produce
many independent relation rows per decomposition event through its multiple
fibres.  The experiment distinguishes raw lifted rows from equations that
remain distinct after pushforward to the P-256 prime-order subgroup.

The P-256 curve is
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`.  The primary
width-164 comparison is `FB1hc72514a2a8d3`; the registered 17-term comparison
`FB1h2f8621cda105` remains in the boundary table.  Neither factor-base artifact
will be rewritten.

## Frozen dependencies

- Round 299 result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_transfer_anchor_round299_20261007/transfer-anchor-result.json`,
  SHA-256 `e6b78d729e785b1831a061cd39e6e57f10c46da8b8e98418321dc72f56752f35`.
- Round 295 selected factor base:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_union_round295_20261007/selected-factor-base.json`,
  SHA-256 `d27516ca40a612ecf3ebabfa8ae04776084c7948110e20f221da438e1e64d8f8`.
- Round 299 toy curve `y^2=x^3-3x+241` over `F_1151`, with
  `#E(F_1151)=1192=8*149`, order-149 generator `(360,129)`, and complete
  cofactor-cleared subgroup.

The runner must reject a dependency hash mismatch.

## Typed correspondence and theorem under test

For a lifted coefficient space `V=F_q^L`, a destination coefficient space
`W=F_q^B`, and the sheet-collapse map

```text
C : V -> W,
C(e_(j,s)) = e_j,
K = ker(C),
```

the exact quotient identity is

```text
dim((span(R)+K)/K) = rank(C(R)).
```

Rows that differ only by choices of sheets over the same destination columns
differ by an element of `K`.  Therefore `r^m` lifted versions of one
`m`-term destination row count as one destination equation, regardless of
their raw lifted-column rank.

The executable bounded correspondence uses multiplication maps
`[k]:E->E`, for `k in {1,2,4,8}`.  These are curve endomorphisms of degree
`k^2`.  The runner enumerates every rational point, every rational kernel
point, and every rational preimage of the selected order-149 subgroup points;
it does not infer fibre sizes from degree.  On the order-149 subgroup the map
is nonzero and injective.  Extra rational sheets, when present in the
cofactor-eight ambient group, differ by rational kernel points.

For P-256, the relevant rational group has prime order `n`.  Any `[k]` with
`gcd(k,n)=1` is a bijection on that group, even though its algebraic degree is
`k^2`.  A general cover/Jacobian claim remains construction-specific, but
rows generated solely by a pushforward fibre kernel are charged according to
their quotient rank, not their raw lift count.

## Frozen factor-base controls

The complete toy experiment uses the following six-column bases imported
from the deterministic Round 299 screen:

| family | folded subgroup labels | uses scalar labels in construction |
|:--|:--|:--:|
| affine-progression | `1,3,5,7,9,11` | yes |
| pair-sum-closure | `3,32,35,38,41,70` | yes |
| low-coordinate-geometry | `8,41,46,49,54,59` | no |
| hash-coordinate-geometry | `4,13,14,22,40,41` | no |

Every base uses arity three.  The destination domain contains all
`2^3*C(6,3)=160` raw signed rows and 80 global-negation-normalized rows.
For a rational fibre size `r`, all `r^3` sheet assignments for every raw row
are enumerated.  No hash sampling or probabilistic branch discard is allowed.

## Measurements

For every base and map `[k]`, record:

- algebraic degree, complete rational kernel and fibre-size histograms;
- raw destination rows, raw lifted rows, unique lifted rows, and source sums;
- target occupancy before and after lift multiplicity;
- rank over `F_149` in lifted columns, the exact deck-kernel rank, union rank,
  quotient rank, collapsed destination rank, and homogeneous versions;
- exact group replay of every lifted row and its `[k]` pushforward;
- false positives, false negatives, missing fibres, duplicate accounting, and
  quotient-identity failures;
- counted source additions, pushforward doublings/additions, modular-rank
  operations, logical bytes, disk bytes, wall time, and peak RSS.

The implementation must also emit a typed transfer assessment with the
source, destination, field, map, subgroup restriction, kernel, exceptional
set, recovery obligation, controls, measured costs, and open construction
obligations.  Measured toy quantities, exact P-256 group facts, analytic lower
bounds, and unimplemented projections remain separate.

## Boundary and promotion gates

The common unit is favourable P-256 group-addition equivalents divided by
`sqrt(n)`, with Pollard rho fixed at ratio 1.  Round 295's free-oracle
independent-log boundary is `13.920747397073491` times rho.  A cover candidate
may improve that number only through increased **quotient** relation rank or a
complete measured reduction in decomposition cost.  Algebraic degree,
rational sheet count, raw row count, and raw lifted-column rank alone receive
no credit.

Promotion requires all of:

1. zero false positives, false negatives, missing fibres, replay failures,
   and quotient-identity failures on every complete cell;
2. at least two independent destination equations from one destination row
   solely because of the cover construction;
3. an explicit P-256 cover/Jacobian construction with fields, dimensions,
   formulas, kernel, exceptional points, inverse recovery, and relation
   pullback;
4. structured residual degree of regularity at most five;
5. projected complete relation collection below `2^120` field/group-equivalent
   operations and below the matched rho cost;
6. projected cost per usable relation below `2^103`;
7. peak projected materialized storage below `2^50` bytes; and
8. no discarded probabilistic branch counted as exhaustive work.

An actual unplanted full-depth P-256 relation is attempted only if all gates
pass.  Otherwise the round publishes a reproducible negative result naming
the exact quotient obstruction.  A successful toy lift is not P-256 evidence
unless the typed construction and scaling obligations are also met.

## Success and stop conditions

The cover-fibre hypothesis succeeds only if complete replay exhibits quotient
rank greater than the collapsed destination rank, or an explicit construction
supplies independent destination relations not generated by the fibre kernel.
The null is accepted within this registered family if every apparent lift
gain lies in `ker(C)` and the quotient identity holds for all buckets.

Stop after all 16 registered base/map cells and their complete row domains
finish, or on the first correctness failure.  Preserve every failure and all
partial counters.  Do not attempt a P-256 relation when the promotion gates
fail.

## Evidence limitation

The transfer assessment profile is available.  Its referenced
`references/methodology.md` and `assets/assessment-template.json` resources
are unavailable in this environment as of preregistration; this limitation
must remain explicit in the assessment artifact.
