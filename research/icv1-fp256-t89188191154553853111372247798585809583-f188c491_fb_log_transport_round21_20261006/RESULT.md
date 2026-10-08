# P-256 FB1 logarithm-transport census, round 21: result

Date run: 2026-10-06

**THE EXACT TRANSPORT CENSUS FINDS ZERO FB1 RELATIONS; PARITY REMAINS BLOCKED
BY ALL 131,458 LOGARITHM DIMENSIONS.**  The positive control succeeds: 257,928
small-multiple equalities on a 4,096-point actual-P-256 doubling orbit replay
exactly and reduce it to one weighted component.  The candidate factor base is
different.  All 8,413,312 images `[a]P_i`, `1 <= a <= 64`, have distinct
abscissae, and all 5,986,344 completed signed Dickson-block images at depths
1--4 have distinct abscissae across blocks and miss the factor base.  The
combined exact relation rank is zero.

This is a reproducible negative result for the two registered transport
families.  It does not prove that no other cross-column invariant exists.  It
does prove that neither small rational multiples through 64 nor complete
depth-1--4 Dickson-block signed sums supply even one of the 131,457 independent
relations required by this route's rho-parity gate.

## One boundary table, one unit

The unit is one P-256 group-addition equivalent and
`S = total additions / sqrt(n)`.  Rows labelled as boundaries omit
implementation overhead and therefore favour the candidate.  They are not
end-to-end attack measurements.

| variant | class | log2 group-addition equivalents | S | candidate / rho | correctness / scope |
|:--|:--|--:|--:|--:|:--|
| Pollard rho | reference | 128.379 | 1.300 | 1.000 | registered reference |
| round-20 direct finite-domain S17 | projected engineering | 147.798 | 911,265.310 | 700,973.316 | exact sampled cells; projected full collection |
| round-20 optimistic finite-domain floor | generic lower boundary | 144.581 | 98,069.141 | 75,437.801 | unimplemented analytic floor |
| round-21 quotient-adjusted two-list floor | generic lower boundary | 137.503 | 725.341 | **557.955** | exact `K=131458`; unimplemented analytic floor |
| round-21 ideal memoryless quotient floor | unimplemented lower boundary | 136.828 | 454.540 | **349.646** | ideal birthday constant; collector absent |

The quotient formulas recover 7.078 bits over round 20's finite-domain image
floor only by granting an unconstrained idealised collector; they do not move
the exact quotient rank or produce an implementation.  Even that stronger
boundary is 8.450 bits above rho.  The balanced two-list boundary is 9.124
bits above rho.  This round is therefore `negative/no-rank-compression`, not
an advance.

The inherited complete projection still governs promotion: direct collection
is `2^151.885` FME, or 31.885 bits above the `2^120` gate; `2^134.811` FME per
usable row, or 31.811 bits above the `2^103` gate; and `2^135.056` materialised
bytes, or 85.056 bits above the `2^50` gate.  The bounded census itself peaks
at 378,599,040 logical bytes and writes zero disk bytes, but that does not
erase the projected full collector's storage.

## Candidate census

The multiplier census sorts complete affine `x` values and retains `y` parity,
column and multiplier.  An equal-`x` pair would prove either
`a*log(P_i)-b*log(P_j)=0` or `a*log(P_i)+b*log(P_j)=0` modulo the registered
P-256 subgroup order.  No truncated key or probabilistic filter is used.

| family | exact images | equal-x candidates | exact relations | independent rank | peak logical bytes |
|:--|--:|--:|--:|--:|--:|
| `[1..64] * FB1` | 8,413,312 | 0 | 0 | 0 | 378,599,040 |
| Dickson depth 1 | 65,724 | 0 | 0 | 0 | 3,154,792 |
| Dickson depth 2 | 147,936 | 0 | 0 | 0 | 6,443,272 |
| Dickson depth 3 | 420,904 | 0 | 0 | 0 | 17,361,992 |
| Dickson depth 4 | 5,351,780 | 0 | 0 | 0 | 214,597,032 |

All registered block depths complete below the frozen 12,000,000-image and
1-GiB caps.  Their maximum selected leaf counts are 2, 4, 8 and 15.  There are
no infinity images, cross-block equalities, factor-base hits, duplicate
relations or replay failures.

Candidate construction charges 14,458,682 projective additions,
72,696,106 normalisation field multiplications and 1,314,580 Dickson field
squarings.  At 17 field multiplications per projective addition, this is
319,808,280 FME (`2^28.253`) for the bounded census.  It is diagnostic work,
not a P-256 discrete-log cost or an amortisable attack precomputation.

## Exact controls

The orbit-positive control uses actual P-256 points `C_i=[2^i]G`.
It materialises 262,144 small-multiple images, finds 257,928 cross-column
equalities, replays every one, detects all 4,095 adjacent doubling edges, and
returns exactly one weighted component with independent rank 4,095.  There
are no inconsistent cycles or replay failures.

The frozen hash-selected 12-column multiplier reference compares the sorted
detector with all 294,528 direct image pairs and returns the same empty
relation set.  That global subset has no multi-leaf blocks through depth 4,
so an additional conservative validation compares Gray enumeration against
direct signed addition on one deterministically selected maximum-width actual
FB1 block per depth.  The four comparisons cover 2, 8, 128 and 16,384 images,
respectively, with zero false positives and zero false negatives.  This added
control changes no candidate count or rank.

The complete result has zero false positives, zero false negatives and zero
failed group replays.  A byte-for-byte second release run reproduced the
canonical JSON.

## Quotient and degree decision

The multiplier graph has 131,458 components, all isolated and of size one.
The block relation matrix is empty.  Hence

```text
K = 131458 - rank(pairwise relations) - rank(block relations in quotient)
  = 131458 - 0 - 0
  = 131458.
```

The frozen structured residual solving-degree maxima remain 3, 3 and 4, and
the local split join degree remains 2.  The degree-at-most-5 gate passes.
Unsplit S17 degree of regularity remains unknown: transport rank neither
establishes nor lowers it.

The quotient-dimension, rho-parity, per-row, total-cost and projected-storage
gates fail.  No full-depth unplanted P-256 relation was attempted.  The next
credible parity experiment must change the factor-base design or prove a
higher-order cross-column action; extending the same generic collision
schedule or quoting local degree cannot close an 8.450-bit ideal-boundary gap.

## Reproduction and custody

Source commit used for the final run:
`8a03b8200f64e39f4c397de5a817ec579af3c829`.

```bash
cargo test --bin p256_fb_log_transport
cargo clippy --bin p256_fb_log_transport -- -D warnings
cargo build --release --bin p256_fb_log_transport
target/release/p256_fb_log_transport \
  --round20 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_batch_collision_round20_20261005/collision-result.json \
  --round6 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_residual_scaling_round6_20261004/degree-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_fb_log_transport_round21_20261006/transport-result.json
```

The final run used Rust 1.98.0 on Linux x86-64.  Wall time and process RSS are
unranked machine telemetry and are not used in the projection; the deterministic
artifact records exact logical storage and zero disk traffic.  The result JSON
is 13,057 bytes with SHA-256
`4aaf5fde72257ed1f262d5674ab2fe81c7b4f90698fa3767615adfed6524e757`.
