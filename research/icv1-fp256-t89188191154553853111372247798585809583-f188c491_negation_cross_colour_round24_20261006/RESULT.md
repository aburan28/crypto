# P-256 negation-folded cross-colour boundary, round 24: result

Date run: 2026-10-06

The negation correction is valid and reaches **ideal collision-constant parity
only for `K=1`**.  It does not reach end-to-end parity.

Both orientations of an equal-abscissa cross-colour hit replay as valid signed
S17 relations:

```text
L =  R = T-B  => T =  L+B
L = -R = B-T  => T = -L+B.
```

The second branch flips every one of the eight left signs without changing its
columns.  Thus the random collision state is the negation quotient of size
`n/2`.  Round 22's `sqrt(pi*K*n)` formula and the projections inherited by
round 23 are superseded.  Their exact zero-transport measurements remain
valid.

## Corrected boundary

With equal left/right sampling and `p_disjoint(B)=C(B-8,9)/C(B,9)`, the exact
expected total sample count to the `K`-th cross-colour collision is

```text
E[T_K] = sqrt(2*n/p_disjoint(B)) * Gamma(K+1/2)/Gamma(K).
```

For `K=1`, this is `sqrt(pi*n/(2*p_disjoint))`; for large `K`, it approaches
`sqrt(2*K*n/p_disjoint)`.  Applying the first-hit constant to `sqrt(K)` was a
second approximation error in rounds 22/23.

| case | legacy S | corrected S | / rho | log2 oracle work | log2 / usable collision | log2 direct work | log2 one-list bytes |
|:--|--:|--:|--:|--:|--:|--:|--:|
| known-log `K=1` target | 1.773 | **1.254** | **0.964** | 128.326 | 128.326 | 131.326 | 133.326 |
| round-23 best width, `FB1h8fc8b5fd8529` | 638.942 | 509.801 | 392.155 | 136.994 | 120.007 | 139.994 | 141.994 |
| selected Dickson, `FB1h2f8621cda105` | 642.817 | 512.893 | 394.533 | 137.003 | 119.998 | 140.003 | 142.003 |
| fixed 138,031 useful rows | 658.692 | 525.559 | 404.277 | 137.038 | 119.963 | 140.038 | 142.038 |

The repository rho reference is `S=1.3`.  The `K=1` random-oracle row is
3.56% below that reference.  This is a lower-bound constant, not an attack:

- direct eight-/nine-term construction averages eight group additions per
  sample and is 7.715 times rho;
- materializing one colour and streaming the other still needs about
  `2^133.326` bytes for `K=1`;
- the exact one-addition sign-Gray controls cover only one fixed support and
  receive no global-column projection credit; and
- a memoryless walk with known scalar labels is a generic rho/claw walk under
  another name, not index calculus.

The actual geometric bases retain `K=B`, so even their optimistic random-oracle
boundaries remain 392--395 times rho.  Fixed-row collection remains
`2^137.038`, its per-row cost remains 16.96 bits above the `2^103` gate, and
one-list storage remains about `2^142` bytes.  Sparse linear algebra is still
negligible beside collection.

## Exact validation

The complete quotient truth table covered 21,365 pairs in cyclic groups of
orders 9, 15, 31, 63, and 127.  Canonical keys matched exactly the union of
`L=R` and `L=-R`: zero false positives and zero false negatives.

Deterministic alternating-colour simulations used 7,680 trials:

| cyclic order | trials | mean S | predicted S | positive / negative hits | FP / FN |
|--:|--:|--:|--:|--:|:--|
| 4,093 | 4,096 | 1.2672 | 1.2533 | 2,093 / 2,003 | 0 / 0 |
| 65,537 | 2,048 | 1.2542 | 1.2533 | 1,031 / 1,017 | 0 / 0 |
| 1,048,583 | 1,024 | 1.2291 | 1.2533 | 495 / 529 | 0 / 0 |
| 16,777,213 | 512 | 1.2311 | 1.2533 | 256 / 256 | 0 / 0 |

The high-order deviations are consistent with the deliberately smaller trial
counts; the exact law, not a fit to these means, supplies the P-256 projection.

Two planted P-256 controls used the same 17 hash-selected, signed, distinct
known-scalar columns.  One forced `L=R`, the other `L=-R`.  Both had equal
abscissae, replayed the appropriate 17-term relation, and independently
recovered and verified the target scalar.  The negative control also verified
that the two group points were exact opposites.

The fixed-support sign-Gray controls exhaustively visited 128 folded eight-sign
states and 512 nine-sign states with one group update per transition and zero
replay failures.  They do not vary all 17 columns and therefore receive no
global selector or cost credit.

## Decision

The ideal `K=1` collision constant passes rho parity, but every complete gate
fails.  Known-log controls are rho-equivalent, not IC; direct generation and
storage fail; geometric quotient dimensions remain full; structured residual
and unsplit S17 regularity remain unknown; and collection remains above the
cost gates.  No full-depth unplanted P-256 relation was attempted.

The next admissible parity target is no longer another coordinate-only factor
base.  It must provide either:

1. exact elliptic-log transport that reduces a geometric base to `K=1` without
   revealing construction logs, plus a global one-addition memoryless selector
   that is demonstrably not generic rho; or
2. a non-generic target invariant that avoids the cross-colour collision
   process entirely.

No screened P-256 base currently has either property.

## Reproduction

```bash
cargo test --bin p256_negation_cross_colour
cargo run --release --bin p256_negation_cross_colour -- \
  --round22 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_factor_base_symmetry_round22_20261006/symmetry-result.json \
  --round23 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_root_orbit_fb_round23_20261006/root-orbit-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_negation_cross_colour_round24_20261006/negation-result.json
```

Canonical artifact: `negation-result.json`, 9,440 bytes, SHA-256
`d1e2fe33e1abbae8a9bf6885ca7cfef8012233363f4cd43e0d30d468a381fc9c`.
