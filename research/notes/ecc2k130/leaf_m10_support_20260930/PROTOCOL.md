# Frozen exact physical-support census on the first degree-263 leaves

Status: preregistered before any leaf slot count from this panel. The source
normal-basis m10 screens count the first slot and carry its physical and
projected-column counts to the rotated slots using source-curve Frobenius.
Coordinate squaring takes a fixed descending leaf to its conjugate curve, so
that shortcut does not justify a fixed-leaf count. This experiment asks what
the **actual per-slot counts** are on two previously certified nonconjugate
degree-263 leaves. It is a structural admission step before an equal-useful-
size PDP comparison, not that comparison itself.

## Fixed inputs and exact census

[`CONFIG.json`](CONFIG.json) freezes the source `E_0: y²+xy=x³+1`, the
normalized leaves for kernel lines `[1,0]` and `[1,4]` from the all-line
certificate, `GF(2^131)` modulus, normal element `β=3`, subgroup order, input
hashes and caps. Use the unchanged rotated m10 x-coordinate domain
`x_i=Σ_{j<d_i} a_ij β^(2^(10j+i))`, with all binary masks, no infinity
selector and both rational sign lifts. The balanced dimensions are ten 13s;
the unequal dimensions are `[14,13,13,13,13,13,13,13,13,13]`. Scan ten
13-dimensional slots and the 14-dimensional slot zero **on each curve**;
reuse a scan only where the two frozen arms literally share a slot/domain.
Do not select a leaf, basis, slot placement, mask or dimension from outcomes.

For every x in each slot, x=0 has the unique point `(0,√b)`. For nonzero x,
the two rational y lifts exist exactly when
`Tr(x + a₂ + b/x²)=0`; otherwise neither exists. Record one canonical row
`mask,x,lift_count,projected_x_or_sentinel\n` as decimal ASCII in binary
Gray-mask traversal order, with SHA-256 for the entire row stream and each
1,024-row chunk, exact lift counts and physical point counts. Sentinel `-2`
means nonliftable and `-1` means `[4]P=O`, including x=0. On a liftable
nonzero x, compute the signed x-class of
`[4]P` through `u=x²+b/x²`, `x([4]P)=u²+b/u²`, with infinity handled when
u=0. Record each slot's distinct nonzero projected signed classes, collision
multiplicities, and the union across all ten slots for each arm. These are
**uncompressed** projected classes on a leaf: do not call leaf coordinate
squaring a cheap same-curve Frobenius action. Count the full physical tuple
product across slots and report `product/q` only as a necessary counting
ceiling, not a hit-rate estimate; the factors can occupy different cofactor
cosets and may collide after projection.

The source is a positive control: every 13-dimensional slot must reproduce
7,977 physical points and the unequal high slot 16,125, matching the frozen
source m10 screens. The leaf counts are unknown at preregistration. A
`PASS_CENSUS` requires 33 complete scans, all row/chunk hashes and summaries
reproduced by a fresh independent arithmetic implementation, and exact group
checks on the first and last two liftable nonzero masks in each unique scan.
The verifier must check the actual leaf coefficient against the locked
all-line certificate, solve sampled y coordinates, check the curve equation,
and compare the closed-form `[4]` x against two full point doublings.
Any mismatch is `FAIL`; a wall/RSS cap or missing scan is `CENSORED`. Preserve
completed rows and failures; do not omit a hard slot or rerun a selected
subset. The producer and verifier have separate 600/1,200-second wall and
1-GiB RSS caps.

## Decision boundary and next gate

Report whether the **single-slot count transfer** fails on either fixed
leaf, and quantify the resulting physical tuple and uncompressed column
differences. A difference is a representation effect, not an easier PDP.
Even equal counts would not restore a free leaf Frobenius action. This census
does not sample natural Q, run a solver, check relation rank or recover a log.
It cannot promote an ECC2K-130 speedup or fill the required four-policy
equal-useful-size/heldout comparison. If a viable fixed leaf domain remains,
the next experiment must match useful projected physical points across
original, transported, leaf-native and exact pullback bases, then measure
natural-target PDP outcomes and complete cost. The separate m10 capacity
[PR #937](https://github.com/aburan28/crypto/pull/937) still requires its
independent exact-head review and one-shot label before its representation
measurement; this census does not satisfy or bypass that gate.

Commit the producer, verifier and their exact source/input SHA lock before
running this panel. Archive all 33 scan receipts, source hashes, resource
measurements, independent replay and a decision in a follow-on outcome PR.
