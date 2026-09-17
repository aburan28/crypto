# Specialized x-only doubling with sparse bridge toward 12 B/s

Status: preregistered before implementation, compilation, collision study or GPU timing.

## Hypothesis

The validated x-only [2]/[3] core is zero-spill but tied selected at about 6.29
B/s because generic rational batch evaluation still performs four field products
per output, plus two products on [3] branches. For doubling,

```
x([2]P) = (x^4+1)/x^2 = (x + 1/x)^2.
```

A denominator-only batch tree can recover `1/x` in three amortized products per
output and form the result with XOR plus a linear square, eliminating the
numerator product. Use doubling whenever `weight(x) mod 32 != 14`; on the sparse
branch apply the already proved x-only bridge `P+sigma^3(P)`. Multipliers 2 and
`1+s^3` generate the full scalar multiplicative group, repairing the [2]/[3]
index-two defect.

The bridge should eventually run in a separate rare-path kernel so its five
extra products and live ranges do not spill the common doubling kernel. This is
only worthwhile if the map retains collision efficiency.

## Frozen collision gate

Run 100 matched planted GF(2^23) DLPs for selected and the proposed map with
identical seeds, DP threshold, canonical Frobenius/negation key and restart1024.
Track scalar coefficients directly, charge every overdue iteration, and verify
every endpoint and recovered scalar. Continue to 2,000 only if mean candidate
work is at most 1.10x selected with no unexplained infinity, cycle or mismatch.

A passing larger study permits an isolated CUDA candidate. First implement the
doubling-only denominator tree and prove `(x+1/x)^2` against independent affine
arithmetic. Then add the sparse bridge queue with exact step accounting,
checkpoint versioning, replay, partial-population, restart and sanitizer gates.
Screen three interleaved equal-work pairs on the same g7.2xlarge / RTX PRO 4500
at 165 W. Both raw and collision-adjusted paired intervals must exceed one;
12 B/s requires a direct complete-iteration measurement. Generic work is
`sqrt(n/262)` times measured collision work, and full-DLP S remains null.

## Preregistered common-kernel diagnostic

The denominator-only doubling kernel passed an independent 2,048-point,
seven-step GF(2^131) affine-x oracle with zero mismatches. Before implementing
the queue, run three interleaved equal-work pairs against selected to measure the
common-kernel ceiling. This diagnostic cannot be promoted because doubling alone
does not implement the validated sparse-bridge map. Preserve all six samples.
Continue queue engineering only if the doubling kernel materially exceeds the
selected rate; no raw result substitutes for the final complete map or its
1.006696 collision-work ratio.

## Preregistered larger-batch common-kernel diagnostic

B16 doubling-only measured 11.622244 B/s median and remains non-promotable.
Before the bridge queue, compile the identical denominator-only kernel at B24
and B32 with split2, four shared X slots and minimum three blocks. Validate at
least one complete checkpoint against affine doubling for each. Run three
interleaved B16/B24/B32 repetitions with approximately equal complete work.
This isolates inversion amortization; none of these rows implements the sparse
bridge or can satisfy the final goal alone. Use the fastest exact geometry as
the queue base.

The first B24 compile stopped at the older x-only source's B16-only guards; no
candidate code was compiled or run. Preserve that log. Apply the exact validated
B16..32 guard ranges from the isolated larger-batch source, leaving runtime code
unchanged, and restart B24 then B32 sequentially.

## Preregistered last-prefix cache diagnostic

The completed geometry screen retains B16: medians were 11.512901 B/s (B16),
11.140197 B/s (B24), and 10.934243 B/s (B32), and every geometry matched its
independent doubling checkpoint. Before adding bridge arithmetic, test the
existing `LAST_SLOT_CACHE=2` storage policy in a new isolated B16 source. Extend
the doubling-only path to keep the final product-before-denominator in its
already declared register, avoiding one global prefix store and load per eight
updates. Do not change the map or any other geometry flag. Require an exact
2,048-point, seven-step affine-x checkpoint, then preserve three interleaved
equal-work pairs against the retained B16 binary. This remains a non-promotable
common-kernel diagnostic and cannot establish the complete-map rate.

## Preregistered outlined unified sparse-bridge candidate

The last-prefix cache candidate was exact but regressed from 11.587505 to
11.300331 B/s median (paired geometric ratio 0.972074, 95% CI
[0.958362, 0.985981]) after increasing local storage from 8 to 24 bytes; reject
it and retain cache mode 1.

Before engineering a compacting queue, test one complete-map implementation
that keeps the denominator-only product tree. For ordinary states use
`denominator=x` and finish with `(x+1/x)^2`. For states with
`weight(normal_x) mod 32 == 14`, call a no-inline helper that returns the
validated bridge numerator and denominator, store the rare numerator in the
otherwise unused x-only y field, and carry the branch flag in bit 3 of the
compact denominator tail after the clean denominator has entered the product.
Clear that metadata bit before field arithmetic. This preserves three
amortized products on the common path and isolates the longer bridge live range
from the walk kernel.

Compile B16/split2/cache1/shared-X4/minblocks3. Reject immediately on a compiler
error, any checkpoint mismatch, or changed common occupancy. If exact, measure
three interleaved equal-work pairs against doubling-only to quantify full-map
bridge overhead, and three against selected for promotion evidence. The result
is promotable only with checkpoint/replay version 4, correct bridge scalar
accounting, complete-map measurement, and the frozen 1.006696 collision-work
ratio. The direct raw complete-map median must reach 12 B/s for the current
performance goal.

## Preregistered split common/rare kernel candidate

The exact outlined unified candidate measured 7.027014 B/s median versus
11.534485 B/s for doubling-only (paired ratio 0.605155, 95% CI
[0.595842, 0.614612]). Its no-inline bridge call forces 128 local bytes per
common thread despite a random-state bridge probability of only
0.006629112049; reject it.

Next isolate the rare arithmetic at global-kernel granularity. For each logical
step, launch a one-step common kernel followed on the same CUDA stream by a
rare bridge kernel. The common kernel performs DP/restart classification on the
original x. It uses denominator x and doubles ordinary states. For bridge
states it enters denominator one into the batch product, leaves x unchanged,
and tags bit 3 of that slot's compact denominator tail after all field
arithmetic. The rare kernel scans those coalesced tail bytes, returns immediately
for ordinary states, and computes the validated sigma^3 bridge only for tagged
states, using an independent field inversion before clearing the step. The
engine repeats this ordered pair for the requested step count and advances
`iterBase` once per pair.

This candidate must retain B16/split2/cache1/shared-X4 and at least three common
blocks per SM. Reject on any stale tag, changed DP timing, checkpoint/replay
version error, scalar mismatch, or affine-x mismatch. If exact, run three
interleaved equal-work pairs against doubling-only. Account every launch and
all rare-kernel work as complete-map wall time; preserve the frozen 1.006696
collision-work ratio and require a direct raw median of at least 12 B/s for the
goal.
