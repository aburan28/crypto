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

## Preregistered block-local bridge prepass

The zero-spill split kernels are exact but measure only 3.691592 B/s median
versus 11.490686 B/s for doubling-only (paired ratio 0.319087, 95% CI
[0.308390, 0.330155]). Flushing shared X and scanning every state in a second
global kernel each step dominates; reject the split launch architecture.

Next keep one persistent walk kernel and isolate bridge arithmetic by phase
inside each block. Add a 2,048-entry shared `uint16_t` event queue plus a
256-byte per-thread bridge mask. At each logical step, all threads perform the
existing DP/restart classification and append only bridge events, then synchronize.
Warp 0 consumes queued events in chunks of 32, computes the sigma^3 rational
numerator and denominator, and stores those rare values. After another block
barrier, all threads execute the existing denominator tree: ordinary states use
x and the doubling finish; queued states use the stored denominator and
numerator. This shares the one batch inversion across both maps and keeps the
last four X slots resident across steps.

The combined 32,388-byte shared allocation must retain three blocks per SM, and
the common kernel must compile without local spills. Require version-4 replay,
built-in solver success, and an independent affine checkpoint before timing.
Then run three interleaved equal-work pairs against doubling-only. Preserve all
complete work and the frozen 1.006696 collision-work ratio; the direct raw
median must reach 12 B/s for the goal.

## Preregistered bridge power-chain follow-up

The validated block-local candidate reaches 9.851270 B/s median, materially
better than both complete-map alternatives but still below 12 B/s. Before
changing its selector or queue, remove one algebraically redundant rare-path
product: after computing `x3=x*x2`, form `x6=square(x3)` instead of
`x6=x2*x4`. Squaring is the existing linear polynomial operation and produces
the identical bridge numerator and denominator. Keep every flag, geometry,
queue rule, and map unchanged. Require zero spills, the full built-in suite,
and the independent checkpoint, then compare three interleaved equal-work
samples against the current block-local binary.

## Preregistered compressed-square screen

Replacing `x6=x2*x4` with `square(x3)` is exact but flat: paired ratio
0.999502 with 95% CI [0.985338, 1.013869]. Retain the original block-local
source.

The doubling-specialized map applies a polynomial square to every ordinary
output, so retest the independently proved compressed square reduction from
the earlier G7 square study in this new workload. Port only its
`ECC_PACKED_SQUARE` implementation into isolated copies of the doubling-only
source. Build mask 1 (compressed even/odd reduction with software spreading)
and mask 3 (the same reduction with native CLMAD interleaving); all map,
batch, cache, launch, and compiler flags remain fixed. Require the existing
2,048-point seven-step doubling oracle and zero spills. Screen three
interleaved equal-work repetitions against the current doubling binary. Only a
qualified variant may then be combined with the complete block-local bridge
source and subjected to its version-4 correctness gates.


## Preregistered shared-mask scan follow-up

Both compressed-square variants are exact and retain 80 registers, 8 local
bytes, 28,032 shared bytes, and three blocks per SM. Mask 1 measured a paired
ratio of 0.997607 (95% CI [0.986259, 1.009085]); mask 3 measured 1.003232
([0.978718, 1.028360]). Neither interval qualifies, so retain the original
multiplication-based square.

The block-local complete map currently uses a shared event queue and one atomic
append for each bridge state. Remove that queue and its counter. Keep only the
per-thread 8-bit bridge masks; after classification, warp 0 scans masks in
lane-strided owner order and computes the same rare rational values for set
bits. This preserves the selector, arithmetic, one block-wide inversion, DP
semantics, state evolution, and scalar accounting while removing atomics and
4,096 bytes of shared queue traffic per block.

Require 80 or fewer registers, unchanged three-block occupancy, the full
built-in suite, and a zero-mismatch 2,048-point seven-step version-4 checkpoint.
If exact, compare three interleaved equal-work repetitions against the current
block-local binary. Retain only a paired improvement whose 95% interval exceeds
one; the direct complete-map median must still reach 12 B/s.


## Preregistered complete-map worker-count screen

The queue-free shared-mask scan compiled with 56-byte stores and 100-byte loads
per thread, versus zero spills for the event-queue reference. It fails the
resource gate and is rejected before correctness or timing.

The 9.851270 B/s complete-map result used 524,288 logical workers. Test whether
the working set and block scheduling leave throughput on the table by running
the unchanged, exact block-local binary at 131,072, 262,144, 524,288, and
1,048,576 logical workers. Use respectively 16, 8, 4, and 2 launches of 1,024
steps so every arm executes exactly 34,359,738,368 complete scalar updates.
Preserve three interleaved repetitions with alternating order on the same GPU
and run ID. No code, map, batch geometry, scalar accounting, or collision ratio
changes. Adopt a worker count only if its paired 95% interval against 524,288
exceeds one and a fresh confirmation agrees; the raw median must reach 12 B/s.


## Preregistered owner-local rare phase

The unchanged complete map was flat across the worker-count range. Relative to
524,288 workers, paired ratios and 95% intervals were 0.979083
([0.958473, 1.000137]) at 131,072, 0.998076 ([0.986017, 1.010283]) at
262,144, and 0.999288 ([0.988336, 1.010361]) at 1,048,576. Retain 524,288.

The failed mask scan assigned eight owner threads to each lane of warp 0, which
extended rare-path loop state and caused spills. Test a simpler phased layout:
classification keeps each thread own 8-bit bridge mask, then after the existing
barrier every active thread computes rare rational data only for its own set
bits. This removes the event counter, all queue atomics, the 4,096-byte queue,
and indirect owner addressing. The rare work remains separated from common
product live ranges by block barriers; bridge masks remain shared for the
backward finish.

Keep the selector, formulas, batch inversion, state format, version-4 scalar
accounting, and all compile flags fixed. Reject on any increase over the
zero-spill reference, occupancy loss, built-in failure, or mismatch in the
2,048-point seven-step affine checkpoint. If exact, measure three interleaved
equal-work pairs against the block-local event-queue reference and require a
paired 95% interval above one.


## Preregistered sigma bridge-1 common-path diagnostic

The exact zero-spill owner-local phase regressed to 8.897410 B/s median versus
9.833647 B/s for the event queue, a paired ratio of 0.910059 (95% CI
[0.898419, 0.921849]). Retain the event queue. Nsight Compute could not collect
hardware counters because this instance denies unprivileged GPU counter access;
the failed read-only attempt is preserved.

Before changing the validated complete map, measure the ceiling of the cheap
identity `x(P+sigma(P)) = x + 1/x`. Create an isolated copy of the exact
x-only doubling source and replace only the final ordinary square with the
unsquared sum. This diagnostic applies scalar multiplier `1+s`, whose generated
subgroup with 2 has index 12, so it cannot be promoted or reported as a complete
rho rate. Build with identical B16/split2/cache/shared-X geometry, require no
resource regression, and validate a deterministic 2,048-point seven-step
checkpoint against independent affine `P+sigma(P)` arithmetic. Then run three
interleaved equal-work pairs against doubling-only. Continue to a separately
preregistered complete-map collision study only if the cheap branch provides
enough measured headroom to offset the current sparse sigma-3 bridge overhead.
