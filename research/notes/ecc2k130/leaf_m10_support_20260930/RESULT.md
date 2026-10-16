# Decision: source m10 slot counts do not transfer to the fixed degree-263 leaves

The one permitted census, Actions run
[36707051160](https://github.com/aburan28/crypto/actions/runs/36707051160),
returned `PASS_CENSUS` on source head
`3813d1879ca03483ac3445ce6d431f0319d25f94`. Do not dispatch that workflow
again. This is a structural representation census. It is not an advance,
engineering gain, relabelling, or accounting correction of an attack cost:
no target was sampled, no solver ran, and no logarithm was recovered.
`PDP_yield`, `full_ECDLP_cost`, and `method_crossover` stay null. The
single-slot count transfer fails on both frozen leaves. A difference in
physical-point counts is a representation effect, not an easier PDP, and it
does not discharge the review-gated m10 capacity attempt in
[PR #937](https://github.com/aburan28/crypto/pull/937).

## What was counted

The preregistration in [`PROTOCOL.md`](PROTOCOL.md) fixes source
`E_0: y² + xy = x³ + 1` over `GF(2^131)` and the normalized degree-263
leaves for kernel lines `[1,0]` and `[1,4]`. Each curve is scanned on ten
13-dimensional low slots and the 14-dimensional high slot of the unequal
arm, using the unchanged beta-3 normal-basis masks. The source curve is the
positive control: every low slot must have 7,977 physical points and the
high slot 16,125, matching the earlier source m10 screen. `GF(2^131)` has
no proper intermediate subfield over `GF(2)`.

Hosted replay and a fresh replay on this Linux x86_64 host
(`Linux-6.12.94+-x86_64-with-glibc2.39`, Python 3.12.8, checkout
`5885d629d8c30261ee8c83157fa51792406545c0`) agree:

| Check | Value |
| --- | --- |
| Status | `PASS` |
| Scans replayed | 33 |
| Raw masks replayed | 294,912 |
| Sample point images checked | 264 |
| Source positive controls | `PASS` |
| `result.json` SHA-256 | `71c6a748a3d65f61e26d139713498aa6731cf64b53c066886ed663f2489ea1d5` |

The producer peak RSS on the hosted run was 45,961,216 bytes. The second-host
receipt was written at 2026-09-30T11:39:25Z, inside the 1,200-second verifier
cap. The GitHub artifact `leaf-m10-support-36707051160` is 8,223,622 bytes,
digest `sha256:cfeb16421b261a5f7c9fc36bfb3a60349cd1f5a6357c1dd49e1874d3303d3344`.
The copied receipts are under
[`evidence/run_36707051160/`](evidence/run_36707051160/).

## Physical points

Low slots are listed in index order `low_0` through `low_9`. Every slot has
within-slot projected-class multiplicity `{1: classes, 2: 0, 3: 0, 4: 0}` and
zero nonzero-x projections to infinity, so the projected signed-class count
equals the liftable nonzero-x count. For every arm, the uncompressed
projected-sign union equals the sum of those slot class counts, so the slot
class sets in that arm are disjoint. The source low slots are the same size;
the leaf low slots are not.

| Curve | Balanced low physical points | High physical points | Balanced union | Unequal union |
| --- | --- | ---: | ---: | ---: |
| source | 7,977 ten times | 16,125 | 39,880 | 43,954 |
| leaf `[1,0]` | 8,145, 8,189, 8,223, 8,233, 8,063, 8,339, 8,121, 8,183, 8,161, 8,103 | 16,259 | 40,875 | 44,932 |
| leaf `[1,4]` | 8,161, 8,095, 8,323, 8,233, 8,243, 8,167, 8,235, 8,297, 8,157, 8,087 | 16,323 | 40,994 | 45,075 |

Leaf `[1,0]` low slots run from 8,063 to 8,339 physical points (86 to 362
above the source control). Leaf `[1,4]` low slots run from 8,087 to 8,323
(110 to 346 above the control). The high slots are 134 and 198 above 16,125.

`N` is the physical tuple product of the arm. `N/q` uses subgroup order
`680564733841876926932320129493409985129` and is a necessary counting
ceiling, not a hit rate. The ten displayed digits match the exact integer
ratio; the source unequal value is the same 3.0987505090 already on the
source m10 screen.

| Curve and arm | Exact physical tuple product `N` | Displayed `N/q` |
| --- | --- | ---: |
| source balanced | 1043268081613366923392350363142874173649 | 1.5329446704 |
| source unequal | 2108900315408742840629516059380574909125 | 3.0987505090 |
| leaf `[1,0]` balanced | 1334234841788902171362588446833961996535 | 1.9604818990 |
| leaf `[1,0]` unequal | 2663391564474617606406915353845719840597 | 3.9135021726 |
| leaf `[1,4]` balanced | 1373563353665728741818015662172737006025 | 2.0182699534 |
| leaf `[1,4]` unequal | 2747295015547811573666887593878885694075 | 4.0367872135 |

## What stays open

Equal or larger leaf counts do not restore a free leaf Frobenius action.
The next comparison, if one is run, has to match useful projected physical
points across original, transported, leaf-native, and exact pullback bases,
then measure natural-target PDP outcomes and complete cost. That comparison
is not this census, and it still depends on PR #937's independent exact-head
review and one-shot label. No `S`, matched-rho ratio, floor ratio, n=83
gate, or ECC2K-130 speedup follows from these slot counts.
