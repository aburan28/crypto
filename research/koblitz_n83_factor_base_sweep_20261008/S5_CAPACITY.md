# Balanced-S5 and compact-orbit transfer gate for N83

This is a source and capacity audit on study-branch revision
`c26ae1e3953823436cb0143d8eb77261cd784783`, after the one-hour pilot.
It performs no N83 index construction, relation search, rank run, or cold
comparison. The primary objective remains the complete single-target cold
runtime on the pinned `a=0` curve.

## Existing implementations and the exact mismatch

`examples/koblitz_s5_sat_instance.rs` is a substantial four-summand
balanced-S5 implementation. It has native XOR, orbit-factorized domains,
pair-sum tries, relative pair support, and group lifting. Its command-line
admission list is `n=7,11,13,17,19,23,37,41,53`; the point base and scalar
rank paths use `u64`, and `pair_sum_trie` is capped at `n<=53`. A degree-83
run therefore needs a full-width point/base/label port, an exact finite
domain, and independent model lifting. Source presence is not an N83 SAT
capacity or throughput receipt. The retained n=53 exceptional-state and
compact-extraction records in `docs/ic/BOUNDARY_TARGETS.md` belong to their
published bases and cannot be treated as this study's N83 measurements.

`examples/koblitz_orbit_dlp_fast.rs` has a `u128` field and group path for
`64<n<=127`, a compact regular-S3 root index, and a recorded n=83 `a=1`
online-after-setup run. That run used 600 orbit columns, the field polynomial
`z^83+z^7+z^4+z^2+1`, a 53-bit subgroup, and a different point-defined-base
JSONL schema. The study objects use `z^83+z^45+z^2+z+1` and the pinned
`n83.factor-base/v1` schema. Although the curves have the same abstract
field size, their polynomial coordinates cannot be interchanged. The wide
runner constructs `KoblitzCurve::new(a,n)` in its selected field basis and
checks a supplied header's modulus against it. It also converts the subgroup
order, Frobenius eigenvalue and rank rows to `u64`. The study's primary
`a=0` subgroup order is `2417851639230796216685689` (81 bits), so its
existing rank and log path cannot run that arm. A diagnostic `a=1` port
would still need an explicit field-basis isomorphism or a pinned-curve
constructor, an object/label importer and independent replay. The existing
runner's historical online timing excludes reusable preparation and is not
the cold objective here.

The current study branch has a separate exact wide relation-rank and
factor-log gate, but the compact-orbit producer is not wired to it. A
transfer should first bind the pinned modulus, generator, subgroup,
Frobenius labels and stored point coordinates; then replay compact-index
witnesses as group equations and wide modular rows on independent small
fixtures. A retained N83 construction under a memory/wall guard and a
natural relation/rank receipt are later gates. A full cold comparison comes
only after all phases, including failed work and I/O, are charged.

## Source-derived index size, not a timing projection

For `K` full signed-Frobenius columns, `build_index128` tries exactly
`83 K^2` ordered `(left,right,relative)` states. Let `S` be the number of
states for which the S3 root solver returns roots. The source retains one
`State128` per such state and allocates a `RootTable128` with
`C=max(16,next_power_of_two(4S))` slots. On this 64-bit build, a direct
`size_of` check gives 48 bytes per state and 32 bytes per table slot.
Thus the two vectors have at least `48S + 32C` bytes of final payload;
`Vec` spare capacity and temporary copies can increase peak allocation.
This excludes temporary arrays, point/label maps, rank rows and the process
runtime. Root-key deduplication does not shrink the allocated table. The
slot-count formula is exact for a realized `S`, while the byte formula is a
lower bound on those two final vectors; `S` is unmeasured on the retained
N83 bases.

The table below is reproducible with integer arithmetic. For each listed
`K`, calculate `T=83*K*K`, `S=T//2`,
`C=max(16,1<<(4*S-1).bit_length())`, and
`GiB=(48*S+32*C)/2**30`. The 48- and 32-byte layout values were checked
with `rustc` and `std::mem::size_of` on the pilot host; a different ABI
requires its own layout check.

| K | Exact states tried, `83K²` | Illustrative GiB for vectors if `S=floor(83K²/2)` |
| ---: | ---: | ---: |
| 64 | 339,968 | 0.039 |
| 256 | 5,439,488 | 0.622 |
| 600 | 29,880,000 | 2.668 |
| 900 | 67,230,000 | 9.503 |
| 1,182 | 115,961,292 | 10.592 |
| 2,048 | 348,127,232 | 39.781 |
| 4,096 | 1,392,508,928 | 159.125 |
| 8,192 | 5,570,035,712 | 636.500 |
| 16,627 | 22,945,941,707 | 2,560.882 |

The half-regular column is an **illustration**, not a measured rate,
probability bound, hardware-capacity result or N83 runtime estimate. The
`K=16,627` four-summand support ceiling in `SIZE_FRONTIER_V2.md` therefore
does not make a complete materialized root index feasible. An implicit or
segmented producer would need its own exact witness-preservation and
resource accounting. The prior n=83 `a=1` record reports about 5.75 GB
peak RSS at K=600, on its different base and field representation; it is
context for the source audit, not an observation for a study candidate.

## Decision for the sweep

The frozen v1/v2 `compact_s3_four_sum` dispositions remain recipes. Do not
rank the stored 54 bases by the historical n=83 online result, by the S5
half-regular illustration, or by construction speed. The smallest useful
next producer gate is a pinned-field, wide-order, replay-bound compact-index
adapter with a guarded K=64 construction; K=1,182 and larger need an
explicit memory strategy before a construction attempt. The balanced-S5
SAT example is an alternative source of four-summand algebra and domain
constraints, not a directly executable N83 backend.
