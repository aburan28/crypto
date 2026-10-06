# P-256 pair-table representation audit, round 34: result

Date run: 2026-10-06

The 33-byte table assumed by round 33 is not directly addable.  Exact P-256
decompression costs about 110 extra group-addition equivalents per segment on
the measured host, raising the candidate-favouring stage projection from
`0.998569` to `1.469555` rho even in the smallest working set.

A 96-byte projective table is the only layout whose cache-resident lower bound
retains the stage margin.  Its 96-KiB row has a median overhead of `0.014566`
and an interval upper endpoint of `0.200739` additions per segment, below the
registered `0.334411` budget.  That result does **not** establish complete-time
or end-to-end parity.  At the largest measured 384-MiB working set the same
layout costs `3.612595` extra additions per segment and projects to `1.014028`
rho (`1.016376` at the interval upper endpoint).  The complete projected table
is 3.318 TB and was not materialized.

Round 34 therefore preserves only a cache-resident projective lower bound as a
stage candidate.  It does not promote the selector and does not attempt a
full-depth unplanted relation.

## Frozen parity budget

The exact round-33 artifact was imported by SHA-256.  For cutoff 219,

```text
round-33 ratio with target correction = 0.998569034150286 rho
mean segment capacity D              = 224.361249307511
allowed extra additions / sample     = 0.001483886469
allowed extra additions / segment    = 0.334410508423
```

The budget is the entire distance from round 33's corrected stage estimate to
rho in the frozen local-oracle unit.  It is only 4.18% of one addition per
lookup across the eight-entry setup.

## Representation results

Every complete layout remains below the original `2^50`-byte storage gate:

| layout | bytes / entry | complete projected bytes | log2 bytes | directly addable |
|:--|--:|--:|--:|:--|
| compressed SEC1 | 33 | 1,140,550,904,196 | 40.052868 | no |
| affine x,y | 64 | 2,211,977,511,168 | 41.008474 | no |
| projective X,Y,Z | 96 | 3,317,966,266,752 | 41.593436 | yes |
| known coefficient | 32 | 1,105,988,755,584 | 40.008474 | no; requires scalar multiplication |

The candidate-favouring row below is the fastest median working set measured
for each layout.  Intervals use sorted ranks 1, 4 and 7 from nine interleaved
release repetitions.

| layout | measured entries | median extra additions / segment | interval upper | median projected / rho | upper / rho | local lower-bound gate |
|:--|--:|--:|--:|--:|--:|:--|
| compressed33 | 1,024 | 110.067343 | 114.071793 | 1.469555 | 1.486690 | fail |
| affine64 | 1,024 | 0.563247 | 0.756132 | 1.000979 | 1.001805 | fail |
| projective96 | 1,024 | 0.014566 | 0.200739 | 0.998631 | 0.999428 | pass locally |
| coefficient32 | 1,024 | 699.745387 | 718.024860 | 3.992828 | 4.071047 | fail |

The corresponding largest measured working sets all fail:

| layout | measured table | median extra additions / segment | interval upper | median projected / rho | upper / rho |
|:--|--:|--:|--:|--:|--:|
| compressed33 | 264 MiB | 116.468055 | 124.504485 | 1.496944 | 1.531332 |
| affine64 | 256 MiB | 3.938756 | 4.083065 | 1.015423 | 1.016041 |
| projective96 | 384 MiB | 3.612595 | 4.161442 | 1.014028 | 1.016376 |
| coefficient32 | 256 MiB | 714.993555 | 746.188409 | 4.058076 | 4.191561 |

The projective cache row is an optimistic lower bound, not a projection that a
3.318-TB random-access table will behave like 96 KiB.  Conversely, this host
measurement is not a hardware-independent lower bound against all possible
batched or locality-preserving implementations.  The concrete unresolved
route is therefore projective storage with a global locality schedule, not
compressed-point lookup.

## Exact replay and deterministic controls

The release run generated 4,096 hash-selected nonzero P-256 scalars and their
compressed, affine, projective and coefficient representations.  It checked
4,096 individual round trips and 4,096 independent eight-entry segments.
Every compressed square root was validated on-curve; every accumulated point
was compared with the independently reconstructed projective reference.

```text
compressed failures      0
affine failures          0
projective failures      0
coefficient failures     0
invalid sqrt/curve cases 0
sizeof(projective point) 96 bytes
```

Deterministic evidence digests:

```text
pool       f9047b0c08436b1271d33cb2efb8d8004fc4227ab14bd44ca79a578d890a5f8d
indices    809a8df8f642037e09482251bc264347e0fe0bda5e0d8a984c3a471cecf45c6e
outputs    91fe172708c1a34b86869ffab8ba3585eeabd65eb766e4e933fcb807e0b8fdfa
semantic   cd34ea31fa8e7eb277a75fc0293fa2569d40df02dc3be3eab08e4a66f37d13e7
raw timing 49cab78929c72b8349fc261ad12bbd803352e52fada9cbda85af2948afe6f950
```

Repeated full release executions reproduced the semantic digest exactly;
host timings and their raw digest were intentionally retained as measured
volatile evidence rather than forced byte-identical.

The host was an AMD EPYC 9V74 running the repository release profile with
Rust 1.98.1.  Timing rows include wall time, process CPU time, execution order,
segment counts and checksums.  These measurements compare layouts on one host;
they are not a fitted subgroup-order exponent.

## Requirement status

| requirement | status | evidence or gap |
|:--|:--|:--|
| exact representation and group replay | verified | 4,096 round trips and 4,096 cross-layout segments; zero failures |
| peak materialized storage below `2^50` | projected for every layout | complete projections range from `2^40.008` to `2^41.593` bytes |
| complete measured time at or below rho | not established | only the 96-KiB projective lower bound passes; every largest measured working set fails and the 3.318-TB table was not built |
| proved P-256 usable-relation probability | not established | round 33 still uses an occurrence upper bound as coverage |
| structured residual degree at most 5 | not established | no new algebraic residual system is solved here |
| collection below `2^120` and per row below `2^103` | not established | no usable full-depth relation or collector |
| demonstrably non-generic end-to-end method | not established | the imported known-log transport remains representation-based |

## Decision

Reject compressed, affine and coefficient pair-table implementations for the
round-33 margin.  Preserve only the projective cache-resident lower bound as a
stage candidate.  The next admissible experiment must construct and price a
batched/locality-preserving projective lookup schedule, including sorting,
scatter state and memory traffic.  No full-depth P-256 relation is authorized.

## Reproduction

```bash
cargo test --bin p256_pair_table_accounting
cargo run --release --bin p256_pair_table_accounting -- \
  --round33 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_biased_restart_round33_20261006/biased-restart-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_pair_table_round34_20261006/pair-table-result.json
```

Canonical result artifact: `pair-table-result.json`, 129,434 bytes, SHA-256
`df68204433a4abf408ac2c39f86d801da7db296056f030360fae8f719bca6a4b`.
