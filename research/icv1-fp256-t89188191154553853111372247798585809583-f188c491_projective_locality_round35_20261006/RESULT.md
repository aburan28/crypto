# P-256 projective pair-table locality, round 35: result

Date run: 2026-10-06

Neither registered locality mechanism preserves round 33's rho margin on the
384-MiB decision row.  Software prefetching is materially useful: the selected
distance-16 row cuts median table overhead from `3.149926` to `0.680998`
P-256-addition equivalents per segment.  The complete admissible budget is
only `0.334411`, so the selected row projects to `1.001483` rho at the median
and `1.005622` at the registered interval upper endpoint.

The best global radix row, using all 65,536 segments in one batch, reaches
`1.006907` rho at the median.  Sorting and scatter state cost more than the
memory latency they remove at the tested batch sizes.

All variants replay exactly.  This is a reproducible negative result for the
registered prefetch and radix implementations, not a claim about every
possible memory system.  No full-depth unplanted relation was attempted.

## Frozen boundary

Round 34 is imported by exact SHA-256.  The experiment preserves:

```text
projective entry                      96 bytes
complete entries                      34,562,148,612
complete projected table              3,317,966,266,752 bytes
round-33 corrected stage ratio        0.998569034150286 rho
mean segment capacity                 224.361249307511
allowed locality overhead / segment   0.334410508423 additions
```

Every timed candidate processes the same 65,536 segments, reads all eight
projective pair entries and performs all eight complete additions.  Index
generation is outside every timed row.  Request filling, counter clearing,
both radix passes, accumulator initialization, sorted reads and scatter
additions are inside each radix timing.

## Decision row: 384 MiB

| variant | parameter | median extra additions / segment | interval upper | median / rho | upper / rho |
|:--|--:|--:|--:|--:|--:|
| random control | - | 3.149926 | 3.869044 | 1.012048 | 1.015125 |
| prefetch | 1 | 0.906343 | 1.503647 | 1.002447 | 1.005003 |
| prefetch | 2 | 0.875896 | 1.722628 | 1.002317 | 1.005940 |
| prefetch | 4 | 0.848525 | 1.456285 | 1.002200 | 1.004801 |
| prefetch | 8 | 0.850325 | 1.607566 | 1.002208 | 1.005448 |
| **prefetch** | **16** | **0.680998** | **1.648164** | **1.001483** | **1.005622** |
| prefetch | 32 | 0.741429 | 1.533477 | 1.001742 | 1.005131 |
| prefetch | 64 | 0.818921 | 2.087435 | 1.002073 | 1.007501 |
| radix global | 256 | 3.109045 | 5.487190 | 1.011873 | 1.022049 |
| radix global | 1,024 | 3.095387 | 4.258264 | 1.011814 | 1.016790 |
| radix global | 4,096 | 2.998680 | 3.567028 | 1.011401 | 1.013833 |
| radix global | 16,384 | 2.300220 | 3.164806 | 1.008412 | 1.012111 |
| **radix global** | **65,536** | **1.948615** | **2.645943** | **1.006907** | **1.009891** |

The selected row is the smallest median ratio on the registered `2^22`-entry
working set.  It fails both the median and interval-upper timing gates.  The
upper interval is deliberately conservative: layout rank 7 is divided by
direct-addition rank 1 across nine deterministically rotated repetitions.

## Cache transition

| working set | selected variant | median extra additions | interval upper | median / rho | upper / rho |
|--:|:--|--:|--:|--:|--:|
| 3 MiB | prefetch-8 | 0.000000 | 0.384037 | 0.998569 | 1.000212 |
| 96 MiB | prefetch-32 | 0.813058 | 1.272054 | 1.002048 | 1.004012 |
| 384 MiB | prefetch-16 | 0.680998 | 1.648164 | 1.001483 | 1.005622 |

Even the 3-MiB row misses at its registered upper endpoint.  The 96- and
384-MiB results are consistent with a DRAM-latency regime; neither supports
projecting a sub-rho 3.318-TB deployment.

## Exact controls

The semantic run uses 4,096 hash-selected P-256 points and 4,096 deterministic
reference segments.  Large radix batches repeat that fixed reference stream
to fill the registered batch while retaining exact per-segment comparisons.

```text
variants checked          13
group mismatches           0
missing requests           0
duplicate requests         0
out-of-range accesses      0
wrong addition counts      0
false positives            0
false negatives            0
```

Deterministic digests:

```text
pool             2426363dfee219f092cf88778966c8e7de5dee45b05a83c52dac824a1bc6dd65
indices          6465aff9b9ca1d60c98c021f86a63645209d8b183ae59cb4095ac6e21d355a90
sorted requests  91e8f8c55af9eef7142d114cd3b095b2f11bc5c248c7782bdd7b5a863256467c
outputs          48f7f05287998eb0598728a4a116b65db44a7f7e39417afa522f87d8f01bea2c
semantic         a69c1c504c8fc2286d10d1409b849f03e92eb213263af597aedc2c8b94090f35
raw timing       06c7bda51d24d9ec397b788d78fe6902fe9b67cd42c212892e42538ce680decd
```

An independent full release replay reproduced semantic digest
`a69c1c504c8fc2286d10d1409b849f03e92eb213263af597aedc2c8b94090f35`
exactly; its expected volatile timing digest was
`76a798590e9d38b80d2721941d3eb44b4e52e04622def7973b57949ec9853b80`.

The measured host is an AMD EPYC 9V74 with Rust 1.98.1.  Results are one-host
implementation measurements, not a subgroup-order exponent or a universal
hardware lower bound.

## Requirement status

| requirement | status | evidence or gap |
|:--|:--|:--|
| exact group and request replay | verified | all 13 variants; zero group, permutation, range or addition-count failures |
| projective locality at or below rho | failed | selected 384-MiB row is `1.001483` median and `1.005622` upper |
| complete measured time at or below rho | failed | locality gate fails before the 3.318-TB table is materialized |
| storage below `2^50` | projected | full projective table is `2^41.593` bytes |
| proved P-256 usable-relation probability | not established | imported occurrence upper bound remains unproved as achieved coverage |
| structured residual degree at most 5 | not established | no algebraic residual solve in this round |
| collection below `2^120` and per row below `2^103` | not established | no usable full-depth relation or collector |
| demonstrably non-generic end-to-end method | not established | imported known-log transport remains representation-based |

## Decision

Reject the registered software-prefetch and two-pass radix locality schedules.
They substantially narrow the memory gap but do not reach rho.  Do not promote
the selector and do not run a full-depth P-256 relation.  A further iteration
must change the representation/access coupling rather than tune prefetch
distance or radix batch size within this design.

## Reproduction

```bash
cargo test --bin p256_projective_locality
cargo run --release --bin p256_projective_locality -- \
  --round34 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_pair_table_round34_20261006/pair-table-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_projective_locality_round35_20261006/projective-locality-result.json
```

Canonical result artifact: `projective-locality-result.json`, 159,053 bytes,
SHA-256
`90704e913a0407fe773ad7550cd5b68c5ead54077ff485978bd51cdf21bcfc89`.
