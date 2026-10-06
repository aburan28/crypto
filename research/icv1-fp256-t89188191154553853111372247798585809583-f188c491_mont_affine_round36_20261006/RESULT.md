# P-256 Montgomery-affine pair table, round 36: result

Date run: 2026-10-06

The registered 64-byte Montgomery-affine representations do not preserve
round 33's rho margin on the 256-MiB decision row.  The selected complete-add
row at prefetch distance 1 costs `1.152984` extra matched additions per
segment, against a complete admissible budget of `0.334411`.  It projects to
`1.003503` rho at the median and `1.005283` at the registered interval upper
endpoint.

The specialized mixed-add formula is exact and removes one of the complete
formula's 17 field multiplications.  Its arithmetic saving is deliberately
present in both the table row and its hot-operand denominator, so it is not
credited as selector-only progress.  Its best 256-MiB row projects to
`1.004379` rho at the median.

All registered representations replay exactly.  The complete projected
affine table falls from 3.318 TB to 2.212 TB and remains below `2^50` bytes,
but the measured access-stage gate fails before full deployment.  This is a
reproducible negative result for the registered affine layouts and prefetch
schedules, not a lower bound for every memory system.  No full-depth
unplanted relation was attempted.

## Frozen boundary

Round 35 is imported by exact SHA-256
`90704e913a0407fe773ad7550cd5b68c5ead54077ff485978bd51cdf21bcfc89`.
The experiment preserves:

```text
complete pair entries                  34,562,148,612
Montgomery-affine entry                64 bytes
complete projected affine table        2,211,977,511,168 bytes (2^41.008)
complete projected projective table    3,317,966,266,752 bytes (2^41.593)
round-33 corrected stage ratio         0.998569034150286 rho
mean segment capacity                  224.361249307511
allowed access overhead / segment      0.334410508423 additions
```

Every timed candidate processes the same 65,536 segments and eight table
lookups.  Complete-affine rows use the existing 17-multiplication complete
formula.  Mixed-affine rows use its exact `z2 = 1` specialization with 16
multiplications.  Each is divided only by an L1-resident direct row using the
same representation and formula.  The projective control retains its own
projective denominator.

## Decision row: 256 MiB affine

| variant | prefetch distance | median extra matched additions | interval upper | median / rho | upper / rho |
|:--|--:|--:|--:|--:|--:|
| **complete affine** | **1** | **1.152984** | **1.569031** | **1.003503** | **1.005283** |
| complete affine | 2 | 1.222227 | 1.720770 | 1.003799 | 1.005932 |
| complete affine | 4 | 1.171828 | 2.013985 | 1.003583 | 1.007187 |
| complete affine | 8 | 1.371238 | 1.900783 | 1.004437 | 1.006703 |
| complete affine | 16 | 1.345401 | 1.897991 | 1.004326 | 1.006691 |
| complete affine | 32 | 1.258444 | 2.139705 | 1.003954 | 1.007725 |
| complete affine | 64 | 1.601352 | 2.077334 | 1.005421 | 1.007458 |
| mixed affine | 1 | 1.564177 | 2.197567 | 1.005262 | 1.007973 |
| mixed affine | 2 | 1.517871 | 2.299295 | 1.005064 | 1.008408 |
| mixed affine | 4 | 1.506912 | 2.309345 | 1.005017 | 1.008451 |
| mixed affine | 8 | 1.480027 | 2.144423 | 1.004902 | 1.007745 |
| **mixed affine** | **16** | **1.357852** | **1.931952** | **1.004379** | **1.006836** |
| mixed affine | 32 | 1.418022 | 1.809179 | 1.004637 | 1.006311 |
| mixed affine | 64 | 1.460019 | 2.268698 | 1.004817 | 1.008277 |

The selected row is the smallest median projected ratio among affine
candidates on the registered `2^22`-entry working set.  It fails both the
median and interval-upper timing gates.  The upper interval is deliberately
conservative: candidate rank 7 is divided by matched-direct rank 1 across
nine deterministically rotated repetitions.

At the same entry count, the projective control occupies 384 MiB and reaches
`1.002302` rho at the median and `1.007242` at the upper endpoint.  It is a
control, not an affine selection candidate.

## Working-set transition

| affine working set | selected affine variant | median extra additions | interval upper | median / rho | upper / rho |
|--:|:--|--:|--:|--:|--:|
| 2 MiB | mixed, distance 2 | 0.000000 | 0.486243 | 0.998569 | 1.000650 |
| 64 MiB | complete, distance 16 | 0.883469 | 1.275452 | 1.002349 | 1.004027 |
| 256 MiB | complete, distance 1 | 1.152984 | 1.569031 | 1.003503 | 1.005283 |

Even the 2-MiB row misses at its registered upper endpoint.  Both larger
working sets are above rho at every reported timing endpoint, so the data do
not support extrapolating a sub-rho 2.212-TB deployment.

## Exact controls

The semantic run uses the same 4,096 hash-selected P-256 points as round 35
and 4,096 deterministic reference segments.  Every affine output is compared
with the projective control.  Identity, doubling, inverse and distinct-point
cases are also checked for both affine formulas.

```text
variants checked          15
special cases checked      8
group mismatches           0
special-case mismatches    0
missing requests           0
duplicate requests         0
out-of-range accesses      0
wrong addition counts      0
curve failures             0
false positives            0
false negatives            0
```

Deterministic digests:

```text
pool          2426363dfee219f092cf88778966c8e7de5dee45b05a83c52dac824a1bc6dd65
indices       6465aff9b9ca1d60c98c021f86a63645209d8b183ae59cb4095ac6e21d355a90
affine table  9556d2a9a961e4dfb859816b89e40f4ebc5bce5344c8cc392a76c2ed528521b0
outputs       48f7f05287998eb0598728a4a116b65db44a7f7e39417afa522f87d8f01bea2c
semantic      b19f2955b94f3e9de2475fc3e054a4fd3d8aeb9cfc39ed3727cedfb64441d562
raw timing    1a38876b8771819a8c828de6b4dcda7be0f07fbb1fe6274e0841fecdfb04cfc9
```

Two pre-publication full release replays reproduced the semantic digest
exactly.  Their volatile timing rows are not substituted for the canonical
run.  The measured host is an AMD EPYC 9V74 with Rust 1.98.1.  Results are
one-host implementation measurements, not a subgroup-order exponent or a
universal hardware lower bound.

## Requirement status

| requirement | status | evidence or gap |
|:--|:--|:--|
| exact group replay | verified | all 15 variants and eight exceptional cases; zero group, range, count, curve, FP or FN failures |
| Montgomery-affine access at or below rho | failed | selected 256-MiB row is `1.003503` median and `1.005283` upper |
| complete measured time at or below rho | failed | access gate fails before the 2.212-TB table is materialized |
| storage below `2^50` | projected | full affine table is `2^41.008` bytes |
| proved P-256 usable-relation probability | not established | imported occurrence upper bound remains unproved as achieved coverage |
| structured residual degree at most 5 | not established | no algebraic residual solve in this round |
| collection below `2^120` and per row below `2^103` | not established | no usable full-depth relation or collector |
| demonstrably non-generic end-to-end method | not established | imported known-log transport remains representation-based |

## Decision

Reject the registered Montgomery-affine complete/mixed representations and
prefetch schedules.  They reduce projected storage by one third but do not
reach rho under matched arithmetic.  Do not promote the selector and do not
run a full-depth P-256 relation.  A further iteration must change the
access-generation coupling rather than only the stored coordinate system.

## Reproduction

```bash
cargo test --bin p256_mont_affine_table
cargo run --release --bin p256_mont_affine_table -- \
  --round35 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_projective_locality_round35_20261006/projective-locality-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_mont_affine_round36_20261006/mont-affine-result.json
```

Canonical result artifact: `mont-affine-result.json`, 213,497 bytes, SHA-256
`f145231c098da83878e890d340f8ae921598519b0f2626bcc7c4cd8f9446d015`.
