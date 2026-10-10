# Aligned S3 root-slot comparison

Replacing the packed 13-byte root slot with a 16-byte aligned full-key slot
preserved every checked S3 relation and scalar on the frozen N41 C800 and N53
C1200 controls. In a paired eight-point N53 C1200 screen, the packed/aligned
one-target online ratio was **0.956×**, with a target-bootstrap 95% interval
**[0.772×, 1.102×]**. The frozen five-percent advancement gate failed, so the
packed layout remains the current candidate.

The aligned slot stores a `u64` canonical key and a `u64` word holding the
same 40-bit witness value. It is explicitly 16-byte aligned. The prefetch,
lookup, hash, probe order, factor base, relation setup, and portable projective
target checks otherwise follow the combined source. Sixteen release tests
passed, including exact root-map and N41/N53 point checks. The aligned source
and its binary hash were committed in `CONTROL_FREEZE.json` before either
control ran.

## Exact controls

Both controls use the frozen T025 public point for their curve, relation seed
`20260928`, 14 preparation workers, one target worker, full requested rank,
and a 32 GiB RSS envelope. `run_control.py` compares the complete raw result
with the existing packed combined candidate, allowing changes only in timing,
peak RSS, root-slot width, and layout name. Both controls passed that semantic
comparison, recovered the known fixture scalar, and stayed under the memory
cap.

| Curve | Columns / usable base points | Aligned online ms | Aligned peak RSS | Same first relation / S3 calls |
| --- | ---: | ---: | ---: | --- |
| N41 | 800 / 65,600 | 0.084208 | 3,989,438,464 B | Yes / 32 |
| N53 | 1,200 / 127,200 | 0.994375 | 12,631,736,320 B | Yes / 3,280 |

The N53 candidate ID is
`IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0hd59cfb6d2ca4`.
The packed reference ID is
`IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0h7d64b2ad00e8`.
Their manifests preserve the exact layout, source, curve, factor base, and
stage configuration.

## Paired N53 screen

`SCREEN_FREEZE.json` fixes public targets T189–T196, alternating AB/BA order,
one fresh process per arm and point, the source and binary hashes, candidate
and workload IDs, 120-second per-run cap, and 32 GiB memory cap. The gate
requires 16 valid runs, a packed/aligned online geometric ratio of at least
1.05×, and a point-bootstrap lower endpoint above 1×. The five exclusive
online phases and all S3 calls are charged to their one public target.

| Point | Packed online ms | Aligned online ms | Packed/aligned |
| --- | ---: | ---: | ---: |
| T189 | 0.790917 | 1.647916 | 0.480× |
| T190 | 0.560167 | 0.600125 | 0.933× |
| T191 | 0.135208 | 0.127167 | 1.063× |
| T192 | 0.087542 | 0.082167 | 1.065× |
| T193 | 1.216708 | 0.983334 | 1.237× |
| T194 | 0.312333 | 0.321333 | 0.972× |
| T195 | 0.591458 | 0.575334 | 1.028× |
| T196 | 0.209500 | 0.188625 | 1.111× |

All 16 native solves reached rank 1,200, recovered the fixture scalar, and
returned the same first relation within each pair. Checked Sage independently
replayed all 16 scalars and verified the generator order. The geometric
packed/aligned S3-call ratio was exactly 1.000×. Peak RSS was
12,647,088,128 bytes. The aligned layout used roughly 0.6 GiB more observed
memory than the packed layout in this screen.

The setup-time ratios to each arm median ranged from 0.845–1.536× for packed
and 0.863–1.497× for aligned, failing the 10% noise gate. There is no
host-wide exclusive CPU partition receipt, so `controlled_speedup` remains
`null`. The entire 78-file archive is under `local/`, including compressed raw
reports, resource counters, stderr, checked Sage runtime information and
replay, failure-preserving ledger, deterministic bootstrap audit, and
`archive_receipt.json`. Reanalysis from the archived copy reproduced
`analysis.json` byte for byte.
