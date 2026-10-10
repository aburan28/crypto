# Combined S3 target prefetch and projective scalar checks

This candidate retains the full-rank N53 C1200 factor base, relation
collector, root table, and one-worker target search. It prepares up to 32
canonical partner keys per eight-state batch, prefetches their exact root
table home slots, and then performs the original witness checks in order.
Target subgroup validation and final scalar replay use the portable
López–Dahab projective multiplier. Both changes are charged within the
single-target online interval.

The source is derived from the verified late-prefetch kernel and the
verified projective scalar-check kernel. Its 16 release tests include exact
projective/affine agreement on N41 and N53 subgroup-boundary scalars. The
T025 control compares the same public point with the frozen C1200 affine,
non-prefetch reference. It requires the same first four-point relation,
factor-base digest, full rank, recovered scalar, and all unchanged report
fields. Prefetch look-ahead may increase the final batch's S3 call count;
the complete charged count remains in the report. A fresh N53 holdout screen
will follow only if the correctness and 32 GiB memory gates pass.

The T025 control passed: all unchanged report fields match the reference,
including its first relation, and checked Sage independently replayed scalar
`429728658711`. The target made 3,280 charged S3 calls versus 3,272 before
look-ahead prefetch. Its exploratory online interval was 0.949 ms versus
1.701 ms for the single earlier reference run. These T025 runs select the
candidate; they do not enter the fresh holdout ratio.

`SCREEN_FREEZE.json` binds previously unused N53 points T189–T196, both
candidate manifests, source and release binary hashes, 14 preparation
workers, one target worker, 32 GiB peak RSS, and alternating AB/BA order.
Every arm receives one point per fresh process with full reusable preparation
outside the online clock. The prespecified advancement gate is all 16 scalar
recoveries verified and a held-out geometric reference/combined online ratio
of at least **2×**. The host-isolation and 10% setup-noise gates still apply
before a controlled wall-time claim.

## Fresh eight-point screening result

The combined candidate is
`IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0h7d64b2ad00e8`.
All 16 one-target solves verified the fixture scalar and native group
equation, reached prepared rank 1200 with 127,200 usable base points, and
stayed below 32 GiB. Checked Sage independently replayed all 16 scalars and
the generator order. Every pair retained the same first relation. The
maximum peak RSS was 12,019,515,392 bytes.

| Point | Reference online ms | Combined online ms | Reference/combined | S3 calls, reference/combined |
| --- | ---: | ---: | ---: | ---: |
| T189 | 1.325 | 0.727 | 1.821× | 2,430 / 2,432 |
| T190 | 1.097 | 0.549 | 1.999× | 1,882 / 1,888 |
| T191 | 0.359 | 0.157 | 2.281× | 298 / 304 |
| T192 | 0.201 | 0.093 | 2.155× | 90 / 96 |
| T193 | 2.079 | 1.096 | 1.897× | 3,940 / 3,952 |
| T194 | 0.695 | 0.360 | 1.931× | 970 / 976 |
| T195 | 1.327 | 0.645 | 2.057× | 1,954 / 1,968 |
| T196 | 0.469 | 0.221 | 2.123× | 568 / 576 |

The frozen local advancement gate passed with a target-level geometric
online ratio of **2.028×**, bootstrap 95% interval **[1.933×, 2.129×]**.
The call-count ratio was 0.985×: extra look-ahead calls in the final batch
are charged rather than omitted. The unchanged setup times ranged from
0.846× to 1.149× their reference-arm median and from 0.891× to 1.153×
their combined-arm median. Both arms fail the 10% host-noise gate, and this
host has no exclusive CPU partition receipt. This result advances the
candidate to independent confirmation; the controlled speedup remains
unknown.

`analysis.json` checks every raw report, hash, candidate/workload/run ID,
phase sum, point, scalar, rank, and first relation. Reanalysis from the
archived `local/paired_screen/` directory reproduced the JSON byte for byte.
`local/archive_receipt.json` binds all 69 archived files by SHA-256. The
separate checked-Sage runtime and replay receipts are also archived there.
