# C3200 projective target checks

This candidate keeps the exact C3200 factor base, S3 root index, relation
collector, rank solve, and target decomposition. It uses the repository's
portable López–Dahab projective scalar engine for the target subgroup check
and the final recovered-scalar replay. Each multiplication performs one
field inversion at its end, replacing the affine double-and-add path that
inverts at every point operation. Both checks remain inside the one-target
online interval and retain the same subgroup and group-equality conditions.

The source includes a cross-engine exact-output test for N41 and N53 on
zero, small, middle, and subgroup-boundary scalars. After it passes, the
frozen T025 control will require the same scalar, first relation witness,
factor-base digest, full rank, root-table layout, and every non-timing report
field as the verified C3200 reference. It will record the five charged target
phases, peak RSS, compressed raw report, and a separate Sage replay. The
validation and recovery phase costs will be compared first; a fresh-target
paired online screen is needed for the owner-set 2× threshold.

After writing the primary T025 result, the control runs 256 balanced
affine/projective check pairs against that same public point and scalar.
Each row records subgroup-check and replay wall time, main-thread CPU time,
and exact outcomes. The repetitions share setup and a known point, so they
measure only this check stage. Their raw ledger is archived separately from
the one-target online result.

## Verified T025 control

The frozen source `d1854989e`, release binary, point T025, and reference
report are bound by `CONTROL_FREEZE.json`. All 17 release tests passed. The
control recovered scalar `429728658711`, matched the reference's first
four-point relation and 498 target S3 calls, reached rank 3200 before target
input, and verified all 339,200 usable base points. Its peak RSS was
30,810,456,064 bytes within the 32 GiB envelope. Checked Sage independently
verified the subgroup order and scalar replay. The canonical candidate is
`IC1N53Ckb1fb339200PDP4rootRCsampleLAgaussTDdirectISO0h0e595bd6143a`.

| T025 target phase | Affine reference ms | Projective ms |
| --- | ---: | ---: |
| Target query, including point validation | 0.632499 | 0.010042 |
| Charged target PDP | 4.775083 | 2.939542 |
| Relation check | 0.043375 | 0.034958 |
| Target descent | 0.038501 | 0.044083 |
| Scalar recovery check | 0.061583 | 0.008000 |
| Complete target online interval | 5.551041 | 3.036625 |

The local complete-interval ratio is 1.828× on this same known point. The
reusable setup times were 173.916 and 80.781 seconds, a strong host-noise
control failure. The host has no exclusive CPU partition receipt, so the
controlled speedup remains unknown. The two T025 online values are retained
as phase diagnostics rather than a one-target speedup estimate.

The 256 balanced check pairs all passed subgroup and scalar replay checks.
The median paired affine/projective wall ratio was **14.017×** for validation
(median 69.855 versus 4.980 microseconds) and **11.380×** for replay (61.188
versus 5.333 microseconds). The corresponding main-thread CPU ratios were
13.287× and 10.886×. These measurements are conditional on one process and
one point. `analyze.py` independently checks every compressed artifact hash,
candidate digest, native result, paired row, and Sage receipt before writing
`analysis.json`.

On the earlier frozen eight-point C3200 screen, deleting validation and
recovery-check time **entirely** would give only a 1.384× geometric online
ratio on those recorded phase totals. The next full-interval iteration needs
to reduce target PDP work; further scalar-check tuning cannot reach the 2×
screen threshold on that panel by itself.
