# Eight-table GF(2) elimination: correctness rejection

The [protocol](PROTOCOL.md) was committed before code or timing. The candidate allowed eight eight-pivot Gray-code tables in one elimination block, with a 4 MiB logical table budget and a four-table fallback on smaller matrices. It was rejected at the exact-output gate before an isolated performance run. The tested implementation and focused test are retained in [correctness_rejected_candidate.patch](correctness_rejected_candidate.patch); the active runtime was restored to the reference source.

## Exact output check

The release [unit-test receipt](CORRECTNESS_CHECK.log) compares four and eight tables on a 90 × 129 dense packed matrix. Both paths return the same rank and textbook fully reduced row echelon form, but their **echelon row 0 differs**. The first differing rows are `four=[1442764906403201025, 7033693675666452550, 0]` and `eight=[9295883729194975233, 12706005379515376712, 0]`. The test therefore fails as required by the frozen raw-output contract. This is expected from grouping more pivots in one block: pivot rows inside a block are reduced against each other, while pivot rows in different blocks are left in echelon form.

A local release [complete F5 output check](F5_OUTPUT_CHECK.json) ran the same binary with `KIC_GF2_TABLES=4/8` on all seven frozen cases at seed `0`, one Rayon thread. It discarded all emitted wall times; this is a correctness check, not an isolated speed benchmark. The [build receipt](BUILD_CHECK.log) and JSON record source/binary hashes, architecture, compiler, exact command, flags, routes and structural counters. The tested patch SHA-256 is `6e6dfcc2ba94a45cec073a36f5e3112d2df5ad768007e28fd1faa0be26b1f0dc`; the unit-test receipt SHA-256 is `732c547516300ce22c35ac00bdc66c22ce5dfee65541e6c4349fc51bc477fc9d`. The changed route was selected only for n20 and n24 degree-4, as designed. All other five cases had equal raw and canonical fingerprints.

| F5 case | Tables | Raw row fingerprint | Canonical row-space fingerprint | Rank | Reduction word ops | Actual table allocation |
| --- | ---: | --- | --- | ---: | ---: | ---: |
| n20 m20 d4 | 4 | `7280d506ac81759a` | `d4a1eeefce856713` | 4,010 | 24,500,060 | 794,624 B |
| n20 m20 d4 | 8 | `6a0afca85dd538d7` | `d4a1eeefce856713` | 4,010 | 27,083,594 | 1,589,248 B |
| n24 m24 d4 | 4 | `3f659516eff553b8` | `ed5234ba018bc079` | 6,924 | 100,213,183 | 1,662,976 B |
| n24 m24 d4 | 8 | `65b535e4d47efb90` | `ed5234ba018bc079` | 6,924 | 109,612,979 | 3,325,952 B |

The n24 candidate changed the raw `Vec<F2BoolPoly>` output and increased counted reduction work by about 9.4%; n20 increased by about 10.5%. Both actual table allocations stayed inside 4 MiB, but the exact-output requirement had already failed. No A/A or A/B wall-time ratio is available or claimed. The experiment yields no accepted complete-call gain and no progress toward a measured 2× result. It is solver-stage evidence only; one-target IC DLP cost remains unknown.

The result distinguishes this table-count change from the prior six-column table-width and column-panel experiments. A future wide-block kernel would need to restore the original 32-pivot echelon rows, including the exact raw fingerprints, before any paired performance comparison. That reconstruction was outside this bounded experiment.
