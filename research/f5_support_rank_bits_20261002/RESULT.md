# Width 6 lowers exact F5 certificate XORs; wall time unproven

The corrected native screen ran all 16 cells on the frozen four n24
degree-4 systems. Every call returned the seven-case suite with the
same exact rank, canonical row space, raw original-row fingerprint,
criterion, row counts, columns, and smaller-case nontiming outputs.
All four widths certified the cut-19 matrix and returned the unchanged
original rows.

| Rank table pattern width | Four-seed counted reduction word XORs | Work / width 8 | Preregistered 90% gate |
| ---: | ---: | ---: | --- |
| 5 | 25,463,537–25,550,188 | 90.57–90.73% | miss |
| **6** | **24,233,029–24,324,342** | **86.15–86.32%** | **pass** |
| 7 | 25,111,954–25,188,124 | 89.27–89.42% | pass |
| 8 | 28,114,847–28,198,391 | 100% | reference |

Width 6 is the best counted-work option. It removes about 3.9 million
64-bit word XORs from the fixed cut-19 rank certificate. The one-shot
Apple ARM64 complete-call and reduction times vary across seeds; these
are nonpromoting diagnostics. A paired complete-call comparison is
registered next. The further-2× claim still requires the unchanged
physical isolated Linux gate against selective echelon.

An initial pilot was **invalid as a reference comparison**: its
selector forced width 8 even for small lower-degree groups whose normal
width is narrower. It is preserved as `PILOT_INVALID.json.gz` with
`PILOT_RAW_SHA256.txt` and `PILOT_ARCHIVE_SHA256.txt`. The corrected
selector caps each group's normal width, so explicit width 8 now
reproduces the fixed cut-19 source's word-XOR counts exactly. The
corrected full receipt is `RECEIPT.json.gz`, with raw and archive
hashes in `RAW_SHA256.txt` and `RECEIPT.sha256`.

The corrected screen driver source SHA-256 is
`b44bced97c8fbdbc5c0ad10eee5620881ab7d5d59e211a7c10de08c3f98d20a7`;
the experimental GF(2) kernel source SHA-256 is
`4b5cacbcc7f1150ed950d5d468c6847678f73fa3a93b736863ebc6f5a0dd2bf9`.
Every run's binary/source hashes, output, phases, and status are in the
receipt. This is a matrix-F5 solver-stage result, not an IC online or
rho speedup.
