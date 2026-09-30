# One-pass matrix-F5 output unpack: rejected at the local screen

The [protocol](PROTOCOL.md) was committed and placed in draft PR #1035 before candidate code or timing. The candidate decoded each packed row once, estimated output capacity from its echelon row index, and grew the vector between 64-bit input words when needed. Its exact tested source and screen are in [measured_candidate.patch](measured_candidate.patch). The opt-in runtime code was removed after the frozen local gate failed. This is a matrix-F5 solver-stage diagnostic, not a one-target IC DLP or rho speedup.

The local Apple ARM64 screen ran one release binary (SHA-256 `a34fa706b762ad1b5cf3aebd8f603189dd4d9c1c4caf30960bdb8f9e1a115249`) through 44 processes: two seeds, two warmups, five reference/reference A/A pairs and five alternating reference/candidate pairs per seed. Every process returned all seven F5 cases. Across both arms and every case, raw and canonical output fingerprints, rank, output terms, column and row counts, criterion/build counts and reduction word operations matched. The candidate route was selected only on n24 degree-4. This was a local `tools/isolated_bench.py busy` screen, not an exclusive CPU reservation or an accepted x86-64 wall-clock claim.

| Seed | Full-call reference/candidate paired median | Exact 95% bootstrap interval | A/A range | Unpack reference/candidate median | Candidate row-capacity growths |
| --- | ---: | ---: | ---: | ---: | ---: |
| `0` | **0.992×** | 0.946–1.069× | 0.932–1.000× | 0.876× | 2,172 |
| `badc0de1` | **0.953×** | 0.925–0.964× | 0.959–1.063× | 0.857× | 2,156 |

For the frozen seed, separate marginal medians were 61.79 ms reference and 62.46 ms candidate for the complete n24 call, and 17.69 versus 20.61 ms for unpacking. On the holdout they were 59.83 versus 62.80 ms complete and 17.82 versus 20.73 ms unpacking. The candidate's maximum allocated row capacity was 8,438 and 8,436 terms respectively; its 2,100-plus growths are a likely source of overhead, but no allocation profile was measured. The smaller cases stayed on the unchanged path; their paired full-call medians were all above their own A/A minima.

Both primary paired medians miss the frozen 1.04× local threshold, and their bootstrap lower bounds miss 1.02×. The holdout primary also falls below its own A/A minimum. No x86-64 promotion or 2× claim follows. The [compressed raw receipt](screen_2026-09-30T051212Z.json.gz) preserves every process output, status, phase timing, source and binary hash, host and route; [SCREEN_MANIFEST.json](SCREEN_MANIFEST.json) records its uncompressed and compressed SHA-256 values and byte counts. The tested source head was `cd8611d79843e5890bf6610a7ae9f3b9f766e0c6`.
