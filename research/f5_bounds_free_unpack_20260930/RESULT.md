# Bounds-free matrix-F5 output unpack: rejected at the local screen

The [protocol](PROTOCOL.md) was committed in draft PR #1037 before candidate code or timing. The candidate validated each packed row's column span once, used unchecked column loads for valid rows, and fell back to the checked decoder for malformed shapes. The [measured candidate patch](measured_candidate.patch) preserves the exact tested implementation and screen. The opt-in runtime path was removed after the local gate failed. This is a matrix-F5 solver-stage result, not an IC one-target DLP or rho speedup.

The Apple ARM64 screen used one release binary, SHA-256 `ce18ab5797246aa8fc85ddc02a44e07aab26a2b0e0a9be3d8b2fbc2df6d60440`, for 44 processes: two seeds, two warmups, five reference/reference A/A pairs and five alternating reference/candidate pairs per seed. Every process returned all seven F5 cases. Across both arms and all cases, raw and canonical fingerprints, rank, output terms, row and column counts, criterion/build counts and reduction word operations matched. The candidate route was selected only on n24 degree-4. Release F5 tests passed 11/11, including sparse, dense, partial-word and long-row fallback checks. The local `tools/isolated_bench.py busy` wrapper is not an exclusive CPU reservation, so its times are exploratory.

| Seed | Complete-call reference/candidate paired median | Exact 95% bootstrap interval | A/A range | Unpack reference/candidate median |
| --- | ---: | ---: | ---: | ---: |
| `0` | **1.004×** | 0.999–1.032× | 0.967–1.087× | 1.006× |
| `badc0de1` | **0.982×** | 0.959–1.178× | 0.986–2.951× | 0.994× |

Separate marginal medians for the frozen seed were 59.45 ms reference and 59.35 ms candidate for the complete n24 call, and 17.90 versus 17.85 ms for unpacking. On the holdout they were 59.82 versus 60.43 ms complete and 17.80 versus 17.90 ms unpacking. The holdout A/A range contains a large timing outlier, so the screen supports no small speed claim. All smaller-case paired full-call medians remained above their own A/A minima.

Both primary medians miss the frozen 1.04× progression gate, and both bootstrap lower bounds miss 1.02×. No x86-64 promotion or 2× claim follows. The [compressed raw receipt](screen_2026-09-30T052850Z.json.gz) preserves every process output, status, phase timing, source/binary hash, host and route. [SCREEN_MANIFEST.json](SCREEN_MANIFEST.json) records its byte counts and SHA-256 values. The tested source head was `8ae66eb764e10e669c1ce7d387bf977a68b23abf`.
