# Direct polynomial-basis sigma map: static result

Decision: **retain the table3 circuit as a GPU-screen candidate; reject the
diagonal and half5 circuits.**  No GPU was used for this result.

The native generator derives every matrix from the repository arithmetic and
proves every generated family on all 131 basis inputs for each `j=3..10`.
Every `L_j` has rank 130, as expected from `I + sigma^j`.  The separately
compiled generated headers also passed 32,768 deterministic dense vectors.

## Static comparison

Apple Clang 17 compiled fixed-`j`, fully inlined wrappers at `-O3`.  The count
is AArch64 host assembly, not PTX or SASS.  Medians are over all eight powers.

| implementation | all-map data | loads / coordinate | median instructions | median text bytes | vs reduced compose |
|---|---:|---:|---:|---:|---:|
| composed, shipping 9-word inverse | existing | 56 sigma-mask loads, shared across the X/Y pair | 513 | 2,044 | |
| composed, reduced-input inverse | existing | 56 sigma-mask loads, shared across the X/Y pair | 498 | 1,984 | reference |
| **direct table3** | **56,320 B** | **220 table words** | **347.5** | **1,382** | **-30.2% instructions, -30.3% text** |
| direct half5 | 56,848 B | 135 table words plus up to 130 conditional XORs | 583 | 2,276 | +17.1% instructions |

The table3 counts by `j=3..10` are 260, 260, 317, 348, 348, 347, 348 and
348.  Reduced-compose counts are 459, 480, 508, 503, 512, 488, 503 and 493.
Shipping-compose counts are 474, 496, 524, 519, 528, 503, 518 and 508.

The deterministic word-data model tells the same ALU story more strongly:
table3 uses 317 shifts/masks/XORs per coordinate, while the existing
reduced-input composition uses 615.  In the actual X/Y forward pair the X
normal-basis conversion remains necessary for Hamming weight.  Applying the
direct map to both coordinates therefore models 742 data-ALU operations
against 1,230, a cut of 488 (39.7%).

## Why this is only a candidate

The direct matrices are dense: 1,724 to 4,383 one-bits.  A diagonal evaluator
needs 228 to 247 of the 261 possible diagonals and 537 to 633 nonzero word
terms, so it is not competitive.

Table3 trades ALU for lookup traffic.  A naive constant-memory implementation
is not credible: the joint `(j, three input bits)` key has 64 values, so a
warp with unrelated walks will serialize many constant-cache addresses.  A
per-block shared layout avoids that serialization and can load each entry as
one aligned 128-bit vector plus one top word.  It would need 88 shared-memory
instructions per coordinate, 176 per X/Y pair, against the current 56 scalar
shared-mask loads per pair.  The 56,320-byte table plus the existing small
allocation permits one block per SM.  A 512-thread block retains 16 resident
warps, matching two current 256-thread blocks.

That shared layout is the only admitted GPU follow-up: batch 16, 512 threads,
one block/SM, identical worker and update counts, matched control geometry,
host replay, corpus identity, A/A noise, then alternating A/B.  The experiment
must be frozen separately before launch.  The static result is not throughput
evidence and does not establish the 26 B/s objective.

Generated artifact hashes from the frozen run:

- table3 header: `868eb934642b118997b5d43e79ac341a4ce001dc987052a47d2b9a594f99b9cd`
- half5 header: `84cf49c4b1df8c18159abd42287b2ce34bc55b48b00c0e33beca021bb78501b1`
- synthesis transcript: `15b690d4181e4d60d6b3f14b4059e61a3fe4a31f588c7d9f8f1effdc0cd9a2f7`

Run `benchmarks/direct-sigma-map/run.sh` from any directory to regenerate and
recheck the native evidence without Python.
