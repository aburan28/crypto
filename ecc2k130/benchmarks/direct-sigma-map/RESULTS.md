# Direct polynomial-basis sigma map: static and GPU result

Decision: **do not promote the direct table.**  The table3 circuit passed the
static screen and every correctness gate, then lost every matched GPU pair by
about 17.3%.  Retain the generator, implementation and negative evidence with
`DIRECT_SIGMA=0` as the default.  The active 26 B/s objective remains unmet.

The native generator derives every matrix from the repository arithmetic and
proves every generated family on all 131 basis inputs for each `j=3..10`.
Every `L_j` has rank 130, as expected from `I + sigma^j`.  The separately
compiled generated headers also passed 32,768 deterministic dense vectors.

## RTX PRO 6000 measurement

One NVIDIA RTX PRO 6000 Blackwell Server Edition, CUDA 13.3.73 and native
`sm_120` code ran the frozen B16/T512/min1 protocol.  Both arms used 385,024
threads, 16 slots and 64 launches of 1,024 steps: **403,726,925,824 complete
scalar updates per timing row**.  Four warmups were excluded, three A/A pairs
measured session noise, and five alternating A/B pairs were ranked.

| arm | median B/s | registers | local bytes | static shared | dynamic shared | result |
|---|---:|---:|---:|---:|---:|---|
| composed shared-sigma control | **15.893066** | 126 | 0 | 1,792 | 0 | reference |
| direct table3 | **13.143106** | 128 | 0 | 0 | 56,320 | reject |

The candidate/control ratios were 0.826654, 0.826692, 0.826986, 0.827032 and
0.826916.  Median 0.826916 is a **17.31% loss**; every pair lost.  A/A maximum
absolute drift was 0.0371%, so noise is two orders of magnitude smaller than
the effect.  Every timing row completed the exact expected work with zero
drops.

Correctness and identity gates passed before timing:

- the dedicated device test checked 1,048 basis cases and 4,096 dense cases;
- both repository arithmetic/storage suites passed;
- each arm replayed 300/300 reports with zero drops;
- both produced 1,710,327 complete v1 records with identical sorted payloads;
- uninterrupted checkpoints and both cross-arm resume directions were byte
  identical; and
- both packed walk kernels had zero stack, spills and local bytes and resident
  one 512-thread block per SM.

The static ALU cut was real but the predicted lookup risk dominated: two maps
perform 176 random shared-memory load instructions per update, with bank
conflicts, in place of 56 shared mask loads.  This falsifies promotion of this
table circuit on this GPU.

The control's 15.893066 B/s is 2.96% above the separately allocated published
B16/T256 fused headline of 15.436677 B/s.  That is not a matched geometry
comparison and therefore not a promotion result.  The earlier fused geometry
panel tested T512 only at batches 32 and 64.  A separate same-allocation
B16/T256/min2 versus B16/T512/min1 panel is justified; it must keep direct
sigma disabled in both arms.

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

## Why the static candidate failed

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

The measured shared layout used batch 16, 512 threads, one block/SM, identical
worker and update counts, matched control geometry, host replay, corpus
identity, A/A noise and alternating A/B.  Its loss shows why the static result
alone was not throughput evidence.

Generated artifact hashes from the frozen run:

- table3 header: `868eb934642b118997b5d43e79ac341a4ce001dc987052a47d2b9a594f99b9cd`
- half5 header: `84cf49c4b1df8c18159abd42287b2ce34bc55b48b00c0e33beca021bb78501b1`
- synthesis transcript: `15b690d4181e4d60d6b3f14b4059e61a3fe4a31f588c7d9f8f1effdc0cd9a2f7`

Run `benchmarks/direct-sigma-map/run.sh` from any directory to regenerate and
recheck the native evidence without Python.

GPU evidence is indexed by `results/gpu-r6-artifact.json`; the independently
reopened terminal audit is `results/gpu-r6-independent-audit.json`.
