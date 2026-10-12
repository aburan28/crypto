# Direct polynomial-basis sigma map: terminal rejection

Decision: **do not promote the direct table.**  The exact table reduced static
ALU counts, passed every arithmetic and walk-equivalence gate, and then lost
every matched RTX PRO 6000 timing pair by about 17.3%.  Current main remains on
the composed shared-sigma network.

This directory preserves terminal evidence only.  The measured native
generator, default-off CUDA implementation, frozen protocol, test kernel and
failed-attempt receipts remain on branch
`codex/ecc2k130-direct-sigma-gpu-20261001`; the measured source commit is
`af2151c5555757e5a605cdcec03f2ce14602bbde`.  The rejected implementation is
not added to the current-main preset by this evidence PR.

## Construction and static screen

For each `j=3..10`, the native C++ generator derived the exact rank-130 map

```
toPolynomial((I + sigma^j) fromPolynomial(p))
```

from repository arithmetic.  It checked all 1,048 map/basis-vector cases and
32,768 deterministic dense host cases.  A three-input-bit table occupied
56,320 bytes for all eight maps.  Fixed-`j` Apple-Clang instruction medians
fell from 498 for the reduced composed path to 347.5 for the table, a 30.2%
static reduction.  Diagonal and half-five-bit alternatives failed their static
screens.

The GPU implementation copied the table once per block into 56,320 bytes of
opt-in dynamic shared memory.  It disabled the 1,792-byte shared sigma-mask
array, so those allocations were not combined.

## Matched GPU result

One RTX PRO 6000 Blackwell Server Edition, CUDA 13.3.73, native `sm_120`,
B16/T512/min1, 385,024 workers, 1,024 steps and 64 launches ran each timing
row.  Every row completed exactly 403,726,925,824 scalar updates with zero
drops.

| arm | median B/s | registers | local bytes | static shared | dynamic shared |
|---|---:|---:|---:|---:|---:|
| composed shared-sigma control | **15.893066** | 126 | 0 | 1,792 | 0 |
| direct three-bit table | **13.143106** | 128 | 0 | 0 | 56,320 |

Candidate/control ratios were 0.826654, 0.826692, 0.826986, 0.827032 and
0.826916.  Median was **0.826916**, a 17.31% loss; every pair lost.  A/A
maximum absolute drift was 0.0371%.

Correctness gates passed before timing:

- 1,048 basis and 4,096 dense direct-map GPU cases;
- packed arithmetic and compact-storage CUDA suites;
- 300/300 replay per arm with zero drops;
- 1,710,327 complete v1 records per arm with identical sorted payloads;
- byte-identical uninterrupted and both cross-arm checkpoints; and
- zero stack, spills and local bytes for both packed walk kernels.

The static ALU cut did not cover the random shared-memory gathers and bank
conflicts.  The result also weakens the proposed sparse-pentanomial compound
route: adding a dense X-to-normal table would raise dynamic shared allocation
to 84,480 bytes and hot shared lookup traffic from 176 to 220 instructions per
update.  That route now needs an end-to-end native/SASS model above the
measured 1.2092x break-even before another GPU run.

The raw archive is retrievable with Modal volume token
`99bd6b3eacaf470ba73e1d1b11777f62`; its SHA-256 is
`2734894090f9885640aea546771ce6d5d042f5abb501bc985d1d7983880b3526`.
`results/artifact.json` indexes hashes and provenance, and
`results/independent-audit.json` records the reopened terminal audit.  This is
benchmark-only evidence; no search, solver, collision recovery or key recovery
ran.  The 26 B/s objective remains unmet.
