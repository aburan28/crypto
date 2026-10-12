# GPU-wide exact-v3 queue: independent post-run audit

Audited source: `7e5281e6b276b89428dbb74a1c05a50ca5b61a02`.

Decision: **PASS evidence; do not promote the GPU-wide queue.** The candidate
preserves the exact table-v3 walk but reduces complete-update throughput from a
7.021655 B/s control median to 5.610974 B/s. The paired geometric mean is
0.7990182843 and all five candidate/control pairs are approximately 0.799.
The active 26 B/s goal remains unmet.

The audit is native C++17 and launches no GPU. It reopens the sealed Modal
archive rather than relying on the producer's comparison receipts.

## Correctness and custody

- Clean source commit, one RTX PRO 6000, native sm_120 and CUDA 13.3.73 are
  bound by the launch, host and source manifests.
- All 24 source-manifest entries, seven retained executables and 79 producer
  artifact-manifest entries verify against their bytes.
- A fresh path/type-safe extraction of the 81,726,971-byte archive is
  byte-identical to the downloaded result tree.
- The native ownership model passed 149,940 cases and the production CUDA
  selector/resolver control passed 49 empty/full/sparse/reset/partial-block
  cases.
- Full and partial production runs replayed 300/300 reports per arm with zero
  drops. The deterministic spread replay matched 300/300 records per arm and
  exercised 299 nonzero trails per arm.

The independent audit scans every v3 record, validates the run-id namespace and
canonical 131-bit key, sorts the complete payload and compares both arms:

| corpus | records per arm | duplicate records | sorted payload SHA-256 |
|---|---:|---:|---|
| full, 96,256 workers | 1,710,384 | 0 | `5cd48f49da20a092a2883cf9c5099ee774712a8bc225d2c91214c72afe79736e` |
| partial, 511 workers | 9,148 | 0 | `7d53f9666ac1809ce24e48494146beb4b83a5c88b2ee92bf8abe98dc52719a58` |
| partial, 513 workers | 9,182 | 0 | `f87932ebe57bfc24a96e865c415f6511f1c6479d07eb20ac7b5bc8c63df52a10` |
| four-launch prefix | 5,242 | 0 | `b92500750dbd87840dcff01570c5d88c6bdc43a366a2683240fe85633766a198` |
| each three-launch continuation | 3,940 | 0 | `76640717c663ab1a1dd23f9963b523f2ec310a302899aa96835dc9cbe7b7cc6f` |

Both prefix checkpoints are byte-identical at iteration 380. All four
control/candidate continuation combinations are byte-identical at iteration
665. Every file has the exact 558,184-byte packed table-v3 checkpoint framing.

The preserved r1 dispatch failure occurred before function spawn, created no
launch receipt and ran zero GPU kernels. It is not a measurement.

## Runtime identity and resources

Every verification, prefix, continuation and timing log has the exact frozen
feature vector, worker count, DP weight, step count, L2 window, launch bounds
and resource tuple. The audit also verifies every timing log against the SHA-256
recorded in `samples.tsv` and requires exactly 201,863,462,912 complete scalar
updates with zero drops.

| arm/stage | registers/thread | local bytes/thread | static shared bytes |
|---|---:|---:|---:|
| block-queue control hot kernel | 128 | 400 | 1,040 |
| GPU-wide candidate hot kernel | 126 | 0 | 0 |
| GPU-wide select kernel | 56 | 0 | 0 |
| GPU-wide resolver kernel | 124 | 400 | 0 |

The candidate removes the block-wide barriers and the hot kernel's local/static
footprint, but charges a counter reset plus select, resolve and hot kernel
launch for every scalar update. That complete scheduling cost dominates the
saved reconvergence work.

## Timing

Five A/A pairs establish maximum symmetric drift of 0.056254%, below the 1%
gate. The five paired candidate/control ratios are:

```
0.7989797513
0.7988113913
0.7989729924
0.7991527520
0.7991745898
```

The native producer result and the independent recomputation agree exactly:

| metric | result |
|---|---:|
| control median | 7.021655 B/s |
| candidate median | 5.610974 B/s |
| paired geometric mean | 0.7990182843 |
| engineering gate | false |
| 26 B/s gate | false |
| terminal decision | `DO_NOT_PROMOTE` |

The result rejects the ordinary per-step three-kernel schedule. A CUDA graph
or another scheduling composition would be a separate experiment; this result
does not predict its rate.

## Immutable artifact

The raw archive has SHA-256
`f677fa729e662a7718ed0beb1d7880473697bb25588d54f19f93c6f39e5caa6f`.
Retrieval metadata and all manifest/audit hashes are recorded in
`results/independent-audit-artifact.json`. The complete machine-readable audit
is `results/independent-audit.json`.
