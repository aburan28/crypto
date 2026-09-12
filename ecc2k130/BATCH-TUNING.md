# Batch depth with native carryless multiplication

The RTX PRO 6000 preset now uses 16 slots per worker and 385,024 workers.
This preserves the previous 6,160,384 logical walks while reducing the state
handled by each worker. It keeps 256 threads per block, minBlocks 2,
256-worker state tiles and the native carryless arithmetic settings.

Run from this directory:

```sh
make bench-rtx-pro6000
make audit-rtx-pro6000
```

The audit accepts an explicit `--batch` and checks the compiled batch,
requested workers, logical walk population and completed scalar count for
every benchmark and collection sample. Its standalone default remains 32;
the Make preset explicitly requests 16. Existing checkpoints require their
original worker and batch counts. These Make targets start fresh workloads.

## Controlled complete-walk results

The [retained comparison](benchmarks/batch-tuning/comparison.json) ran on
one NVIDIA RTX PRO 6000 Blackwell Server Edition using CUDA 13.3.73, native
`sm_120` code and disabled PTX JIT. Every timed row completed
**201,863,462,912 scalar updates**: 6,160,384 walks, 1,024 steps per launch,
32 launches. The final timer synchronizes pending reseeds.

| Workload | Batch 32 median B/s | Batch 16 median B/s | Ratio-of-medians gain |
|---|---:|---:|---:|
| Complete scalar benchmark | 8.673447 | **13.206088** | **52.2588%** |
| DP34 collection | 8.508196 | **12.936060** | **52.0423%** |

Three alternating pairs per workload favored batch 16. Benchmark ranges
were 8.673304–8.675858 B/s for the control and 13.191508–13.208166 B/s for
the candidate. Collection ranges were 8.507229–8.509269 and
12.915063–12.948114 B/s, respectively. Five warmups and the initial screening
rows are excluded from these medians.

All six collections produced 5,149 records, 164,768 bytes and zero drops.
Their sorted record multisets share SHA-256
`ab237b6352380547fd37fcdc9e2aa83f6ae1a718b842590b9a19e4224e5d336b`.
The retained artifact records those hashes and sizes; transient corpus and
checkpoint payloads were not retained.

The benchmark median is still below the 15 B/s objective. These paired
measurements establish the selected configuration's gain on the measured
workload.

## Public-command audit

The separate [public audit](benchmarks/batch-tuning/native-audit.json) ran
`make audit-rtx-pro6000` from commit
`8ff25658aef2c213624ce14b8838d8a8a40f7be8` on a fresh allocation:

| Workload | Median B scalar updates/s | Range across three repetitions |
|---|---:|---:|
| Complete walk benchmark | **13.283756** | 13.151469–13.434216 |
| DP34 collection | **12.813626** | 12.800253–12.815682 |

All six samples completed 201,863,462,912 scalar updates with the requested
batch 16, 385,024 workers and native carryless mode. Each collection recorded
5,149 points and zero drops. GPU arithmetic and full client integration
passed. The [artifact audit](benchmarks/batch-tuning/native-audit-review.json)
recomputes the source identity, batch/worker geometry, completed counts and
summary statistics against the frozen public checkout.

This result validates the published command on a separate allocation. It
retains corpus counts and sizes; content equality is supplied by the paired
comparison. Its binary identity is reported before validation, and its
artifact does not contain a post-run hash or full native-code capture.
The raw public result SHA-256 is
`cf91f9ae23742aaff5bd5505bb3b0fbf1712e0ea336460eb77a0068788555587`.

## Screening and state tradeoff

All modes used the same source, compiler, arithmetic and logical population.
The predefined rule required the winning candidate to exceed the faster of
two bracketing controls by 0.5%, followed by the paired confirmations above.

| Batch | minBlocks | Workers | Registers/thread | Screen B/s |
|---:|---:|---:|---:|---:|
| 32 | 2 | 192,512 | 122 | 8.673328 / 8.677384 |
| 16 | 2 | 385,024 | 122 | **13.228833** |
| 32 | 4 | 192,512 | 64 | 7.833332 |
| 16 | 4 | 385,024 | 64 | 8.748415 |
| 8 | 4 | 770,048 | 64 | 12.575303 |

The nonwinning rows have screening measurements only. All selected kernels
and their field helpers compiled without stack or register spills; runtime
reports showed zero local and shared bytes. The CUDA occupancy API reported
two or four resident blocks per SM as expected. These predictions are not
measurements of achieved occupancy or cache hit rates.

Smaller batches amortize inversion over fewer walks, increasing field
products per complete scalar update from `5 + 5/32 = 5.15625` to
`5 + 5/16 = 5.3125`. At the same two-block occupancy, the source inventory
of hot state falls from about 244.91 MiB to 121.54 MiB. This suggests a cache
capacity explanation, but the experiment did not measure cache residency.
The result demonstrates that less work per worker can outweigh additional
arithmetic for this native implementation.

## Correctness and provenance

Before timing, all five binaries passed GPU arithmetic tests covering
3,120 Frobenius vectors, 2,526 reductions, 18,194 ordinary products,
18,194 paired products and 1,157 squares. All five passed full client
replay, restart and resume checks. Common logical states matched across
geometries at 64 and 16,384 walks.

The checkpoint gate executed 63 children: 42 successful operations and
21 expected geometry rejections. It checked split/resume at 8 and 257
workers, same-batch resume across residency settings, and preserved bytes
when batch or worker headers were incompatible. Equal checkpoint payload
length does not make different geometries compatible. Cross-batch state
equality was checked on the seeded valid workloads; arbitrary zero
denominators can affect different groups under the existing batch inversion.

The [independent artifact audit](benchmarks/batch-tuning/comparison-review.json)
checks all 23 timed rows, exact counts, source and binary hashes, complete
native code, geometry gates, sample order and paired statistics. The
source is the native implementation at commit
`087ac51aa84e8574cadb281fc6f18fe92788c9a8`; all 157 source files are identical
across modes. This change promotes compile and worker settings without
editing the CUDA arithmetic or walk formula.

The raw comparison SHA-256 is
`19f1f247587d3f08766d2efbf79961f1d4e6e5426e7e048c4f6b8826bd634311`.
The audit SHA-256 is
`9a89491145351b62864d018aa36c16fdbb32133845cea30f6d0992926f4deb3f`.
