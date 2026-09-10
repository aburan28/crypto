# RTX PRO 6000 benchmark and audit preset

Run from the `ecc2k130` directory with Modal installed and authenticated:

```bash
make bench-rtx-pro6000
make audit-rtx-pro6000
```

The first target runs three complete walk benchmarks. The second performs
GPU arithmetic and integration checks, three benchmarks, and three DP
collection runs; it writes `build/rtx-pro6000-audit.json`. Both allocate one
RTX PRO 6000 Blackwell Server Edition through Modal. Set `MODAL=/path/to/modal`
when using a specific client environment.

The preset selects CUDA 13.0.0, the packed backend, batch 32, 128 threads per
block, minBlocks 4, 1024 steps per launch and 32 launches. It enables the single
polynomial product, denominator cache, by-value operands, both Frobenius
networks, polynomial chains, explicit inversion schedule and paired products.
The preset also enables [polynomial coordinate storage](POLYNOMIAL-STATE.md).
These settings retain the existing iteration, DP report and packed-checkpoint
semantics. The linked comparison validates normal-to-polynomial resume and
the reverse direction, and measures both modes on one GPU.

CUDA 13.0 requires a compatible driver; the tested driver is 580.95.05.
NVIDIA lists 580.65.06 for the Linux CUDA 13.0 GA toolkit in its
[release notes](https://docs.nvidia.com/cuda/archive/13.0.0/cuda-toolkit-release-notes/index.html#cuda-driver).

## Historical standalone CUDA 13 audit with normal storage

The [standalone audit](benchmarks/cuda13/native-audit.json) used the normal
`packed_audit.py` launch path with the earlier normal-storage preset parameters,
on the CUDA 13.0.0 image:

| Workload | Median B scalar iterations/s | Range across three repetitions |
|---|---:|---:|
| Complete walk benchmark | **6.317452** | 6.301125–6.317701 |
| DP cutoff 34 collection | **6.185531** | 6.185358–6.193231 |

This audit used a different GPU allocation from the compiler comparison
below. Its absolute rates should not be used to estimate the compiler's
isolated gain. Both sets of results retain their full environment identity.

## Compiler comparison on one GPU

The [complete comparison](benchmarks/cuda13/compiler-comparison.json) compiled
the same source with nvcc 12.8.93 and 13.0.48, then ran both binaries on the same
RTX PRO 6000 allocation. Both used native sm120 code and identical build knobs.
Three repetitions alternated CUDA 12.8 then CUDA 13.0 for each workload.

| Workload | CUDA 12.8 median B iterations/s | CUDA 13.0 median B iterations/s | Change |
|---|---:|---:|---:|
| Complete walk benchmark | 6.205207 | **6.268586** | +1.02% |
| DP cutoff 34 collection | 6.075874 | **6.159745** | +1.38% |

Benchmark ranges were 6.171362–6.209920 B/s for CUDA 12.8 and
6.245650–6.309846 B/s for CUDA 13.0. Collection ranges were 6.075587–6.079389
and 6.158061–6.162090 B/s, respectively. Rates drifted downward during the
comparison; the CUDA 13 binary remained faster in every paired repetition.
The compiler gain in this comparison was modest.

Every timed run completed 100,931,731,456 scalar walk iterations. Each
collection repetition wrote 2,633 records, dropped zero, and matched its
corpus-file size. Walk synchronization, prior-launch restarts and host report processing were
timed; initial setup was excluded. These historical binaries did not wait for
any remaining final asynchronous reseed before stopping the timer. Current
binaries include that final wait; see the polynomial-state comparison for
measurements of the corrected timing boundary. CPU trail replay is disabled during timing
after separate correctness checks pass.

Both binaries passed 3,120 GPU Frobenius vectors, 1,261 reductions, 18,185
individual polynomial products and 18,185 paired cases checking both outputs.
Integration checks cover report replay across launches and restarts, exact
resume, scalar accounting, preservation of incompatible checkpoints and the
overdue-walk guard. Checkpoints from both compilers match byte-for-byte at
common worker and step counts.

The measured source is commit `4bd7f5d1142896f53be7ed4fe73249dbe3fa4f25`,
with source digest `69912aa7bd5afcef721bc97030d5741515c5eb40cb4820d7e59210a63587d95b`.
The original preset added command aliases around those build and runtime
options; the current preset also selects polynomial coordinate storage. Compiler versions, binary hashes, GPU identity and complete raw
output are included in the comparison artifact.

These measurements count complete scalar walk iterations. They do not
establish 60 B iterations/s, performance on other GPUs, or live clock traces;
the GPU-state metadata is a snapshot taken before validation.
