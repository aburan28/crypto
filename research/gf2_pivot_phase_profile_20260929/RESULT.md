# Matrix-F5 pivot loop: strip clearing is the dominant subphase

The frozen [protocol](PROTOCOL.md) completed in [CI run 36623102489](https://github.com/aburan28/crypto/actions/runs/36623102489) from experimental head `35f789bb1b6404813d4dd3c2a5a3f337f129cf7f` (the Actions checkout reported synthetic merge SHA `271f364987cc706ab5a42186fb9cf439d70d85d0`). All 88 processes succeeded: four seeds, one warmup per mode, five broad/broad A/A pairs and five alternating broad/fine pairs per seed. All seven F5 cases matched raw and canonical row fingerprints, rank, output terms, criterion/build counts, and reduction word operations. The shared-kernel and F5 correctness tests passed before timing.

The runner was an AMD EPYC 7763, Linux x86-64, Rust 1.98.1, with AVX2/BMI2, one pinned CPU and one Rayon thread. The [complete compressed receipt](runs/36623102489-t1.json.gz) has SHA-256 `c49bb401621193ceceeb7d5b89ac8754d38a1bd67464696814cd0ead27fdda74`; the uncompressed JSON has SHA-256 `eebcd7b7ee0b45ee5aa9cf5af829c6774a4f5de64c42994629053c2efc9e5823` and 1,884,222 bytes. It retains every process, status, host load, affinity, source/binary hash, exact output and phase timer. No failure, timeout or OOM row was dropped.

## One-thread n24 degree-4 pivot breakdown

Each entry is the median of five fine-profiled paired calls, in milliseconds. The four fine phases are exclusive within the broad pivot bucket. Medians of separate columns need not add exactly.

| Seed | Pivot total | Strip clear | Strip scan | Selected-row XOR | Earlier-pivots XOR | Timer/other residual |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `0` | 28.983 | **18.936** | 4.100 | 2.446 | 2.511 | 1.034 |
| `badc0de1` | 29.066 | **18.952** | 4.128 | 2.487 | 2.514 | 1.029 |
| `5eed2026` | 29.029 | **18.895** | 4.125 | 2.498 | 2.533 | 1.016 |
| `f5c02a28` | 29.089 | **18.945** | 4.123 | 2.483 | 2.510 | 1.024 |

On frozen `0`, strip clearing is about 65% of profiled pivot work and 31% of the 61.084 ms profiled elimination. The matrix has 6,924 rows, 12,951 columns, rank 6,924 and 248 pivot blocks. Its other broad phases are nonpivot row clearing 21.789 ms, Gray-code table construction 8.135 ms and strip loading 2.143 ms. The phase ordering repeats on every holdout.

The fine timers perturb the kernel: frozen complete-call marginal medians were 119.645 ms broad-only and 122.069 ms broad-plus-fine; reduction medians were 59.057 and 61.714 ms. The paired broad/fine complete-call ratio was 0.990× (A/A range 0.996–1.016), and the reduction ratio 0.978× (A/A range 0.985–1.020). Across holdouts, paired broad/fine complete-call ratios were 0.989–0.991× and reduction ratios 0.973–0.980×. The measured ordering is stable despite this overhead; these are diagnostic costs of the instrumented kernel, not a speed claim for the original.

The next implementation target is the strip clear loop: for each pivot it conditionally XORs one 64-bit strip value into every unpivoted row. Try a branchless update that can vectorize, preserving the exact strip, rank, row output and operation counts. Test the full F5 call on frozen and holdout inputs; a faster strip loop alone cannot establish the requested 2× complete-call gain. The exact [instrumentation patch](rejected_instrumentation.patch), [workflow](WORKFLOW.yml) and [paired script](paired_f5.py) are preserved for replay. Runtime profiling was removed from production after measurement. This is a solver-stage diagnostic, not an IC online-time or DLP speedup.
