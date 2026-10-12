# Matrix-F5 elimination: pivot work is the largest measured phase

The opt-in phase profile in [PROTOCOL.md](PROTOCOL.md) completed in
[CI run 36621118791](https://github.com/aburan28/crypto/actions/runs/36621118791)
from PR head `68f96c147565446884776cf320ce90fa39318c33` (the Actions
checkout reported synthetic merge SHA
`12b04c6d1f66dfdd33a3952742586df42f51896d`). All 88 processes
succeeded: four seeds, one warmup per mode, five profile-off/off A/A
pairs and five alternating off/on pairs per seed. Every F5 case
matched raw and canonical row fingerprints, rank, output terms,
criterion/build counts, and reduction word operations. The x86
release-mode shared-kernel and F5 tests and benchmark job passed.

The runner was an AMD EPYC 7763, Linux x86-64, Rust 1.98.1, with
AVX2/BMI2, one pinned CPU and one Rayon thread. The
[complete compressed receipt](runs/36621118791-t1.json.gz) has
SHA-256 `88174073bab96a1b16033e1e098bc55ca5d8a83fc1ebd4f5db01b576f0f32da1`;
the uncompressed JSON has SHA-256
`13ff10fd30af0fcc8bfbc3882d586b7761044c5871fc035523b8a83a9baaf75c`
and 1,001,325 bytes. It retains every process, failure status, host
load, affinity, source/binary hashes, exact output and phase timer.
No failure, timeout or OOM row was dropped.

## One-thread n24 degree-4 elimination

Each entry is the median of five profile-on paired calls, in
milliseconds. The totals are inside the shared eliminator and
exclude F5 row build and polynomial unpack. Each call's exclusive
phases sum to its profiled total up to the small residual; medians
of separate columns need not sum exactly.

| Seed | Total | Strip load | Pivot work | Table build | Row clear | Residual |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `0` | 59.348 | 2.037 | **27.382** | 8.005 | 21.959 | 0.043 |
| `badc0de1` | 58.472 | 1.851 | **27.190** | 7.969 | 21.535 | 0.045 |
| `5eed2026` | 59.341 | 1.960 | **27.339** | 8.101 | 21.940 | 0.043 |
| `f5c02a28` | 59.357 | 1.920 | **27.235** | 8.045 | 22.008 | 0.046 |

On frozen `0`, pivot work is about **46%** of elimination, clearing
nonpivot rows 37%, Gray-code table construction 13.5%, and strip
loading 3.4%. There were 248 pivot blocks and rank 6,924. The
phase ordering repeats on every holdout. The above-pivot reverse
phase is zero for selective-echelon output. The current `pivot work`
bucket still includes searching the strip, reducing the selected
pivot row by earlier pivots, reducing earlier pivots by it, and
clearing the pivot bit from the strip; it needs a second split before
choosing an implementation change.

Instrumentation was not a speed candidate. On frozen `0`, complete
F5-call marginal medians were 121.954 ms profile-off and 122.868 ms
profile-on; reduction medians were 59.225 and 60.013 ms. The
complete-call off/on paired median was 0.993× (A/A 0.990–1.010),
and the reduction off/on median 0.984× (A/A 0.989–1.009).
Thus the timer overhead is visible, but much smaller than the
roughly 5 ms gap between pivot work and row clear. The profiled
phase totals are diagnostic estimates, not timings of the
uninstrumented kernel. The 2× complete-call goal remains open.

The exact [instrumentation patch](rejected_instrumentation.patch),
[workflow](WORKFLOW.yml), and [paired script](paired_f5.py) are
preserved for replay. The runtime profiling option was removed from
production after measurement. This is a matrix-F5 solver-stage
diagnostic, not an IC online-time or DLP speedup.
