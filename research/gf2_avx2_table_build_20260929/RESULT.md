# AVX2 Gray-code table construction: reduction gain, full-call gate missed

The one-thread x86 paired experiment in [PROTOCOL.md](PROTOCOL.md)
completed in [CI run 36611401578](https://github.com/aburan28/crypto/actions/runs/36611401578)
at measured head `7f9db38fca6e4f0ddcf1c0a15840061cbced77c1`. All 88
processes succeeded. The benchmark verified that the AVX2 table builder
was selected only in the new primary arm, and every case on every seed
matched raw and canonical row fingerprints, rank, output term count,
criterion/build counts, and reduction word XORs. The x86 exactness tests,
release build and benchmark job passed.

The runner reported AMD EPYC 7763 with AVX2/BMI2, Linux x86-64,
Rust 1.98.1, one pinned CPU and one Rayon thread. The
[compressed complete receipt](runs/36611401578-t1.json.gz) has SHA-256
`876270ff25f64ee459f4c6204cf7afdf9de60b9f0b0938c97add7ca38532c7f2`;
the uncompressed JSON SHA-256 is
`1909ecdb8fe08e473e4a58fb3fa37927cad1ed5cfb97768bc172f549417f5681`.
It retains every call, host load, affinity, source/binary hashes, exact
outputs and phase timings. No failure, timeout or OOM row was removed.

## Frozen and holdout complete-call comparison

Each ratio is the median of five paired prior/new ratios; the intervals
are exact 3,125-resample bootstrap 95% intervals over those five pairs.
`wall_ms` covers criterion through full row unpack. The A/A range is
five prior/prior pairs from the same workload. Values above one favor
AVX2 table construction.

| Seed | Full call prior/new (95% interval) | A/A range | Reduction prior/new (95% interval) |
| --- | ---: | ---: | ---: |
| Frozen `0` | **1.012×** (0.996–1.044) | 0.982–1.001 | 1.036× (1.014–1.045) |
| Holdout `badc0de1` | 1.013× (0.991–1.016) | 0.992–1.006 | 1.021× (0.981–1.038) |
| Holdout `5eed2026` | 1.018× (1.006–1.019) | 0.988–1.007 | 1.033× (1.026–1.040) |
| Holdout `f5c02a28` | 1.013× (1.003–1.020) | 0.997–1.103 | 1.031× (1.012–1.037) |

On the frozen primary, marginal medians were 125.146 ms prior and
123.051 ms new for the complete call, 59.842 and 57.568 ms for
reduction, and 56.587 and 57.024 ms for unpack. The reduction-phase
improvement does not reach the requested further 2× complete-call
target. The frozen complete-call median and its lower interval bound
both miss the preregistered 1.03 incremental gate. One smaller case,
`f5_n20_m20_d3` on holdout `f5c02a28`, had a 0.984× median against
its 0.997 A/A minimum. The tested AVX2 build option therefore is not
promoted or retained in production.

The exact implementation and tests are preserved in
[rejected_candidate.patch](rejected_candidate.patch), and the run
configuration in [WORKFLOW.yml](WORKFLOW.yml) and
[paired_f5.py](paired_f5.py). Apply the patch to the frozen source to
replay. This is a matrix-F5 solver-stage measurement, not an IC
online-time or DLP speedup. The further 2× complete-call goal remains
open.
