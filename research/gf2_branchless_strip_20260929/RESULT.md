# Branchless AVX2 strip clearing: opt-in gate passed

The frozen [protocol](PROTOCOL.md) completed its opt-in comparison in [CI run 36626590158](https://github.com/aburan28/crypto/actions/runs/36626590158) from experimental head `d4ad64b901d797856f05fdb1eb9f97938b29530b` (Actions checked out synthetic merge SHA `9ac23bc378823c754be16d3b888d0be6c905c090`). All 88 processes succeeded: four seeds, one warmup per mode, five reference/reference A/A pairs and five alternating reference/candidate pairs per seed. All seven F5 cases matched raw and canonical row fingerprints, rank, output terms, criterion/build counts and reduction word operations. The shared GF(2) kernel and F5 correctness tests passed before timing.

The host was an AMD EPYC 7763, Linux x86-64, Rust 1.98.1, AVX2/BMI2, with one pinned CPU and one Rayon thread. The [full compressed receipt](runs/36626590158-t1.json.gz) has SHA-256 `504bc6e731912fbaedab057d4ba036d8da56a8b75198e578d3b444e5a96bbeb1`; the uncompressed JSON has SHA-256 `decade2e9fd4f093f004e4368340e6c83557dc0ea62519b18cdaf3a21537e07d` and 807,986 bytes. It preserves each process and status, host load, affinity, source/binary hashes, outputs and timings. No failure, timeout or OOM was dropped. The measured source hash was `521bc0dc51e0c41b7da241d186bfabb5897dc4c951da95a4e41fcc011277abd2` for `gf2_elim.rs`.

| Seed | Paired complete-call reference/candidate median (95% bootstrap interval) | A/A range | Paired elimination ratio |
| --- | ---: | ---: | ---: |
| `0` | **1.159× (1.143–1.164×)** | 0.990–1.003× | 1.375× |
| `badc0de1` | 1.145× (1.131–1.161×) | 0.997–1.007× | 1.357× |
| `5eed2026` | 1.168× (1.159–1.170×) | 0.995–1.013× | 1.384× |
| `f5c02a28` | 1.160× (1.152–1.170×) | 0.987–1.001× | 1.377× |

On the frozen primary n24 degree-4 case, marginal median complete-call times were 129.539 ms reference and 111.841 ms candidate; elimination was 65.098 versus 47.403 ms. Build stayed 7.833 versus 7.806 ms and unpack 55.730 versus 55.734 ms. These marginal medians are separate from the paired ratios. Every smaller-case complete-call median improved on every seed; the minimum ratio was 1.065×. The preregistered incremental gate passes, but the requested further 2× gate does not. The exact [measured candidate patch](measured_candidate.patch) and [opt-in workflow](WORKFLOW_T1.yml) are preserved.

The remaining work before enabling by default is the separately frozen one-thread and four-thread default-on confirmation. This is a matrix-F5 solver-stage gain. It establishes no IC one-target online DLP or rho speedup.
