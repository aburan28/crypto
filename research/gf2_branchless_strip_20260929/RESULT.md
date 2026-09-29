# Branchless AVX2 strip clearing: default-on controls passed

The frozen [protocol](PROTOCOL.md) completed its opt-in comparison in [CI run 36626590158](https://github.com/aburan28/crypto/actions/runs/36626590158) from experimental head `d4ad64b901d797856f05fdb1eb9f97938b29530b` (Actions checked out synthetic merge SHA `9ac23bc378823c754be16d3b888d0be6c905c090`). All 88 processes succeeded: four seeds, one warmup per mode, five reference/reference A/A pairs and five alternating reference/candidate pairs per seed. All seven F5 cases matched raw and canonical row fingerprints, rank, output terms, criterion/build counts and reduction word operations. The shared GF(2) kernel and F5 correctness tests passed before timing.

The host was an AMD EPYC 7763, Linux x86-64, Rust 1.98.1, AVX2/BMI2, with one pinned CPU and one Rayon thread. The [full compressed receipt](runs/36626590158-t1.json.gz) has SHA-256 `504bc6e731912fbaedab057d4ba036d8da56a8b75198e578d3b444e5a96bbeb1`; the uncompressed JSON has SHA-256 `decade2e9fd4f093f004e4368340e6c83557dc0ea62519b18cdaf3a21537e07d` and 807,986 bytes. It preserves each process and status, host load, affinity, source/binary hashes, outputs and timings. No failure, timeout or OOM was dropped. The measured source hash was `521bc0dc51e0c41b7da241d186bfabb5897dc4c951da95a4e41fcc011277abd2` for `gf2_elim.rs`.

| Seed | Paired complete-call reference/candidate median (95% bootstrap interval) | A/A range | Paired elimination ratio |
| --- | ---: | ---: | ---: |
| `0` | **1.159× (1.143–1.164×)** | 0.990–1.003× | 1.375× |
| `badc0de1` | 1.145× (1.131–1.161×) | 0.997–1.007× | 1.357× |
| `5eed2026` | 1.168× (1.159–1.170×) | 0.995–1.013× | 1.384× |
| `f5c02a28` | 1.160× (1.152–1.170×) | 0.987–1.001× | 1.377× |

On the frozen primary n24 degree-4 case, marginal median complete-call times were 129.539 ms reference and 111.841 ms candidate; elimination was 65.098 versus 47.403 ms. Build stayed 7.833 versus 7.806 ms and unpack 55.730 versus 55.734 ms. These marginal medians are separate from the paired ratios. Every smaller-case complete-call median improved on every seed; the minimum ratio was 1.065×. The preregistered incremental gate passes, but the requested further 2× gate does not. The exact [measured candidate patch](measured_candidate.patch) and [opt-in workflow](WORKFLOW_T1.yml) are preserved.

## Default-on confirmation

The separately frozen [CI run 36628083311](https://github.com/aburan28/crypto/actions/runs/36628083311) passed both one-thread and four-thread jobs at source head `c906854b5452e29158a314c9c576223b7953cfe8` (Actions checkout synthetic merge SHA `3b12dfeaaec6239511a3a1cd340cab0148e6b3f2`). Mode 0 explicitly sets `KIC_GF2_BRANCHLESS_STRIP=0`; mode 1 leaves it unset and exercises the new AVX2 default. Each job has 88/88 successful calls. Every case matches raw/canonical fingerprints, rank, terms, counts and word operations; the one- and four-thread output signatures also match each other. Both jobs passed the release-mode shared-kernel and F5 correctness tests.

The one-thread job used one pinned CPU of an AMD EPYC 7763. The four-thread job used four pinned CPUs of a distinct AMD EPYC 9V45. Ratios below are paired **within each job**; no absolute time or ratio is carried between the two hosts. Values above one favor the default-on path.

| Seed | One-thread complete call (95% bootstrap interval) | One-thread A/A range | Four-thread complete call (95% bootstrap interval) | Four-thread A/A range |
| --- | ---: | ---: | ---: | ---: |
| `0` | **1.252× (1.229–1.262×)** | 0.992–1.002× | 1.336× (1.283–1.394×) | 0.930–1.026× |
| `badc0de1` | 1.230× (1.219–1.258×) | 0.983–1.012× | 1.348× (1.346–1.447×) | 0.995–1.086× |
| `5eed2026` | 1.239× (1.233–1.244×) | 0.989–1.006× | 1.369× (1.335–1.492×) | 0.937–0.992× |
| `f5c02a28` | 1.237× (1.226–1.240×) | 0.992–1.007× | 1.330× (1.287–1.403×) | 0.975–1.029× |

On the one-thread frozen primary, marginal complete-call medians were 132.194 ms with the explicit old path and 105.730 ms with the AVX2 default; elimination was 68.409 and 42.909 ms. Unpacking remained 54.617 and 54.580 ms. All six smaller cases improved on every seed: the smallest one-thread complete-call ratio was 1.100×, and the smallest four-thread ratio was 1.131×. The one-thread gain clears the preregistered 1.05× median and 1.02× lower-bound gates; the four-thread nonregression gate also clears. No twofold complete-call gain is established.

Both full receipts are archived: [one thread](runs/36628083311-t1.json.gz), compressed SHA-256 `5faa47d454c920ea78188c998c0142ed240f10e7328da52d7a4388734f45590b`, uncompressed SHA-256 `c36993d1a530de567b43a3c9ae03dc88079aed7bf74106033756aedebb90bdf8` (807,098 bytes); and [four threads](runs/36628083311-t4.json.gz), compressed SHA-256 `a439b9d4e7fd8cf5ce122b00358782bd5a9a1dfbff026b7f2afde2c87dc42df4`, uncompressed SHA-256 `29eb59024fd6c351b7f57461745ae3d86362f63c793349c3c310c9c18e578226` (806,967 bytes). They retain every call, status, host, source/binary hash, workload, exact output and timing. The [confirmation workflow](WORKFLOW_CONFIRM.yml) is preserved and removed from the active CI set after the run. The default is enabled only where AVX2 is available; `KIC_GF2_BRANCHLESS_STRIP=0` is the rollback, while other CPUs keep the old default loop and the scalar branchless path remains opt-in.

This is a matrix-F5 solver-stage gain. It establishes no IC one-target online DLP or rho speedup. The requested further 2× complete-call goal remains open, with unpacking still about 55 ms of the 106 ms measured default-on call.
