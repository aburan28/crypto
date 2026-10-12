# Matrix-F5 byte-subset unpack: rejected

The [protocol](PROTOCOL.md) was frozen before candidate code or timing. This is a solver-stage experiment with exact `Vec<F2BoolPoly>` output; it is not a single-target index-calculus DLP comparison. The candidate precomputed all 256 ordered subsets of each eight-column group and copied their terms into output rows. The measured implementation is retained as [measured_candidate.patch](measured_candidate.patch), and the [candidate workflow](WORKFLOW_CANDIDATE.yml) is archived. The runtime candidate and one-off workflow were removed after this result.

## Unchanged-source baseline

The [CI baseline run 36657129349](https://github.com/aburan28/crypto/actions/runs/36657129349) completed on AMD EPYC 7763, Linux x86-64, Rust 1.98.1, with AVX2/BMI2 and one Rayon thread. Its ten-second preflight saw 0.19 outside CPU seconds; the privileged reservation moved 152 threads, left zero user threads eligible on the pinned CPU, and recorded zero contended samples. All seven F5 cases returned exact output. The n24 degree-4 full call took 106.119 ms, including 43.309 ms elimination and 54.528 ms output unpacking, with 6,924 rows and 13,734,979 terms.

The unmodified binary SHA-256 is `604a81550ed7fc1caa232fffa2aee5850658a2092e6f3bd5dad8d26c43abd293`. Its compressed binary and hash are in [BINARY_ARCHIVE.json](runs/36657129349-attempt1/BINARY_ARCHIVE.json); the [baseline manifest](runs/36657129349-attempt1/baseline-manifest.json) and [workflow](WORKFLOW_BASELINE.yml) preserve source, host, exact output and isolation. This standalone baseline time is never used as the candidate denominator.

## Same-binary paired candidate

[CI run 36657747432](https://github.com/aburan28/crypto/actions/runs/36657747432) used binary SHA-256 `5b74cc523423edb3bb3420a70780ed3808f0278c0ae8bb17ff3e4ca80447ca38` for both explicit `KIC_F5_UNPACK_BYTE_TABLE=0/1` arms. Each thread-count job completed all 88 calls, seven cases per call, with matching raw and canonical row fingerprints, rank, terms, columns, builder and criterion counts, and reduction word operations. Mode 1 selected the table on n20 and n24 degree-4 cases; the other cases used the reference path. The largest table allocation was 14,927,180 bytes, below the frozen 16 MiB cap.

Neither wall-clock job qualified under the frozen isolation rule: the one-thread reservation recorded four contended samples and the two-thread reservation recorded three. Both left zero user threads eligible on reserved CPUs. The ratios below are **exploratory, not accepted speedup measurements**. They compare paired reference/candidate calls on the same host and binary; a ratio below 1 means the candidate was slower.

| Frozen n24 degree-4 case | One thread | Two threads |
| --- | ---: | ---: |
| Paired full-call reference/candidate median | 0.645× | 0.667× |
| Exact bootstrap 95% interval | 0.569–0.681× | 0.649–0.704× |
| Reference full-call median | 106.34 ms | 96.93 ms |
| Candidate full-call median | 164.83 ms | 145.22 ms |
| Reference unpack median, including mode setup | 54.89 ms | 37.53 ms |
| Candidate unpack median, including table build | 110.59 ms | 84.85 ms |
| Frozen A/A full-call ratio range | 0.952–1.036× | 0.986–1.049× |

All four n24 degree-4 seeds regressed: one-thread paired medians were 0.645–0.660× and two-thread paired medians were 0.667–0.685×. The n20 degree-4 case also regressed, with frozen-seed full-call ratios 0.704× and 0.763× at one and two threads. These deviations are much larger than the respective A/A ranges. The table path missed the preset 1.10× full-call target and 1.05× lower-bound gate by a wide margin; no default-on confirmation is warranted. Its 14.9 MB table footprint likely displaced useful column and output data from cache, but that explanation is an inference, not a measured cache-miss attribution.

The complete raw receipts, including every call and contention sample, are deterministically compressed with SHA-256 manifests for [one thread](runs/36657747432-attempt1-t1/MANIFEST.json) and [two threads](runs/36657747432-attempt1-t2/MANIFEST.json). [RUN_CANDIDATE.json](RUN_CANDIDATE.json) preserves GitHub's run metadata. The failed CI qualification jobs are retained rather than rerun for a favorable time. This experiment supplies no 2× result.
