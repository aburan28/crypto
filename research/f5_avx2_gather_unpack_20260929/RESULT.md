# Matrix-F5 compact AVX2 gather unpack: rejected

The [protocol](PROTOCOL.md) was frozen before candidate code or timing. The candidate decoded four set-bit column indices at a time, gathered four values from a capped `u32` column map with AVX2, widened them to `u64`, and wrote four exact `F2BoolMono` terms. The measured source is retained as [measured_candidate.patch](measured_candidate.patch) and the [candidate workflow](WORKFLOW_CANDIDATE.yml) is archived. The opt-in runtime option and active workflow were removed after the frozen complete-call gate failed. This is a matrix-F5 solver-stage result, not a one-target IC DLP or rho speedup.

## Same-binary paired result

[CI run 36668480279](https://github.com/aburan28/crypto/actions/runs/36668480279) passed the shared GF(2) and F5 exactness suites and completed both paired jobs on Linux x86-64 with Rust 1.98. The tested PR head was `bc3e2cdc8ae012f392616481e7a058ce6f17bfdd`; GitHub's pull-request merge checkout was `1947a37e5e094457521723788ecf38c430c4af8f`. Both jobs used binary SHA-256 `ecab7aaedeb9f3d6587ef31f107e821b4cebce2931ee9ecdf39e845d95eaa329` for explicit `KIC_F5_UNPACK_AVX2_GATHER=0/1` arms. The **one-thread job ran on Intel Xeon Platinum 8573C** and the **two-thread job on AMD EPYC 9V74**. Each ratio compares only processes on its own host; absolute times across these hosts are not compared.

For each seed and thread count, the first isolated, exact block was selected before reading its timing. Each selected block contains two warmups, five A/A reference pairs, and five alternating A/B pairs: 22 processes and all seven F5 cases per process. All 176 selected processes matched raw and canonical row fingerprints, rank, terms, columns, builder and criterion counts, and reduction word operations. Every selected reservation left zero eligible user threads on its reserved CPUs and had zero contended samples. The one-thread job retained five preflight refusals with zero benchmark calls; the two-thread job selected every first attempt. Mode 1 used the gather path only on n20 and n24 degree-4, with actual map allocations of 24,784 and 51,804 bytes, respectively, below the frozen 256 KiB cap. Source/binary hashes, CPU flags, affinity, route and complete process output are in the receipts.

The table reports the n24 degree-4 **complete F5 call**, including map conversion and allocation. Ratios are paired reference/candidate medians; separate arm medians need not divide to that paired value. Intervals are exact 3,125-resample bootstrap 95% intervals over the five A/B pairs for each seed. A ratio below one means the candidate was slower.

| Rayon threads and host | Seed | Reference median, ms | Candidate median, ms | Paired full-call ratio | 95% interval | A/A ratio range | Unpack ratio |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1, Intel | frozen | 94.438 | 97.325 | **0.9652×** | 0.9630–0.9874× | 0.9848–1.0398× | 0.9488× |
| 1, Intel | holdout A | 91.647 | 95.756 | 0.9711× | 0.9199–0.9804× | 0.9757–1.1002× | 0.9528× |
| 1, Intel | holdout B | 90.536 | 92.683 | 0.9868× | 0.9282–0.9921× | 0.9727–1.0378× | 0.9589× |
| 1, Intel | holdout C | 92.712 | 92.443 | 0.9946× | 0.9745–1.0292× | 0.9217–1.0018× | 0.9911× |
| 2, AMD | frozen | 73.140 | 80.023 | **0.9140×** | 0.8887–0.9394× | 0.9408–1.0435× | 0.8108× |
| 2, AMD | holdout A | 75.202 | 81.452 | 0.9223× | 0.9134–0.9890× | 0.9680–1.0561× | 0.8287× |
| 2, AMD | holdout B | 75.565 | 80.319 | 0.9389× | 0.9200–0.9441× | 0.9791–1.0350× | 0.8547× |
| 2, AMD | holdout C | 75.748 | 80.508 | 0.9354× | 0.9120–0.9634× | 0.9710–0.9999× | 0.8545× |

The one-thread frozen full-call median and lower bound miss the preregistered 1.05× and 1.02× gates. The frozen one-thread n24 result also falls below its own A/A minimum, as do all four two-thread n24 results. Unpacking itself regressed on both hardware classes; moving four lookups into one AVX2 gather did not compensate for gather and index-preparation cost on these hosts. That causal explanation is an inference from the measured phase ratios, not a hardware-counter result. This candidate supplies no accepted incremental gain or 2× result. IC end-to-end cost remains unknown.

## Raw receipts and reproduction

- [One-thread manifest](runs/36668480279/MANIFEST-t1.json) and `segments-t1.tar.gz`: all nine attempts, including five preflight refusals; archive SHA-256 `ea2ce8086034143a2ddc342a7a6b2ba45aa5987d34751d409feb4c74efc3525d`.
- [Two-thread manifest](runs/36668480279/MANIFEST-t2.json) and `segments-t2.tar.gz`: four first-clean attempts; archive SHA-256 `b26b473905d11beabdc3690ee149a556b376d7b004132d38b0dbc866cf6fd351`.
- [Run metadata](RUN_CANDIDATE.json) records both successful CI jobs. Each manifest lists every file's byte count and SHA-256, selected attempt, binary/source hashes, exactness, ratios and rejected attempts. Extract an archive with `tar -xzf segments-t1.tar.gz` or `tar -xzf segments-t2.tar.gz` from this result directory.

To reproduce the measured candidate, apply `measured_candidate.patch` to the protocol's pinned source, restore `WORKFLOW_CANDIDATE.yml` to `.github/workflows/f5-avx2-gather-segments.yml`, and run the workflow with both modes in one release binary. Preserve the first-clean selection rule. Local ARM correctness tests and an x86-64 Linux type check passed before CI; the supported x86-64 CI tests executed the AVX2 path.
