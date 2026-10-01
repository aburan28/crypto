# Matrix-F5 compact paired unpack: rejected

The [protocol](PROTOCOL.md) was frozen before candidate code or timing. The tested implementation converted the F5 column masks to a capped `u32` lookup and decoded two set bits at a time into the exact `Vec<F2BoolPoly>` output. Its source is retained as [measured_candidate.patch](measured_candidate.patch); the [candidate workflow](WORKFLOW_CANDIDATE.yml) is archived. The runtime option and active workflow were removed after the frozen gate failed. This is a matrix-F5 solver-stage result, not a measured one-target IC DLP or rho speedup.

## Same-binary paired result

[CI run 36664115975](https://github.com/aburan28/crypto/actions/runs/36664115975) completed on Linux x86-64, AMD EPYC 7763, Rust 1.98. The tested PR head was `a37f8b6ff814a4568dfadb820df5236ff75df056`; the pull-request merge checkout was `150b356b2ad0ff91db5836cc8b5bd7b4a020e8d1`. Both thread-count jobs used binary SHA-256 `56d1a8596b86a081d9818707827f4fb5b097e9892c707db21fa63fc1d3ade982` for explicit `KIC_F5_UNPACK_COMPACT_UNROLL=0/1` arms. The benchmark and GF(2) kernel SHA-256 values, CPU flags, affinity, source hashes, and all process outputs are in the receipts.

For each seed and thread count, the first isolated, exact block was selected before inspecting timing. Each selected block has two warmups, five A/A reference pairs, and five alternating A/B pairs: 22 processes and all seven F5 cases per process. All 176 selected processes matched raw and canonical row fingerprints, rank, terms, column and builder counts, and reduction word operations. Both modes used the direct pack and direct unpack routes; mode 1 selected compact paired unpack only on the n20 and n24 degree-4 cases. The largest actual compact-map allocation was 51,804 bytes, below the frozen 256 KiB cap. Every selected reservation left zero eligible user threads on the reserved CPUs and had zero contended samples. The one-thread job also retained four preflight refusals with zero benchmark calls; the two-thread job selected every first attempt.

The table reports the n24 degree-4 **complete F5 call**, including map conversion and allocation. Ratios are paired reference/candidate medians, so independent arm medians need not divide to the displayed ratio. Intervals are exact 3,125-resample bootstrap 95% intervals from each seed's five A/B pairs. A ratio below one means the candidate was slower.

| Rayon threads | Seed | Reference median, ms | Candidate median, ms | Paired full-call ratio | 95% interval | A/A ratio range | Unpack ratio |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | frozen | 107.594 | 107.615 | **0.9998×** | 0.9986–1.0045× | 0.9969–1.0096× | 0.9999× |
| 1 | holdout A | 107.037 | 106.618 | 1.0047× | 0.9961–1.0104× | 0.9848–1.0158× | 1.0053× |
| 1 | holdout B | 106.364 | 105.184 | 1.0090× | 0.9891–1.0178× | 0.9880–1.1235× | 1.0064× |
| 1 | holdout C | 106.881 | 108.203 | 1.0006× | 0.9840–1.0118× | 0.9883–1.0090× | 0.9970× |
| 2 | frozen | 95.165 | 93.038 | 1.0240× | 1.0107–1.0550× | 0.9824–1.0252× | 1.0735× |
| 2 | holdout A | 96.403 | 93.843 | 1.0257× | 1.0094–1.0410× | 0.9879–1.0231× | 1.0738× |
| 2 | holdout B | 95.793 | 93.392 | 1.0247× | 1.0079–1.0406× | 0.9793–1.0272× | 1.0927× |
| 2 | holdout C | 96.071 | 93.472 | 1.0278× | 0.9152–1.0359× | 0.9861–1.0022× | 1.0884× |

The one-thread frozen complete-call ratio and its lower bound miss the preregistered 1.05× and 1.02× gates. Its 1.000× unpack ratio shows the paired decode did not accelerate the dominant output phase on that host. The n16 degree-3 full-call median also fell below its own A/A minimum on the frozen one-thread seed and holdout-A two-thread seed. The two-thread gain does not repair the one-thread failure. This candidate supplies no accepted incremental gain or 2× result. The observed result is engineering evidence about a solver stage only; IC end-to-end cost remains unknown.

## Raw receipts and reproduction

- [One-thread manifest](runs/36664115975/MANIFEST-t1.json) and `segments-t1.tar.gz`: all eight attempts, including four preflight refusals; archive SHA-256 `63ebaee42057918586ac56e5cb70a203aedde74471f9d63134b992fb47376590`.
- [Two-thread manifest](runs/36664115975/MANIFEST-t2.json) and `segments-t2.tar.gz`: four first-clean attempts; archive SHA-256 `1856d8c61b03b29b9006f17b2a247afa08a756cc970d53ebca3c07fc1b7c2daf`.
- [Run metadata](RUN_CANDIDATE.json) records both successful CI jobs. Each manifest lists every file's byte count and SHA-256, selected attempt, binary/source hashes, exactness, ratios, and rejected attempts. Extract an archive with `tar -xzf segments-t1.tar.gz` or `tar -xzf segments-t2.tar.gz` from this result directory.

To reproduce the measured candidate, apply `measured_candidate.patch` to the protocol's pinned source, restore `WORKFLOW_CANDIDATE.yml` to `.github/workflows/f5-compact-unrolled-segments.yml`, and run the workflow at the frozen PR head. Preserve the paired baseline/candidate processes and the first-clean selection rule; absolute times from a different host are not a denominator for these ratios.
