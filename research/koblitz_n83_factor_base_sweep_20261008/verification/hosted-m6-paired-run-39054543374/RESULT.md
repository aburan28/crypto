# Hosted six-summand K600 Gray/hash pair

[GitHub Actions run 39054543374](https://github.com/aburan28/crypto/actions/runs/39054543374) completed on 2026-10-17 at exact source and workflow commit `17f03c518377468eeb49331ecfd8c3dedc2d72d0`. The uploaded artifact is ID `12111425305`, named `n83-hosted-m6-39054543374-1`. The [frozen protocol](../hosted-m6-paired-protocol-20261017/PROTOCOL.md) fixes the same public a=0 N83 target, first-trial coefficients, sequential K600 base, three S3-backed archive identities, and 4+2 six-summand search for both arms. The only changed factor-base input is the K600 residual pool: Gray-prefix versus hash-derived.

Both arms imported successfully and returned **`UNKNOWN_search_cap`**, each after one SAT call, one identity pattern and **20,000 conflicts with zero models**. Neither pattern was refuted. A conflict cap is a censored search result; the zero observed models cannot be turned into a zero-yield estimate, an UNSAT result, or a base ranking.

| Measurement | Gray-prefix pool | Hash-derived pool |
| --- | ---: | ---: |
| Sequential-base signed Frobenius orbits | 600 | 600 |
| Shared signed orbits | 503 | 0 |
| Residual signed orbits | 97 | 600 |
| Residual points | 16,102 | 99,600 |
| Domain clauses | 237,280 | 323,796 |
| Peak variables / clauses / XOR rows | 118,684 / 547,615 / 6,059 | 120,696 / 634,131 / 6,059 |
| Peak cgroup bytes | 380,858,368 | 639,315,968 |
| Search container wall, local observation | 1,497.377 s | 1,555.412 s |
| Solver calls / conflicts / models | 1 / 20,000 / 0 | 1 / 20,000 / 0 |
| Raw status | `UNKNOWN_search_cap` | `UNKNOWN_search_cap` |

Both no-network ARM64 containers used one CPU, 4 GiB hard memory and swap limits, 64 PIDs and a read-only root. They exited 0 without OOM and were removed. The GitHub job ran from 17:56:24 to 18:56:01 UTC, 3,577 seconds including pilot regeneration, build, imports, searches and upload. The paired supervisor reported 3,059.359 seconds. The accounting charges the **full 7,500-second hosted envelope**, in addition to earlier separately charged work. The observed wall difference is not a complete-runtime comparison.

The [independent audit](audit.json), implemented by [replay.py](replay.py), checks the plan and source commit, public input, three regenerated archive hashes, manifest, built worker, pinned image, container commands, limits, mounts, imports, search counters, exit/cleanup states, and outer classification. It passed as `PASS_MATCHED_4_PLUS_2_ONE_QUERY_REPLAY`. The three input archives were freshly downloaded from the requested S3 factor-base prefix and passed SHA-256, BLAKE3, length and single-part ETag checks in the [S3 receipt](factor-bases-s3-check.json). A separate [inventory](pilot-s3-inventory.json) confirms that all 54 manifest-listed pilot factor-base objects were present at the expected sizes; that inventory did not rehash all 54 contents.

The [custody manifest](custody.json) hashes every one of the 65 raw artifact files. A deterministic 9,473,710-byte local archive was created and byte-replayed against those files; its SHA-256 is `221e47216c6315764a6381dc95aa94f850ed6ca5ef6b06cff1843b5466d3905e`. The raw GitHub artifact remains downloadable from the run, and the local archive is at `/Volumes/SSD990/llm/tmp/n83-hosted-m6-run-39054543374.tar.gz`. Automatic export review rejected upload of this broader archive to S3 because it also contains a worker binary, logs and configuration; it is **not** represented as an S3 round trip. The three factor-base archives themselves are stored and freshly verified in S3.

For independent replay, download the GitHub artifact and the three `s3_uri` objects in the S3 receipt, then run:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/hosted-m6-paired-run-39054543374/replay.py \
  --artifact ARTIFACT_DIR --s3-bases S3_BASE_DOWNLOAD_DIR
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/hosted-m6-paired-run-39054543374/custody.py verify \
  --artifact ARTIFACT_DIR --archive LOCAL_ARCHIVE --receipt \
  research/koblitz_n83_factor_base_sweep_20261008/verification/hosted-m6-paired-run-39054543374/custody.json
```

The first command requires the `blake3` Python package. A fresh archive can be generated with `custody.py create` before the second command. The separate tests in [test_evidence_tools.py](test_evidence_tools.py) reject a changed search status and a changed archive member.

This is one public-target query per arm. Natural partial-relation rate, graph closure, useful matrix rank, linear algebra, individual-log extraction, verification and matched rho remain unmeasured for this comparison. `selected_best_total_runtime` and `total_index_calculus_runtime_ms` remain null.
