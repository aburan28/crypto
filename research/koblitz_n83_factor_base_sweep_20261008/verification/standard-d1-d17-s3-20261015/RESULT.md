# N83 standard-subspace dimension frontier in S3

The committed exporter was run once at every supported dimension 1–17 from source commit `604c2631fcd13393ab346ea6be142854111b4a83`. The exact native build and reconstruction passed for dimensions 2–17. Dimension 1 was rejected with `cofactor projection contains no nonidentity point`; an [independent group check](d1-independent.json) confirmed that its three source points all project to the identity. It has no stored archive. The prior [dimension-18 archive](../standard-d18-s3-20261015/RESULT.md) completes the valid standard-subspace frontier through the exporter's dimension-18 cap.

| Dimension | Projected subgroup points | Sign-folded columns | Compressed bytes |
| ---: | ---: | ---: | ---: |
| 2 | 2 | 1 | 636 |
| 3 | 6 | 3 | 763 |
| 4 | 12 | 6 | 945 |
| 5 | 24 | 12 | 1,304 |
| 6 | 52 | 26 | 2,161 |
| 7 | 122 | 61 | 4,242 |
| 8 | 258 | 129 | 8,249 |
| 9 | 524 | 262 | 15,909 |
| 10 | 1,048 | 524 | 30,779 |
| 11 | 2,036 | 1,018 | 58,590 |
| 12 | 4,054 | 2,027 | 115,139 |
| 13 | 8,086 | 4,043 | 228,033 |
| 14 | 16,122 | 8,061 | 453,526 |
| 15 | 32,654 | 16,327 | 914,231 |
| 16 | 64,904 | 32,452 | 1,810,844 |
| 17 | 130,604 | 65,302 | 3,635,635 |
| 18 | 261,444 | 130,722 | 7,263,179 |

The 16 new compressed objects total 7,280,986 bytes. The exporter binary SHA-256 is `56003eed67a803eda956fe6003c5f64ccdfda96345eac2b55d13bc0b58497cb8`, the same binary pinned in the dimension-18 construction receipt; the committed exporter source and `src` tree are checked in each manifest replay. Each valid dimension has a content-addressed object at `s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/standard-subspace/a0/objects/<compressed_blake3>.jsonl.gz` and a manifest under the same prefix at `d<dimension>/manifest.json`. [s3-replay.json](s3-replay.json) records every exact key, SHA-256, write-once `If-None-Match: *` response, AES-256 S3 HEAD checksum, fetched-byte comparison, and native reconstruction result. Together with the 57 previously retained signed-Frobenius archives and the dimension-18 standard archive, this branch accounts for 74 stored factor-base archives.

The separate [independent audit of the S3 downloads](audit-s3.json) checked the full compressed and uncompressed hashes, sorted unique point rows, every sign pair, 1,028 sampled curve equations and 73 independent subgroup multiples across dimensions 2–17. It also proved exact nesting of the projected point sets from dimension 2 through dimension 18; the dimension-18 comparison is byte-pinned to its previously audited archive. The [one-byte dimension-17 corruption control](tamper-control.json) and a [dimension-2 sign-label corruption with all byte hashes recomputed](semantic-tamper-control.json) were rejected by both the independent auditor and the native verifier. [construction.json](construction.json) and the append-only [construction](construction.jsonl) and [upload](upload-progress.jsonl) journals retain the dimension-1 failure and each successful step. Construction wall values are unisolated L0 observations on a saturated host.

The new archives expand the factor-base size frontier but do not rank complete index-calculus runtime. Cofactor projection generally moves points outside the original x-subspace, so a source-point-to-projected-column adapter is still needed before an algebraic subspace solver can consume these archives. No dimension in this frontier has a measured natural relation yield, useful rank, graph closure, sparse linear algebra, fresh individual-log runtime, or matched rho comparison. `selected_best_total_runtime` remains null. The approved first mixed-orbit producer gate remains separately unspent after its [busy-host preflight](../mixed-orbit-first-gate-preflight-20261015/RESULT.md).
