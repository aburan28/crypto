# Published-Q development control, before the new evaluation freeze

The unmodified main-source release binary and one candidate release binary (W64 `KIC_S3_PREFILTER=off` and `blocked`) ran on the already published S3-batch n41/K255 and n53/K440 point-only Q. This was a correctness/control pass, not an isolated timing panel or a held-out speed claim. The [independent receipt](DEVELOPMENT_RECEIPT.json) replays all three base-log and rank traces and all 6,144 recovered target logs. In each n, all three modes have identical base, rank transitions, first accepted four-point witness, probe count, recovered scalar and group verification for all 1,024 Q. The original S3 and scalar recovery tests plus a new exhaustive GF(2^5) no-false-negative test pass.

| n | rank + target root keys | definite table misses skipped | false positives | true table hits | filter bytes | recovered logs |
|:--|--:|--:|--:|--:|--:|--:|
| 41 | 3,220,797 | 3,187,392 | 32,126 | 1,279 | 4,194,304 | 3×1,024 |
| 53 | 27,632,196 | 27,379,072 | 251,660 | 1,464 | 16,777,216 | 3×1,024 |

The skipped-miss rate is about 99% in both sizes. It does **not** show a complete-process improvement: the filter adds construction, hashing and memory, and these development runs were on an unisolated Arm64 laptop. Their exploratory in-process times are retained in the raw data but do not enter the pre-registered decision. The full frozen Linux off/filter/rho comparison on fresh Q is the next gate. The filter has no route to improve the K+L attempt floor or native m≥3 PDP yield.

Reproduce the receipt by extracting [the 4.5 MiB raw archive](development_raw.tar.gz) into a scratch directory and running `python3 verify_development.py --raw <scratch>/s3filter-dev-20260930 --out <new-receipt.json>`. The archive SHA-256 is `833423d385a01f24d0cdbb098979dcf070d801f53e60d96761a559eeffcd59f8`, the original source SHA-256 is `e998120842faf393910a6b4a9740b3af9880eb626578fa88c575a32e5254f2c7`, and the candidate source SHA-256 is `702a0a05709bc14bc10bafbb08edb2a2f5f86794e970d84c317da3bc65cdcf38`. The original public-Q and fixture hashes are in the previous [freeze](../compact_s3_batch_20260929/FROZEN.json); every raw file hash is in the receipt.
