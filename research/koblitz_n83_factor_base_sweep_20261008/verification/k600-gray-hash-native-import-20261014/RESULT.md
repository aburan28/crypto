# Frozen K600 Gray/hash native import check

Status: `PASS_NATIVE_IMPORT_ONLY`. The native macOS worker imported both retained
K600 pairs from the S3-backed archives, and the independent Python audit
replayed their point domains and signed-orbit overlap. No SAT relation search
was launched in this check. Graph rows, useful matrix rank, target-log work,
and complete index-calculus runtime remain unmeasured, so this check selects
no runtime-winning factor base.

The frozen input plan is `../k600-gray-hash-paired-protocol-20261014/plan.json`
(SHA-256 `fbd3fb943bc0f181797d106c41ee9961b8d70768b5c65a5584002e3c1b5e488f`).
Both arms use the 600-column sequential base object `3888e259…cab1a0` from
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/a0/objects/`.
The Gray pool is `966a0493…66dc47`; the hash pool is `b4dda45e…9a0363e`.
The plan and `result.json` carry full object hashes, S3 URIs, compressed SHA-256
digests, source hashes, and exact commands.

| Stage | Outcome | Wall time, L0 | Max RSS, L0 | Import/domain result |
| --- | --- | ---: | ---: | --- |
| Clean offline release build | exit 0 | 863.081 s | 3,996,385,280 B | Rust 1.98.1; binary SHA-256 `c5b9af9930f293533de819c0a673a0b7149ad2418800fa3fee7509d772d207e7` |
| Gray pool | `PASS` | 31.482 s | 150,241,280 B | 600 + 600 columns; 503 shared signed orbits; 97 residual orbits / 16,102 residual points |
| Hash pool | `PASS` | 20.817 s | 167,804,928 B | 600 + 600 columns; 0 shared signed orbits; 600 residual orbits / 99,600 residual points |

The build log reports a `rust-objcopy` debug-info stripping warning because
that tool could not load `libLLVM.dylib`. Cargo nevertheless exited 0, and the
built executable passed both imports. The raw warning is retained unchanged.

The two worker receipts report `solver_search_executed=false`. The runner
checked each worker JSON receipt with `receipt_audit.audit_verify` and then
separately parsed the compressed point archives with `receipt_audit.domains`.
After copying the raw streams into this directory, their SHA-256 values were
recomputed and matched against `result.json`. Source commit
`ec2b665457c9c4e76b48481af35c2318510a4c00` was clean during the build and
both imports; the pinned Rust `src` tree was
`f5dac3ceb6d264227ef74b197e49f3cf01d1164e`. The frozen source and archive
checks were repeated before each worker call.

A separate `replay.py` invocation then checked the committed source bytes at
that commit, the executable hash, every compressed archive SHA-256, all six
retained raw-stream hashes, both worker JSON receipts, and both parsed point
domains. Its `replay.json` status is
`PASS_FROZEN_SOURCE_BINARY_ARCHIVES_AND_DOMAINS`. The replay used the same
local archive bytes audited against S3 in the frozen paired protocol; it did
not perform a fresh S3 download.

The Gray and hash calls ran in that order on one host. Their import timings
have cache and order effects and omit solver work. The source-only query floors
(33,975 Gray; 889 hash, conditional on the stated half-K row model) are capacity
bounds, not observed relation throughput. Neither these import timings nor the
source bounds rank total attack runtime.

The first attempt stopped at a toolchain-version assertion before any build or
import; see `PRELAUNCH_ATTEMPT.md`. The corrected run above is preserved
separately. Raw files: `build.stdout`, `build.stderr`, `gray_pool.stdout`,
`gray_pool.stderr`, `hash_pool.stdout`, `hash_pool.stderr`, and `result.json`.
The executable was copied to
`/private/tmp/n83-native-k600-import-20261014/koblitz_n83_mixed_orbit_worker.preserved`
and rehashed to the same SHA-256. Only the task-owned Cargo `target/`
intermediates (478,112 KiB before removal) were deleted to recover system-disk
headroom. `cleanup.json` records the exact path and hash; the independent
replay was run again against the preserved executable and passed. Its
`post_cleanup_replay.json` is byte-identical to the original `replay.json`.
The executable itself is not committed to Git and can be rebuilt from the
pinned source if this scratch path is lost.
