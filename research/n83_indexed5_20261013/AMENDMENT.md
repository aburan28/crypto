# N83 five-byte index: real-field parity gate

The indexed five-byte table returned the same first S3-root witness as the
existing packed table for every inserted key from six archived K512 base
representatives, and agreed on 2,047 additional canonicalized probes. The
test constructs the actual `F_2^83` polynomial basis, all 83 Frobenius shifts
of each representative, and the corresponding S3 roots. A ten-bit fingerprint
match always triggers comparison with the full root reconstructed from its
retained state. This establishes exact lookup behavior on real N83 roots in
addition to the synthetic collision and token-boundary tests.

| Evidence | Value |
| --- | --- |
| Final test source SHA-256 | `64171575e58b97c4ad9041b4b49e5f2b5697d26f71be73861fc36ea879db9161` |
| Command | `TMPDIR=/Volumes/SSD990/llm/tmp CARGO_TARGET_DIR=/Volumes/SSD990/llm/tmp/n83-indexed5-cargo-target cargo test --release --locked --offline --example koblitz_orbit_dlp_fast root_bloom_table_tests -- --nocapture` |
| Result | 6 passed, 0 failed, 9 filtered |
| Raw test log | `example_tests_v2.log`, SHA-256 `a8d645364c9f7081857dc1548062e4089b94c318fabc93f10269b190b51a9536` |
| Source formatting | `rustfmt --check examples/koblitz_orbit_dlp_fast.rs` passed |

The first frozen full-pipeline K512 control used the earlier source and binary
hashes in `PROTOCOL.md`. It reached its 900-second cap with exit 124, peak RSS
844,267,520 bytes and no summary, base dump, or target row. The independent
analysis records `TIMEOUT` and leaves witness parity unknown for that full
run. Its raw files and hashes are unchanged. The test amendment is a new code
revision and cannot retroactively turn that run into a verified target result.

The five-byte layout reduces an allocated `2^30`-slot root table from 8 GiB
to 5 GiB; a `2^31`-slot table is 10 GiB. These are allocation calculations,
not end-to-end timing measurements. The K1798 and K2500 gates remain unrun
under this revision. The local macOS data volume had 3.4 GiB free when the
parity gate completed, below the frozen 10 GiB K1798 preflight. A new
full-pipeline run requires a new frozen source/binary receipt and its own
write-once run directory.
