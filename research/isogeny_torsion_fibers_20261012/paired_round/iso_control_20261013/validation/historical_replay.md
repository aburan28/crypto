# Native replay of the retained failed runs

Date: 2026-10-13 local. V1 and v2 source files were reconstructed by
reversing the separately versioned v3 and v2 runner changes. The
[v1 reconstructed source](../runner-v1-reconstructed.rs) SHA-256 is
`e74f5e1664da2162d26c504adf7ce889ce07d4b37363040aa9800a84af6c8a27`,
**identical** to the source hash recorded before v1 execution. The
[v2 reconstructed source](../runner-v2-reconstructed.rs) SHA-256 is
`9e20e7dc0c7cfeef4194adea3220912600e6cd8b231bf88c750288ff5d3c0a6a`;
no pre-run v2 source hash was recorded, so byte identity to the original
v2 source cannot be certified. Both source files are native Rust and
retain the corresponding protocol parameters.

Both were rebuilt with `rustc --edition 2021 -O -L
dependency=target/release/deps --extern
crypto_lib=target/release/deps/libcrypto_lib-2f648d1af79d9c46.rlib`.
The library artifact SHA-256 was
`160187f29fb1d6db80445cd19f93fe227c7d309ee2e01a0c0c813e9bd5a35961`.
The reconstructed v1 and v2 binaries had SHA-256
`ef75307745ab1b2a9877e9e8f14480c8754e347a48eba1ff2bcbdd99194d0b31`
and
`a0ddd8e2f78ae53bfe6fd50782895045eb8860b2c269b9bd8400d15137b98ee3`.
These direct-rustc binary hashes differ from the original Cargo v1 binary
hash and do not identify the original v2 binary.

Each rebuilt binary exited 2, as the original failed run did. `cmp` passed
byte for byte for both pairs of stdout and stderr:

| run | raw TSV | run log | original run status |
|:--|:--|:--|:--|
| v1 | [raw-v1.tsv](../raw-v1.tsv) = [replay-v1.tsv](../replay-v1.tsv) | [run-v1.log](../run-v1.log) = [replay-v1.log](../replay-v1.log) | six-class cap at F13 target 4; later targets not reached |
| v2 | [raw-v2.tsv](../raw-v2.tsv) = [replay-v2.tsv](../replay-v2.tsv) | [run-v2.log](../run-v2.log) = [replay-v2.log](../replay-v2.log) | six-class caps at F13 targets 4, 5, 6 |

The original raw/log evidence was not rewritten. This replay validates
deterministic behavior under the reconstructed source and current native
library artifact. It does not turn either failed run into an accepted
five-arm performance comparison or provide independent-host validation.
