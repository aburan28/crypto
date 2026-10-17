# Software verification receipts

All records here are software/build checks, not additional AES scientific
trials or performance measurements. The original 41-command suite remains in
the linked `cryptanalysis` archive.

| Check | Result | Record |
| --- | --- | --- |
| `cargo +stable fmt --manifest-path tools/aes-discover/Cargo.toml -- --check` | pass | native-fmt.log |
| Release locked all-target Clippy with warnings denied | pass | native-clippy.log |
| Release locked native tests | 12 library + 2 integration pass; 24 OpenSSL comparisons inside integration test | native-tests.log |
| Eleven-file companion parity | pass | source-parity.log, native-source.sha256 |
| Changed-reference software fixture | expected mismatch rejected | source-parity-mutant.log |
| Nine blocked backlog cards | structural validation pass | backlog-validation-with-base-problem.log |
| Required root release library suite | compilation fails with 26 errors before tests | crypto-lib-attempt-2.log |

`native-environment.log` records compiler, OpenSSL, source base and native
executable hash. `base-library-identity.log` binds the failing unchanged base
inputs. The initial root attempt's missing sparse-checkout dependency is
preserved separately; the dependency and compile-time fixtures were materialized
before the substantive second attempt. See the
[baseline problem](../../../problems/aes-integration-root-library-baseline.md).

The backlog CLI was compiled from its unchanged source and unchanged SHA-1/
SHA-256 modules using the retained isolated `backlog-validator.Cargo.toml`,
lockfile and hash-module adapter. This avoids linking the broken unrelated
library for a card-structure check; it does not make the root suite pass.
`backlog-validator-source.sha256` binds those original modules. The first
isolated build omitted its existing libc dependency; the corrected build and
successful validation are preserved separately. No checker assertion was
changed. To reproduce, copy the retained Cargo files and adapter into a temporary
directory, name them `Cargo.toml`, `Cargo.lock`, `hashes.rs`, build, and run
`backlog --root /workspace/crypto validate` with the nine `XP-20261017-aes-*` IDs.
Adjust the original checkout prefix if using a different location.
