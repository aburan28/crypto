# Dedicated executable validation, 2026-10-09

The CLI and release integration were implemented in commit `5bcc591f3`. The underlying
research sources and canonical curve records remain unchanged. Native release builds
use Rust/Cargo 1.93.1; local validation used an Apple M4 Pro running macOS 26.6.

| Requirement | Evidence and result |
| --- | --- |
| One installed executable without a checkout or build tools at runtime | Copied only `isogeny` into `/private/tmp/isogeny-installed.Q5Mcu2`; ran all commands from that directory |
| Dedicated package tests | 25 tests passed, including extension arithmetic, field arithmetic, exact map identity and coefficient-mutation rejection, hash vectors, and CLI input validation |
| P224 degree-1471 construction | PASS; output byte-identical to the sealed study map, SHA256 `ecf1109dbff2a72000fd84d9939127db68e2287e55e4ec10d143229874c230c8` |
| P224 default search includes independent replay | PASS; kernel, codomain, exact rational map, standard subgroup and 20 scalar-transport checks |
| P192 degree-10453 construction | PASS; output byte-identical to the sealed study map, SHA256 `aff763af2e9adb95765a7f29223aed7dc46c81f488301a3ba06db10de5595f29` |
| P192 degree-73 full enumeration | PASS; two independently certified maps, coverage COMPLETE |
| Screening from 1010 through 65537 on both source curves | Candidate records equal the sealed `small-orders-65537.json` after normalizing object-key order |
| Independent verification does not trust the producer's hash alone | Changed a kernel coefficient and updated the receipt hash; `verify` rejected the result with status FAIL and exit 1 |
| Workflow syntax and shell checks | Actionlint passed; YAML parsed successfully |
| Frozen research preservation | Continuation manifest verified all 209 bound files after implementation |
| Root release-library gate | Rerun with `cargo test --release --lib --locked --offline`; still fails with 643 pre-existing errors |

The independent high-degree P192 certificate remains the sealed study certificate.
This packaging validation repeated its construction and checked byte identity; it did
not repeat the approximately 11-minute P192 independent replay. The P224 replay and
P192 degree-73 pair were independently verified through the installed executable.
No new performance comparison is claimed.

The source-snapshot branch referenced by the original research receipts is published
as `codex/p192-p224-large-degree-isogenies-build-snapshots`. The research and release
implementation are in [PR #1606](https://github.com/aburan28/crypto/pull/1606).

The release workflow validates installed binaries on Linux glibc, Linux musl,
Apple silicon and Windows, and packages standalone archives with README and BUILD.json.
Cross-platform results remain pending until those CI jobs finish; local macOS checks
are not cross-platform validation.

ReleaseMe is `dev-build-deploy/release-me` v1.1.0, pinned to
`cede87fd8a4829ef571dff6bde807ebbfcaeb3e9`. Explicit stable release versions use the
crate's major/minor and `crate patch + workflow run number`, since ReleaseMe marks
versions with a prerelease suffix as prereleases. The tag is created at the exact
built commit before ReleaseMe creates the release, preventing branch advancement
during a build from changing the tag's source identity. Publication remains restricted
to main after all release jobs pass. PR dry runs publish no GitHub release.
