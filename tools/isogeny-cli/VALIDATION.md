# Dedicated executable validation, 2026-10-09

The latest native setup/hashing continuation is recorded in
[`research/isogeny_catalogue_setup_20261009/RESULTS.md`](../../research/isogeny_catalogue_setup_20261009/RESULTS.md).
Against commit `5cbc3a2dd` with the same 349-model input, five matched pairs per panel
used 67.9589% fewer instructions for screening, 0.1478% fewer for P224/1471 verification,
and 5.2924% fewer for the complete P192/73 search, including the constructor child,
receipt/executable hashing, output and independent replay. All 15 paired output
comparisons passed, including full catalogue JSON and byte-identical map output.
These are ARM Linux Callgrind instruction counts; native macOS throughput remains
unmeasured. Two external Docker interruptions are preserved with their lost counter
costs marked unmeasured; no measured campaign-total gain is claimed.

The candidate passed 60 CLI tests, 3 field-reference tests, 114 standalone algorithm
tests and the forced-software receipt-hash reference check. Release CI checks the
software fallback on all target platforms. SHA-256 values match the unchanged
repository implementation, and the compact index preserves exact curve/subgroup data.
Root compilation retains the same 643 errors; hosted CI and publication remain gated.

## Earlier field optimization and coverage round

The optimization and coverage continuation is recorded in
[`research/isogeny_coverage_optimization_20261009/RESULTS.md`](../../research/isogeny_coverage_optimization_20261009/RESULTS.md).
The final canonical catalogue contains 349 models: 204 prime, 72 binary, 60 Koblitz,
4 subfield and 9 extension models. Of these, 174 support short-Weierstrass prime-field
construction through 640 bits and 146 have complete public generator data for
independent replay. Montgomery/Edwards, binary and extension map adapters remain
unresolved; their known-order screening and catalogue records are available.

| Continuation check | Evidence and result |
| --- | --- |
| Dedicated CLI tests | 58 passed |
| Independent field reference tests | 3 passed, including unchanged arithmetic and BigUint references through P521 |
| Existing standalone algorithm suite | 114 passed |
| Wider-field full enumeration | P384/19 and P521/7 each produced two independently certified maps with 20 scalar-transport checks per map |
| High-degree replay through the optimized field | P224/1471 passed |
| Canonical model certificates | 349 verified, 0 invalid, 0 unsupported |
| Matched instruction-count panels | Five rounds per panel; screening used 96.814% fewer instructions, P224 verification used 6.307% fewer, and the full P192/73 search used 1.185% more |

These measurements use native ARM Linux inside the local Docker VM with the same
compiler, release flags and frozen 345-model inputs. They do not establish native
macOS wall-time gains. Candidate lists, certificates and full-search map bytes match
the baseline. The report manifest binds the study, source and canonical artifacts.
The root library still has the previously recorded 643 compiler errors; hosted CI
and release publication remain pending.

## Initial packaging validation

The CLI and release integration were implemented in commit `5bcc591f3`. The underlying
research sources and canonical curve records remained unchanged in that packaging phase. Native release builds
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
