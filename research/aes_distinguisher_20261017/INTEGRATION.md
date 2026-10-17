# Native AES laboratory integration

Classification: software engineering. The current user authorized this task to
proceed in a standalone checkout after reviewing the prepared integration patch.
The Conductor service remains unreachable; no reservation, task lease, or change
to the repository's standing coordination rules is claimed.

## Delivered scope

`tools/aes-discover` is a standalone native Rust crate. Its sources, tests,
report-rendering example and Cargo files mirror
`cryptanalysis@2fe7393e06e58c1f3504bdbf57255a1fb36017d0` byte for byte. The engine
and original protocol were frozen before scientific execution at
`cryptanalysis@8a3ee1d85cbb32b0fb78aa8b92cc363e3e0caffd`.
The integration base is `crypto@c6d64436a33359320ad9d5b42c4acebee2968b9a`.

The CLI discovers finite differential endpoint predicates and validates a frozen
plan on fresh fixed keys using disjoint pairs, a matched lazy PRP, conservative
bounds and family correction. It also provides exact S-box tables, bounded
trail proposals, integrals, component quartet ablations, schedule/boomerang
diagnostics and key-disjoint datasets. It uses table-based reference AES and
public synthetic fixtures. See [the tool guide](../../tools/aes-discover/README.md).

The nine advisory cards in [hypotheses.json](hypotheses.json) retain their
original implementation-status boundaries. Their source identifiers resolve
through [sources.json](sources.json). Each also has a corresponding blocked
[backlog proposal](BACKLOG.md) with controls, falsification, missing capabilities,
proposed native interfaces and explicit admission requirements. The complete literature assessment,
explanatory diagram and PDF remain in the linked companion study. No scientific
trial, claim promotion or held official campaign arm is performed by this
integration. Complete quartet attacks, SAT certificates, exhaustive toy hulls,
model training and GPU work remain separate implementation obligations.

## Verification

The eleven mirrored implementation/Cargo files compare byte for byte; a
deliberately altered temporary reference is rejected. Native formatting,
release Clippy with warnings denied, twelve library and two integration tests
pass, including all three FIPS vectors and 24 generated OpenSSL comparisons.
No cryptanalytic trial is run by these checks. Receipts are in [checks](checks/).
The parity script compares both file sets and file contents; documentation and
integration-only scripts may differ by repo.
CI checks this mirror against the exact companion commit, then runs formatting,
release Clippy, tests and shell syntax checks.

All nine blocked backlog cards pass the repository's native `backlog validate`.
For this structural check, the unchanged `src/bin/backlog.rs`, its card/cairn
modules and the library's unchanged SHA-1/SHA-256 modules were compiled in an
isolated temporary Cargo build, avoiding a second full-library build. The card
validator and cryptographic hash source files were not rewritten. No cairn
command was invoked. The required full `cargo test --release --locked --lib`
does not reach tests because unchanged base ecbench/Gaudry sources fail to
compile. The [named baseline problem](../../problems/aes-integration-root-library-baseline.md)
links the full error log, source identities and exact acceptance conditions.
The first attempt's missing sparse-checkout dependency was materialized before
the substantive second attempt. Root regression success is not claimed.

| Requirement | Status and evidence |
| --- | --- |
| Native AES binary and component/endpoint commands | Implemented; fourteen native tests, formatting and release Clippy pass |
| Companion implementation/Cargo parity | Verified across eleven files; altered-reference rejection also passes |
| Nine repository backlog proposals | Structurally validated; all blocked on explicit implementation/protocol/checker obligations |
| Required root-library regression suite | Blocked by unchanged base compile errors; named problem and receipts linked above |
| Complete quartet/SAT/toy-hull/learned/GPU extensions | Proposed obligations, not implemented by this integration |
| New scientific trials or official campaign transitions here | Not performed |
| PR merge | Requires all applicable current-head CI and review gates; no bypass authorized |

No ECDLP measurement, throughput improvement, curve identity or accepted frontier
changes here; canonical curve graphs, IC scoreboards and performance
leaderboards require no update. The original suite's measured counts, graphs
and PDF remain immutable in the companion evidence tree.
