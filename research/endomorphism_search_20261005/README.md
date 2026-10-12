# Endomorphism decision-tree integration

The [native command](../../tools/endomorphism-search/README.md) implements exact
order search, concrete bounded map discovery, subgroup checks and scalar
decomposition. Its machine-readable decisions retain unsupported cases and
unmeasured gains as unresolved. This PR is tooling and correctness validation.

The immutable [original archive](historical/python-source-and-results.zip)
preserves the delivered Python source, tests, raw pilot timings, scope notes
and SHA-256 manifest. These are historical Python measurements; this native
integration does not promote their two exploratory gain labels to accepted
native performance results. No whole-pipeline ECDLP cost or speedup is supplied.

The [protocol](PROTOCOL.md) was written before native replay. The final replay
reproduces all 500 order summaries, all four synthetic curve workloads and all
15 maps, and verifies 64,589 accelerated scalar results. The final Cargo run
passes 21 search/CLI tests plus the 10 shared SHA-256/SHA-224 known-answer tests.
Formatting and strict Clippy checks pass with Rust 1.98. Tampered evidence and
an altered manifest are both rejected without a success result.

`results/` retains the first replay's wrong extraction-path failure, the first
successful native validation, and the final validation after charging return
isomorphism constants and adding malformed-lattice rejection. Earlier receipts
are preserved and are not results for the final source. Final replay receipts
bind the native module hashes. Cargo log durations are test diagnostics, not
benchmarks. `native-manifest.sha256` binds the source, protocol, archive and
receipts; verify it from the repository root with:

```sh
sha256sum --check research/endomorphism_search_20261005/native-manifest.sha256
```

The canonical IC scoreboard receives no new performance figure: this work
supplies no IC/rho run, reference re-match, IC relation model or end-to-end
measurement. Curve registration and catalog/UI publication are separate from
the bounded correctness replay. The tool reports generated identities and
their unregistered status rather than mutating catalog data.
