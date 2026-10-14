# PR #1680 integration provenance

The integration preserves current Eigen-Séta KEM/PKE implementation and maintained reports from parent `2e506d3472c4bf28aad15e7ddb3fae5a40220fbd`. PR #1680's incoming source is `74bec88bdc6fb43724493d452d72df43ca3ea735`.

The distinct paper, KAT/robustness tests, public-filter example and `run-20261011T1955Z-verify` records are retained. The five conflicting incoming claim/report versions are archived byte-for-byte under [archives/pr-1680-74bec88bd](archives/pr-1680-74bec88bd/ARCHIVE.md), with hashes and original Git blobs in its [manifest](archives/pr-1680-74bec88bd/MANIFEST.json). Maintained reports and all later PKE/recovery run records remain in place.

| Historical evidence | Source binding | Interpretation after integration |
|---|---|---|
| Benchmark and first Toy16 instruction-count run | `19c641d12` | Existing historical measurements, unchanged |
| Second Toy16/L1 verification run | `d504af28e0a262cd97cdca0a58c8168e90ce6f60` | Source-bound correctness and instruction-count evidence for that build |
| Public degree-filter rate | `72051ac67` | Public-filter-only diagnostic; excludes later key-generation rejection stages |
| Paper and conflicting incoming summaries | `74bec88bdc6fb43724493d452d72df43ca3ea735` | Historical research manuscript and claims; see [paper provenance](../../papers/eigen-seta/INTEGRATION_PROVENANCE.md) |

The later scalar-recovery and public-input changes have not repeated those historical instruction-count measurements. Passing current correctness tests does not establish current constant-time behavior, independent cryptanalysis, or the paper's security conjectures. The maintained `SIDE_CHANNELS.md` limitations remain applicable, including its warning that temporary buffers are not all zeroised.

The imported live Rust source is formatted, and current `tests/common/mod.rs` adds `pub mod testrng;` for its deterministic known-answer tests. Expected Toy16, Toy32 and L1 KAT fingerprints are preserved. Original incoming test/example bytes remain in the archive. Validation logs for the combined branch are recorded separately; this note does not relabel historical measurements as measurements of the integration.
