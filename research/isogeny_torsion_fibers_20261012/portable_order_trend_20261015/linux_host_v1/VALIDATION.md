# Independent-host result validation

- GitHub Actions [run 38831567671](https://github.com/aburan28/crypto/actions/runs/38831567671)
  completed its independent Linux job successfully. Each of its three
  Rust replay receipts says `pass`, with 84 child records and 1,785
  correct scalar outputs per cell.
- The native `different_j_host_performance` analyzer checked 252/252
  Linux-versus-Darwin arithmetic vectors, 5,355 scalar outputs and
  37,061 relation witnesses before writing `COMPARISON_V1.json`. Its
  final source SHA-256 is
  `d3e87b05c90508f31873e5c58de7521450700d3f554c0bca297b42833cd5686c`;
  the local Darwin release binary SHA-256 was
  `cef30532655c11d7cf4ecb024f7245d665778c3ea5ede9b99b4316ad7fcef824`.
  Adding graph round ranges did not change the comparison JSON; a
  byte-for-byte `cmp` of the pre/post-renderer output passed.
- The raw archive was extracted and all 1,533 fresh relative-path SHA-256
  entries passed `sha256sum -c`. The original hosted manifest's
  self-reference issue is described in [README.md](README.md).
- `cargo check --release --example different_j_order1201_host_performance`
  and `cargo test --release --example different_j_order1201_host_performance`
  passed. `cargo test --release --example different_j_host_performance`
  passed after the graph update; these example test binaries contain
  zero embedded unit tests. Functional validation is the independent
  checker and arithmetic comparison above.
- `actionlint` passed on both replay workflows. The repository's existing
  PR workflow-budget checker reported 0 violations across 364 workflows
  after both jobs were set to 25 minutes. `bash -n` passed on the new
  order-1,201 Linux script. Typst compiled the report and results PDFs;
  both rendered as one legible A4 page, and the SVG-derived graph PNG
  rendered at 1,406 × 703 pixels (150 ppi).
- The repository release-library suite was run earlier in this PR and
  remains a red receipt due sandbox socket denial and shared-host
  timeouts, as described in [the parent validation record](../validation/VALIDATION_V1.md).
  No `src/` library file changed in this follow-on.
