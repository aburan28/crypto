# Five-order sensitivity validation

- `rustfmt --check examples/different_j_price_sensitivity_five.rs` and
  `cargo test --release --example different_j_price_sensitivity_five`
  exited zero. The example test binary contains zero embedded unit tests;
  the functional check is the frozen native derivation and independent
  DLP replay receipts it requires.
- Before execution, source, binary and five-order trend SHA-256 values
  matched [`RUN_SPEC_V1.json`](RUN_SPEC_V1.json). The executable exited
  zero and emitted 60 unique matched rows and five N64 summaries.
- The first 48 row objects and four earlier N64 summary objects were
  byte-equal after `jq -S` canonical JSON serialization to
  [`../sensitivity/SENSITIVITY_V1.json`](../sensitivity/SENSITIVITY_V1.json).
  The fifth cell's three raw `(G_full,G_rho,rho_steps)` triples and
  thresholds are listed in [RUN.md](RUN.md).
- Typst compiled the one-page A4 [run](RUN.pdf) and [results](RESULTS.pdf)
  PDFs. Both pages were rendered and visually inspected. The graph SVG
  and 150 ppi PNG were generated from the same native SVG plot source;
  their benchmark copies have matching SHA-256 values.
- The earlier full `cargo test --release --lib` receipt remains red for
  sandbox socket denial and shared-host timeouts; no library source is
  changed by this deterministic follow-on.
