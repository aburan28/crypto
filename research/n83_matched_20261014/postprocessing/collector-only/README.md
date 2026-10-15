# Separate post-hoc collect-only extraction

The original frozen controller's `collect` command failed after the full V3
independent replay passed. Its exact error was `resume raw records do not
reproduce the frozen observation`: it reconstructs observations after parsing
and reserializing floats, then requires exact Value equality. The original
controller, failure logs, observations and measured files remain unchanged.

This separate native Rust program performs collection only. It never launches
solvers, reconstructs `observe`, recalculates a clock, or changes an existing
verification flag from false to true. Its only subprocess is `/usr/bin/sha256sum`.
It preserves the frozen collector's schema, row splitting, seed/rank fields and
metadata. Additional provenance explicitly identifies this extraction and the
failed original collection attempt.

Before creating output, it verifies independently pinned receipt, protocol,
manifest, checker source/executable/helper and original runner identities. The
V3 receipt must report a complete PASS:218 scalar replays,104 IC witnesses,
five A/A semantic pairs,32 complete processes and zero invalid/missing answers.
All336 unique file snapshot hashes are rechecked, including every observation,
execution, original controller source, binary and frozen input. The whole
campaign tree must contain only sealed files. Extraction reads and parses the
same cached bytes that passed their SHA check.

Every process and target identity is checked against the prospective32-process
matrix and independent audit rows. Only already-verified observations are
eligible; independent verification is then added. The resulting input has
32 process rows,192 cohort target rows and26 auxiliary singleton/A/A rows.
`serde_json` is pinned with `arbitrary_precision` to preserve original decimal
number tokens when cloning sealed observations, avoiding another f64 round trip.
Outputs and the extraction receipt must be new files outside the sealed campaign.

The checked config is [config.json](config.json). Its SHA pins include the
retained failed original command's stdout, stderr and exit2 files. It is a local
execution record; reproductions use new absolute paths and separately reviewed
pins without changing the frozen record.

Build and test, using the session's pinned Rust toolchain:

```sh
CARGO_TARGET_DIR=/workspace/n83_matched/collector_only_build CARGO_BUILD_JOBS=1 cargo test --release --locked --offline --manifest-path research/n83_matched_20261014/postprocessing/collector-only/Cargo.toml
CARGO_TARGET_DIR=/workspace/n83_matched/collector_only_build CARGO_BUILD_JOBS=1 cargo build --release --locked --offline --manifest-path research/n83_matched_20261014/postprocessing/collector-only/Cargo.toml
```

After source/config review, the authorized execution is:

```sh
/workspace/n83_matched/collector_only_build/release/n83-matched-collect-only --config research/n83_matched_20261014/postprocessing/collector-only/config.json --out /workspace/n83_matched/postprocess-20261015T0018Z/statistics_input.json --receipt /workspace/n83_matched/postprocess-20261015T0018Z/collector_extraction_receipt.json
```

The extraction receipt binds source, manifest/lock, executable, config, audit,
every copied observation and normalized output by SHA256. The output explicitly
states `original_frozen_collect_succeeded:false`. This program does not establish
an additional scalar proof, timing reproduction, host comparison or ledger
promotion; it consumes the fully verified V3 receipt and its immutable file seal.
