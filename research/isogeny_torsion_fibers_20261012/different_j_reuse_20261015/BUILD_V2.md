# V2 native build and public input freeze

The [V1 protocol](PROTOCOL_V1.md) was committed as `80992a205` before
implementation. Its scalar-ordering leak was corrected in
`077d34bc2` before constructing batch inputs or running a batch arm.
The [V1B correction](TARGET_ORDER_CORRECTION_V1B.md) preserves the
original failed rule for audit.

The source was formatted with `rustfmt examples/nonisomorphic_extension_isogeny.rs`.
`cargo check --release --example nonisomorphic_extension_isogeny` passed
in 1m 28s; its exact stdout/stderr are retained in
[BUILD_V1.stdout](BUILD_V1.stdout) and [BUILD_V1.stderr](BUILD_V1.stderr).
`cargo build --release --example nonisomorphic_extension_isogeny` then
passed; see [BUILD_V2.stdout](BUILD_V2.stdout) and
[BUILD_V2.stderr](BUILD_V2.stderr). No Python research path ran.

The input-only command was:

```sh
target/release/examples/nonisomorphic_extension_isogeny freeze-batch \
  research/isogeny_torsion_fibers_20261012/different_j_full_image_20261015/PUBLIC_V2.json \
  research/isogeny_torsion_fibers_20261012/different_j_reuse_20261015/PUBLIC_BATCH_V2.json \
  research/isogeny_torsion_fibers_20261012/different_j_reuse_20261015/CHECKER_BATCH_V2.json
```

Its [receipt](FREEZE_V2.jsonl) reports 60 distinct source and 60
distinct mapped points. The [public solver input](PUBLIC_BATCH_V2.json)
contains no answer or scalar fields; the separate
[checker file](CHECKER_BATCH_V2.json) contains the known answers for
later replay. `RUN_SPEC_V2.json` pins source, binary, input, checker,
script and freeze-receipt hashes before any batch measurement.
