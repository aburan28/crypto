# Native source and build custody before the relation run

Host: Darwin arm64, macOS 26.6, shared L0. The frozen
[protocol](PROTOCOL_V1.md) was committed as `c8fabd428` before any new
relation run. The V2 solver input was derived by the native command:

```sh
target/release/examples/nonisomorphic_extension_isogeny freeze-full-image \
  research/isogeny_torsion_fibers_20261012/different_j_20261015/PUBLIC_V1.json \
  research/isogeny_torsion_fibers_20261012/different_j_full_image_20261015/PUBLIC_V2.json
```

Its [freeze receipt](FREEZE_V2.jsonl) records 4,575 source points, 1,525
image points, three preimages for every image point, 21 selection
candidates and eight columns. The original public fixture SHA-256 is
`fa41958038b0fe01ec1664ca69acf3c44cc67af9feb38a7c8dbf6bf6a3041b65`;
the new public fixture SHA-256 is
`75fded656d5fba32418e8a6231fb81e46412b8ae088b496dc4a1427e9fd14967`.

The first `cargo build --release --example nonisomorphic_extension_isogeny`
attempt exited 101 at compile time: `serde_json::json!` reached the
default macro recursion limit while expanding the enlarged relation
record. The source was not frozen for that failed build and it produced
no updated executable or experiment output. Adding the crate-level recursion limit
of 256 resolved it. The corrected build exited zero. `rustfmt --check`,
`zsh -n` for the run script, and `jq empty` for both specs exited zero.
`cargo test --release --example nonisomorphic_extension_isogeny` exited
zero and compiled the example test harness; it contains zero built-in
tests. The frozen relation runs and separate arithmetic replay supply
the substantive correctness evidence.

Final pre-run source SHA-256:
`3e7452d3d7c0cbb421c60da6d90b2ba42cddeea7e1cb90d51423511d9f33f28e`.
Final pre-run executable SHA-256:
`54254777a50c0d86ab1a87a25b6941426595d28e5872b5797f0a13364ac5c294`.
The exact [V2 run spec](RUN_SPEC_V2.json) pins these bytes and the
[run script](run_v2.zsh) before relation measurement.
