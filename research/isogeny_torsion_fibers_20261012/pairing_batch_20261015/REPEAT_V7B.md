# V7B same-input repeat after controlled release suite exits

Frozen before the repeat. Use the exact [V7 spec](RUN_SPEC_V7.json),
source SHA-256 `55e23e6c7c361ce295ea7c3313f82ad4734ef900eb7a71bdafc9ac78dbbcea6e`
and binary SHA-256 `d4eb7405211a1a433f570663334c413f8abdce6195a6f43c7c3f31fa087a67ab`.
Keep the 66 public targets, H selection, nine-entry source baby table,
rounds 0–2, five prefixes, rotated arm order, 144-step rho cap and 32
retries unchanged. Start only after the controlled `cargo test --release
--lib` process exits; preserve its final stdout/stderr and status first.

Hypothesis: direct source BSGS remains the lower complete cold CPU arm at
`N=64` and `N=66` without that test-suite process. Reference: the lower
of same-target direct source BSGS and direct prime-field rho. A valid
round requires all 66 scalar answers from all three arms, charged
input/setup/target/verification work, primitive and RSS receipts, and
independent replay of all `3*3*66=594` outputs. Any failed arm or scalar
mismatch censors its round. Stop after three rounds and one replay. Write
new `*_V7B` files; never overwrite V7. This remains shared-host L0, so
other host activity and a portable counted-`S` conversion stay unresolved.

From the repository root, execute rounds sequentially, writing each to
`COMPACT_BSGS_401_R<round>_V7B.jsonl` and `.stderr`, concatenate the three
JSONL rows to `COMPACT_BSGS_401_V7B.jsonl`, and run the separate checker:

```sh
target/release/examples/torsion_pairing_batch compact-bsgs-run \
  research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/INPUTS_V4.jsonl 2 \
  "$round" research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/PUBLIC_401.json \
  research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/SELECTION_401_V1.json
target/release/examples/torsion_pairing_batch compact-verify \
  research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/INPUTS_V4.jsonl 2 \
  research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/PUBLIC_401.json \
  research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/ANSWERS_401.json \
  research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/COMPACT_BSGS_401_V7B.jsonl
```

Cell index two and all named files were checked against V7 before the
repeat. The answer file is read only by `compact-verify`.
