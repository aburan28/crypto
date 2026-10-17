# X25519 fixed-radix correctness validation

This is a correctness record for the next arithmetic implementation in the
[ECC portfolio](../PLAN.md). It does not measure time or establish a speedup.
The public X25519 interface now uses five radix-2^51 field limbs, the
[RFC 7748](https://www.rfc-editor.org/rfc/rfc7748.html) Montgomery ladder,
branchless conditional swaps in source, and a fixed inversion addition chain.
The existing `BigUint` field and Ed25519 implementation are separate paths.

## Frozen inputs and outputs

- Repository HEAD before this uncommitted change:
  `9f588c4f455d98072aa678314e50a85b673a0f94`.
- Original `src/ecc/x25519.rs` SHA-256:
  `0d3faadadcc0e29a770641d7fadf530e62e4e58f5e531eee9be057d2b6eef44e`.
- New `src/ecc/f25519.rs` SHA-256:
  `b5a6169fff9fcb8c3f62becc43f6dc83a30c8d536e947a138562978507406a8d`.
- New `src/ecc/x25519.rs` SHA-256:
  `f945190011682bd6e282d63617a18e670d9abfa14f06e1653e23ce5fa1928867`.
- [Field rows](field_rows.txt) SHA-256:
  `cd975d3bcc1488522e55aa173bccdd601486bde22e78b9998015b4cb38c30932`.
- [X25519 rows](x25519_rows.txt) SHA-256:
  `8f518fdab646592e84064179a2524bff4ff09f945c1b301d09ea1fdb286f40cf`.

The machine was Apple ARM64 T6041 on Darwin 25.6.0; Rust was
`rustc 1.93.1 (01f6ddf75 2026-02-11) (Homebrew)`, and Python was 3.13.1.
The working tree was dirty with the source change and this study directory.
No wall-time values are reported from the contended host.

From the repository root, the exact validation commands were:

```sh
rustc --edition=2021 -C opt-level=2 research/ecc_arithmetic_20261013/x25519_validation/ecc_f25519_check.rs -o /private/tmp/ecc_f25519_check
/private/tmp/ecc_f25519_check > research/ecc_arithmetic_20261013/x25519_validation/field_rows.txt
python3 research/ecc_arithmetic_20261013/x25519_validation/ecc_f25519_check.py
rustc --edition=2021 -C opt-level=2 research/ecc_arithmetic_20261013/x25519_validation/ecc_x25519_check.rs -o /private/tmp/ecc_x25519_check
/private/tmp/ecc_x25519_check > research/ecc_arithmetic_20261013/x25519_validation/x25519_rows.txt
python3 research/ecc_arithmetic_20261013/x25519_validation/ecc_x25519_check.py
```

The field checker passed 1,000 deterministic operand pairs, including zero,
`p-1`, `p`, `p+1`, a high-bit input and all-ones inputs, for addition,
subtraction, multiplication and inversion against Python integer arithmetic
modulo `2^255 - 19`. The second checker passed 100 deterministic full X25519
outputs against an independently written integer Montgomery ladder. Its Rust
driver also passed the two RFC 7748 section 5.2 vectors before emitting the
rows. The source-level Rust tests additionally include a differential check
against the original `BigUint` ladder. A focused release Cargo build was
[stopped under host load](CARGO_RELEASE_ATTEMPT.md) before it emitted a test
result, so this record does not call the full-library test passed.
The smaller [source test harness](../source_tests/README.md) then compiled
these same modules and passed 13 focused tests, including the source-level
field and full-ladder differential tests. It does not replace the full
repository library test.

An [ARM64 assembly check](ASM_AUDIT.md) found two limb-dependent branches in
the first canonical encoder and removed them with an opaque selection mask.
The retained pre-fix and current assembly, exact source reconstruction and
branch scanner make that narrower observation reproducible. An x86-64 Linux
cross-compile of the same source found zero conditional jumps in field multiply
and nine fixed-counter jumps in X25519; that object was not executed here.

## Scope and next gates

The checkers cover arithmetic equality for the generated inputs, including
masked high bits and noncanonical coordinates. They do not constitute a
constant-time machine-code audit or exhaustive proof. The next gates are the
repository's release library tests, compiled-code side-channel review, a
matched X25519 baseline/candidate panel with an A/A control on isolated
ARM64 and x86-64 hosts, and complete protocol regression tests.
Ed25519, the NIST prime curves, ECC2K point/walk operations, Koblitz scalar
recoding, and Metal/CUDA/HIP remain separate requested tracks.
