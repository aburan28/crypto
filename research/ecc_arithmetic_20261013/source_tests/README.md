# Focused ECC source test harness

This small offline Cargo crate imports the repository's actual `f25519.rs`,
`x25519.rs`, `curve25519.rs`, and `p521_field.rs` by path. It compiles their
`#[cfg(test)]` modules without building the repository's large library test
binary. The `ecc` re-export in `src/lib.rs` preserves the original
`crate::ecc::curve25519` reference used by the X25519 differential test.

From the repository root:

```sh
CARGO_BUILD_JOBS=1 CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/ecc-arith-source-tests cargo test --offline --locked --release --manifest-path research/ecc_arithmetic_20261013/source_tests/Cargo.toml -- --nocapture
```

The first attempt compiled dependencies but failed with Rust E0282 in the new
P-521 test-only modulus expression. Its [stderr](cargo_attempt1.stderr) and
[stdout](cargo_attempt1.stdout) are retained. Adding an explicit `BigUint`
type and `usize` shift fixed that compile error; the production field
operations were unchanged. The second attempt completed with exit code 0:
[stdout](cargo_attempt2.stdout) reports 13 passed tests, zero failed, and
[stderr](cargo_attempt2.stderr) records the optimized build and test binary.

The frozen lockfile has SHA-256
`b6f076779e65baa75ae036ae8aa7c20a7d75070034825f01ce4dd4dc9a76d8c0`;
`src/lib.rs` has SHA-256
`d1fec8356bf5bc8283695accd16b5e8edb3c67f3a6e5fb29c9442532e599e65c`.
The successful stdout has SHA-256
`d70480304a36ba246cc4d99e0efd1d9bed0f804bdf40648720eec72ce755801f`.
This focused harness verifies the imported source modules and local unit
tests. It does not replace the full library test, check downstream point
integration, or measure performance.
