# P-521 fixed-width field validation

This record covers the new nine-word arithmetic modulo `2^521 - 1` in
[`src/ecc/p521_field.rs`](../../../src/ecc/p521_field.rs). It is a correctness
check of canonicalization, addition, subtraction, multiplication and a fixed
inversion chain. It contains no timing or point-arithmetic result.

## Frozen source and environment

- Repository base before this change: `9f588c4f455d98072aa678314e50a85b673a0f94`.
- Field source SHA-256: `a5f0e5f163f2a0aef53bb4f2d3f6adab0b04db279f8b3fe7f8a4d6e2b2416380`.
- Rust driver SHA-256: `0959e00d5e8750bc9be50b1c1116b91c9a879ed50a3e054567a50e86e2e3ace1`.
- Python checker SHA-256: `3edafd4e14246e89dc91ba235cb8fad6c041d1ad8fa491c537427b9f7d932ce4`.
- [Compressed raw rows](p521_rows.txt.gz) SHA-256: `6fb5e7512443266211059f8a9f7c1e578762be6699282c14f037a1191f02c0a4`.
- Compiler stderr SHA-256: `ac30ca74a93e30fbfb633e7cbe737f21de5e2556050de04b3098f88aced53a85`.
- [x86-64 Linux object](p521_field_x86_64.o) SHA-256:
  `bf09265607057fcb0727095d5679dbe62eda815987a26f47a48fe104e4a25da9`.
- Host: Apple ARM64, Darwin 25.6.0; Homebrew `rustc 1.93.1
  (01f6ddf75 2026-02-11)`; Python 3.13.1. The working tree was dirty with
  this implementation and the X25519 work.

From the repository root, reproduce with:

```sh
rustc --edition=2021 -C opt-level=2 research/ecc_arithmetic_20261013/p521_validation/ecc_p521_check.rs -o /private/tmp/ecc_p521_check 2> research/ecc_arithmetic_20261013/p521_validation/rustc.stderr
/private/tmp/ecc_p521_check > research/ecc_arithmetic_20261013/p521_validation/p521_rows.txt
gzip -n research/ecc_arithmetic_20261013/p521_validation/p521_rows.txt
python3 research/ecc_arithmetic_20261013/p521_validation/ecc_p521_check.py
```

The seeded driver emitted 1,000 input/output rows. The Python checker uses
independent arbitrary-precision integers and passed all 1,000 for each
operation. Rows include zero, all-one input bytes and the exact field modulus;
the rest use a fixed xorshift seed `0x41b76d8e579d1531`. The first checker
invocation preceded compression and failed solely because its expected `.gz`
file did not exist; it was rerun after compression and passed. No correctness
conclusion is drawn from that failed invocation.

The inverse uses `f[k] = a^(2^k - 1)`, builds `f[512]` and `f[7]`, then
computes `a^(2^521 - 3)` with a fixed sequence. Inversion of zero maps to
zero; callers must reject zero when a mathematical inverse is required.
The source also cross-compiled to an x86-64 Linux object with Rust 1.98.0:

```sh
rustup run stable rustc --edition=2021 --target x86_64-unknown-linux-gnu --crate-type lib -C opt-level=3 --emit=obj src/ecc/p521_field.rs -o research/ecc_arithmetic_20261013/p521_validation/p521_field_x86_64.o
```

That object was not executed on x86-64.
The [source test harness](../source_tests/README.md) separately compiled this
module's `#[cfg(test)]` code and passed its BigUint differential test as part
of 13 focused tests. The first harness build exposed and retained a test-only
type inference error, fixed before the successful rerun.

## Scope and next gates

This sample check is not an exhaustive proof or machine-code side-channel
certificate. The module has not yet been selected by a P-521 point path.
Next gates are the repository release test, compiled-code review, point and
scalar differential tests, and isolated matched benchmarks with an A/A
control. The [diagram](p521_flow.svg) and [PDF report](REPORT.pdf) show the
implemented and observed validation path.
