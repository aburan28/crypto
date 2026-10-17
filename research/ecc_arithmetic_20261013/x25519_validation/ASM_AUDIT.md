# ARM64 and x86-64 assembly checks for X25519

This is an instruction-flow observation on one compiler and target, not a
general constant-time certificate. Source, assembly and checker artifacts
have `internal` provenance.

The first fixed-radix source used a Rust bitmask to select between a field
value and that value minus `p` in `Fe::to_bytes`. On the Apple ARM64 T6041
host, Rust 1.93.1 with `-C opt-level=3 -C target-cpu=native` turned that
selection into two `b.lo` branches whose predicate depends on field limbs.
The exact prior field source is
[`f25519_before_barrier.rs`](f25519_before_barrier.rs), SHA-256
`960c2a11e03b81d3b9da3563e611d8b1a55118cd555f73bbc312d401e5dd52c0`.
The original tool output was not retained as a file; the
[`reconstruct_pre_barrier.py`](reconstruct_pre_barrier.py) script reproduced
that exact source hash and generated the
[`pre-fix assembly`](arm64_x25519_before_barrier.s), SHA-256
`c469eb42ebb1a42d5afef353112cd65548e2d9defc6014b62aafdf3fd9f3ebd5`.

The current source makes the all-zero/all-one selection mask opaque with
`core::hint::black_box`. Its
[`assembly`](arm64_x25519.s), SHA-256
`7a84f5057baccba102afda214f2a3dd6a3d2a7d5fa77860993ca46a5c0b0b328`,
has no conditional branch in `Fe::mul` and nine `b.ne` branches in the
X25519 function. Each of those nine is a fixed-count loop branch in the
inversion addition chain; the two limb-dependent `b.lo` instructions are
absent. [`arm64_branch_audit.py`](arm64_branch_audit.py) lists and checks
those counts from the retained assembly. The checker result was:

```text
pre-fix Fe::mul: 0; X25519: 9 b.ne plus 2 b.lo
current Fe::mul: 0; X25519: 9 b.ne
```

Both builds used the same public X25519 driver and the same target options:

```sh
rustc --edition=2021 -C opt-level=3 -C target-cpu=native --emit=asm research/ecc_arithmetic_20261013/x25519_validation/ecc_x25519_before_check.rs -o research/ecc_arithmetic_20261013/x25519_validation/arm64_x25519_before_barrier.s
rustc --edition=2021 -C opt-level=3 -C target-cpu=native --emit=asm research/ecc_arithmetic_20261013/x25519_validation/ecc_x25519_check.rs -o research/ecc_arithmetic_20261013/x25519_validation/arm64_x25519.s
python3 research/ecc_arithmetic_20261013/x25519_validation/arm64_branch_audit.py
```

The source change did not alter the retained 1,000 field or 100 X25519
correctness outputs.

The current source was also cross-compiled to x86-64 Linux assembly with
Rust 1.98.0. The retained [`x86-64 assembly`](x86_64_x25519.s) has SHA-256
`32fbebcaa70f818aa268d03fa69266e957363d3a63e99fe61a4a714adce15c22`.
[`x86_64_branch_audit.py`](x86_64_branch_audit.py) finds zero conditional
jumps in `Fe::mul` and nine in `x25519`: one unsigned `jb` on the fixed
255-step ladder counter, and eight `jne` on fixed inversion-chain counters.
No secret-derived branch was identified in those function bodies for this
build. The exact commands were:

```sh
rustup run stable rustc --edition=2021 --target x86_64-unknown-linux-gnu -C opt-level=3 --emit=asm research/ecc_arithmetic_20261013/x25519_validation/ecc_x25519_check.rs -o research/ecc_arithmetic_20261013/x25519_validation/x86_64_x25519.s
python3 research/ecc_arithmetic_20261013/x25519_validation/x86_64_branch_audit.py
```

The x86-64 object compiled but did not execute locally. These scans are
specific to the standalone Rust builds and compiler options. They do not
establish constant-time behavior for the full library, memory/cache behavior,
other compilers or CPUs, or Ed25519. Those checks remain open.
