# Configured scalar worker validation

The requested change is a separately configured exact legacy N83 worker. It
keeps the production library, earlier worker examples, and frozen binaries
unchanged. Both scalar kernels use the existing Fold (`Gf2_128`) field and
production `StrongRhoG` walk, in two monomorphizations of one executable.

The strict configuration requires `backend: "fold"`,
`scalar_kernel: "reference" | "barrett"`, `lanes: 64 | 128 | 256`,
`dp_bits: 6 | 8 | 10`, `escape_fruitless: false`, and
`step_cap_factor: 2000`. Unknown, missing, duplicate, and out-of-grid fields
are rejected. Input is exactly one point array containing two decimal coordinate
strings. The worker neither accepts nor opens a known-answer sidecar.

The reference scalar is the production `u64` implementation. At this 53-bit
subgroup its multiplication uses exact `u128` remainder, without the library's
small-modulus floating approximation. The candidate's hot multiplication uses
an exact integer reciprocal and the upper half of an exact 256-bit product.
Its hot addition uses one overflow-aware subtraction on canonical residues.
Subtraction, inversion, and infrequent modular multiplication delegate to the
production exact routines. Uniform sampling delegates every method to
`UniformInt<u64>`, including both one-shot methods used by `gen_range`.

The native Rust checks cover:

- Exact high-product arithmetic against `BigUint`, on limb edges and 16,384
  deterministic random full-width operand pairs.
- Exact reciprocal, quotient error at most one, residual below twice the
  modulus, corrected product, overflow-aware addition, and delegated rare
  operations, across 276 moduli and edge plus random canonical operands. The
  moduli include the exact subgroup, powers of two, and `u64` boundary values.
  Both a correction and an additive carry are explicitly exercised.
- Scalar conversions and inversion delegation.
- Constructed and one-shot inclusive/exclusive uniform samplers, with matching
  RNG positions after every sequence of draws.
- An exact N83 solve on a correctness-only target `[12345]G`, outside the timed
  corpora. Jump points, jump scalars, sample canonicalizations, final scalar,
  all outcome counters, charges, and remaining RNG position match reference.
  Exact primitive arithmetic and identical randomness supply the deterministic
  trajectory argument; the tests additionally observe equal full-solve counts.

The first compile failed because two static `UniformInt` calls needed explicit
`<u64>` type arguments. The failure is retained in
`new_worker_initial_compile_failure.log`; the corrected release suite passes
all seven tests. These correctness checks make no speed claim.

Reproduction from repository root:

```sh
export RUSTUP_HOME=/workspace/n83_evaluation/toolchain/rustup
export CARGO_HOME=/workspace/n83_evaluation/toolchain/cargo
export CARGO_TARGET_DIR=/workspace/n83_evaluation/build
export CARGO_BUILD_JOBS=4
export PATH=/workspace/n83_evaluation/toolchain/rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin:/usr/bin:/bin
cargo test --release --locked --example n83_rho_speed_eval
cargo build --release --locked --example n83_rho_speed_eval
```

No additional compiler flags or CPU-specific `RUSTFLAGS` are used.
Conductor was attempted before edits and is unavailable (`conductor: command
not found`); no `.conductor/WORKFLOW.md` was present.

The worker source SHA-256 at build is
`bd2e5b2ed4ca4229f187b3f1c5c613633d9cb773200c961fdd8cbe8fb1d65241`.
The release executable SHA-256 is
`997f6102564cd0ebc72ee46e17e834f3d5970ed129c86d9c56b9274d5405a2ec`
(3,029,536 bytes). The campaign must bind its copied executable hash rather
than assume this build path remains immutable.

The reciprocal proof and independent audit are maintained separately in the
study's `audit/` directory. Holdout measurements and selection rules are owned
by the campaign, not inferred from these tests.
