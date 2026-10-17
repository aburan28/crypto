#set page(paper: "a4", margin: 17mm)
#set text(size: 9pt)
= P-521 fixed-width field validation

*Status:* local arithmetic correctness observation. No timing, point,
scalar, or side-channel conclusion follows. The source hashes, exact commands,
compiler stderr, and compressed raw rows are in
`research/ecc_arithmetic_20261013/p521_validation/README.md`.

== Implementation boundary

The new field uses nine little-endian 64-bit words, a 2^521 = 1 end-around
reduction, 81 schoolbook word products per multiply, and a fixed inversion
chain for 2^521 - 3. It is exposed as a module but no P-521 point path calls
it yet. The source was uncommitted relative to `9f588c4f4` during the run.

== Observation

On Apple ARM64 with Rust 1.93.1, the Rust source emitted 1,000 deterministic
rows for canonicalization, add, subtract, multiply, and invert. A separately
written Python integer checker matched every row modulo 2^521 - 1. The raw
rows are retained compressed with a SHA-256 digest. The checker invocation
before compression failed because the `.gz` file was missing; its repeat
after compression passed.
The same source cross-compiled to an x86-64 Linux object with Rust 1.98.0;
that object was retained but did not execute on an x86-64 host.
The source-bound Cargo harness passed 13 focused tests, including this
module's BigUint differential test, after a retained test-only compile fix.

#image("p521_flow.svg", width: 100%)

== Remaining gates

The repository release test, compiled-code side-channel review, P-521 point
integration and scalar correctness, and isolated A/A plus matched A/B timing
panels remain. This diagram represents the observed field check only.
