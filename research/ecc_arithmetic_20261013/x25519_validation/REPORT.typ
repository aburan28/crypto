#set page(paper: "a4", margin: 17mm)
#set text(size: 9pt)
= X25519 fixed-radix correctness validation

*Status:* observations from a local correctness run. No performance or
constant-time conclusion is drawn here. The exact commands, source hashes,
drivers and raw rows are in
`research/ecc_arithmetic_20261013/x25519_validation/README.md`.

== Implementation boundary

The public X25519 entry point now uses five radix-2^51 limbs, a dedicated
square, a fixed inversion chain, and the RFC 7748 ladder with source-level
branchless swaps. The original `BigUint` field remains a separate reference
path; Ed25519 is unchanged. The source was uncommitted relative to repository
HEAD `9f588c4f4` during this run.

== Observations

On Apple ARM64 T6041 with Rust 1.93.1 and Python 3.13.1, the actual Rust
field source emitted 1,000 deterministic operand-pair rows. The separate
Python integer checker matched addition, subtraction, multiplication and
inversion modulo `2^255 - 19` in all rows. The actual Rust X25519 source
passed two RFC 7748 section 5.2 vectors, then emitted 100 deterministic
full-operation outputs. A separately written Python integer ladder matched
all 100. Inputs and outputs are retained as `field_rows.txt` and
`x25519_rows.txt`.

#image("validation_flow.svg", width: 100%)

An ARM64 assembly scan found two field-limb-dependent `b.lo` branches in the
first encoder. After the opaque-mask fix, the retained assembly has zero
conditional branches in field multiply and nine fixed-loop `b.ne` branches
in the X25519 function, with no `b.lo` there. The exact pre-fix and current
assembly and scanner are recorded in `ASM_AUDIT.md`.
An x86-64 Linux cross-compile with Rust 1.98.0 found zero conditional jumps
in field multiply and nine fixed-counter jumps in X25519. It did not run on
an x86-64 host; the retained assembly and scanner are also in `ASM_AUDIT.md`.

== Missing checks

The release library Cargo test waited behind the shared build lock, then was
stopped after more than twenty minutes of compilation as host load reached
256; it emitted no test result. The host was heavily contended, so no
wall-time comparison was taken. A Modal x86-64 validation launch was
rejected by automatic approval review before source export because the new
eight-file bundle lacks specific authorization. Machine-code timing behavior,
isolated ARM/x86 A/A and A/B comparisons, and complete protocol regression
remain open. The diagram is a flow of observed checks, not a speedup graph.
The smaller source-bound Cargo harness subsequently passed 13 focused tests
across the imported 25519 and P-521 modules; the full library test remains open.
