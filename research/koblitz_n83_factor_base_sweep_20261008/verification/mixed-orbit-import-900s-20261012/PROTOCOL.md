# Bounded N83 mixed-orbit import gate

The user authorized **up to 900 additional local seconds** after the earlier
11,700-second allowance was fully charged. This round uses all 900 seconds as
a conservative cap: 600 for a Linux offline build (450-second outer wall),
240 for import of two exact retained S3 objects (240-second outer wall), and
60 for receipt and cleanup work. There is one launch per stage and no retry.
The exact authorization record binds a fresh output directory and the purpose
`mixed_orbit_build_import_900_gate`.

The compiled worker is `examples/koblitz_n83_mixed_orbit_worker.rs` from a
clean commit. The pinned image is
`sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31`,
Linux arm64. Both stages run with no network, an immutable image and a
container cleanup audit. Build has four CPUs and 8 GiB; import has one CPU,
4 GiB and zero swap. The small K=1,182 and envelope K=2,048 compressed objects
must match SHA-256 values `eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380`
and `aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a`.
The object URIs and BLAKE3 hashes are in the parent
[2,400-second search protocol](../mixed-orbit-supervisor-source-20261012/PROTOCOL.md).
The current local copies originated in the earlier verified S3-backed base
panels and are rehashed before launch.

The import worker validates the exact manifests, prefix of small orbits in the
envelope, curve and subgroup metadata, representative and point labels, and
returns a `n83.mixed-orbit-import/v1` receipt. The outer supervisor validates
that receipt and freezes source and input hashes. Keep raw stage output,
container states, elapsed walls, final receipt, budget receipt, and any
failure. `PASS_import_only` means the wide factor-base adapter accepted the
two retained objects. No SAT trial, partial relation, graph row, matrix rank,
individual logarithm, or complete cold runtime is measured in this gate.
The existing 2,400-second search protocol remains separate. Its authorization
record was frozen against the earlier 11,700-second baseline; update that
accounting to 12,600 seconds and obtain a new run-specific grant before any
later search launch.
