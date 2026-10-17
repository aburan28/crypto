# Native quintic roots in implicit binary domains

Source-only checkpoint. Benchmark reports and raw performance evidence remain local pending explicit publication approval. No8x or fastest-tool claim is made.

The monic-quintic extractor computes gcd(H,L_D) without expanding the implicit domain polynomial, caches characteristic-two quotient squares, and partitions affine subspaces using lower annihilators. It returns every root within the declared domain, with exact multiplicities. It does not solve the remaining multivariate coefficient-family search.

Run cargo run --release --manifest-path research/quintic_domain_roots_20261017_v1/Cargo.toml -- controls from the repository root. The independent tiny-domain oracle uses Hasse derivatives rather than shared polynomial division. The wide component is binary-only; prime/odd domains and the full matched F4/F5/F6 gate remain separate.

All mathematical and solver code is Rust. The included earlier native field source is an exact frozen dependency and provides the prime/binary/odd-extension predicate prototype; its diagnostics are not a qualified comparison. No report payload or private inputs are included.
