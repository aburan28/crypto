# ISO-1 polynomial-record replay

This package builds the repository's existing SHA-256 and identity-certificate
modules directly from their source paths. It runs the same report checker as
the top-level example without compiling unrelated cryptanalysis modules.
It performs polynomial-record replay and mutation controls; it does not
certify census point counts or execute a formal theorem prover.

```sh
cargo run --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --bin iso1_report_check -- research/iso1_weak_classes_20261007/identity_certificates.json
cargo test --release --locked --manifest-path research/iso1_weak_classes_20261007/report_check_runtime/Cargo.toml --bin iso1_report_check
gp -q -f research/iso1_weak_classes_20261007/algebra_identity_check.gp
```

To regenerate the records, omit the final JSON argument and save standard
output to a new file, then compare it with the retained records. The issuer
uses a fixed statement and seed. Replay adds 64 domain-separated points to
the 32 recorded evaluations; a constant perturbation of either identity
must be refused by the mutation test. PARI/GP separately expands both
polynomials exactly over the integers.
