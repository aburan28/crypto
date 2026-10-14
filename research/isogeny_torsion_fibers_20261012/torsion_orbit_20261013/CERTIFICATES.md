# Cubic equivariance polynomial certificates

These records support the written derivation in [REPORT.md](REPORT.md).
The [native issuer/checker](../../../examples/isogeny_torsion_identity.rs)
uses `identity.certificate/v1`. It certifies two polynomial equalities
over the schema's field `F_(2^61−1)`. The statements do not encode the
elliptic-curve group law or the three projective exceptional points;
those obligations are explicit in the report and finite verifier.

| Record | Statement | Sealed id | Sealed + extra points | Checker result |
| :--- | :--- | :--- | ---: | :--- |
| [Group-law reduction](group-law.identity.certificate.json) | `x³(y+1)(y−2)−(y−1)³(y+1)−x³(y−3)+(y−1)²(y²−x³−1)=0` | `IDC1h41adc1f4e1c5522c` | 64 + 64 | Accept; `log₂` false-accept bound −7510.79 |
| [Möbius conjugacy](mobius-conjugacy.identity.certificate.json) | `((y−3)+(y+1))((y−1)−(y+1))=((y−3)−(y+1))(y−1)` | `IDC1ha0c3d83775e32754` | 64 + 64 | Accept; `log₂` false-accept bound −7680 |

The issuer source SHA-256 is
`696985397490407d24b5881deaa80f9cd69cc80dfdc4e85b863c93bcbd0f9d8c`.
The two JSON SHA-256 values, in table order, are
`7fc79d51fc464dfed11dd9b591b0b105ea8684de0d7ca2de74c298ac6917b31e`
and `70a3ed318d1e0fdb43f8caf60a9af7fad293c47ed9b1c4ea1e48adbebd1c4f79`.
The records' `proof` fields are null; the partial Lean declarations and
their limits are in [PROOF_MAP.md](PROOF_MAP.md).

The initial issuance used a direct native `rustc` build from the unchanged
repository release library while the full library test build held Cargo's
target lock. Its binary SHA-256 was
`93aba9f6b4c582a133d098a713dd25acfb840fb46b29891a106e7e41ba74b792`.
The prebuilt release-library archive SHA-256 was
`776cbf9cf27a903db2a824ffa7b8d2d832e2d9c16c372969b365c48dad39bcd5`.
The direct build and issuance commands were:

```sh
rustc --edition=2021 -O examples/isogeny_torsion_identity.rs -L dependency=target/release/deps --extern crypto_lib=target/release/deps/libcrypto_lib-2f648d1af79d9c46.rlib --extern serde_json=target/release/deps/libserde_json-aa1f2095037ebc36.rlib -o /private/tmp/isogeny-torsion-identity-20261013
/private/tmp/isogeny-torsion-identity-20261013 issue group_law > research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/group-law.identity.certificate.json
/private/tmp/isogeny-torsion-identity-20261013 issue conjugacy > research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/mobius-conjugacy.identity.certificate.json
/private/tmp/isogeny-torsion-identity-20261013 verify research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/group-law.identity.certificate.json
/private/tmp/isogeny-torsion-identity-20261013 verify research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/mobius-conjugacy.identity.certificate.json
```

The direct `rustc --test` execution and
`cargo test --release --example isogeny_torsion_identity` each passed
both issuer/replay and mutation-control tests. The standard Cargo
regeneration also matched both sealed JSON files byte for byte (`cmp`
exit zero):

```sh
cargo run --release --example isogeny_torsion_identity -- issue group_law > /private/tmp/isogeny_group_law_cargo.json
cargo run --release --example isogeny_torsion_identity -- issue conjugacy > /private/tmp/isogeny_mobius_cargo.json
cmp group-law.identity.certificate.json /private/tmp/isogeny_group_law_cargo.json
cmp mobius-conjugacy.identity.certificate.json /private/tmp/isogeny_mobius_cargo.json
```

The two `cmp` commands are shown from this study directory. The Cargo
commands ran from the repository root.
