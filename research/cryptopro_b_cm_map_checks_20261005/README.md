# CryptoPro-B public CM map correctness checks

The Rust example replays the public degree-155 elliptic CM endomorphism,
composed of degree-5 and degree-31 maps followed by an isomorphism. Its
independent affine group arithmetic checks agreement with the supplied
subgroup scalar on 39 known-input cases, including 32 frozen random scalars
and seven boundary cases. Some boundary cases generate the same point.
It also checks 16 additivity pairs, negation, and
`omega^2 - omega + [155] = 0` on each input.

From the repository root, run:

```sh
cargo test --release --example cryptopro_b_cm_map_checks
cargo run --release --example cryptopro_b_cm_map_checks > /tmp/cryptopro_b_cm_native.json
```

The verifier also checks exact component polynomial curve identities,
source/target chaining, kernel-square denominators, final isomorphism,
generator order, the CM discriminant identity, the conductor factor
product, the scalar CM polynomial, and the proposed polarization matrix
determinant restricted to the rational subgroup. That last check alone
does not audit the endomorphism ring or prove geometric indecomposability.

The input `endomorphism155.json` is unchanged from the supplied
`cryptopro_b_evidence.zip` and has SHA-256
`51ad7edaf9c931364be2513a1163a95474858381152a715ad9140648fad6a042`.
The example checks this hash and parses large integer tokens exactly.

`evidence/legacy/` preserves the original Python source snapshot, README,
and recorded result without alteration. The legacy result reports 39
known-input cases, 16 additivity pairs, and zero failed checks. Its scalars
were generated with Python RNG seed `202610059`; the native replay reads
those frozen scalars instead of substituting another random stream.
The Python snapshot is historical evidence. The executable verifier and
CI replay use Rust and the repository's existing Cargo dependencies.
CI uploads its native JSON report as `cryptopro-b-cm-native-replay`.

These checks cover the elliptic forward map. No explicit genus-two curve
equation or transfer formulas are supplied, and no genus-two transfer is
verified. Relation collection, relation solving, and total computational
cost remain unmeasured. The report leaves costs and speedup null.

See [PROTOCOL.md](PROTOCOL.md) for the replay contract and acceptance gate.
