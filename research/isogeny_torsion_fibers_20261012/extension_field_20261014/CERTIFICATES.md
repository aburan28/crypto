# Scaling polynomial identity certificate

The exact identity

`(u²+3)(u⁴−3u²+9)=u⁶+27`

is the algebraic step used with each checked witness `u²=−3` to obtain
`u⁶=−27`. The [native issuer](../../../examples/isogeny_torsion_extension.rs)
emits the [identity.certificate/v1 record](SCALING.identity.certificate.json)
over the schema field `F_(2^61−1)`. The separately invoked
[native checker](../../../examples/isogeny_torsion_identity.rs) accepted its
64 id-bound evaluation points and 64 extra checker-seed points in
[SCALING_VERIFY.stdout](SCALING_VERIFY.stdout):

`IDC1hc0fabee63dedc307 Accept { points: 64, extra_points: 64, false_accept_log2: -7477.124799907692 }`

The certificate SHA-256 is
`2d6592a361345d6c5e842b17e95be4082006edc535694b5b61a3446f2d81d836`.
The issuer source SHA-256 at issuance is
`e03af656ab5f50544ed373ca78ded0f12d2fb1bbe571ef66d615047c481b7918`.
The certificate has a null `proof` field. A separate universally quantified
integer theorem, [`scaling_sixth_power_factorization`](lean/ScalingWitness.lean),
is proved by the pinned Lean 4 build. Its `#print axioms` reports only
Lean's standard `propext` and `Quot.sound`; no research-claim axiom or
unverified native evaluation was used. The four exact modular coefficient
witnesses report no axioms. See the [proof map](PROOF_MAP.md).

Run from the repository root:

```sh
cargo run --release --example isogeny_torsion_extension -- issue-scaling-identity > research/isogeny_torsion_fibers_20261012/extension_field_20261014/SCALING.identity.certificate.json
cargo run --release --example isogeny_torsion_identity -- verify research/isogeny_torsion_fibers_20261012/extension_field_20261014/SCALING.identity.certificate.json
```

This polynomial identity and its checks establish the coefficient relation;
the curve equation, inverse scaling and group-law scope are stated separately
in the [results report](REPORT.md). The certificate is not a proof of the
elliptic-curve group-isomorphism statement.
