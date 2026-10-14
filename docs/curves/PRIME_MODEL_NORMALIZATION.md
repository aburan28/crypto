# Explicit prime model normalization

The catalogue records short Weierstrass coefficients. Historical source
records may use another coordinate model. The native `registry_models`
tool consumes the declared form and preserves the original record, source
hash and coordinate map in `sources/prime-models.json`. The Python registry
harvester consumes those native field results and rejects stale or missing
normalization metadata.

Run from the repository root:

```sh
cargo run --release --locked --manifest-path tools/isogeny-cli/Cargo.toml --bin registry_models
cargo run --release --locked --manifest-path tools/isogeny-cli/Cargo.toml --bin registry_models -- --check
python3 scripts/test_prime_model_normalization.py
```

For the declared model `y² = x³ + a x² + b x`, the native translation is
`u = x + a/3`, `v = y`. Its short coefficients are
`A = b − a²/3`, `B = 2a³/27 − ab/3`, in the recorded prime field.
The tool checks the source's original ICV1 and slug after translation.
This adapter supports fields up to 64 bits and checks primality with the
existing deterministic native Miller–Rabin test; larger fields are rejected.
Native tests check the coordinate bijection and both counts on small fields,
and count both actual recorded AGM cases in their original and short forms.

During integration, two catalogue rows had interpreted the AGM source's
quadratic coefficients as short coefficients. Their erroneous fingerprints
were `56b6783a1485` and `7ee99c50b352`. Those rows denoted different models
and are not mathematical aliases for the measured curves. The corrected
entries use the frozen producer's original identities,
`icv1-fp21-tm700-d62ba3a6` and `icv1-fp21-tm700-07d957b2`, and preserve
their recorded order of 1,049,284. `sources/prime-model-corrections.json`
retains the deprecated extraction labels and their correction provenance.

The historical AGM JSON, producer and reports remain unchanged. The three
larger round4 prime models keep their recorded orders with an unknown
certification status; this correction does not promote those claims.
