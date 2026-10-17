# L1 plaintext-encryption round trip

Protocol frozen on 2026-10-11 at 19:56 UTC, before execution. Source:
`eigen_seta/` imported from repository commit `cd47963f9` on base
`7f7672d01`, with scalar-recovery correction, public-input validation and
the new ChaCha20-Poly1305 KEM–DEM PKE layer. Command from `eigen_seta/`:
`cargo run --release --example pke_run -- L1`. The example encrypts a fixed
39-byte message with fixed associated data using `OsRng` for KeyGen,
encapsulation and the AEAD nonce. Random seed not captured.

Acceptance: exit 0 and print `correct true` after decrypting the message.
Record KeyGen, Encrypt and Decrypt wall seconds, public-key and ciphertext
sizes; retain failure. One sample only on an AMD EPYC 7763 shared cloud VM
with rustc 1.99.0, no isolation or performance comparison. Toy16 tampering
and associated-data rejection are covered separately by the integration
test.

Result: exit 0, `correct true`. Measured wall-clock samples: KeyGen 25.044 s,
Encrypt 0.3897 s, Decrypt 11.208 s. Public key 1344 bytes; complete PKE
ciphertext 1411 bytes for a 39-byte message. One sample per operation, no
spread or interval, no isolation or performance comparison. The exact output
is in `records/run_output.txt`; `timings.svg` is embedded in `report.pdf`
from `report.typ`. This shows one authenticated plaintext round trip, not a
security or constant-time proof.
