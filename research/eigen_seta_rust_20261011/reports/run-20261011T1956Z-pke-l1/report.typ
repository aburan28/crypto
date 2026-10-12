#set page(paper: "a4", margin: (x: 1.8cm, y: 2cm))
#set text(size: 10pt)
#set table(stroke: 0.4pt, inset: 5pt)
= Eigen-Séta v2 L1 plaintext-encryption record

Run `run-20261011T1956Z-pke-l1`, protocol frozen 2026-10-11 19:56 UTC.
Source: same-repository commit `cd47963f9` imported into base `7f7672d01`,
with the scalar-recovery fix, public-input validation, and ChaCha20-Poly1305
KEM–DEM PKE layer.

Command from `eigen_seta/`: `cargo run --release --example pke_run -- L1`.
Rust 1.99.0, release build, AMD EPYC 7763 shared cloud VM, x86-64 and five
logical CPUs. No core pinning, isolation or frequency control. `OsRng`
seed not captured. Requested and completed: one key and one 39-byte message
encryption/decryption with fixed associated data. Exit 0, plaintexts equal.
No failures.

#table(columns: (1fr, 1fr, 1fr),
  [*KeyGen (s)*], [*Encrypt (s)*], [*Decrypt (s)*],
  [25.044], [0.3897], [11.208],
)

Each value is a measured wall-clock sample (n=1); no spread or uncertainty
interval was recorded. Public key: 1344 bytes. Complete PKE ciphertext:
1411 bytes, including the KEM ciphertext, nonce and authentication tag.
Exact output: `records/run_output.txt`.

#figure(image("timings.svg", width: 100%), caption: [Wall times, one sample per operation, log seconds. Run run-20261011T1956Z-pke-l1.])

This result demonstrates one authenticated plaintext round trip, not a
security, constant-time or comparative-performance conclusion.
