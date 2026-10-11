# L3/L5 plaintext-encryption round trips

Protocol frozen on 2026-10-11 at 19:58 UTC, before execution. Source:
`eigen_seta/` imported from repository commit `cd47963f9` on base
`7f7672d01`, with scalar-recovery correction, bounded and validated public
input handling, and the ChaCha20-Poly1305 KEM–DEM layer. Commands from
`eigen_seta/`, run concurrently:
`cargo run --release --example pke_run -- L3` and
`cargo run --release --example pke_run -- L5`. Each run encrypts and decrypts
the same fixed 39-byte message with fixed associated data; KeyGen,
encapsulation and nonce use `OsRng` (seeds not captured).

Acceptance per set: exit 0 and print `correct true`. Record the example's
KeyGen, Encrypt and Decrypt wall seconds, sizes and failures. One sample per
set on a shared AMD EPYC 7763 cloud VM with rustc 1.99.0, no isolation or
performance comparison. Toy16 tampering and associated-data rejection are
covered by an integration test.

Both runs exited 0 with `correct true`. On the shared VM, the measured
single-sample wall times and sizes were:

| Set | KeyGen (s) | Encrypt (s) | Decrypt (s) | pk (bytes) | PKE ct (bytes) |
|---|---:|---:|---:|---:|---:|
| L3 | 186.427 | 1.1051 | 47.138 | 1992 | 2059 |
| L5 | 298.906 | 2.2990 | 104.223 | 2632 | 2699 |

Each size used a 39-byte message. No spread or uncertainty interval was
recorded. The two runs were concurrent, and other local validation was also
active; no isolated performance comparison is claimed. Exact output and
commands are in `records/run_output.txt`; `timings.svg` is embedded in
`report.pdf` from `report.typ`. This is correctness evidence, not a security
or constant-time proof.
