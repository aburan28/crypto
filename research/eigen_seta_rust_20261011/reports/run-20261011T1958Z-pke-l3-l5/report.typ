#set page(paper: "a4", margin: (x: 1.8cm, y: 2cm))
#set text(size: 10pt)
#set table(stroke: 0.4pt, inset: 5pt)
= Eigen-Séta v2 L3/L5 plaintext-encryption record

Run `run-20261011T1958Z-pke-l3-l5`, protocol frozen 2026-10-11 19:58 UTC.
Source: same-repository commit `cd47963f9` imported into base `7f7672d01`,
with scalar-recovery fix, public-input validation and the ChaCha20-Poly1305
KEM–DEM PKE layer.

Commands from `eigen_seta/`, run concurrently:
`cargo run --release --example pke_run -- L3` and
`cargo run --release --example pke_run -- L5`. Rust 1.99.0, release build;
AMD EPYC 7763 shared cloud VM, x86-64 and five logical CPUs. No core
pinning, isolation or frequency control. Other local validation ran during
the same period. `OsRng` seeds not captured.

#table(columns: (auto, auto, auto, auto, auto),
  [*Set*], [*KeyGen (s)*], [*Encrypt (s)*], [*Decrypt (s)*], [*Correct*],
  [L3], [186.427], [1.1051], [47.138], [true],
  [L5], [298.906], [2.2990], [104.223], [true],
)

Requested and completed: one 39-byte authenticated plaintext round trip
per set, with fixed associated data. Both processes exited 0; no failures.
Each time is a measured wall-clock sample (n=1); no spread or uncertainty
interval was recorded. Public key and complete PKE ciphertext sizes were
1992 / 2059 bytes at L3 and 2632 / 2699 bytes at L5. Exact output:
`records/run_output.txt`.

#figure(image("timings.svg", width: 100%), caption: [Wall times, one sample per operation and set, log seconds. Run run-20261011T1958Z-pke-l3-l5.])

These results demonstrate functional authenticated plaintext round trips.
Concurrent work on a shared VM prevents an isolated performance comparison.
They do not establish security or constant-time behavior.
