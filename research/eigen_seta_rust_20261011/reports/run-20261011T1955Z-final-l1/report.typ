#set page(paper: "a4", margin: (x: 1.8cm, y: 2cm))
#set text(size: 10pt)
#set table(stroke: 0.4pt, inset: 5pt)
= Eigen-Séta v2 L1 KEM round-trip record

Run `run-20261011T1955Z-final-l1`, protocol frozen 2026-10-11 19:55 UTC.
Source: repository commit `cd47963f9` imported into base `7f7672d01`, with
the scalar-recovery and bounded public-input validation edits. This run
preceded the later explicit 3-power membership check and PKE layer.

Command from `eigen_seta/`: `cargo run --release --example run -- L1 1`.
Rust 1.99.0, release build, AMD EPYC 7763 shared cloud VM, x86-64 and five
logical CPUs. No core pinning, isolation or frequency control. `OsRng`
seed not captured. Requested and completed: one key and one KEM round trip.
Exit 0; shared keys equal. No failures.

#table(columns: (1fr, 1fr, 1fr),
  [*KeyGen (s)*], [*Encaps (s)*], [*Decaps (s)*],
  [31.050], [0.3032], [11.464],
)

Each value is a measured wall-clock sample (n=1); no spread or uncertainty
interval was recorded. Public key and KEM ciphertext: 1344 bytes each.
Sampled secret degree: 277 bits. Exact output: `records/run_output.txt`.

#figure(image("timings.svg", width: 100%), caption: [Wall times, one sample per operation, log seconds. Run run-20261011T1955Z-final-l1.])

This result demonstrates one functional round trip, not security,
constant-time behavior or a comparative performance conclusion.
