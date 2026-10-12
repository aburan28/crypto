#set page(paper: "a4", margin: (x: 1.8cm, y: 2cm))
#set text(size: 10pt)
#set table(stroke: 0.4pt, inset: 5pt)
= Eigen-Séta v2 L1 round-trip record

Run `run-20261011T1944Z-l1-fix`, protocol frozen 2026-10-11 19:44 UTC.
Source: repository commit `cd47963f9` imported into base `7f7672d01`, with
the local constant-time scalar-recovery edge-case correction. This run does
not cover the subsequent public-input validation edit.

Command from `eigen_seta/`: `cargo run --release --example run -- L1 1`.
Rust 1.99.0, release build, AMD EPYC 7763 shared cloud VM (x86-64, five
logical CPUs). No core pinning, isolation, or frequency control. `OsRng`
seed not captured. Requested: one key, one encapsulation and decapsulation;
completed: one each. Exit 0, shared keys equal; no failures.

#table(columns: (1fr, 1fr, 1fr),
  [*KeyGen (s)*], [*Encaps (s)*], [*Decaps (s)*],
  [25.788], [0.2242], [11.406],
)

Each value is a measured wall-clock sample (n=1); no spread or uncertainty
interval was recorded. Public key and ciphertext: 1344 bytes each. Sampled
secret degree: 278 bits. Source data and exact output:
`records/run_output.txt`.

#figure(image("timings.svg", width: 100%), caption: [Wall times, one sample per operation; log seconds. Run run-20261011T1944Z-l1-fix.])

The result is a functional round trip, not a security claim, timing bound,
or constant-time validation. Historical runs in the study benchmark report
were performed on another host and are not pooled with this sample.
