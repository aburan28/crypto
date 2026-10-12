#set page(paper: "a4", margin: (x: 1.8cm, y: 2cm))
#set text(size: 10pt)
#set table(stroke: 0.4pt, inset: 5pt)
= Eigen-Séta v2 L3/L5 round-trip record

Run `run-20261011T1946Z-l3-l5-fix`, protocol frozen 2026-10-11 19:46 UTC.
Source: repository commit `cd47963f9` imported into base `7f7672d01`, with
the local constant-time scalar-recovery edge-case correction. The runs do
not cover the subsequent public-input validation edit.

Commands from `eigen_seta/`, run concurrently:
`cargo run --release --example run -- L3 1` and
`cargo run --release --example run -- L5 1`. Rust 1.99.0, release build;
AMD EPYC 7763 shared cloud VM, x86-64, five logical CPUs. No core pinning,
isolation or frequency control. `OsRng` seeds not captured.

#table(columns: (auto, auto, auto, auto, auto),
  [*Set*], [*KeyGen (s)*], [*Encaps (s)*], [*Decaps (s)*], [*Correct*],
  [L3], [112.612], [1.3059], [47.275], [1/1],
  [L5], [259.034], [1.9862], [107.221], [1/1],
)

Requested and completed: one key and one encapsulation–decapsulation per
set. Both processes exited 0; no failures. Each time is one measured
wall-clock sample (n=1); no spread or uncertainty interval was recorded.
Public key and ciphertext were 1992 bytes for L3 and 2632 bytes for L5.
Sampled secret-degree bit lengths: 409 and 536. Exact output:
`records/run_output.txt`.

#figure(image("timings.svg", width: 100%), caption: [Wall times, one sample per operation and set, log seconds. Run run-20261011T1946Z-l3-l5-fix.])

This demonstrates functional round trips. Concurrent execution on a shared
VM prevents an isolated performance comparison. These samples do not
establish security or constant-time behavior.
