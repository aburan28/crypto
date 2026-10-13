#set page(paper: "us-letter", margin: 0.62in)
#set text(size: 9pt)
= Bounded isogeny control run v1
2026-10-13 16:30 UTC · shared Darwin arm64 · L0 · exit 2

The source, degree-three destination, scaling-isomorphism destination, and
negation rho arms verified the first five F13 targets. The six-class rho arm
failed at target index 4 after reaching the 100,000-operation walk cap.
The executable stopped before F13 indices 5–7 and all F199 workloads.
The raw file has 25 rows: 24 verified, one failed. The chart distinguishes
failure from targets never reached.

#image("verification-v1.svg", width: 100%)

Raw: `raw-v1.tsv`, SHA-256 `43b95449…613a8b87ace`.
Recorded source hash: `e74f5e16…84af6c8a27`; binary hash:
`fd995586…e04d2c`. The source was reconstructed with an identical hash;
a native rebuild reproduced raw and log bytes. The original Cargo binary
was not retained. Commands and counters are in `README.md`,
`validation/historical_replay.md`, and the raw TSV. This run does not
support a complete-panel comparison.
