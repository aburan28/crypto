#set page(width: 21cm, height: 10cm, margin: 1.5cm)
#set text(size: 9pt)
#align(center)[
  #text(size: 14pt, weight: "bold")[P-256 trapdoor search, round 3 — run report]
  #v(0.3em)
  Run id run-20261012T2330Z · 2^20 deterministic X9.62 candidates · `cargo run --release --example p256_trapdoor_round3 -- 1048576`
]
#v(1em)
#table(
  columns: (auto, auto, auto, auto),
  inset: 5pt,
  [*Probe*], [*Published*], [*Honest (n = 1,048,170)*], [*corrected p*],
  [QR(b)], [true], [0.500], [1.000],
  [cubic residuosity of b], [false], [0.332], [1.000],
  [Hamming weight of b], [134], [128.00 ± 8.00], [1.000],
  [CF max quotient of b/p], [1275 (11 bits)], [≥ published: 0.5092; max seen 30 bits], [1.000],
  [Hamming weight of r], [125], [127.50 ± 7.99], [1.000],
  [QR(r)], [true], [0.500], [1.000],
)
#v(1em)
#align(center)[#image("graph.svg", width: 92%)]
#v(0.5em)
Recorded values only. No outlier flagged; seed-preimage phrase probe: 1000 phrases, no hit (expected).
