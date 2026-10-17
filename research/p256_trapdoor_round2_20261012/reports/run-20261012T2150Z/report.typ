#set page(width: 21cm, height: 12cm, margin: 1.5cm)
#set text(size: 9pt)

#align(center)[
  #text(size: 14pt, weight: "bold")[P-256 trapdoor search, round 2 — run report]
  #v(0.3em)
  Run id: run-20261012T2150Z · deterministic probes · commit: see PR · source: RESEARCH_P256_BACKDOOR_SPECULATION.md §4 ranks 1–3
]

#v(1em)

#table(
  columns: (auto, auto, auto, auto),
  inset: 5pt,
  align: horizon,
  [*No*], [*Hypothesis*], [*Probe*], [*Measured*],
  [1], [H1 coefficients], [ord(b) = ord(r) = ord(−3) mod p], [253 bits, exact; contains the 160-bit prime Q of p−1],
  [2], [H1], [ord(j) mod p], [252 bits, exact; contains Q],
  [3], [H1], [residuosity], [b, r QR non-cubes; j QR cube; −b flips the QR bit],
  [4], [H1], [Hamming weights], [b z=+0.81, r z=−0.31, j z=−0.12, n z=+4.88 (n order-derived)],
  [5], [H1], [perfect powers / nearest c·2^k], [none; gaps ≥ 2^49 bits],
  [6], [H1], [continued fractions (64 terms)], [max partial quotients 1275, 62, 107],
  [7], [H2 prime], [p−1 factorization], [2·3·5²·11·101·281·1531·65537·490719·6602113·Q₁₆₀ (complete)],
  [8], [H2], [p+1 factorization], [2^96·7·273921· 140-bit composite (partial)],
  [9], [H2], [ord_p(2)], [(p−1)/2 exactly (2 is a QR; p ≡ 7 mod 8)],
  [10], [H2], [twist embedding degrees ≥2^20 subgroups], [none for k ≤ 100],
  [11], [3 endo/isog], [GLS gap 4p−t²], [258 bits],
  [12], [H3], [rational roots of x³−3x+b], [1 (⇒ 1 rational 2-isogeny)],
)

#v(1em)

#grid(
  columns: (1fr,),
  align: center,
  image("graph.svg", width: 88%),
)

#v(0.5em)

Recorded values only; every number copied from records.json. No anomaly flags (threshold ±5σ).
