#set page(width: 8.5in, height: 11in, margin: 0.7in)
#set text(font: "New Computer Modern", size: 10pt)
= Run: cubic torsion-equivariant coordinate

2026-10-14 05:51 UTC. Protocol commit `b9f3259739`. The native Rust
verifier checked the proposed `chi(P+T)=1/(1−chi(P))` identity on all
rational points of four source curves and checked the order-three
projective action. All assertions passed.
An initial replay failed when the temporary binary was absent. After
recompilation, raw and log files matched the original byte-for-byte.

#table(
  columns: (1fr, 1fr, 1.5fr),
  align: center,
  table.header([*p*], [*r*], [*Rational points checked*]),
  [617], [103], [618],
  [2,633], [439], [2,634],
  [9,677], [1,613], [9,678],
  [38,561], [6,427], [38,562],
)

#v(0.14in)
#image("../../equivariant.svg", width: 100%)

Source, binary, raw hashes and exact commands are in `README.md`.
The coordinate does not factor through the isogeny because it takes
three distinct values on the three kernel points. This is an exact
coordinate check, without a decomposition oracle or complete DLP cost.
