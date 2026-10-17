#set page(width: 8.5in, height: 11in, margin: 0.7in)
#set text(font: "New Computer Modern", size: 10pt)
= Rational three-torsion orbit screen

2026-10-13 local / 2026-10-14 UTC. Conductor T-2611. Source-linked full
report: `REPORT.md`; frozen protocol: `PROTOCOL.md`; native counts:
`raw-v1.tsv`; code: `screen.rs`; proof map: `PROOF_MAP.md`.

The tested question was whether three rational lifts of each of twelve
factor-base points add new two-sum relation targets after torsion correction
and projection to the prime-order subgroup. They did not: support was
identical to the twelve-point control in all four groups.

#table(
  columns: (1fr, 1fr, 1.2fr, 1.2fr, 1.2fr),
  align: center,
  table.header([*Subgroup r*], [*12 points*], [*36 lifts*], [*36 distinct*], [*Lifts / 12*]),
  [103], [64], [64], [103], [1.000],
  [439], [74], [74], [351], [1.000],
  [1,613], [75], [75], [559], [1.000],
  [6,427], [78], [78], [642], [1.000],
)

#v(0.15in)
#image("graph-support.svg", width: 100%)

All values are exact distinct subgroup support counts. `screen.rs` checked
torsion order, point distinctness, source and destination equations,
isogeny and cofactor invariance, subgroup order, and point-addition against
a complete scalar index. The conditional quotient-support implication built
in Lean 4.34.1 with no axioms.

#pagebreak()

== A nonfactorable equivariant coordinate

On `E:y²=x³+1`, let `T=(0,1)` and `chi=(y+1)/(y−1)`. The proposed
action is `chi(P+T)=f(chi(P))` with `f(z)=1/(1−z)` on the projective line.
The three kernel points map to `1, infinity, 0`, so this coordinate does
not factor through the degree-three quotient. The action has order three.

#image("equivariant.svg", width: 100%)

The native Rust verifier checked both the `+T` and `+2T` identities on
all 51,492 rational source points across the four fields, plus `f³=id`
on all 51,492 projective-line points. Its inputs and counts are in
`EQUIVARIANT_PROTOCOL.md` and `equivariant-v1.tsv`. This is a coordinate
identity result; a complete extension-field oracle comparison is pending.

For generic `P`, the addition slope `(y−1)/x` gives
`y(P+T)=y−2−(y−1)³/x³=(y−3)/(y+1)` on `y²=x³+1`.
Conjugating this rational map by `chi` gives `f`. The two
denominator-cleared polynomial steps have native identity certificates
`IDC1h41adc1f4e1c5522c` and `IDC1ha0c3d83775e32754` beside the run.
The points `O,T,2T` use direct projective checks. The Lean proof map
contains a compiled conditional orbit-support theorem and the order-three
linear-action consequence; its curve-level bridge is still open.

The earlier four-size matched rho and full cold CPU rows are in
`../full_cost_20261013/REPORT.md`. This run made no new complete-solve
measurement, so it has no speed or `S` ratio. The next experiment is the
same-target extension-field torsion-symmetric decomposition oracle in
`NEXT_GATE.md`, compared with direct source, source-orbit and invariant
quotient controls. It must price construction, relation collection,
linear algebra, recovery, verification and matched rho on the same target.
