#set page(paper: "us-letter", margin: 0.58in)
#set text(size: 9pt)

= Isomorphism-only transport and six-class rho controls
2026-10-13 local · T-2372 · v3: 80/80 verified source-target answers

The same 16 public source targets were solved through source relations,
degree-three transport, scaling-isomorphism transport, negation rho, and
six-class rho. The destination pairs have subgroup orders 19 and 211.
The isogeny is φ(x,y) = ((x³+4b)/x², y(x³−8b)/x³); the isomorphism
is ψ_u(x,y) = (u²x,u³y), with u = 2 or 14. The frozen inputs and
run parameters are in `PROTOCOL.md`, `REMEDIATION_V2.md`, and
`RESTART_V3.md`; the original target set is `../raw-v4.tsv`.

The chart uses S = (group additions + doublings)/√r, including cold
relation setup, collection, descent, and verification. Each mean and
observed min–max range is over eight targets. The destination arms
include source-point cross-verification. Native map, lookup, row, and
memory work is recorded separately, so this is a counted-operation
diagnostic.

#image("relation_group_ops-v3.svg", width: 100%)

#table(
  columns: (0.55fr, 1.75fr, 1fr, 1fr),
  inset: 3pt,
  stroke: 0.3pt,
  [r], [arm], [mean S], [verified],
  [19], [rho negation], [10.725], [8/8],
  [19], [rho six classes], [16,320.893], [8/8],
  [19], [source relations], [15.371], [8/8],
  [19], [degree-three relations], [16.145], [8/8],
  [19], [scaling relations], [19.816], [8/8],
  [211], [rho negation], [5.938], [8/8],
  [211], [rho six classes], [6.342], [8/8],
  [211], [source relations], [50.186], [8/8],
  [211], [degree-three relations], [41.917], [8/8],
  [211], [scaling relations], [43.423], [8/8],
)

#pagebreak()
= Relation collection and charged rho restarts

At order 211, degree-three and scaling-isomorphism relations used 359
and 358 probes over the same eight targets, with 3,102 and 3,107
collection group operations. Their full counted totals were 4,871 and
5,046, so scaling/degree-three = 1.035927 in this unit. At order 19,
the totals were 563 and 691, ratio 1.227353. Both maps reach the same
destination model, while φ = ψ_u ◦ [m] on these rational groups
(m = 10 mod 19 and 29 mod 211). The fixed scalar streams therefore
need not visit the same destination points. These exact observations
do not isolate a degree-three-specific relation effect.

#image("rho_group_ops-v3.svg", width: 100%)

The first two run versions are retained: v1 had 24 verified rows and
one six-class cap failure before stopping; v2 had 77 verified rows and
three six-class cap failures. V3 followed the precommitted 16-seed
restart schedule. Three F13 targets required two attempts; all failed
first attempts are charged, giving six-class rho 569,129 counted group
operations across eight F13 targets versus 374 for negation rho. At
order 211 the respective totals are 737 and 690. These are fixed
implementation controls, not a best-rho calibration.

The v3 TSV and log replay byte for byte. An independent Rust checker
verified all 80 scalar answers, frozen target coordinates, map images,
phase sums, seed/restart accounting, and six-class point multipliers.
Its receipt is `independent_v3.log`. Lean checked the conditional
coefficient update, scalar-value scaling, and negation lemmas in
`SixClassCoefficients.lean`; the actual curve group-action bridge remains
open. The host was shared Darwin arm64,
L0; no admitted wall-time claim is made. Curve-specific map discovery
and exhaustive preflight checks are out of band for these target rows.
Full native-work conversion,
a base-field non-isomorphic pair, larger orders, independent-host replay,
and a multi-point rational torsion fiber remain open. Full evidence and
requirement status are in `REPORT.md`.

#pagebreak()
= Explicit transport diagram

The source and destination identities, both map degrees, and the checked
relation of the maps on each rational subgroup are shown below. The
degree-three kernel has three geometric points and one rational point.

#image("transport-v3.svg", width: 100%)
