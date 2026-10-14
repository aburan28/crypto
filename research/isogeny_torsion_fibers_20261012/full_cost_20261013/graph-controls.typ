#set page(width: 10.2in, height: 5.2in, margin: 0.48in)
#set text(font: "New Computer Modern", size: 10pt)
#let fail = rgb("f9d7cc")
#let pass = rgb("dcefdc")
#let cell(value, ok) = box(width: 1.45in, inset: 6pt, fill: if ok { pass } else { fail })[#value]

= Duplicate rho control across the retained run versions

Each cell is the full-work duplicate-rho / matched-rho ratio. Both copies
performed the same counted group operations on the same target and seed.
V1 uses elapsed wall time; V2 and V3 use current-thread CPU time.
V2/V3's frozen admissible A/A band is 0.8–1.25.

#v(0.20in)
#grid(columns: (1.3in, 2.4in, 2.4in, 2.4in), gutter: 0.08in,
  [*subgroup r*], [*V1 wall*], [*V2 CPU*], [*V3 CPU*],
  [103], [#cell([420.987×], false)], [#cell([0.827×], true)], [#cell([0.810×], true)],
  [439], [#cell([0.771×], false)], [#cell([0.924×], true)], [#cell([0.918×], true)],
  [1,613], [#cell([2.071×], false)], [#cell([0.754×], false)], [#cell([0.819×], true)],
  [6,427], [#cell([5.690×], false)], [#cell([0.817×], true)], [#cell([0.824×], true)],
)

#v(0.25in)
Pink marks a failed or visibly unstable control. Green marks the frozen
CPU ratio gate only. The residual 8–19% A/A CPU imbalance limits the
interpretation of the smaller source-versus-destination differences.
Source: `summary-v1.tsv`, `summary-v2.tsv`, `summary-v3.tsv` and retained
raw records beside them.
