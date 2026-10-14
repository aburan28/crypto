#set page(width: 10in, height: 5.8in, margin: 0.38in)
#set text(font: "New Computer Modern", size: 10pt)
#let blue = rgb("2875a8")
#let orange = rgb("d66a3a")
#let bar(w, color) = rect(width: w, height: 0.13in, fill: color, radius: 1pt)
#let row(r, a, b, wa, wb) = [
  #grid(columns: (0.65in, 1.35in, 3.8in, 0.8in), gutter: 0.07in,
    [*#r*], [12 subgroup], [#bar(wa, blue)], [#a],
    [], [36 torsion lifts], [#bar(wa, blue)], [#a],
    [], [36 distinct], [#bar(wb, orange)], [#b])
  #v(0.12in)
]

= Exact two-sum subgroup support

Four rational three-torsion, non-isomorphic degree-three pairs. Bar length is
linear in the number of distinct reachable subgroup elements; full scale
642 elements = 3.5 inches. Values come directly from `raw-v1.tsv`.

#v(0.04in)
#row([103], [64], [103], 0.349in, 0.561in)
#row([439], [74], [351], 0.403in, 1.913in)
#row([1,613], [75], [559], 0.409in, 3.047in)
#row([6,427], [78], [642], 0.425in, 3.500in)

#line(length: 7.8in, stroke: 0.6pt + gray)
#v(0.08in)
#text(fill: blue)[■] orbit-equivalent support; #text(fill: orange)[■]
independent subgroup points. Each orbit arm has 36 rational representatives
but only twelve quotient classes. Counts are exact on these finite groups;
this graph is a relation-support diagnostic, not whole-method speed.
