#set page(width: 10.5in, height: 6.0in, margin: 0.45in)
#set text(font: "New Computer Modern", size: 10pt)
#let blue = rgb("2875a8")
#let orange = rgb("d66a3a")
#let bar(w, color) = rect(width: w, height: 0.16in, fill: color, radius: 2pt)
#let row(size, srcw, dstw, src, dst) = [
  #grid(columns: (0.9in, 1.0in, 4.2in, 1.0in), gutter: 0.08in,
    [*#size*], [source], [#bar(srcw, blue)], [#src×],
    [], [transport], [#bar(dstw, orange)], [#dst×])
  #v(0.20in)
]

= Complete cold thread-CPU work versus matched negation rho

Four non-isomorphic degree-three pairs; eight public targets and five measured
rounds per size. Bars use a log₂ axis: one bit of ratio is 0.61 inches.
Exact values and paired target intervals are in `summary-v3.tsv`.

#v(0.10in)
#row([103], 1.05in, 1.09in, [3.307], [3.449])
#row([439], 1.98in, 2.01in, [9.524], [9.789])
#row([1,613], 2.80in, 2.85in, [24.068], [25.459])
#row([6,427], 3.90in, 3.88in, [83.930], [81.820])

#line(length: 7.6in, stroke: 0.6pt + gray)
#v(0.10in)
#text(fill: blue)[■] source relations · #text(fill: orange)[■] destination relations.
All 320 V3 relation records recovered and independently checked their
scalars. Rho is 1× on this axis. A/A thread-CPU ratios were 0.810–0.918;
shared-host wall timing failed its noise control. This is an L0 local CPU
observation and is not an admitted wall-clock speed claim.
