#set page(width: 10.0in, height: 7.0in, margin: 0.45in)
#set text(font: "New Computer Modern", size: 11pt)
#let navy = rgb("1d4e70")
#let ink = rgb("1f2937")
#let bar(w) = rect(width: w, height: 0.24in, fill: navy, radius: 2pt)

= Cold pairing cost / same-target rho

*V2 frozen pairing transfer · 2026-10-15 05:55 UTC · shared Darwin arm64, L0*

Each bar is the median of three *paired full current-thread CPU ratios*.
The pairing total includes construction, candidate search, both pairings,
the multiplicative-field solve, separate source-scalar replay, and the
post-solve pairing checks. Blue bar length: 1 inch = 1,000× rho.

#v(0.24in)
#grid(columns: (1.2in, 5.15in, 2.45in), gutter: 0.08in,
  [*Field / order*], [*Median ratio*], [*Three-round min–max*],
  [$p=61$, $r=19$], [#bar(2.78181in) #h(0.08in) *2,781.81×*], [1,741.21–4,344.32×],
  [$p=271$, $r=61$], [#bar(0.69815in) #h(0.08in) *698.15×*], [697.24–1,046.10×],
)

#v(0.16in)
#grid(columns: (1.28in, 5in), [], [#line(length: 5in, stroke: 0.8pt + ink)])
#grid(columns: (1.28in, 5in), [], [#grid(columns: (1in, 1in, 1in, 1in, 1in, 1in),
  [0], [1,000], [2,000], [3,000], [4,000], [5,000×])])

#v(0.2in)
*Source:* `SUMMARY_V2.json` derived from 30 native raw commands in
`runs_v2/`, commit `b4e2bff9cb`. Rho A/A and all six independent scalar
replays passed. The min–max interval is a three-run range, not a confidence
interval. CPU ratios on this shared host are descriptive.
