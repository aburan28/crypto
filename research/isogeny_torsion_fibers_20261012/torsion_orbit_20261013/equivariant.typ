#set page(width: 10in, height: 4.6in, margin: 0.42in)
#set text(font: "New Computer Modern", size: 11pt)
#let blue = rgb("2875a8")
#let orange = rgb("d66a3a")
#let box(content, color) = rect(width: 1.2in, height: 0.48in,
  fill: color.lighten(88%), stroke: 0.8pt + color, radius: 3pt,
  align(center + horizon, content))

= A cubic coordinate that distinguishes the three kernel points

For each `p` in 617, 2633, 9677, 38561: source `E/F_p:y²=x³+1`,
destination `E'/F_p:y²=x³−27`, `T=(0,1)`, and the verified degree-three
quotient `phi:E→E'`.

#v(0.14in)
#grid(columns: (1.4in, 1.3in, 0.55in, 1.3in, 0.55in, 1.3in, 0.55in, 1.3in),
  gutter: 0.05in,
  [*source orbit*], box([$O$], blue), [→ +T],
  box([$T$], blue), [→ +T],
  box([$2T$], blue), [→ +T], box([$O$], blue))

#v(0.18in)
#grid(columns: (1.4in, 1.3in, 0.55in, 1.3in, 0.55in, 1.3in, 0.55in, 1.3in),
  gutter: 0.05in,
  [*chi values*], box([$1$], orange), [→ f],
  box([$infinity$], orange), [→ f],
  box([$0$], orange), [→ f], box([$1$], orange))

#v(0.18in)
`chi=(y+1)/(y−1)` and `f(z)=1/(1−z)`; the three distinct values show
that `chi` does not factor through `phi`. In contrast,
`phi(O)=phi(T)=phi(2T)=O` on the destination. The finite verifier checked
`chi(P+T)=f(chi(P))` on all 51,492 rational source points across the four
fields and checked `f³=id` on all 51,492 projective-line points.
