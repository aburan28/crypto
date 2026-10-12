#set page(margin: 1.5cm)
#set text(size: 9.5pt)
= Run report: m = 3 curve-vs-presentation at n = 31 (PinX3)
2026-10-11 UTC 23:57 · pre-registered design `m3-presentation-vs-curve-n31-design-20261012.md` · 4-vCPU shared container. 128/128 cells in 872 s, all BSGS-verified, 0 inconsistent relations; 20661 relations / 256000 probes.

#table(columns: 3,
  [*metric*], [*ICC_curve*], [*95 % CI*],
  [yield per |F|³], [0.000], [[0, 0.116]],
  [ln µs per call], [0.048], [[0, 0.216]],
  [work per call], [0.228], [[0.032, 0.403]],
  [ln projected solve], [0.067], [[0, 0.219]],
  [_diag:_ |F|], [0.239], [[0.025, 0.437]],
  [_diag:_ pins per call], [0.230], [[0.035, 0.406]],
  [_diag:_ roots per pin], [0.000], [[0, 0.138]],
  [_diag:_ work per call ÷ |F|], [0.086], [[0, 0.297]],
)
#image("icc_n31_m3.svg", width: 95%)
