#set page(margin: 1.5cm)
#set text(size: 9.5pt)
= Run report: geometric-V cost, m = 3 curve-vs-presentation (n = 19), PinX3
2026-10-11 from UTC 21:50 · base 4b6e09b78 + this PR · 4-vCPU shared container (wall clock indicative). Records under `experiments/koblitz_m3/` and `experiments/koblitz_presentation/`. No failures; all solves BSGS-verified.

== Per-call cost (measured ms)
#table(columns: 6,
  [*n, l*], [*V*], [*e₃ serial*], [*e₃ clmul*], [*e₃ values/call*], [*roots/call*],
  [17, 5], [mono], [19.39], [9.61], [8192], [1019.2],
  [17, 5], [geom], [28.71], [9.29], [7790], [966.0],
  [19, 6], [mono], [211.24], [106.47], [65536], [16409.6],
  [19, 6], [geom], [503.19], [122.77], [65536], [16345.6],
)
PinX3, geometric V, 40 probes: 0.052 ms (n = 17) and 0.154 ms (n = 19), against 9.70 and 123.60 ms for e₃ on the same probes; identical relations (2/40, 1/40).
#image("per_call.svg", width: 80%)

== m = 3 curve-vs-presentation, n = 19 (128/128 cells; 4180 relations / 51200 probes)
#table(columns: 3,
  [*metric*], [*ICC_curve*], [*95 % CI*],
  [yield per |F|³], [0.098], [[0, 0.224]],
  [ln µs per call], [0.011], [[0, 0.137]],
  [work per call], [0.000], [[0, 0.051]],
  [ln projected solve], [0.000], [[0, 0.069]],
)
#image("icc_n19_m3.svg", width: 75%)

== PinX3 full solves (measured)
#table(columns: 5,
  [*n, a₂, l*], [*V*], [*seconds*], [*probes*], [*unknowns*],
  [23, 1, 8], [mono], [0.5], [418], [118],
  [23, 1, 8], [geom], [0.7], [435], [118],
  [31, 0, 10], [mono], [33.5], [5612], [525],
  [31, 0, 10], [geom], [45.9], [5876], [507],
)
#image("full_solves_pin.svg", width: 70%)
