#set page(margin: 1.6cm)
#set text(size: 10pt)
= Run report: m = 3 symmetric-S₄ solvers
2026-10-11, UTC 19:55 · base 4b6e09b78 · 4-vCPU shared cloud container (wall clock indicative)

== Runs
+ Full solves: e₃ enumeration, n = 17 (l = 5) and n = 19 (l = 6), 4 curves × {monomial, geometric g1}, 50 yield probes then closure. Records: `experiments/koblitz_m3/full_solves_e3.jsonl`.
+ Solver comparison at n = 17, a₆ = 1, monomial V, 4 probes, l = 3, 4, 5. Records: `experiments/koblitz_m3/solver_compare_n17.jsonl`.

== Per-call seconds (measured), n = 17
#table(columns: 4,
  [*l*], [*e₃ enumeration*], [*symmetric Gröbner*], [*S₃ chain*],
  [3], [0.000116], [0.0106], [1.61],
  [4], [0.00141], [4.10], [9.21],
  [5], [0.0216], [timeout (>300 s / 4 probes)], [45.2],
)
#image("per_call_n17.svg", width: 85%)

== Full solves (measured wall seconds; all 16 verified against BSGS)
#table(columns: 3,
  [], [*monomial V*], [*geometric V (g1)*],
  [n = 17, l = 5], [3.5, 5.3, 7.9, 3.9], [3.9, 9.7, 11.0, 14.5],
  [n = 19, l = 6], [93.6, 104.2, 73.7, 78.8], [287.3, 260.9, 168.3, 124.8],
)
#image("full_solves.svg", width: 85%)

== Failures
One cell timed out (symmetric Gröbner, l = 5): a cost, not a mathematical result.
