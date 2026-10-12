# Exact disclosed-query audit: SAT missed two existing n19 relations

The [preregistered](PROTOCOL.md) exact group oracle classified all **137
distinct ordinary queries** recorded by the [20-job recovery
pilot](../generic-backend-recovery-pilot/RESULTS.md). It replayed all 16
source-bound worker admissions before classification and retained the four
no-report timeouts as unknown. The [complete per-query result](RESULT.json)
contains each scalar, original PDP verdict, exact existence label and an
independently re-added witness when one exists. Its SHA-256 is
`a7119332e6b50dfd9bceb3b15dcf6d62bf579ff065bb5026e49ec0a5db09329c`.
No worker, target or query was rerun or generated.

The exact oracle uses each cell's verified **geometric** factor base before
cofactor projection. It indexes all unordered point pairs by their group
sum, subtracts every third base point from `[a]G`, and re-adds the first
found triple with independent Python field and curve arithmetic. That
exhausts every three-summand decomposition, including repeated points. The
geometric set can include points that project to the identity; the `fb<B>`
candidate label remains the distinct nonidentity subgroup-usable count.

| Cell | Geometric points | Usable `B` | Folded columns |
| --- | ---: | ---: | ---: |
| n17a1 | 63 | 62 | 29 |
| n19a0 | 65 | 62 | 27 |
| n23a0 | 75 | 72 | 33 |
| n23a1 | 53 | 52 | 23 |
| n31a0 | 69 | 66 | 27 |

“Feasible” is exact group-level relation existence on the **reported**
ordinary queries. “Missed” means that a feasible query ended in the worker's
bounded-incomplete verdict; it is neither proved UNSAT nor a verified
relation from that worker. Timed-out jobs have no report and therefore no
auditable query count or yield.

| Cell | Solver | Reported queries | Worker witnesses | Exact feasible | Feasible but incomplete | Exact infeasible | Original rank |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| n17a1 | F4 | — | — | — | — | — | — |
| n17a1 | F5 | 104 | 29 | 33 | 4 | 71 | 29/29 |
| n17a1 | SAT native-XOR | — | — | — | — | — | — |
| n17a1 | SAT CNF | — | — | — | — | — | — |
| n19a0 | F4 | — | — | — | — | — | — |
| n19a0 | F5 | 24 | 2 | 2 | 0 | 22 | 2/27 |
| n19a0 | SAT native-XOR | 24 | 0 | 2 | 2 | 22 | 0/27 |
| n19a0 | SAT CNF | 24 | 0 | 2 | 2 | 22 | 0/27 |
| n23a0 | F4 | 4 | 0 | 0 | 0 | 4 | 0/33 |
| n23a0 | F5 | 4 | 0 | 0 | 0 | 4 | 0/33 |
| n23a0 | SAT native-XOR | 4 | 0 | 0 | 0 | 4 | 0/33 |
| n23a0 | SAT CNF | 4 | 0 | 0 | 0 | 4 | 0/33 |
| n23a1 | F4 | 4 | 0 | 0 | 0 | 4 | 0/23 |
| n23a1 | F5 | 4 | 0 | 0 | 0 | 4 | 0/23 |
| n23a1 | SAT native-XOR | 4 | 0 | 0 | 0 | 4 | 0/23 |
| n23a1 | SAT CNF | 4 | 0 | 0 | 0 | 4 | 0/23 |
| n31a0 | F4 | 1 | 0 | 0 | 0 | 1 | 0/27 |
| n31a0 | F5 | 1 | 0 | 0 | 0 | 1 | 0/27 |
| n31a0 | SAT native-XOR | 1 | 0 | 0 | 0 | 1 | 0/27 |
| n31a0 | SAT CNF | 1 | 0 | 0 | 0 | 1 | 0/27 |

On n19a0 the F5, native-XOR SAT and CNF SAT reports have the **same 24
ordered natural queries**. Trials 3 and 10 (probe scalars `81459` and
`104177`) are the only two that admit a three-point relation. F5 returned
verified witnesses at both. Each SAT arm used its full 100,000-conflict cap
on each query and reported `incomplete`, including those two. The exact
oracle independently found the same triples as F5, in index sets
`{4,22,45}` and `{14,46,63}`. Thus the observed zero SAT yield in this cell
is a **solver-budget miss of two existing relations**, not evidence that
these 24 queries had no relations. These two disclosed queries are suitable
future correctness controls, but choosing them after seeing the result means
they cannot estimate natural relation yield in a new experiment.

On n17a1, 33 of F5's 104 ordinary queries were exactly feasible. F5
witnessed 29 and reached full rank; four feasible queries ended
bounded-incomplete at trials 4, 67, 71 and 75. The verified one-target F5
recovery and its 15.236 s online interval remain as reported in the parent
PR. The exact oracle does not turn missed PDP attempts into worker relations,
does not revise the measured rank or timing, and does not qualify SAT.

As a descriptive uncertainty check, nominal 95% Wilson intervals under an
independent-query model are 33/104 = 31.7% (23.6–41.2%) for exact n17
feasibility, 29/104 = 27.9% (20.2–37.2%) for F5 witnessed yield, and
2/24 = 8.3% (2.3–25.8%) for exact n19 feasibility. For the n19 SAT arms,
0/24 witnessed yield has a nominal upper bound of 13.8%. The zero-exact-yield
n23 panels have only four queries per arm (upper 49.0%); n31 has one (upper
79.3%). These intervals **do not** account for one fixed seed, common queries,
the full-rank stopping rule, or variation between targets. The n19 SAT and
F5 arms are paired on the same query stream, so their observed two-versus-zero
witness difference is exact for this stream; there are no independent
replications from which to estimate a paired cost distribution. The group
oracle is a correctness diagnostic, not an executable IC candidate or a
performance baseline.

The next SAT experiment should freeze a changed encoding or backend, first
use the two n19 relation-positive queries as *correctness controls*, then
measure natural yield on newly fixed ordinary queries and require a complete,
source-bound one-target solve on a disclosed point. Only after that gate can
fresh paired IC/incumbent/rho targets address the requested speed comparison.
No claim is drawn from the four timeouts or from zero-yield cells with one or
four samples.
