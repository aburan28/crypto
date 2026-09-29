# Decision: native three-summand semantic and representation gate

**PASS_SEMANTIC_AND_CAPACITY for the preregistered expanded controls.** This
is a construction and solver-interface result, not a native-leaf PDP yield or
ECDLP result. The expanded [protocol](PROTOCOL.md) was committed before its
candidate revision `d068fcffb2725a06b2f148c046983a84511859d2`, then the
20-file [source/input lock](FROZEN.json) was committed and pushed before the
[final held run](evidence/final/receipt.json). The first run and its narrower
x=0 SAT witness remain intact in [`evidence/initial`](evidence/initial/README.md).

On GF(2^5), the three implicit factor sets each had three physical points.
The independent point oracle enumerated all 27 ordered triples and found 27
distinct supported points among 44 rational curve points. For every triple,
the full native Boolean chain accepted the exact intermediate and sum and
rejected a different rational target with the other inputs held fixed. The raw
27-row stream has SHA-256
`744ac9a63965ae2e7046c40dbc65658264982b9f7961deab79c28f0d8815443a`.
It includes 24 generic/generic paths and three inverse/copy paths. A separate
shift/reduce Koblitz implementation recomputed and replayed every row.

| Fixed solver control | Result | Verified interpretation | SAT child wall |
| --- | --- | --- | ---: |
| First supported non-O target | SAT | Complete CNF model, T+T+T=T, inverse/copy | 0.0075 s |
| First target whose every witness is generic/generic | SAT | Complete CNF model, exact generic/generic point sum | 0.0066 s |
| First unsupported rational target | UNSAT | Independently exhaustive oracle confirms absence; no external UNSAT proof | 0.0063 s |

All three toy DIMACS files had 1,197 variables and 4,145 clauses; the
independent verifier checked every raw clause and both positive models, then
lifted every factor point, intermediate point and target against the separate
group law. The toy solver is a disclosed correctness control on 27 possible
factor triples, **not** a natural-target runtime or checked general UNSAT
engine.

| Actual degree-263 leaf | Full chain DAG nodes | Primary Boolean inputs | Exact target-O DIMACS bytes | Producer wall | Peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: |
| `[1,0]` | 650,376 | 1,312 | 48,285,979 | 1.879 s | 202,506,240 B |
| `[1,4]` | 650,391 | 1,312 | 48,287,148 | 1.943 s | 204,800,000 B |

Both native chains used their literal archived `b` coefficient and the same
131-dimensional normal-basis split `[44,44,43]`. The independent verifier
regenerated the basis by bit-polynomial squaring, checked its rank, rebuilt
both circuits and recomputed the DIMACS byte arithmetic. All producer and
verifier children exited 0 inside their frozen wall/RSS caps. A fresh
`ci_replay.py --evidence .../evidence/final` returned
`ARCHIVE_REPLAY_PASS` for all three arms.

The n131 figures establish only that a complete native three-summand circuit
and its deterministic CNF representation fit the selected caps. No n131 CNF
was solved, no useful physical factor count or relation yield was measured,
and no matched native/pullback target stream, rank or full logarithm was
computed. The next gate must use a separately frozen natural-target stream,
check every positive model, distinguish externally proved UNSAT from solver
UNKNOWN/resource exits, and establish physical-support/rank admission before
charging all stages against original, transported and exact pullback arms.
Keep full-DLP cost, n131 transfer assumptions and rho crossover **unset**.
