# Result: native XOR cuts storage, but the public n83 PDP remains indeterminate

The frozen native-XOR representation passed all exporter, solver-convention,
planted-model, and numerical replay controls.  On the exact public n83 target,
CryptoMiniSat 5.16.0 returned `INDETERMINATE` at the preregistered 120-second
cap and produced no model.  The decision is
`BOUNDED_NATURAL_TARGET_TIMEOUT`: this is a verified representation reduction,
not a relation, discrete logarithm, or faster-than-rho algorithm.

## Correctness controls

All seven native module tests passed, including exhaustive type-II-ONB n=3
field laws, exhaustive complete addition over every valid n=3 point pair,
exact n83 rho-scalar replay, and both complete and inductive planted circuits.
The toy extended-DIMACS test accepted the true XOR/AND assignment and rejected
independently corrupted XOR and AND outputs.

The pinned CryptoMiniSat binary also interpreted the committed parity syntax
as intended: the exact `z=a xor b` control exited 10/SAT with model
`v 1 2 -3 0`, while its one-bit-corrupted twin exited 20/UNSAT.

## Planted native-XOR certificate

The deterministic planted chain retained the frozen inductive DAG and emitted:

| Quantity | Planted native-XOR chain |
| --- | ---: |
| Primary inputs | 2,996 |
| XOR gates / parity rows | 495,192 |
| AND gates | 195,682 |
| DAG nodes / variables | 693,872 |
| Remaining CNF clauses | 587,049 |
| Total constraint rows | 1,082,241 |
| Extended-DIMACS bytes | 22,859,716 |

The exporter streamed every parity row and CNF clause against the complete
model.  A fresh process then replayed every DAG gate, verified all ten curve
points, and recomputed their group sum.  Both checks passed.  The planted
extended-DIMACS BLAKE3 is
`91b95c8c73ff777cb1466624b41ff8443a1c3ddab71034793e233ac331d1865a`
and the complete-model SHA-256 is
`d86e724e638bf4bf794a77ea6330541a3c2769af6a6fb60bb93b2729f3842053`.

## Frozen public-target attempt

The exact public Q reproduced the parent inductive DAG prefix
`e8715d805517fc35cb70f6b96c5d387a28579c0800e259f69e642c91e5636f88`
and produced:

| Quantity | Ordinary CNF in PR #1346 | Native XOR here | Change |
| --- | ---: | ---: | ---: |
| DAG nodes / variables | 694,818 | 694,818 | 0.00% |
| XOR gates / parity rows | 495,808 gates | 495,808 rows | representation change |
| AND gates | 196,012 | 196,012 | 0.00% |
| CNF/total constraint rows | 2,571,271 | 1,083,847 | -57.85% |
| Serialized bytes | 57,642,729 | 22,893,789 | -60.28% |
| Peak solver RSS | 778.93 MiB | 412.18 MiB | -47.08% |
| Terminal status | `UNKNOWN` | `INDETERMINATE` | no relation |

The public extended-DIMACS BLAKE3 is
`6cc3d18c15b6cbd22dea85b03f03f42380bf20acdbdb3ff6125185db397a21f6`.
CryptoMiniSat parsed 588,039 CNF clauses and 495,808 XOR rows in 0.27 seconds.
At the stop it reported 22,075 conflicts, 48,666 decisions, approximately 562
million propagations, 422,068 kB peak RSS, and 120.30 seconds in its worker
thread.  The harness observed 120.45 seconds and exit 15.

The anticipated Gaussian limitation occurred exactly as disclosed.  After
preprocessing, the solver found one connected residual parity component with
481,127 rows and 666,687 columns.  This exceeded its default 100,000-row and
100,000-column limits, so it enabled zero Gauss-Jordan matrices.  Raising those
limits after seeing the result would attempt a prohibitively large dense
matrix and would violate the frozen protocol; it was not done.

## Decision and next admissible step

Native parity rows are a substantial storage and memory improvement, but
gate-level XOR preservation alone does not make this public PDP practical.
The solver still sees one enormous connected parity component and cannot apply
its Gaussian engine to it.

A technically distinct follow-up would quotient out definitional XOR cones:
represent each XOR node as an affine form over primary inputs and AND outputs,
materialize only parity forms that cross nonlinear AND boundaries, and choose
explicit cut points that keep every parity block below a preregistered matrix
limit.  That approach must first prove equivalence on exhaustive toy DAGs and
replay the planted n83 model; it must not merely increase this run's timeout or
matrix limits.

The only complete solve remains the strong signed-Frobenius rho receipt with
201,733,439,488 charged walk iterations (`2^37.553659...`) and recovered scalar
`467066815623456506232910`.  This experiment has no public relation, relation
rank, factor logarithms, target recovery, or matched cold-cost crossover, so
it makes no speed claim.  Degree 51 remains a separate composite-extension
case (`51=3*17`) and cannot establish the prime-degree n83 claim.

## Reproduction and evidence storage

The frozen source is commit `9deefe38109b5745e974fbbb6808f332cd4c65eb`
and the preregistration head is
`2ee9c33b10cf941aaff0844415e4b8b382237c01` in PR #1352.  `FROZEN.json`
records commands, caps, source digests, official solver-asset identity, and the
pre-outcome state.  `MANIFEST.json` records generated-file sizes and hashes.
The deterministic 22.9-MB instances and 5.3-MB planted model are hash-retained
but not committed; the small receipts, full solver log, control receipt,
UNKNOWN model, and replay summary are committed.
