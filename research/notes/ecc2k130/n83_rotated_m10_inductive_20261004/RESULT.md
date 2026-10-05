# Result: inductive validity shrinks the n83 gate, but the public PDP still times out

The separately frozen inductive-validity circuit passed every proof and replay
control, and it is materially smaller than the complete-per-edge circuit.  It
did not find a relation for the exact public n83 target within the frozen
120-second CaDiCaL window.  The decision is therefore
`BOUNDED_NATURAL_TARGET_TIMEOUT`: this is a verified representation
improvement, not a natural relation, discrete logarithm, or faster-than-rho
algorithm.

## Proof and drift controls

Before the public run, the validity-free addition relation was exhaustively
checked on the type-II-ONB n=3 curve for every valid input pair, candidate
result, and slope.  It accepted exactly the numerical group sum and admitted a
slope for every valid pair.  The n83 deterministic planted chain then replayed
all Boolean gates, all ten factor curve equations, and the independently
recomputed group sum in a fresh process.

The parent complete public circuit was also rebuilt before the inductive
attempt.  Its DAG prefix BLAKE3 remained
`689697d223df16ad1cb65edfa6f44f4ecbb3f0032ca55f8b31587da8696525cf`
and its ordinary-CNF BLAKE3 remained
`bd4c5148fa97d7af749c5736bd389b77fe938a2b978aed901d893b8bd463b892`.
Thus the comparison is against the frozen parent representation rather than a
drifting reconstruction.

## Verified planted certificate

The ten slot dimensions remain `[9,9,9,8,8,8,8,8,8,8]` and partition all 83
normal-basis coordinates.  The inductive planted chain produced:

| Quantity | Planted chain |
| --- | ---: |
| Primary inputs | 2,996 |
| XOR gates | 495,192 |
| AND gates | 195,682 |
| Total DAG nodes / CNF variables | 693,872 |
| CNF clauses | 2,567,817 |
| CNF bytes | 57,564,571 |

The full-model SHA-256 is
`d86e724e638bf4bf794a77ea6330541a3c2769af6a6fb60bb93b2729f3842053`.
Fresh-process replay returned `PASS_MODEL_REPLAY`: every gate matched, each
decoded point was on `y^2+xy=x^3+1`, and the ten-point sum equalled the fixed
planted target.

## Paired public-target comparison

Both circuits use the same curve, public Q, factor slots, ordinary Tseitin
encoding, CaDiCaL 3.0.1 binary, one thread, and 120-second real-time cap.

| Quantity | Complete per edge | Inductive validity | Change |
| --- | ---: | ---: | ---: |
| Primary inputs | 2,996 | 2,996 | 0.00% |
| XOR gates | 716,288 | 495,808 | -30.78% |
| AND gates | 279,052 | 196,012 | -29.76% |
| DAG nodes / CNF variables | 998,338 | 694,818 | -30.40% |
| CNF clauses | 3,702,311 | 2,571,271 | -30.55% |
| CNF bytes | 83,446,635 | 57,642,729 | -30.92% |
| Solver peak RSS | 1,340.11 MB | 778.93 MB | -41.88% |
| Conflicts at stop | 255,395 | 311,472 | +21.96% |
| Decisions at stop | 270,483 | 367,314 | +35.80% |
| Propagations at stop | 388,730,839 | 560,249,682 | +44.12% |
| Terminal model | `c UNKNOWN` | `c UNKNOWN` | no relation |

The public inductive DAG prefix is
`e8715d805517fc35cb70f6b96c5d387a28579c0800e259f69e642c91e5636f88`
and its CNF BLAKE3 is
`052b4e522a4b8d8f188e4983bb155fc491af80a01dcc546a86a94a497c6ded30`.
CaDiCaL used 120.00 seconds of real time and 778.93 MB peak RSS, exited zero at
the cap, wrote `c UNKNOWN`, and supplied no assignment to replay.

## Decision and implication

The representation-win criterion passes, but the solver criterion does not.
Removing redundant curve equations saves roughly 30% of the Boolean system
and 42% of peak memory; it does not turn the direct ordinary-CNF search into a
working public-target PDP.  The higher conflict and propagation counts in the
smaller circuit also show that raw clause volume is not the only obstruction.
Extending the same timeout would be post-outcome retuning and is not done.

A defensible next experiment must change the solving representation, not only
trim more copies of the same equations.  The circuit contains 495,808 explicit
XOR gates whose parity structure is flattened into ordinary Tseitin clauses.
A separately frozen native-XOR export with an XOR-aware solver is therefore a
more informative next gate.  It still must pass planted-model validation and a
bounded public target before any relation-yield panel.

The strong signed-Frobenius rho receipt remains the only complete n83 solve:
201,733,439,488 charged walk iterations (`2^37.553659...`) recover
`467066815623456506232910`.  This attempt has no relation, rank, factor
logarithms, target recovery, or matched cold-cost crossover, so it makes no
speed claim.  Degree 51 remains a separate composite-extension case
(`51=3*17`); a subfield-assisted n51 result would not establish a result for
prime extension degree 83.

## Reproduction and retained evidence

The frozen source is commit `203eb86df1d0cc5ee28d99e0de39f3bbe5f7d511`
and the preregistration head is
`9cc5ec0c62fa662c6a91b6e7da9cd5d1afc993df` in PR #1346.  Commands and caps
are fixed in `FROZEN.json`.  `MANIFEST.json` records every generated file's
size and SHA-256.  Large deterministic CNFs and the planted model are retained
locally and hash-pinned but not committed; the small receipts, replay summary,
UNKNOWN model, and complete solver logs are committed.
