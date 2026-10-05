# Result: the native n83 relation is exact, but the public PDP hit the cap

The frozen native implementation passed its representation controls and did
not solve the exact public n83 target within the preregistered 120-second
CaDiCaL window.  The decision is therefore
`BOUNDED_NATURAL_TARGET_TIMEOUT`: this complete-point encoding is a working,
replayable relation gate, but it supplies no natural relation, DLP, or
faster-than-rho result.

## Exact cross-repository identity

The identity stage reproduced the cryptanalysis repository's
`EC1N83Ckb1h876c2921cb64` instance in the crypto repository's native code:

- G and Q are on `y^2+xy=x^3+1` in the frozen type-II optimal normal basis;
- `[2417851639230796216685689]G=O`; and
- `[467066815623456506232910]G=Q`.

The last equality independently replays the scalar recovered by the strong
signed-Frobenius rho run.  That run charged 201,733,439,488 walk iterations
(`2^37.553659...`).  Its receipt does not provide a reconstructible matched
single-worker wall time, so this result does not invent one.

## Complete planted certificate

The ten frozen slot dimensions are `[9,9,9,8,8,8,8,8,8,8]`; their coordinate
positions partition all 83 normal-basis coordinates exactly once.  The
complete planted chain built in 0.209 s and produced:

| Quantity | Planted chain |
| --- | ---: |
| Primary inputs | 2,996 |
| XOR gates | 715,672 |
| AND gates | 278,722 |
| Total DAG nodes / CNF variables | 997,392 |
| CNF clauses | 3,698,857 |
| CNF bytes | 83,368,485 |

The generated full Boolean model was replayed in a fresh process.  Every DAG
gate matched, all ten decoded factors satisfied the curve equation, and their
independently recomputed group sum equalled the planted target.  The model
SHA-256 is
`f6956638c2bb75d3befa58ba29a1d356dca2ebc2e5c76b29092861ef1ec05a1a`.

## Frozen public-target attempt

The exact rho target produced a closely matched but independently fixed
circuit:

| Quantity | Public target |
| --- | ---: |
| Build wall time | 0.218 s |
| Primary inputs | 2,996 |
| XOR gates | 716,288 |
| AND gates | 279,052 |
| Total DAG nodes / CNF variables | 998,338 |
| CNF clauses | 3,702,311 |
| CNF bytes | 83,446,635 |
| Solver real time | 120.03 s |
| Peak solver RSS | 1,340.11 MB |
| Conflicts | 255,395 |
| Decisions | 270,483 |
| Propagations | 388,730,839 |
| Terminal model | `c UNKNOWN` |

CaDiCaL exited 0 at its real-time limit, wrote no assignment, and neither
proved UNSAT nor found SAT.  The solver stdout SHA-256 is
`7bf253783de5a01b1ca363dfd9a24219e00e7513c5d46faed4b0323fac4c2a49`.
The ordinary CNF's BLAKE3 digest recorded by the streaming exporter is
`bd4c5148fa97d7af749c5736bd389b77fe938a2b978aed901d893b8bd463b892`;
its independent SHA-256 is
`22db3390b0f398f3f68e80ea3186be36b8cb17ee7722c967dafa541daf935176`.

## Decision and next admissible step

This run rejects any claim that the present direct Tseitin/CaDiCaL route is
already a practical natural-target PDP solver.  It does not prove that Q is
outside the tuple domain, and it is not a no-go theorem for high-arity
descent.  It does show that merely reaching adequate counting capacity is not
enough: one bounded PDP query consumed the whole solver budget and produced
zero usable relations.

A follow-up may preregister a structurally smaller native circuit—most
plausibly a Karatsuba-style normal-basis multiplication schedule plus removal
of redundant intermediate curve-validity equations, backed by exhaustive toy
proofs.  It must use a fresh target/cap gate and may not relabel this timeout.
Even a verified SAT relation would only open relation-yield testing.  A speed
claim still requires full relation collection, rank, factor logarithms,
held-out target recovery, and complete cold cost below matched strong rho.

The 51-degree case remains separate: 51 is composite (`3*17`) and admits
proper subfields, so an n51 subfield descent cannot establish an n83
prime-extension result.

## Reproduction and evidence storage

The frozen source is commit `422753df47b3420031dbe78b3e84fdf964745c21`
and the preregistration head is `2428ae0ceac3fc1bdb6920e4db295892de522ce5`
in PR #1343.  Commands and caps are in `FROZEN.json`.  `MANIFEST.json` records
the byte count and SHA-256 of every generated file.  The two 83-MB CNFs and
7.7-MB planted model are deterministic and hash-retained but intentionally
not committed to avoid permanent repository bloat; the small identity,
decoded planted receipt/replay, public result, UNKNOWN model, and complete
solver logs are committed.

