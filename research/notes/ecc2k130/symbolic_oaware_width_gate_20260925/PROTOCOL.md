# Symbolic O-aware rotated-PDP width and normal-form gate (pre-outcome)

## Question and fixed scope

Can the current Rust SAT/Semaev stack express a **symbolic, multiword,
full-point** rotated PDP for the n131 unequal m9/m10 arms without enumerating
prefix states? The bounded answer here is an interface/width admission audit,
not a solver benchmark. The candidate relation is the complete affine
full-point relation in [DESIGN.md](DESIGN.md); it is a specification and is
**not yet a Boolean circuit**. The separate Kohel μ4 point-map check asks
whether that normal form has an exact full-point conversion on the n13
`y²+xy=x³+1` toy curve. It does not test any addition circuit.

The frozen parent is the repository commit in `FROZEN.json`. Inputs are
`INPUT.json`: n13 modulus `0x201b`, the K0 source equation, m9 slots
`[15,15,15,15,15,14,14,14,14]`, and m10 slots
`[14,13,13,13,13,13,13,13,13,13]`. The n131 source field/curve is used
only for width arithmetic. #781's complete n13-m5 O-aware one-hot exporter
and #785's 32-target SAT panel are the semantic comparators, not input to a
new solver run. Existing `coordinate_search.rs` and its exotic-coordinate
note already discuss Kohel μ4; this gate claims no discovery of that form.

## Frozen toy map and controls

For K0 and characteristic 2, test the explicit map to the split μ4 model

`(X0+X2)^2=X1*X3`, `(X1+X3)^2=X0*X2`:

`(x,y) -> (x², x²+y, 1, x²+x+y)`, and `O -> (1,1,0,1)`.

For every source point, including O, check both quadrics, nonzero projective
coordinates, scaled-coordinate invariance, inverse roundtrip, and uniqueness.
Independently enumerate all normalized `X2=1` μ4 solutions by iterating
`X0`, taking its unique square root `x`, and solving the first quadric for
`X1`; compare complete point sets. For `X2=0`, prove and test the unique
projective point `(1:1:0:1)`. The inverse is
`x=(X1+X3)/X2`, `y=(X0+X1)/X2` when `X2!=0`.

Fixed negative controls: corrupt the image of source `(0,1)` by flipping
`X1` and require quadric rejection; corrupt the `X2=0` chart with
`(1,0,0,1)` and require rejection; try affine inversion of genuine O and
require a guarded chart exception. Check source 4-torsion chain
`(1,0) -> (0,1) -> (1,1) -> O` under a separate bit-serial group law.
No SAT, relation, or target enumeration is run here.

## Admission decision and caps

`audit.py --preflight` only validates frozen file hashes and input schema.
The draft PR and hash-only CI must exist and pass before `--run` reads any
toy point-set or admission outcome. Then run exactly one cold
`python3 run.py --output-dir evidence` from this directory. The wrapper
runs `audit.py` and a separate `verify.py` sequentially under a total
60-second child wall cap and 256 MiB child-RSS acceptance cap. It archives
both exact child commands, UTC start/end, exit, stdout/stderr, wall, CPU,
RSS upper samples, all file hashes/bytes, and the run receipt SHA/bytes.
The independent verifier uses polynomial-product/reduction arithmetic,
Euclid inverses, and a complete Artin-Schreier-root table; it does not call
the producer's field arithmetic, trace, half-trace, or point enumeration.
If the first attempt fails, retain all files and its failed receipt under
`evidence/failure_0/`; never silently overwrite it. `ci_replay.py` checks
the committed transcripts/hashes and recomputes the independent toy
verifier from the archived producer receipt. Any hash drift, point mismatch,
failed control, timeout, or missing receipt fails the gate.

The n131 admission checker is deliberately fail-closed. It computes the
raw affine-chain floor and inspects the pinned `u64` source contracts.
It must return `NOT_ADMITTED` until a multiword symbolic field, unbounded
Boolean monomials/model lifting, a complete O-aware circuit, all 32 toy
Q+T outcomes with independent model replay and negative certificates,
and bounded n131 export/cost evidence are actually provided and checked.
The required artifacts are named in `INPUT.json` but absent here. A
successful width audit means the current path was correctly refused;
it is not an n131 solver-feasibility or attack-speed result.

Kohel's μ4 paper gives several bidegree-(2,2) addition-law charts, each
with an exceptional divisor; selecting among them is a separate circuit
obligation. Binary Edwards complete-formula literature is a different
potential route, but would still need a verified full-point conversion,
including O and torsion, for this exact curve. Neither route is admitted
by the toy point-map test alone.
