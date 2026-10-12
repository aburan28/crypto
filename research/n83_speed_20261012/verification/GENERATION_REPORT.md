# Fresh degree-83 targets: generation and independent replay

Observed native checks, 2026-10-12 UTC. Eight fresh points passed generation
and a separately invoked fixture replay. No solver timing is reported here.
The final 139-process answer panel has a separate receipt when it completes.

| Requirement | Status and evidence |
| --- | --- |
| Exact legacy degree-83 field, curve, generator and subgroup | Verified by independent arithmetic in `native_helpers.hpp` |
| Three fresh training points and five fresh holdout points | Verified; eight frozen point-only files |
| Disjoint from original Q, all K600 points, all eight earlier optimization targets and each other | Verified; 609 initially excluded signed orbits |
| Point membership, subgroup membership and scalar replay | All eight pass `on(Q)`, `rQ=O`, `dG=Q` |
| Deterministic sampling and exclusion stream | Reproduced from SplitMix64 seed 202610126201; eight draws, zero rejections |
| Separation of solver points and validation scalars | Point-only inputs and a separate checker-only sidecar |
| All 139 solver answers, including losing configurations | Outside the generation receipt; full-panel checker is implemented |

## Inputs and observed checks

The field polynomial is `z⁸³+z⁷+z⁴+z²+1`, with low terms `[0,2,4,7]`.
The curve is `y²+xy=x³+x²+1`, `a=b=1`, and its subgroup has prime order
`r=8569786107849059`. The generator is
`G=(6614244651398184860292015,8868260410029249391539529)`.
The field degree 83 is different from the approximately 53-bit subgroup size.

The native checker verifies the Rabin irreducibility conditions for the prime
degree 83, deterministic 64-bit Miller–Rabin primality for r, the Frobenius trace
recurrence, curve order and cofactor, exact G and original Q membership, and
`[lambda]G=tau(G)` for `lambda=2998004823950457`. It checks `lambda^83=1`,
`lambda!=1` and its characteristic polynomial. It replays the 600 base
representatives and all 99,600 signed Frobenius coordinates and labels.
These checks use a standalone shift/XOR field multiplier and Euclidean inverse,
not the optimized solver's field code.

It also reproduces all eight previous targets from the earlier seed
`202610128301`, checks their scalar sidecars only in the verifier, and binds
their point files to the frozen earlier receipt. Thus previous exclusion inputs
are not trusted merely because their metadata says they passed.

## Signed-orbit exclusion invariant

For each nonidentity point, the implementation uses the key
`min{x(P),x(tau(P)),...,x(tau^82(P))}`. Negation on this binary curve is
`−(x,y)=(x,x+y)`, so the two signs have the same x coordinate. For a fixed
nonzero x, the curve equation has exactly the paired roots y and x+y whenever
one root exists. Thus matching a Frobenius x coordinate means matching a
point up to sign; equality of minimum keys detects equality of signed orbits.
An odd prime-order subgroup has no nonidentity point with x=0, since such a
point equals its negation and has order two.

The verified lambda has exact order 83, because 83 is prime and lambda!=1.
Its powers cannot equal −1: their group has odd order. Every nonidentity
subgroup point therefore has 83 distinct Frobenius points and 166 distinct
signed Frobenius points. This is the algebraic reason that checking the
minimum x key excludes the entire signed orbit, rather than only one point.
The native checks and source support this argument; a separate proof-assistant
formalization of the orbit argument is not asserted.

Initially excluded: 600 K600 representative orbits, the original public Q orbit,
and eight prior optimization target orbits, a total of 609 distinct keys.
Every accepted new point adds one key before the next proposal is considered.
All eight new keys are disjoint from that initial set and from one another.

## Sampling, binding and cost scope

SplitMix64 produces reproducible 64-bit pseudorandom proposals. For bound
`r−1=8569786107849058`, the threshold is `2^64 mod (r−1)=4564369618378800`.
A draw below that threshold is rejected; otherwise `d=draw mod (r−1)+1`.
This conversion avoids modulo bias for uniform 64-bit proposals. The fixed-seed
pseudorandom stream and signed-orbit exclusion condition are explicit parts
of the workload; no claim of cryptographic randomness is made.

Observed: eight draws, zero range rejections, zero excluded-orbit rejections.
The first three accepted points are training_1 through training_3; the next
five are holdout_1 through holdout_5. Each file contains one coordinate pair,
with no scalar. The sidecar is only for independent validation after an answer.

Generator receipt SHA-256:
`2b2a4987a5f93e751f5a28eb7856d31f8833a5b4d55f9baa92eb2c5ec427e0ae`.
Five-target combined point-file SHA-256:
`463d2cc915a37cf095ec81dfbdcded89ad2f32729068e26ff355b386f80e9a9c`.
Fixture replay receipt SHA-256:
`7e8919532fd8332d56f749835e13570502b17db791ff045eb3be12e9381f20be`.

Generation and replay are verification work, excluded from both timed solver
arms. They do not measure solver speed, independently reproduce timing, prove
every intermediate solver state or supply the missing admission/second-host
gates. Full commands, source and input hashes are in the adjacent README and
receipts. See [generation.svg](generation.svg) for the exclusion and validation
flow; it is a diagram because this step contains no comparable timing values.
