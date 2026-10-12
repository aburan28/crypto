# Fixed-addition rho screen toward 12 B/s on G7

Status: complete negative result. The fixed-add kernel improved raw throughput
to 7.601382 B/s median, but its 1.399081x measured collision-work penalty reduced
effective throughput to 5.433126 B/s. Retain the selected walk. See
[RESULTS.md](RESULTS.md).

## Hypothesis

The selected walk applies `R <- R + sigma^j(R)`. Its forward path converts X
and Y from polynomial to normal basis, routes two Frobenius permutations, XORs
the coordinate deltas, and converts both deltas back to polynomial basis.
Instead use the same eight-way, normal-Hamming-weight selector and add one of
eight fixed subgroup points `sigma^(3..10)(P)`. The point-addition denominator
and numerator are then polynomial XORs with constant coordinates. This retains
the selected polynomial product, weighted prefix, block inversion, backward
point formula, compact state, B16/split2 geometry, DP predicate, canonical
Frobenius/negation key, restart accounting, and full seed replay.

For branch j the replay coefficient changes additively by `s^j`, rather than
multiplicatively by `1+s^j`. A collision between endpoints
`A = epsilon sigma^c(B)` therefore solves
`k=(a_A-epsilon*s^c*a_B)/(epsilon*s^c-1) mod ell`. Constants are generated from
curve parameters and copied in polynomial basis; no runtime conversion of the
addends is allowed. The candidate must reject incompatible configurations.

This is an algorithm change, so raw updates/s alone cannot qualify it. First run
2,000 planted GF(2^23) solves for selected and fixed-add walks, using identical
start seeds, DP threshold, eight-way partition, orbit key, and verified recovered
scalars. Report mean complete iterations and the fixed/selected collision-work
ratio. A material collision regression closes the candidate even if raw speed is
higher.

## GPU gates

Freeze an isolated copy of the selected CUDA 13.3-front-end / 13.4-assembler
source and selected executable. Add a default-off `ECC_PACKED_FIXED_ADD` mode.
Before timing require an independent fixed-point coordinate/scalar table proof,
all built-ins, exact candidate GPU state and DP output against an independent
CPU fixed-add reference, seed replay and verified collision recovery, partial
populations, checkpoint/resume, restart boundaries, zero-denominator fixtures,
actual-walk storage checks, and memcheck/initcheck/synccheck. The disabled mode
must match all selected native entries.

The initial speed screen uses three interleaved selected/candidate repetitions on
the same g7.2xlarge / RTX PRO 4500 at 165 W, each with 34,359,738,368 complete
scalar updates. Report raw rate, paired 95% log-speedup interval, measured
collision-work ratio, and collision-adjusted rate. Promotion requires both raw
and adjusted paired intervals wholly above 1, then fresh benchmark and DP34
confirmation with exact outputs. Generic work remains `sqrt(n/262)` times the
measured collision-work ratio. Full-DLP S remains null until an end-to-end
recovery measurement exists.
