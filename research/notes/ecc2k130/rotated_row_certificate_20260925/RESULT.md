# Frozen Frobenius-compressed relation-row certificate

**Decision: PASS for row semantics on the frozen witnesses.** At the source
commit `ad2e447a15c4e57b2d8aabf7938646b83c1e9a2d`, a producer and a
separately implemented field/group-law verifier reconstructed every one of
the 17,917 positive full-point `Q+T` witnesses saved by [PR #762](https://github.com/aburan28/crypto/pull/762).
The identity

`Σ_i [λ^i][4]τ^(n−i)(P_i) = [4]Q`

held for each saved tuple in all four rotated n=13 arms. The verifier rebuilt
the factor bases and complete 2,003-target × four-torsion-coset roster per arm
from the #762 archive, checked every saved source witness and every miss,
recomputed signs and aggregated coefficients, and evaluated the final row by
its own group law. The full archived data passed an independent archive-only
replay after packaging. This certifies the row representation for these
finite inputs; it does not certify row independence, a linear-system solve,
new PDP yield, n=131 factor-base coverage, or a discrete logarithm.

| Frozen rotated arm | Positive `Q+T` rows | Distinct canonical columns appearing in terms / nonzero rows | Row weights 0 / 1 / 2 / 3 | Zero projected terms | Repeated-column terms | Producer arm wall |
|---|---:|---:|---:|---:|---:|---:|
| n13, beta 3, m5 | 1,799 | 2 / 2 | 1 / 226 / 1,572 / 0 | 1,673 | 3,952 | 1.759 s |
| n13, beta 3, m6 | 6,307 | 2 / 2 | 1 / 470 / 5,836 / 0 | 7,288 | 18,412 | 7.026 s |
| n13, beta 7, m5 | 1,799 | 2 / 2 | 1 / 226 / 1,572 / 0 | 1,679 | 3,946 | 1.812 s |
| n13, beta 7, m6 | 8,012 | 3 / 3 | 1 / 261 / 3,464 / 4,286 | 13,326 | 14,699 | 9.523 s |

The four zero-weight natural rows remain in the archive, and no canonical
column cancelled to zero after coefficient aggregation in these saved rows.
Small column counts reflect the tiny n=13 factors of #762; they are not a
matrix-rank measurement. For each beta/m, the separate all-`(0,1)` control
produced a zero row, while `P,τ(−P),(0,1),…` produced one repeated
canonical column and a nonzero row. Its coefficients are `(1,−89)` for beta 3
and `(−1,+89)` for beta 7, because the canonical representative flips sign.
Both source sums and cofactor-four `Q+T` splits passed in both implementations.

The optional n=131 check used **eight public synthetic planted tuples**, four
each for `(m,d)=(5,25)` and `(6,21)` with beta 3 and fixed SHA-256 masks. It
made no PDP search. The m5 four rows had weights `1,4,5,5` and 15 distinct
canonical columns in their terms; m6 had weights `1,5,6,6` and 18 columns.
Each tuple's raw source sum, rational torsion class, [4] projection, inverse
Frobenius transport, `λ^i` action and signed aggregation passed independent
full-point replay. The planted repeated-mask row compressed to one column for
each arity; the x=0 tuple had exactly one projected-zero term. These checks
show the row algebra works on the fixed n=131 polynomial model, with no claim
about the Certicom target or normal-base relation generation.

The frozen run occupied 12:17:01–12:18:16 UTC on the local macOS host. The
producer used 21.238 s wall and 20.199 s child user CPU; its four n=13 arms
were each below 120 s and 512 MiB and its two n=131 cells took 0.490 and
0.585 s, each below 30 s. Independent verification used 52.914 s wall and
51.042 s child user CPU. The launcher reported a 177,078,272-byte high-water
child RSS across the run. Full operation counts, stage costs, controls and
per-row witnesses are in the raw archive. Source/archive setup and
constructor costs are charged in the wrapper and frozen construction receipt.
Wall times are local-host diagnostics, not matched attack timings.

The next falsifiable gate is to compute exact modular-q rank on these frozen
rows, then test point-only target-log recovery on a separately frozen holdout
whose labels are hidden until validation. The 17,917 correct row identities
are **not** 17,917 independent relations; each n=13 arm has only two or three
canonical columns. Any n=131 admission then needs a compressed-base
cardinality/rank-cost model, a separately budgeted PDP solver, relation yield,
memory bound and same-workload rho comparison. The separate projected-PDP
corpus and solver prerequisites are tracked outside this certificate.
