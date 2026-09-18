# Point-halving mixed rho

Question: can replacing some affine-addition iterations by point halvings beat
the measured **17.952 B updates/s** table walk on one RTX PRO 6000?

Answer: **not with this representation.** The exact two-product halving
primitive reaches **16.053 B/s**, 0.894 of the matched reference, before a
rho-producing addition branch or mixed-walk control is charged. A mixture of
an operation slower than the reference and the reference cannot exceed the
reference without unpriced parallelism. The mixed walk was therefore rejected
before changing the campaign protocol.

## Boundary and target

Unit: billions of complete point operations per second on the same GPU. The
reference is the current table-walk preset (`PACKED_INLINE_POLY=3`, pair ILP,
L2 persistence, slot unroll 2), freshly rerun three times. The success target,
set by the request, is strictly above its **17.952 B/s** median with correct
subgroup points.

Point halving alone is a permutation and is not a rho iteration. A valid mixed
walk must also pay for a non-permutation branch. A point-dependent branch
diverges across CUDA warps and executes both paths; a warp-uniform alternating
phase avoids divergence only by making the phase part of the state.

## Single table

| variant | products / point op | measured B/s | / 17.952 reference | correct | class |
|---|---:|---:|---:|---|---|
| table-add reference | 5.3125 | **17.952** | 1.000 | 300/300 reports, 0 dropped | reference |
| halving, first phase-map build (superseded) | 2 | 12.125 | 0.675 | 512/512 | engineering |
| **halving, direct square-root map + trace-only subgroup selection** | **2** | **16.053** | **0.894** | **512/512** | **engineering, rejected** |
| 50/50 uniform alternation, raw | 3.65625 average | 16.950 | 0.944 | extrapolation from measured rows | model |
| 50/50 uniform alternation, phase-priced effective | 3.65625 average | 11.985 | 0.668 | phase doubles the rho state | model |
| success boundary | — | **>17.952** | **>1.000** | required | boundary |

The alternation rows are extrapolations, not kernel measurements. Their raw
rate is the harmonic mean of the two measured operation rates. The effective
row divides by `sqrt(2)` because an autonomous alternating map walks on
`(point, phase)`. They are included to prevent the raw-operation count from
being reported as rho progress.

## Arithmetic and correctness

For `P=(X,Y)` on `y²+xy=x³+1`, a half `Q=(x,y)` is recovered from

```
lambda² + lambda = X
x = sqrt(Y + X + lambda X)
y = x(lambda + x)
```

The two rational halves differ by 2-torsion. Choosing either root is not enough:
252 of the first 512 test points selected the wrong coset. The implementation
selects the unique odd-subgroup half with a second 2-descent bit. If
`lambda2²+lambda2=x`, then four-divisibility is tested without another field
product:

```
Tr(x_quarter) = Tr(x(lambda+x)) + Tr(x) + Tr(lambda2 x).
```

In the self-dual ONB each product trace is just parity of coordinate-wise AND.
The final primitive therefore remains two field products. The CUDA result was
checked on `[i]P`, `i=1..512`, against scalar `[2^-1 mod ell]`; all coordinates
matched, every result doubled back to its input, and all generated quadratic
roots and square roots passed.

The remaining cost is linear algebra, not multiplication count. The optimized
kernel contains 1,744 static SASS instructions and six `CLMAD` instructions
per halving. Directly generating the inverse-Frobenius permutation raised the
first 12.125 B/s build to 16.053 B/s, but the two quadratic solves and subgroup
selection still make halving slower than the highly amortised affine-add path.

Frozen measurements and commands are in
[`benchmarks/point-halving/result.json`](benchmarks/point-halving/result.json).
