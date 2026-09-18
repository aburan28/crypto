# Point-halving mixed rho

Question: can replacing some affine-addition iterations by point halvings beat
the measured table walk on one RTX PRO 6000?

Answer: **not with either representation measured.** Polynomial-state halving
improves the exact two-product primitive from 16.053 to **17.921 B/s**, but is
still only 0.890 of the current 20.134 B/s table walk before a rho-producing
addition branch or mixed-walk control is charged.

## Boundary and target

Unit: billions of complete point operations per second on the same GPU. The
reference is the current verified B17 table-walk preset at **20.134 B/s**.

Point halving alone is a permutation and is not a rho iteration. A valid mixed
walk must also pay for a non-permutation branch. A point-dependent branch
diverges across CUDA warps and executes both paths; a warp-uniform alternating
phase avoids divergence only by making the phase part of the state.

## Single table

| variant | products / point op | measured B/s | / 20.134 reference | correct | class |
|---|---:|---:|---:|---|---|
| selected B17 table-add reference | 5.2941 | **20.134** | 1.000 | 300/300 reports, 0 dropped | reference |
| halving, first phase-map build (superseded) | 2 | 12.125 | 0.602 | 512/512 | engineering |
| halving, normal-basis state | 2 | 16.053 | 0.797 | 512/512 | engineering, superseded |
| **halving, polynomial state, 512 threads** | **2** | **17.921** | **0.890** | **512/512** | **engineering, rejected for mixed rho** |
| 50/50 uniform alternation, raw | 3.6471 average | 18.963 | 0.942 | extrapolation from measured rows | model |
| 50/50 uniform alternation, phase-priced effective | 3.6471 average | 13.409 | 0.666 | phase doubles the rho state | model |
| 25 B/s target | — | **25.000** | **1.242** | required | above both primitives |

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

The remaining cost is linear algebra, not multiplication count. Directly
generating the inverse-Frobenius permutation raised the first 12.125 B/s build
to 16.053 B/s. Keeping state and both products polynomial raises it again to
17.921 B/s, but the quadratic solves and subgroup selection still make halving
slower than the highly amortised affine-add path.

Frozen measurements and commands are in
[`benchmarks/point-halving/result.json`](benchmarks/point-halving/result.json).
The polynomial-state follow-up is in
[`benchmarks/point-halving-poly`](benchmarks/point-halving-poly/README.md).
