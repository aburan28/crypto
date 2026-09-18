# Point-halving mixed rho

Question: can replacing some affine-addition iterations by point halvings beat
the measured table walk on one RTX PRO 6000?

Answer: lambda-affine polynomial state produces a genuinely cheaper primitive:
**28.953 B/s with one field product**, 1.438 times the current 20.134 B/s table
walk. It still does not produce a 25 B/s rho walk: halving alone is a
permutation, and even an ideal free-dispatch 50/50 mixture has raw harmonic
mean 23.752 B/s before representation changes or collision constants.

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
| **halving, polynomial λ-state, selected** | **1** | **28.953** | **1.438** | **512/512** | **cheaper primitive; not a rho map alone** |
| ideal 50/50 λ-halving/table-add mixture | 3.1471 average | 23.752 | 1.180 | free-dispatch extrapolation | upper-bound model |
| 30 B/s primitive target | — | **30.000** | **1.490** | required | not met |

The mixture row is an extrapolation, not a kernel measurement. It is the
harmonic mean with free branch dispatch and omits the λ-to-affine conversion
needed by an addition branch, so it is an upper bound rather than a claim.

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
The one-product λ-state progression is in
[`THROUGHPUT-30B-HALVING.md`](THROUGHPUT-30B-HALVING.md).
