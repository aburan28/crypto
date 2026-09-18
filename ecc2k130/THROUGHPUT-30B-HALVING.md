# 30 B/s cheaper-iteration attempt

Target: **30 billion complete point halvings/s** on one RTX PRO 6000.
Result: **28.953 B/s**, target not met.

## Boundary

Unit: complete subgroup point halvings per second. The reference rho iteration
is the verified 20.134 B/s table-add walk. A halving primitive is useful only
inside a mixed walk because `[2^-1]` is a permutation of the prime subgroup.
Consequently its raw rate is not an ECDLP speedup by itself.

Success required median above 30 and exact agreement with scalar
`[2^-1 mod ell]` on every test point. A mixed-rho claim additionally requires
the addition branch, representation changes, dispatch and collision constant.

## Single table

| variant | products | measured B/s | / 30 | / table-add | correct | class |
|---|---:|---:|---:|---:|---|---|
| selected B17 table-add rho step | 5.2941 | 20.134 | 0.671 | 1.000 | 300/300 | reference |
| λ-state halving, first one-product build | 1 | 19.695 | 0.657 | 0.978 | 512/512 | engineering |
| + closed second-descent trace | 1 | 23.344 | 0.778 | 1.159 | 512/512 | engineering |
| + select before square root | 1 | 27.725 | 0.924 | 1.377 | 512/512 | engineering |
| + product/reducer CLMAD balance | 1 | 27.888 | 0.930 | 1.385 | 512/512 | engineering |
| **+ CLMAD quadratic-prefix scan** | **1** | **28.953** | **0.965** | **1.438** | **512/512** | **engineering; target not met** |
| dual independent chains/thread | 1 | 28.917 | 0.964 | 1.436 | 512/512 | engineering, neutral |
| ideal free-dispatch 50/50 mixed upper model | 3.1471 avg | 23.752 | 0.792 | 1.180 | extrapolation | model, not result |
| 30 B/s target | — | 30.000 | 1.000 | 1.490 | required | boundary |

## Algebraic reductions

In lambda-affine state `(u,lambda_Q)`, a half has

```text
lambda_P^2 + lambda_P = u
x_P^2 = u (lambda_Q + u + lambda_P + 1).
```

This uses one product. For the Koblitz curve's trace-zero `a`, selecting the
odd-subgroup root initially required evaluating a second half-trace on
`x_P`. Two identities move the decision before square root:

```text
Tr(x_P lambda_P) = Tr(t (lambda_P + u))
Tr(x_P H(x_P))   = Tr(t H(t)),       t = x_P^2.
```

If the candidate is wrong, toggle `t` by `u` and `lambda_P` by one, then take
one square root. In phase order, the remaining quadratic form satisfies

```text
Tr(t H(t)) = w + binom(w,2) mod 2
```

for the Hamming weight excluding phase bit zero. Finally, solving the
quadratic uses a prefix XOR; two `CLMAD.lo` instructions prefix the low and
high 64-bit chunks, with parity carried between them.

These are exact accounting changes, not dropped work. The generated code
checks the closed form on all basis vectors and all basis-vector pairs; the
device differential test checks full coordinates on 512 subgroup points.

## Rho verdict

The primitive crosses 25 B/s but does not establish a 25 or 30 B/s rho walk.
Halving alone has no birthday collisions. The optimistic harmonic mean of
equal measured halving and addition steps is only 23.752 B/s even before GPU
branch divergence, lambda-to-affine conversion and the mixed-walk collision
constant. A practical mixed implementation must be measured end to end.

Frozen evidence:
[`benchmarks/throughput-30b-halving`](benchmarks/throughput-30b-halving/README.md).
