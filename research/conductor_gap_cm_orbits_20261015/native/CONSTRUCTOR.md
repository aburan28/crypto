# Codomain-only constructor: arithmetic to review

This describes the exact identities used by `sqrt_velu.cpp` and its batched
variant. It is an explanatory derivation, not a machine-checked proof or a
point-map certificate. The two implementations have separate binary receipts.

Let K have odd prime order ell, n=(ell-1)/2, choose B>0, and put
G=floor(n/(2B)). Baby points have odd indices 1,3,...,2B-1. Giant points have
indices 2B,6B,...,(4G-2)B. For each giant, adding and subtracting the B babies
covers the odd indices in one consecutive block. Together they cover
1,3,...,4BG-1, exactly 2BG representatives modulo sign. Direct addition covers
the remaining odd indices 4BG+1,...,ell-2. No representative is omitted or
counted twice. Baby and giant abscissae differ, so the inversions below exist.

For a giant point with abscissa X, write F=X^3+aX+b and A=3X^2+a. For a baby
abscissa x, put z=1/(X-x). The addition law gives two abscissae u ± sqrt(delta),
where

```
u     = X - A*z + 2*F*z^2
delta = -4*F*z + 12*X*F*z^2 - 4*F*A*z^3 + 4*F^2*z^4
```

Their first three power sums are 2u, 2(u^2+delta), and
2(u^3+3u*delta). Thus only the first six reciprocal power sums
S_j=sum_babies 1/(X-x)^j, together with S_0=B, are needed.
The baby product H(T)=product_babies(T-x), its first six derivatives and the
logarithmic-derivative/Newton identities supply these S_j. FLINT evaluates
the derivatives at the giant abscissae using product/remainder trees.
Batching changes the evaluation schedule, not the covered index sets or sums.

The separately versioned `sqrt_velu_batched_v4.cpp` expands these three moment
polynomials directly instead of multiplying arrays with many zero entries.
Its 4,096 coefficient checks compare every coefficient against the original
array products over both a 53-bit and a 96-bit field. It also delays reduction
inside short integer sums; GMP accumulation is exact and the final reduction
preserves the field result. Both forced-small-batch and full-batch reference
campaigns reproduce all 100 known codomains. These checks do not replace the
independent certificate required for a newly constructed upper-limit endpoint.

The input twist is normalized so that a0=d^2*a and b0=d^3*b. Scaling back
multiplies its kernel abscissae by d, so the kth kernel moment is multiplied
by d^k. For the resulting half-kernel moments M1,M2,M3, Vélu's formulas are

```
T  = 3*M2 + a0*n
W  = 5*M3 + 3*a0*M1 + 2*b0*n
a1 = a0 - 10*T
b1 = b0 - 14*W
```

All field operations use exact GMP/FLINT arithmetic. The floating-point square
root chooses a window size only; the integer G and direct remainder determine
exact coverage. Independent direct-reference equality and PARI endpoint
certificates check different parts of this construction. A matching trace,
ordinary discriminant D*ell^2 with prime ell, and j different from the class-
number-one maximal-order j identify the floor conductor in these families.
The CM classification implication is not claimed as formalized in Lean.

The output is a chosen-trace curve model. Full isogeny evaluation on arbitrary
points, transported generators and dual-map certificates are separate work.
