# One-product lambda-state point halving

Store a point as `(x, lambda=x+y/x)` in the polynomial basis. If `Q=2P`,

```text
lambda_P^2 + lambda_P = x_Q
x_P^2 = x_Q (lambda_Q + x_Q + lambda_P + 1)
```

so repeated subgroup halving needs one polynomial product. The odd-subgroup
root is selected with the exact second-descent trace test.

On one RTX PRO 6000 the verified median is **23.904 B halvings/s**, 1.187×
the 20.134 B/s table-add walk but only 0.956 of 25 B/s. All 512 subgroup
points matched scalar multiplication by `2^-1 mod ell`.

Halving alone is a permutation, not a rho walk. Even a free-dispatch 50/50
mixture with the table-add rate has raw harmonic mean 21.858 B/s before
collision-constant accounting. This primitive therefore does not establish a
25 B/s effective rho walk.

`result.log` contains the build, differential check and timing output;
`result.json` freezes the boundary ratios.
