# 30 B/s lambda-halving attempt

The verified one-product lambda-state primitive reaches **28.953 B
halvings/s**, 0.9651 of the 30 B/s target and 1.158 of 25 B/s.

Selected changes:

- polynomial `(x, lambda)` state and one nonlinear product;
- evaluate the second-descent bit on `t=x_half^2` before square root, removing
  the second inverse-Frobenius permutation;
- closed-form `Tr(x H(x))` from Hamming weight;
- two-CLMAD phase-order prefix scan;
- CLMAD top correction and narrow reducer rebalance;
- one 768-thread block per SM.

All 512 test points match scalar multiplication by `2^-1 mod ell`. The
primitive is not itself a rho walk: halving is a permutation. Even a
free-dispatch 50/50 mixture with the 20.134 B/s table-add path has raw harmonic
mean 23.752 B/s before representation and collision-constant costs.

`result.log` preserves the compiler, differential test and three samples.
`result.json` freezes the configuration and boundary ratios.
