# Polynomial-state point halving

The original normal-basis primitive measured 16.053 B halvings/s. Keeping
state and both nonlinear products in the polynomial basis raises the verified
median to **17.921 B/s**, a 1.116× engineering gain.

The implementation converts only around the ONB-linear half-trace, square-root
and trace-pairing operations. Both products use `mulPolynomial131`.

This does not produce a faster rho walk. Halving alone is a permutation, and
the primitive remains only 0.890 of the selected 20.134 B/s table-add walk
before a non-permutation mixing branch is charged.

`result.log` preserves the build, 512-point CUDA differential check and three
timing samples. `result.json` freezes the ratios and configuration.
