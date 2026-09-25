# Degree-73 descent and Boolean stage diagnostic at (2^{37})

This experiment constructs a rational, separable degree-73 isogeny from

```text
E: y^2 + x*y = x^3 + 1 over GF(2^37)
```

and verifies a matched two-/three-summand Boolean workload on its codomain.
The field exponent is deliberately 37, not one of the eight-bit controls.

The exact structural certificate in
`research/koblitz_structural_fields_20260924/` gives conductor
`194399 = 73 * 2663` for the Frobenius order and maximal source endomorphism
order of discriminant `-7`. Since `(-7/73) = -1`, 73 is inert in
`Q(sqrt(-7))`. The rational 73-isogenies at this top vertex are descending;
the target order has discriminant `-7 * 73^2`. Sage's finite-field isogeny
enumerator must return all 74 base-field-defined cyclic kernels, and the run
records one map, its kernel polynomial, rational functions, target model, and
checks. The kernel subgroup is Galois-stable; its points are not assumed to be
individually rational over the base field. The squarefree degree-36 kernel
polynomial is checked to divide the degree-2664 73-division polynomial.

## Frozen protocol

`contract.json` was fixed before measurement. It records two seeds (one
development split and one holdout), four factor-base representatives from a
prime-order subgroup of order `230603167`, eight planted targets for each of
two and three summands, and ordered/canonical encodings. All support points and
targets on the codomain are the exact images under the same isogeny. Source and
codomain ground truth is checked by signed point addition. The executable also
checks the codomain equation after binary Weierstrass normalization, source and
target cardinality, the division-polynomial kernel certificate, and sampled
homomorphism pairs.

The Boolean systems use binary Semaev `S3` and `S4` summation equations in the
standardized model `y^2 + x*y = x^3 + A*x^2 + b`. Every solver root is enumerated and
compared with the independent point-addition workload; extra algebraic roots
are retained in the raw output. The Boolean F4/F5 reference has six variables
at most. The reported XORs count its matrix stage and criterion only.

The registered pass condition is at least a 10% F5 matrix-stage XOR reduction
for the isogenous model in every seed/summand/encoding cell, with no paired
increase in completion degree. The comparisons are diagnostic measurements,
not an attack: no full DLP is solved, no field-operation conversion is
calibrated, and full-DLP cost, `S`, rho/floor ratios and speedup remain null.
The two deterministic repetitions check counter replay; they are not
independent samples.

## Reproduction

Requires SageMath 10.6 and Python from Sage:

```sh
sage -python research/koblitz_isogeny_descent_37_20260925/experiment.py \
  --output /tmp/koblitz-37-run
python3 research/koblitz_isogeny_descent_37_20260925/verify.py /tmp/koblitz-37-run
```

The workflow runs the frozen contract in the official `sagemath/sagemath:10.6`
container and retains the compressed raw output as a workflow artifact.
