# Oriented binary Vélu transport

This package implements a full point map for odd cyclic kernels with rational
abscissae on ordinary binary curves `y²+xy=x³+a*x²+b`. It follows the structural
preflight merged in [PR #703](https://github.com/aburan28/crypto/pull/703).
The image includes the correct y coordinate; downstream relation collection
no longer needs to choose a sign by lifting an x coordinate.

This is a correctness prerequisite for paired factor-base experiments. It
measures no PDP yield, relation rank, complete ECDLP cost or speedup.

## Interface and limits

```python
from research.ecc2k130_oriented_transport_20260924.oriented_velu import BinaryVeluMap

phi = BinaryVeluMap.from_generator(source, kernel_curve, generator, degree)
image = phi(point)
assert phi.codomain.on_curve(image)
```

The curves use the existing `FastGF2m`/`Koblitz` interface in
`research/ecc2k130_relations`. The generator has exact odd order `degree >= 3`
and is rational on `kernel_curve`. That curve must have the same field and b
as `source`; its a may differ, allowing the degree-263 twist construction.
The constructor enumerates the cyclic kernel and rejects incorrect orders,
including non-exact composite orders. The kernel's abscissae form a frozen set
at `phi.kernel_abscissae`. General kernels with nonrational abscissae are out
of scope. The caller supplies a valid field implementation.

`None` is infinity. The map accepts canonical `(x, y)` tuples on its source,
maps infinity and valid kernel points to infinity, and rejects off-curve or
malformed points before considering an exceptional denominator. The source's
2-torsion point is handled by the ordinary formula. Input validation and a
codomain equation check remain enabled. This is variable-time research code.

`phi.then(psi)` composes maps with identical intermediate normalized curve
models. It does not infer a missing isomorphism or claim to construct a dual.
The validation constructs reverse quotients on tiny fields and records when
their composition equals `[degree]` or `[-degree]` after normalization.

## Formula and orientation

Let S contain one abscissa u for each pair of nonzero kernel points ±Q,
`t=sum(u)`, and `c=u/(x+u)`. The normalized quotient coefficients are
`a'=a`, `b'=b+t+t²`, and

```
X = x + sum(c+c²)
Y = y + sum(u*(x²+y)/(x+u)² + c³+c).
```

The formulas use field addition (XOR). They apply away from the kernel;
the exceptional kernel and infinity images are supplied explicitly.

To derive Y, pair Q and -Q in the full Vélu sum and then normalize the raw
codomain by `y -> y+t`. Writing `L=(y+v)/(x+u)`, the pair contributes
`c*L²+c²*L+c*(u+a)+c³+c`. Using the source equation and the equation of Q
reduces this to `u*(x²+y)/(x+u)²+c³+c`. For a twist lift `v+s*u`, the additional
terms cancel through `s²+s=a+kernel_a`; consequently the final expression
requires only u. It is independent of a choice of twist lift and determines
one oriented homomorphism, rather than an arbitrary lift of X.

The independent reference is the standard direct full Vélu sum implemented
with bit-polynomial arithmetic and a formal quadratic extension in
`research/ecc2k130_direction_review_20260924/redteam_velu_replay.py`. Tiny
rational-kernel checks also evaluate the direct group-law sums without using
the simplified Y expression.

## Evidence

The protocol was committed before execution in `a37f70768cd29da7e8ad19fdb5cb57a95582d623`.
Additional generic-interface controls were registered before execution in
`ecd1ba6873883bc25ce353884d2464e18e3b2935`.

| Control | Verified workload |
|---|---:|
| Fixed small fields n=3,5,7 | 10 quotient maps |
| Exhaustive small-field additivity | 72,292 full-point equalities |
| Non-unit b=7, composite degree 9 | 1,296 full-point equalities |
| Composite-degree direct rational Vélu oracle | 27 full-coordinate matches |
| Normalized reverse compositions | 8 full-group controls |
| Exact public target, both saved seeds | 8 representative maps |
| Independent exact-target direct Vélu oracle | 40 full-coordinate matches |
| Signed scalar transport | 56 checks, including 0 and -1 |
| Public P,Q image subgroup order | 16 checks |
| Actual degree-263 rational kernel inputs on twist | 2,096 |
| Malformed point and incompatible-parameter rejection | 122 checks |

All passed with zero disagreements. The final local validation took 12.67
seconds; this is a replay timing, not comparative performance evidence.
`validation.json` stores hashes, complete row counts, oriented target images,
command and host metadata. `validation_initial.json` preserves the earlier
pass before the additional composite-kernel controls; its source hashes
correspond to commit `30d892601eb8145cbacf213b9ddaafd338990d97` and the original
protocol. It is not the current validation receipt.

The generic formula is supported by the paired-sum derivation; finite tests
alone do not prove an arbitrary-input theorem. The exact-target replay checks
eight representatives, not every full-point map in the 264-line neighborhood.
A general dual-isogeny constructor is not implemented. The tiny composition
controls do not claim a degree-263 dual. No challenge logarithm is recovered.

## Replay

```sh
python3 research/ecc2k130_oriented_transport_20260924/validate.py \
  --out /tmp/oriented-velu-validation.json
```

The output path must not exist. Failures are retained as FAIL receipts and a
nonzero exit status. The focused GitHub workflow reruns the validation and
compares deterministic rows and all current source/input hashes to the saved
receipt. A downstream cost experiment must charge kernel setup, transport,
lifting, relation verification and all subsequent ECDLP phases separately.
