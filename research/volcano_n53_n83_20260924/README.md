# Characteristic-two Koblitz volcano audit at n=53 and n=83

Classification: **structural accounting**. This extends the bounded
prime-field metadata validation in PR #614 to the standard binary Koblitz
models; the previous Rust volcano backend does **not** support these fields.
The calculations are exact for the four named curves, while edge counts are
theoretical kernel counts, not a computational census. No relation collection,
Gröbner regularity, transported factor base, or ECDLP speed was measured.

We use `K_a: y² + xy = x³ + ax² + 1` for `a = 0, 1` over `F_(2^n)`.
Direct enumeration over `F_2` gives traces `t_1=-1` for `K_0` and `t_1=+1`
for `K_1`. The Frobenius recurrence is `t_0=2`,
`t_n=t_1*t_(n-1)-2*t_(n-2)`. For the auxiliary Lucas sequence `u_0=0,
u_1=1`, `u_n=t_1*u_(n-1)-2*u_(n-2)`, the identity
`t_n²-4*2^n = -7*u_n²` proves that the CM field has fundamental
discriminant `-7` and `f_pi=|u_n|`. Already at `n=1`, the Frobenius
order has discriminant `-7`, so `End(K_a)=O_(-7)` is maximal for both
curves and remains so over extensions. The base curve is on the surface;
another isogenous curve with the same trace need not be there.

| Model | Trace | Rational group order | `f_pi` | Prime factors and maximum local depth |
|---|---:|---:|---:|---|
| `K_0`, 53 | -56,619,371 | 9,007,199,311,360,364 | 68,476,319 | `68,476,319`, depth 1 |
| `K_1`, 53 | +56,619,371 | 9,007,199,198,121,622 | 68,476,319 | `68,476,319`, depth 1 |
| `K_0`, 83 | -6,151,469,093,347 | 9,671,406,556,923,184,866,742,756 | 347,450,761,417 | `6,473 × 53,676,929`, depth 1 at each prime |
| `K_1`, 83 | +6,151,469,093,347 | 9,671,406,556,910,881,928,556,062 | 347,450,761,417 | `6,473 × 53,676,929`, depth 1 at each prime |

Trial division through the integer square root certifies each stated prime
factor. For every other odd prime ell, `v_ell(f_pi)=0`, so every ordinary
curve in this isogeny class is locally on the ell-surface. At the listed
primes, endomorphism conductors of isogenous curves can *in principle*
divide `f_pi`; actual codomain orders have not been certified.

| Field size | Candidate degree ell | (-7/ell) | Horizontal kernels on maximal base | Descending kernels on maximal base | ell divides #K_a(F_(2^n))? |
|---:|---:|---:|---:|---:|---|
| 53 | 68,476,319 | +1 | 2 | 68,476,318 | No, for both `a` |
| 83 | 6,473 | -1 | 0 | 6,474 | No, for both `a` |
| 83 | 53,676,929 | -1 | 0 | 53,676,930 | No, for both `a` |

Here horizontal plus descending equals `ell+1` cyclic kernels *over the
algebraic closure*. Because `ell | f_pi`, Frobenius is scalar on
`E[ell]` and every such subgroup is Frobenius-stable, giving an isogeny
defined over the base field. Yet `ell ∤ #E(F_(2^n))` in every row, so
there is **no pointwise-rational generator of order ell**. A rational-point
kernel sampler would return empty and incorrectly suggest there are no
isogenies. This is especially relevant for `ell=6,473` at 83; exhaustive
enumeration of all 6,474 kernels is still a separate task. Because the
class number of `-7` is one, the two split horizontal kernels at 53 are
self-isogeny edges at the level of j-invariants; kernel multiplicity matters.

For the repository's existing `K_0/F_(2^53)` benchmark, the documented
subgroup order is `21,044,858,204,113`; the computed full order divided by
it is exactly `428`. This is an independent compatibility check against
`docs/ic/runs/koblitz-n53-one-unit-20260921.json`. No corresponding 83-bit
benchmark or frozen factor-base specification was found in the inspected
repository snapshot. The groups for `K_0` and `K_1` at the same n have
*different* orders, so they are not in the same rational isogeny class.

## What a fair index-calculus comparison would require

The mathematical transport invariant is conditional: for a separable
`F_(2^n)`-isogeny `phi:E->E'`, a target and an ordered factor base in a
subgroup on which `phi` is injective have exactly matching decomposition
tuple counts after transporting every point. The degree candidates above
are coprime to their respective full rational group orders, so a map of
one of these degrees, once actually constructed, would be injective on
rational points. This statement concerns point decompositions. Reusing the
same **x-coordinate rule** on the codomain defines a different factor base,
and Gröbner degrees depend on the resulting polynomial encoding.

An empirical study must construct and verify the actual binary-field maps
(including codomain, dual/composition certificates), distinguish transported
and native supports, hold subgroup/target/seed/cardinality fixed, and count
zero-yield instances and end-to-end cost against the same rho boundary.
Factor-base generation, relation collection, solver degree of regularity,
and scalar recovery are unmeasured for the neighbors in this report.
The n=53 repo workflow is a baseline on one surface curve, not a measurement
of floor-curve index calculus. The n=83 data here are structural only.

## Reproduce and limits

Run `python3 research/volcano_n53_n83_20260924/structural.py`.
The script writes `results.json`; `run.log` records the four rows. As
independent low-degree controls, it enumerates **every (x,y)** on both
binary models for `n=1,3,5,7` (eight model/field cases) and verifies
point counts against the recurrence. For `n=53/83`, the order is derived
from the Frobenius recurrence, not by point enumeration. The predicted
isogeny edge numbers use the standard CM volcano theorem and have not been
verified by constructing those edges. No conclusion about easier ECDLP or
solver regularity follows from conductor factorization alone.

The general volcano kernel-count result appears in A. V. Sutherland,
*Isogeny volcanoes*, §2.7–2.8, <https://arxiv.org/html/1208.5370v3>.
