# Decision: all rational degree-263 line duals admitted for transport controls

**PASS_ALL_LINES, with a scoped claim.** The [frozen producer](RESULT.json)
completed 264/264 canonical projective kernel directions on the first saved
Certicom twist-torsion basis and all six preregistered directions on the
independent second basis. The [fresh-process replay](VERIFY.json) returned
`PASS`. The protocol was pushed at `c3d1172a`, transparently amended before
the campaign after a developmental timing pilot, the final source pushed at
`f41cd690`, and the exact eleven-file [SHA-256 lock](LOCK.json) pushed at
`0292ceba` before the result. No random stream or selected post-result line
was used.

| Evidence | Producer | Independent replay |
| --- | ---: | ---: |
| Tested line directions | 264 + 6 | 270 |
| Nonzero forward/reverse twist-kernel inputs | 141,480 | 141,480 |
| Invalid kernel-abscissa controls | 540 | 540 |
| Fixed direct bit-polynomial reference lines | 12 | 12 |
| Other negative input controls | 20 | 20 |
| Distinct exact first-basis kernel sets | 264 | 264 |

For every tested line, the forward map and sign-corrected dual satisfy both
full-coordinate `[263]` identities on the public `P` and `Q`, and infinity
maps to infinity. The raw reverse quotient composed as `[-263]` on all 270
lines, so the interface consistently applies negation. All 262 nonzero
forward and 262 nonzero reverse twist-kernel points per line map to infinity,
while off-curve inputs with those same abscissae reject before the kernel
shortcut. The six second-basis exact kernel sets all occur among the first
basis's 264, independently checking the basis-line enumeration. The twelve
fixed controls reconstruct both generators and forward/reverse full
coordinates with the separate bit-polynomial group law and direct Vélu sums;
the independent replay repeats the complete all-line and exceptional-input
checks from the locked inputs.

The full first basis gives 263 distinct normalized codomain models: the two
horizontal directions `[1,39]` and `[1,64]` both return to source `b=1`,
and the other 262 are descending. This agrees with the prior
[endomorphism-order certificate](../../../ecc2k130_endo_ring_263_20260925/README.md):
source order discriminant `−7`, first descending leaf conductor `263` and
discriminant `−484183`. The rational endomorphism algebra is unchanged. The
changed order and lost cheap geometric Frobenius on a fixed leaf still give
no measured native PDP advantage. Transported and pullback factor bases
remain same-support/rank controls; native selection can differ in cost only
after its construction, action, solver, rank, and recovery are charged.

The producer took 598.296 seconds and peaked at 28,278,784 RSS bytes, below
the preregistered 2,100-second and 2 GiB caps. Its first-basis map setup
from the **supplied** torsion basis totalled 321.163 seconds across 264
lines; all first-basis construction and controls totalled 523.402 seconds.
The verifier took 327.334 seconds. These are certificate/replay timings,
not a cold factor-base construction or ECDLP cost. The producer receipt SHA-256
is `0506e506ddf4cc1580bca18821f872a65508f9ced7c53f5a648301af5e86262a`;
both processes agree on deterministic digest
`1d61faaa8501b80d0abbeb263d5ba4514bdc35f4b8c6b9c2dd4772b2d650718a`.

This closes the rational degree-263 forward/dual correctness and finite
exceptional-input gap for the frozen public curve and bases. It does **not**
support arbitrary nonrational extension-field kernels, certify a cold torsion
discovery pipeline, solve a degree-263 PDP instance, measure natural target
yield, establish relation rank, or recover an ECC2K-130 scalar. The next
high-arity exact-target gate remains the review- and label-gated
[m10 representation-capacity PR #937](https://github.com/aburan28/crypto/pull/937).
Only a passing capacity gate should be followed by a separately frozen
natural-target solver/yield and equal-useful-size native/original/transported/
pullback comparison. Preserve the matched automorphism-aware rho and fully
charged n37/n41/n53 log panels as the end-to-end reference; no n131 speedup
is set here.
