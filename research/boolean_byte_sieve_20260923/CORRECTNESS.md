# Producer argument for projection, recovery and compiled schedules

This argument and its tests cover generated Boolean systems of the declared
bounded sizes. They are not external review or a machine-checked proof, and do
not establish a production index-calculus or rho result.

## Projections are necessary conditions only

Packing equation values into bits gives a u32 word F(x), with F(x)=0 exactly when
every equation vanishes. The first filter evaluates bits 0–7 and the second bits
8–15. A nonzero projected word proves F(x) is nonzero. A zero projection proves
nothing about the other equations, so every surviving point receives a complete
evaluation before acceptance. Tests include a contradiction only in bit 31,
forcing both filters to pass every point while the complete check rejects it.

Projection commutes with GF(2) addition and multiplication by a Boolean scalar.
Thus the same exact quadratic-update schedule can maintain either full words,
bytes or transposed bit planes. No independence assumption about equations is
needed for correctness. Survival frequencies are measured, not inferred from
the number of selected equation bits.

## Coefficient updates and second-stage packing

Write the six low coordinates as y and the remaining coordinates as h. Then

    F(y,h) = C(h) + sum(i<6) L_i(h)y_i + Q_low(y).

C is quadratic in h, each L_i is affine in h, and Q_low is fixed. In reflected
Gray order, a transition flips j=ctz(step). Its high-only difference equals its
linear coefficient, its fixed lower-neighbor quadratic coefficient for j>0,
and the current higher-coordinate contributions. Flipping j updates only the
scheduled differences for indices below j, and changes L_i by Q_(i,j).

The second-stage implementation stores six eight-bit coefficients in disjoint
byte positions of one u64. Updating that packed word by XOR updates every
coefficient independently, without carries. Selecting byte positions by the low
assignment and folding by XOR computes the corresponding affine contribution.
The universal selection masks contain only coordinate patterns, not instance
solutions or preprocessing. All instance coefficient images are built inside
the charged context.

The discovery-only deferred implementation applies the XOR difference between
the current and last checked high assignments. Since each L_i is affine, even
numbers of intervening flips cancel. Direct tests compare this identity after
arbitrary gaps. The selected two-stage implementation does not depend on that
optimization; its discovery result remains separately archived.

## Assignment order and full-domain coverage

The retained method uses four low binary coordinates and a high Gray index t.
Write t=4T+r, with 0<=r<4. Its next two coordinates equal

    (r XOR (r>>1)) XOR (2*(T mod 2)),

and the remaining coordinates equal Gray(T). Visiting the four low subgroups in
that order preserves the first satisfying assignment exactly when moving from
16-point to 64-point blocks. Tests check the identity for all high steps through
the n24 bound and compare complete models with the retained enumeration.

The lowest four high-coordinate flip indices repeat a fixed fifteen-transition
sequence between group boundaries. Literal-index updates implement the same
XOR recurrence as the general loop. Group boundaries still use the general Gray
transition. Every individual block checks its cap before advancing; an incomplete
budget returns UNKNOWN. Tests compare compiled and looped models, work and traces
at every small cap and every unique-model placement in a ten-variable system.

## Representations and controls

NEON and SSE2 byte operations preserve the eight equation bits independently.
Each load and store addresses a complete 16-byte group. Zero-lane masks are checked
for all 64 positions, affine images spanning all constant byte values, and explicit
values with the sign bit set.
The bit-plane representation transposes the same 512 Boolean values: OR of its
eight truth planes marks every failing point, and complement gives the survivors.
It uses native NEON on AArch64 and a portable word implementation elsewhere.

Full-word controls keep all equation planes, and quiet controls remove enumeration
checksums. Quiet traces are serialized as null, not as zero-valued certificates.
Instrumented test modes still compare diagnostic traces. Retained methods keep
their previous traces; no trace-equality claim is made across different policies.

The same-order first models are checked on original equations. Partial scans in
search leaves use the retained prefix and scatter local coordinates back through
the same labels. Tests and evidence checks preserve prefix work and returned
models. Leaf quiet modes retain prefix diagnostics; only enumeration checksums
are removed.

## Measurement and claim boundary

Each random method order is followed by its reverse. The verifier independently
reconstructs every permutation and binds its seed to fixture identity. Equivalent
dispatcher and component paths must have identical models and semantic work.
Timing differences between them may reflect code layout or measurement context;
they are not evidence of different mathematical work.

Every computed lane, failed filter, survivor check, context construction and model
validation is charged. Partial and full evaluations are separate counters.
u32 coefficient-update counts and packed-u64 update counts remain separate and
are not a calibrated operation unit. Passing correctness checks or improving a
stage does not establish a cryptanalytic crossover.
