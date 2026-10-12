# Frozen protocol: n131 rotated arity and factor-base cardinality screen

Status: preregistered before any new arity's exact rational F0 or projected
column outcome. The source model, normal beta=3 and polynomial encoding are
those independently checked in [#762](https://github.com/aburan28/crypto/pull/762)
and the literal public-point import [#771](https://github.com/aburan28/crypto/pull/771).
The separate beta-3 m6,d21 exact census in draft
[#773](https://github.com/aburan28/crypto/pull/773) is a **reference to compare
only after its independent full replay is accepted**, not an input or fallback
for this experiment. Freeze this protocol, source, inputs, caps and hashes in a
draft PR with green hash-only CI before reading any new-arity count.

## Hypothesis and fixed arms

The fixed curve is `y²+xy=x³+1` over
`F2[z]/(z^131+z^13+z²+z+1)`, with prime subgroup order
`q=680564733841876926932320129493409985129`. For each arm use the first
`d` vectors `beta^(2^(m*j))`, j=0..d-1, for F0 and its m Frobenius-rotated
slots. The preselected arms are `(m,d)=(7,18),(8,16),(9,14),(10,13)`;
all have m*d≤131. Check that the full 131 beta conjugates have rank 131,
all m slot bases have rank d and combined rank m*d. These are four different
factor bases; do not substitute an easier beta or dimensional variant.

For each arm enumerate all `2^d` coefficient masks exactly once. `L` is the
number of nonzero F0 x with a rational y lift. The trace criterion is
`Tr(x+1/x)=0`, equivalently to the curve's Artin–Schreier lift condition.
The zero mask contributes the unique `(0,1)` point and projected O; the
normal-basis coordinate proof excludes x=1 from every F0. The physical F0
size is exactly `F=1+2L`. For each liftable nonzero x, compute
`u=(x+1/x)²`, `v=(u+1/u)²`; nonzero v identifies one canonical ±[4]
projected F0 column. Let C be the number of distinct v, excluding O.
Record the multiplicity of x per v (at most four), C, F, the exact
`N=F^m` ordered physical tuple count, and the **necessary** uniform-target
coverage ceiling `min(1,N/q)`. This counts no target PDP or rank. A column
can be further dependent under Frobenius and relations; C is not matrix rank.

The two-axis *follow-up screen* is `N>=q` **and** `C<=16384`. An arm passing
is a candidate for a branch-complete implicit PDP feasibility test, not an
attack or a relation-yield claim. Among passing arms choose the smallest m;
if none passes, stop this parameter family. The threshold corresponds to a
factor-base unknown count at most 1/64 of the approximate one-million-column
m6 design while retaining a nonvacuous 100% counting ceiling. Passing does
not bound solver time, actual support, relation rank, linear algebra, or rho
cost. Report all four arms even if one fails or caps out; a censored arm is
unknown and cannot be selected.

## Independent replay, caps and evidence

The producer uses #762's polynomial-reduction field and Euclidean inverse,
Gray-order x updates and a precomputed trace mask. It writes one canonical
`mask,x,lift,v` SHA-256 row digest, 8192-row chunk digests, counters, exact
column histogram, and cold setup/enumeration timing; it does not retain all
mask rows. The verifier reconstructs each x directly from its coefficient
bits using independent bit-serial multiplication and polynomial-division
Euclid inversion. It independently derives the trace mask, recomputes every
row and digest, and checks every exact count. For 64 SHA-selected literal
mask ordinals per arm, independently reconstruct rational y lifts by half
trace, check the full curve equation and two group-law doublings against v.
This is an arithmetic semantic sample atop the complete x/column replay.

Run arms in listed order, one producer then one verifier per arm, in separate
cold children; each producer has a 600-second wall cap, each verifier 1200 s,
both 256 MiB RSS acceptance caps. The runner also has a 30-second process
grace, preserves UTC intervals, command/exit, stdout/stderr, result hashes,
partial files and failure receipts. If a cap or check fails, preserve it and
do not replace that arm. Native operation counters and wall/RSS are stage
costs, never end-to-end S. Only after all four full independent replays pass
may a Pareto/follow-up decision be stated. Compare to #773's m6 full verified
count only if #773 has merged, and label any unverified number provisional.
Neither `F^m/q` nor a toy support rate transfers as actual n131 PDP yield.
