#set page(paper: "us-letter", margin: 0.65in)
#set text(size: 10pt)

= Isogeny and torsion fibers: a corrected cofactor statement
2026-10-12 · Structural correction · T-2070

The geometric kernel of a separable isogeny has size equal to its degree.
Only the rational kernel controls rational fibers. The previous fiber
report said the full degree d divides the cofactor h when the target prime
r does not divide d. That implication is false. The corrected statement is:
if K is the rational kernel and the curve order is h r, then |K| divides h.

Proof: |K| divides h r by Lagrange and divides d by the kernel order.
Since r is prime and does not divide d, |K| is coprime to r. Euclid's
lemma gives |K| divides h. Hence cofactor multiplication kills every
rational kernel ambiguity in a lifted relation. The arithmetic step is
checked by `lean/FiberCofactor.lean`; the geometric bridges are explicit
hypotheses of that Lean declaration and remain to be formalized.

#image("fiber.svg", width: 100%)

On E/F₁₃: y² = x³ + 2, the curve order is 19 and h = 1. A Frobenius-stable
order-3 kernel at x = 0 has only the identity over F₁₃. The standard
quotient formula maps to E': y² = x³ + 11. A native check verified all 19
images and every one of 361 rational ordered-pair sums. This is an exact
bounded control; full geometric formalization remains open.

#image("census.svg", width: 100%)

The frozen scan found 349 prime-order witnesses among 1,014 curves at 21
base primes. Two independent native point counts agreed for every curve.
These are selected-family counts, not a density estimate or DLP cost.

For a speed question, a destination coordinate solver, transport costs,
full relation and algebra costs, target recovery, and matched rho still
need a same-target comparison. The next experiment is recorded in
`problems/nonrational-isogeny-fiber-dlp-cost-gap.md`. Prior work on
small rational torsion and isogeny-class transfer remains separate:
Faugère et al. (EUROCRYPT 2014); Galbraith, Hess and Smart (EUROCRYPT
2002); Voloch, "The discrete logarithm problem on elliptic curves and
descents." Exact URLs and source links are in `REPORT.md`. The full
release-library test was initially unverified under shared-host
contention; see `validation/README.md` for the retained attempts and
the later controlled pass recorded in the paired follow-on.

#pagebreak()
= 2026-10-13 addendum: non-isomorphic pair

The next bounded control used `E/F617: y²=x³+1` and
`E′/F617: y²=x³+590`. These order-618 models are not isomorphic over
the base field: 590 is not a sixth power. The explicit degree-three
map has rational kernel `{O,(0,1),(0,616)}`. Its 206-point rational
image has three preimages per point, and the independent checker
verified the three-point fiber for each of eight public targets. The
target subgroup has order 103 and cofactor 6; both relation factor
bases were projected into it.

#image("nonisomorphic_pair_20261013/group_ops-v2.svg", width: 100%)

All 24 source-relation, destination-relation and matched-rho arm records
verified their recovered scalars. Counted group-operation totals over
the eight targets were rho 664, source relations 3,312, destination
relations 3,794; destination/source = 1.145531. The destination arm
used fewer operations on three targets and more on five. The first
timing wrapper failed after the native child finished, and its record
is retained; the direct v2 run and replay both exited zero and matched
byte for byte. The result PDF, exact identities, diagram, phase table,
independent checker and limits are in
`nonisomorphic_pair_20261013/REPORT.md`.

These are complete counted group operations for this small control,
with unpriced field, lookup, row, memory and construction work. An
end-to-end native-cost and scaling result remains open. Earlier
F13/F199 graphs remain frozen because their curves and targets differ.
