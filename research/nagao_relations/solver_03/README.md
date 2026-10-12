# Hybrid root retention and fair coefficient conditioning

This round tests two successors to solver_02 without modifying its evidence.
The frozen contract and source snapshot precede execution. Targets, subspace
bases, target strata, seed, three-summand domain and three-second search budget
are unchanged. Chained S3 and RR norm are rerun as contemporaneous controls.

## Hypotheses

* H5: retain one explicit root and constrain the residual quadratic to split
  in V; this may reduce the cost of imposing support compared with the fully
  root-free cubic. A lack of paired completion improvement over norm rejects
  this implementation on the tested panel.
* H6: give each of the four coefficient branches a share of the remaining
  time instead of allowing the first branch to consume everything. A lack of
  completion improvement rejects this schedule on the tested panel.

Neither hypothesis predicts more mathematical relations: fixed factor base
and arity fix genuine target support. The complete-cost goal remains 20%
improvement at three sizes against the strongest matched Semaev baseline.
These diagnostic SAT trials cannot meet that gate; full-DLP S, ratio to rho
and calibrated solver operation counts remain unmeasured.

## Hybrid equivalence proof

Use the monic residual norm cubic H and coefficient chart of solver_02. For
an explicit z in V, perform monic synthetic division:

```
H(X) = X³ + h2*X² + h1*X + h0
u = h2 + z
v = h1 + z*u
Q(X) = X² + u*X + v
H(X) = (X+z)*Q(X) + H(z).
```

Impose `H(z)=0`, `Q divides L_V`, and `Q(z)!=0`, alongside the previous
`b!=0`, `H(0)!=0`, `H(r)!=0`. Because L_V has simple roots exactly V,
Q divides L_V precisely when Q has two distinct roots in V. The additional
Q(z) condition prevents z from duplicating either of those roots. Thus H
has exactly three distinct roots in V. The H(0) and H(r) conditions remove
the prohibited abscissas, and the norm/divisor argument from solver_02
recovers a valid signed relation.

Conversely, for any admissible cubic H dividing L_V, choose any of its three
roots as z. Its quotient Q divides L_V and does not vanish at z, so every
restricted relation is represented. This encoding has three choices of z
per admissible coefficient pair; it does not canonicalize that symmetry.
Future symmetry constraints must preserve at least one of these choices.

Modulo Q in characteristic two, `X² = u*X+v`. Consequently the square of
`A+B*X` reduces to `(A²+B²*v)+(B²*u)*X`. Repeated application computes
`X^(2^j) mod Q`; combining these remainders with the fixed linearized
coefficients of L_V gives its exact remainder. This is the implemented
quadratic recurrence, with two variable multiplications per squaring step.

The validation compares all six output field vectors (two remainder
coefficients, H(z), H(0), H(r), Q(z)) against independent scalar polynomial
arithmetic for all 7,936 choices `(a,b!=0,z in V)` at n=5,d=3 on one target.
It also repeats the previous 992 cubic support checks. This is finite
implementation validation, not a formal proof of the field/CNF compiler.

## Fair scheduling semantics

One persistent solver handles the four low-two-bit assignments of b in
order 0,1,2,3. Each branch receives remaining wall time divided by remaining
branches, sharing one deadline. An early UNSAT leaves more time for later
branches. An unknown branch is retained as unknown while later branches run.
Only four UNSAT branch results establish global UNSAT; any verified SAT
establishes SAT. A solver return or final certificate check can overrun the
soft deadline. Construction/loading and solve/verification times are recorded
separately; no uncharged restart budget is supplied.

## Evidence and reproduction

See contract.json, run.py, support.py, raw.jsonl and summary.json in this
directory. Every emitted SAT model is checked against CNF, curve arithmetic
and the exhaustive independent projected-relation oracle. Restricted UNSAT
must agree with the oracle. Timeouts remain unknown and all failed cases
are retained. Known-decomposable targets are reported separately from
uniform targets; overlap between strata is possible.

To reproduce without overwriting evidence, copy the three source/contract
files into a new sibling directory, commit the snapshot, and run its run.py
with pycryptosat installed. No new dependencies beyond solver_02 are required.

## Next algebraic route: condition the quadratic coefficient

There is an exact linear-algebra replacement for quadratic divisibility once
u is fixed. Define the F2-linear map `T_u(w)=w²+u*w` on V. For a nonzero
`u in V`, its kernel is exactly `{0,u}`, because `w²+u*w=w*(w+u)`.
Thus its image `W_u=T_u(V)` has dimension `d-1`.

It follows that `X²+u*X+v` has two distinct roots in V if and only if
`u in V\{0}` and `v in W_u`. For the forward direction, the sum of the
two roots is u, and either root maps to v. Conversely, the preimage of v
under T_u consists of a pair `{w,w+u}` in V, giving exactly the desired
roots. This proof holds for any F2-subspace V in the field; V need not
be closed under multiplication or squaring.

For each fixed u, Gaussian elimination on the images of a basis of V
therefore gives `n-d+1` independent binary linear membership equations
on v. This could replace the nonlinear modular-squaring circuit. In the
hybrid norm system impose `h2+z=u`, retain `H(z)=0` and all exclusions,
and substitute these linear membership equations. The tradeoff is up to
`2^d-1` coefficient branches and their setup cost. No improvement follows
from the proof alone, and this solver variant has not been run here.

The next experiment should freeze a shared budget across u branches,
charge image-space construction and every attempted branch, validate on
all toy coefficients, and compare against norm on the identical targets.
Discard it if branching erases the saving in the membership constraints.
This is a derived candidate, not a claim of novelty or a measured speedup.

## Completed results and verdicts

All 96 trials completed with zero validation errors; 992 cubic and 7,936 hybrid vector checks passed. Counts below are verified decisions / attempts, with uniform (U) and known-decomposable (D) strata separate.

| Variant | n5 U | n5 D | n9 U | n9 D | n11 U | n11 D | Class | Full-DLP S / rho |
|---|---:|---:|---:|---:|---:|---:|---|---|
| chained-s3 | 4/4 | 4/4 | 4/4 | 4/4 | 0/4 | 2/4 | Engineering | Unmeasured |
| rr-norm | 4/4 | 4/4 | 4/4 | 4/4 | 3/4 | 3/4 | Engineering | Unmeasured |
| support-fair | 4/4 | 4/4 | 0/4 | 0/4 | 0/4 | 0/4 | Engineering | Unmeasured |
| hybrid | 4/4 | 4/4 | 1/4 | 1/4 | 0/4 | 0/4 | Engineering | Unmeasured |

H5 is rejected against the norm control on this panel. Hybrid does recover
one uniform and one supported target at n=9 where the fully root-free
variants failed, but resolves neither stratum at n=11. This is a limited
improvement over root-free support, not an improvement over norm. H6 is
rejected on the tested panel: all 16 larger-field fair-schedule trials
made four calls across four branches and still timed out. Starvation was
fixed, but fixing it did not improve completion.

At n=5 all variants report two SAT and two UNSAT in the uniform stratum.
At n=9 both controls report three SAT and one UNSAT there; hybrid reports
one SAT and three unknowns. At n=11 norm's uniform results are two SAT,
one UNSAT and one unknown. All supported-stratum decisions are SAT.
Strata can overlap; the two hybrid successes are instance slots, not a
claim of two distinct uniformly sampled relations.

On the first uniform n=9 target, CNF variable counts are 1,357 (S3),
1,797 (norm), 4,287 (support-fair), and 3,130 (hybrid). Hybrid reduced
circuit size relative to cubic support without overtaking the controls.
Circuit size is not total solver work; no cost or asymptotic claim follows.
The contemporaneous controls reproduced the previous completion counts.
All recorded source hashes and all 96 aggregates were checked against raw
records. No hypothesis was revised after observing results.
