# Compiled support constraints: proof and attempted falsification

The root-free formulation is exact in the stated domain, but this SAT encoding
lost to the explicit-root norm formulation. This round contains 96 matched
solver trials and 992 exhaustive coefficient checks. No validation error was
recorded. It establishes no attack improvement.

## Four hypotheses and their verdicts

| ID | Idea | Evidence | Verdict |
|---|---|---|---|
| H1 | Encode support with `H divides L_V`, retaining only function coefficients | Algebraic argument below; all 992 coefficient pairs on one five-bit target matched the independent scalar remainder and exclusion checks; zero-b CNF rejected | Mathematical equivalence in the restricted domain; finite implementation validation passed, not a formal circuit proof |
| H2 | Fewer unknown field elements improve solver completion | Both support variants resolved 0/8 cases at each of n=9 and n=11; norm resolved 8/8 and 6/8 | Rejected for this circuit, solver, panel and budget; root elimination in general remains open |
| H3 | Condition two bits of b to help root-free search | Same completion counts; every n=9/n=11 conditioned trial exhausted its budget in branch zero | No benefit under this policy; other schedules and conditioning choices remain open |
| H4 | An advantage survives increasing size and total-cost accounting | Norm has a promising first-relation completion result at n=11, but strongest Semaev comparison and calibrated costs are absent | Unestablished; the frozen 20% all-phase improvement target is not met |

The contract was committed before execution. See [contract.json](contract.json),
[support.py](support.py), [run.py](run.py), [raw.jsonl](raw.jsonl), and
[summary.json](summary.json). The publication record maps the local source
snapshot recorded by the runner to its published equivalent.

## Matched results

All entries below are **verified decisions / four attempted instances**. A
decision is either a verified first relation or oracle-confirmed restricted
UNSAT. Uniform and known-decomposable strata remain separate; the latter is
not an estimate of natural yield. Targets may overlap across strata.

| Variant | n=5,d=3 uniform | n=5,d=3 supported | n=9,d=4 uniform | n=9,d=4 supported | n=11,d=5 uniform | n=11,d=5 supported | Classification | Full-DLP S / rho ratio |
|---|---:|---:|---:|---:|---:|---:|---|---|
| Chained S3 control | 4/4 | 4/4 | 4/4 | 4/4 | 0/4 | 2/4 | Engineering control | Unmeasured |
| RR norm | 4/4 | 4/4 | 4/4 | 4/4 | 3/4 | 3/4 | Engineering observation | Unmeasured |
| Root-free support | 4/4 | 4/4 | 0/4 | 0/4 | 0/4 | 0/4 | Engineering negative result | Unmeasured |
| Support, two-bit conditioning | 4/4 | 4/4 | 0/4 | 0/4 | 0/4 | 0/4 | Engineering negative result | Unmeasured |

At n=5 the uniform stratum contains two SAT and two UNSAT targets; at n=9,
three SAT and one UNSAT. At n=11 norm's uniform decisions are two SAT and one
UNSAT, with one timeout. All supported-stratum decisions are SAT. Every
timeout is unknown, never mathematical evidence of nonexistence.

The solver/search budget was three seconds per instance. Construction and
solver loading were recorded separately, not subtracted from that budget;
the raw diagnostic times include all three phases. Verification and root
extraction occur within the measured search phase, though a final verification
or solver return can overrun the deadline. Oracle setup is recorded separately
as validation infrastructure. No wall-time ratio is a full-cost speedup.

Removing field-level roots did not remove Boolean work: on the first uniform
n=11 target the source CNF has 1,831 variables for S3, 2,512 for norm and 7,972
for support. This is consistent with expensive remainder constraints; it does
not by itself establish the cause of solver difficulty. The conditioned
variant reuses one solver and one deadline; its first branch can starve the
remaining branches, as happened in every larger-field trial.

This experiment stops at the first verified relation. Earlier follow-up 01
used Hamming-weight bases and full projected enumeration. Those timings and
completion counts cannot be used as a matched comparison with this round.

## Restricted algebraic equivalence

Work on the nonsingular binary curve `y² + xy = x³ + 1`, with finite target
`R=(r,s)`, and require three distinct nonzero factor abscissas in V, all
different from r. Put

```
f = x² + a*x + c + b*y,           b != 0,
c = r² + a*r + b*(r+s).
```

Then f vanishes at -R. Its norm under `(x,y) -> (x,y+x)` is

```
N(X) = X⁴ + (b+b²)X³ + (a²+a*b)X² + b*c*X + c²+b².
H(X) = N(X)/(X+r).
```

The division is exact because R lies on the curve. Let
`L_V(X) = product over v in V of (X+v)`. This polynomial has simple roots
exactly V. Therefore a monic cubic H divides L_V if and only if it has three
distinct roots in V. The extra conditions `H(0)!=0` and `H(r)!=0` exclude
zero and the target abscissa. For each remaining root x, the equation f=0
recovers a unique y because b is nonzero. Substituting this y into the curve
equation gives N(x)/b², hence gives zero.

The pole of f at infinity has order four: x² has pole order four, y order
three, and x order two. Its four distinct zeros are the recovered points
and -R, so its divisor is `P1+P2+P3+(-R)-4O`. The elliptic-curve divisor
criterion then gives `P1+P2+P3=R`.

Conversely, a relation in this domain makes that degree-zero divisor
principal. Its function belongs to `L(4O)=span{1,x,y,x²}`. The x² coefficient
cannot vanish, since a nonzero function with pole order at most three cannot
have four distinct zeros. Normalize it to one. The y coefficient cannot
vanish either: a quadratic in x cannot vanish at four distinct abscissas.
Thus this relation has the stated coefficient representation, and its H is
the product of its three factor-abscissa terms, so H divides L_V.

Finally, the circuit computes successive `X^(2^j) mod H` by squaring and
reducing terms of degrees four and three using the monic cubic. Induction on
j shows its three output field elements equal `L_V mod H`. This proves the
recurrence, assuming correct field primitives and their CNF translation.
The 992 checks validate that implementation only at n=5,d=3 on one target;
they do not replace a proof of every compiler primitive at every field size.

These are standard norm, splitting and divisor facts applied to this model,
not a novelty claim. Primary-source context is in the parent study note.
Repeated abscissas, infinity summands and other charts are outside this claim.

## Ideas to pursue next

1. **Hybrid root retention.** Keep one explicit factor abscissa and impose
   subspace support on the residual quadratic. The current negative result
   motivates this middle ground. First prove soundness and completeness with
   the same exclusions, exhaust toy coefficients against the oracle, then
   compare first-relation completion on these exact targets. Reject this
   implementation if it gives no completion improvement over norm under the
   same budget. This experiment has not yet been run.
2. **Fair conditioning schedule.** Split a fixed budget across the four b
   branches, or use bounded conflicts before switching branches. Current
   branch starvation disproves no such schedule. Freeze the policy before
   running and charge every branch and restart; retain it only if paired
   completion improves. This is a successor to H3, not a reinterpretation
   of the failed policy.
3. **Strongest-baseline gate.** Put the repository's symmetrized S4 method
   behind the same subspace, target, arity and verification interface. Count
   construction, solving, recovery and rejected candidates in calibrated
   units. A 20% all-phase cost improvement at three increasing sizes is the
   existing goal; until this gate is met, norm is a promising control, not
   an established improvement over the best Semaev approach.

At fixed base and arity, exact encodings have the same set of genuine
relations. Changing equations alone cannot improve the counting yield.
Full attack costs additionally require relation dependence, matrix solving
and recovery, with rho measured in the same operation unit. No exponent fit,
new counting-boundary ratio, or full-DLP S is justified by these trials.

To rerun, copy contract.json, support.py and run.py into a new sibling folder
under research/nagao_relations, commit that snapshot, and run its run.py with
pycryptosat installed. The runner refuses to overwrite raw or summary files.
