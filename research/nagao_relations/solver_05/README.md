# Direct quadratic function solver

The solver replaces a SAT solve inside each branch with two exact
Artin-Schreier solves. It resolved every target/mode slot in its three-size
panel, including all eight eleven-bit first decisions and all eight
complete projected enumerations. The maximum recorded all-phase time was
1.092205 seconds under the three-second budget. No correctness failure was
recorded. This is a small-panel engineering result, not an ECDLP breakthrough.

## Algebraic derivation

Write the residual cubic as H=X³+h2 X²+h1 X+h0 with

```
h2 = b²+b+r
h1 = a²+a*b+r*h2
h0 = b*c+r*h1,    c=r²+a*r+b*(r+s).
```

If H divides L_V, h2 is the sum of its three roots, hence belongs to V.
Enumerate h2 in V and solve `b²+b=h2+r`. Over an odd-degree binary field,
trace zero gives the two roots by half-trace, and trace one gives none.
Exclude b=0.

Choose one candidate root z in V, excluding 0 and r. Also exclude z=h2:
then the other two roots would have sum zero and coincide. Expanding H(z)
cancels the mixed a*b*r terms and gives

```
H(z) = (r+z)*a² + b*z*a + K,
K = H(z) evaluated with a=0.
```

Here r+z and b*z are nonzero. Put `t=b*z/(r+z)` and
`D=K/((r+z)*t²)`. Then `a=t*w` and `w²+w=D`. Another trace/half-trace
solve produces at most two a values. Thus neither a nor b is scanned over
the whole field. For each resulting function, check the exact cubic
condition `H divides L_V` and the previous exclusions before extracting
points and verifying their signed sum.

Every restricted relation is found: its h2 occurs in V, its nonzero b is
one of the first quadratic's roots, each factor abscissa is an eligible z,
and its a is one of the second quadratic's roots. Conversely, only functions
passing full support, exclusion and direct curve verification are admitted.
Function duplicates are retained in work counters and skipped before
repeated extraction; projected x-tuples are counted once per target.

For the first-d-coordinate normal-basis subspace, trace is nonzero on V.
Exactly half the h2 values yield two b solutions before excluding zero.
There are at most 2^d b branches, at most 2^d z branches each, and at most
two a candidates per branch. This upper bound describes candidate search,
not complete bit complexity: half-traces, inversions, support tests and
extraction still cost work. It supplies no improved generic counting yield.
This hybrid conditions on one root abscissa; it is not the original pure
root-free SAT formulation. No novelty claim is made.

## Validation and accounting

Complete projected oracle equality holds on all 43 affine five-bit targets,
covering both positive and negative cases. The 48 benchmark slots use the
exact solver_04 targets and modes, with no favorable reselection; every
complete result matches the exhaustive oracle. These targets were fresh
relative to solver_02/03, but became development data in solver_04.

Counted field API operations include search and support/extraction/verification,
with inversions expanded into primitives. Raw results include duplicate
functions and all target costs. Python control, allocation and coordinate
conversion are not field API operations; SAT work has no comparable count.
The combined report therefore treats wall time as preliminary and leaves
full-DLP S and rho ratios unmeasured. The generator checks its soft deadline
between yielded candidates; actual elapsed time must be checked for overruns.
None exceeded three seconds in this campaign.

The combined report ../goal_round_20260913.md compares chained S3, norm,
conditioned SAT, and both matched S4 controls. S4 wins at smaller sizes;
the three-size 20% success criterion remains unmet.

## Next falsifiable improvements

Move the a=0 residual norm out of the inner z loop; build L_V once per
instance; reuse inverses where their inputs are unchanged. Freeze this
optimized implementation before testing a new held-out panel. Charge all
precomputation and compare identical amortization policies for the SAT
controls. Then extend to larger supported odd field degrees and multiple
subspace dimensions. Reject a purported gain that disappears with setup,
failed targets, duplicate models or verification included.
