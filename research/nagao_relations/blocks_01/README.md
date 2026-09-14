# Certified block rejection for the support equation

This experiment tests whether rejecting whole coefficient-search blocks saves
work against the identical bilinear solver with pruning disabled. The strongest
implemented Semaev comparator, the direct S3 pair table, remains in the matched
suite alongside the filtered hybrid, chained S3 SAT and symmetric S4 SAT.
The frozen [contract](contract.json) precedes timing. Results belong in
`RESULTS.md`; a correct rejection rule alone is not a performance result.

## Exact circuit and its scope

On y²+xy=x³+A x²+B, write the finite target as (r,s), with r nonzero.
The general norm construction and normalization are proved in
[the predecessor note](../bounds_01/NEXT.md). Put

    h = b²+b+r in V,   u=h+z,   t=r+z,   v=w²+u*w,
    b != 0, z,w in V, z not in {0,r,h}.

Multiplying the normalized support equation by b²r² gives

    F_b(z,w) = t²v² + b²*r*z*v + r²*z²*u² + b⁴*B = 0.       (1)

For fixed b and hence fixed h, expansion has just a constant, two binary-linear
maps and a binary-bilinear map:

    C       = b⁴*B,
    Lz(z)   = r²*h²*z² + r²*z⁴,
    Lw(w)   = r²*w⁴ + r²*h²*w²,
    M(z,w)  = z²*w⁴ + ((r²+h²)*z²+z⁴+b²*r*z)*w²
                         + b²*r*(h*z+z²)*w.

Every Frobenius power here is F2-linear. Therefore F=C+Lz+Lw+M is bilinear
in the binary coordinates of z,w, despite its higher degree over K. We keep
those maps as columns and update them by XOR when fixing z coordinates.
No summation polynomial or expanded subspace polynomial is constructed.
This is a coefficient/support hybrid: w is an auxiliary support witness,
so this prototype does not realize the stricter root-free objective.

At an exact leaf solve the binary linear equation in w, including its full
kernel. Recover

    gamma = t/(b*r), eta = z*u/b, D=(r+s)/r,
    a = gamma*v+r+b*D+eta.

The norm identity then yields roots z,w,w+u. Distinctness, nonzero roots,
exclusion of r, curve membership and the exact signed sum are all checked.
Every admissible decomposition has b nonzero: otherwise its normalized
pole-four function is a quadratic in x vanishing on four distinct abscissas.
Enumerating h in V and both roots of b²+b=h+r therefore covers every b.
The r=0 chart uses the unchanged predecessor, explicitly identified in output.

## Why a block rejection is sound

Fix the high coordinates of z to a prefix and leave k low bits free. Write
the resulting equation as

    C' + sum_i Z_i*z_i + sum_j W'_j*w_j + sum_ij M_ij*z_i*w_j = 0.

Relax each product z_i*w_j to an independent binary variable. This enlarges
the solution set. If C' lies outside the span of all nonconstant columns,
even the relaxation is inconsistent. Gaussian elimination supplies a binary
linear functional lambda with lambda(C')=1 and lambda(column)=0 for every
column. Applying it to the displayed equation proves 1=0 for every member
of the block. This certificate is sufficient for rejection; consistency is
never sufficient for acceptance. Full-rank relaxed spans cannot reject.

The implementation partitions z at fixed b, not both coordinates together.
The no-pruning control uses the same partition updates and exact leaves.
Every rejected block records b, h, prefix, free-bit count and lambda. The
independent checker derives columns by finite differences of unexpanded (1),
not by reusing the circuit expansion. It checks every certificate, including
ones from timed-out searches. Complete searches must partition all admissible
branches into visited and certified-rejected leaves without gaps.

## Boundaries, accounting and falsification

The existing signed-triple yield ceiling is unchanged:
Pr[decomposable uniform affine target] <= min(1,T_valid/(#E(K)-1)), with
T_valid excluding infinity and target-abscissa collisions. Exact masses for
these spaces are in [bounds_01](../bounds_01/RESULTS.md). Changing the solver
does not move that counting boundary. The reference cost is the measured
direct S3 table on identical targets and support. The old hybrid's seven
multiplications per branch is NOT a bound on this new circuit.

Circuit setup costs O(d²) field primitives per b. The complete no-pruning
control still visits Theta(2^(2d)) leaves when the Artin--Schreier fiber has
Theta(2^d) admissible b values. Block pruning has no proved asymptotic gain:
the number and rank of surviving relaxations must be measured. There is no
new universal ECDLP bound, fitted exponent, or rho comparison here.

Count additions, multiplications and squarings separately; inversion internals
are expanded, with inversion calls retained as diagnostics rather than added
twice. The sum called `fieldOperations` weights those three API calls equally,
not by a measured machine-cost conversion. Binary preimage/kernel/separator
XORs and AND-popcount-parity operations are separate. Internal field vectors
at these n<=30 sizes fit a 64-bit word. Allocation, indexing, bit tests and
control are uncalibrated; timing charges them but API counters do not. SAT
internals have no calibrated common counter. Summing partial vectors cannot
produce an honest total-operation speedup, full-DLP S or rho ratio.

All setup, failed blocks/leaves, extraction and signed verification are inside
the cold deadline. Cooperative deadline overruns remain marked. Independent
oracle construction and certificate replay are harness-only, never solver
inputs. First verified relation and complete enumeration are separate trials.
The 16 fixed inputs comprise eight predecessor inputs and eight fresh targets;
192 trials include all six variants, both modes, one repetition and three
seconds per trial. Uniform and supported targets remain separate. This is
a family-specific exploratory screen, not the broad WDSat regression or
full-DLP promotion gate described in the repository's AGENTS.md.

Any false rejection or incorrect complete result invalidates correctness.
Continue to a matched batch only with at least four complete equal-output
enumeration pairs, at least 20% fewer field API operations in aggregate and
no increase in separately counted binary word operations. Otherwise retain
the negative evidence and stop this candidate. Passing that screen does not
meet the original calibrated 20% saving across three increasing sizes.

## Reproduction

Use Python with pycryptosat installed. From the repository root:

    python research/nagao_relations/blocks_01/experiment.py --validate-only

The timing run uses `experiment.py` without arguments, in a separate checkout
of the frozen source revision recorded in `raw.jsonl` and `publication.json`.
Both the target-freezing command (`--freeze`) and the raw writer refuse to
overwrite evidence. `analyze.py` independently replays saved certificates,
checks matched coverage and hashes, and generates the summary and scoreboard
fragment without rerunning timing.
