# Yield and work bounds for the frozen decomposition chart

These bounds describe uniform finite targets in the full curve group and the
specific three-distinct-abscissa chart. They do not describe arbitrary prime
subgroup targets, unrestricted decomposition, or every index-calculus method.
No novelty claim is made for the counting, trace, or meet-in-the-middle bounds.

## Exact group orders

Enumerating the F4 curve gives six points, including infinity. Its Frobenius
trace is -1. With s0=2 and s1=-1, the recurrence

    s_j = -s_(j-1) - 4*s_(j-2),    N_(2j) = 4^j + 1 - s_j

gives the orders of its extensions. This follows by taking powers of the two
roots of the characteristic polynomial T^2+T+4. The point-count and separability
identities are standard; see Sutherland's [point-counting notes](https://math.mit.edu/classes/18.783/2021/LectureNotes7.pdf)
and [Frobenius characteristic equation](https://math.mit.edu/classes/18.783/2022/LectureSlides8.pdf).
The script checks the F4 and F64 counts by enumeration, and records actual
ambient orders. Field degree is not prime-subgroup security size.

## Remove the excluded triples exactly

Let K be the number of admissible nonzero abscissas and F the corresponding
M=2K signed points. There are T0=8*C(K,3) unordered signed triples with
distinct abscissas. Define A_P(Q) as the number of unordered signed factor
pairs summing to Q, with neither abscissa equal to x(P). Pair multiplicities
must be retained: distinct sign choices can have the same sum.

    Z = (1/3) * sum_(P in F) A_P(-P)
    X = sum_(P in F) A_P(-2P)
    T = T0 - Z - X

Z counts triples summing to infinity, each seen at each of its three factors.
For a finite sum R with x(R)=x(P), either R=P or R=-P. R=P would force the
other two factors to be inverses, contradicting distinct abscissas. Thus only
R=-P is possible, giving the equation Q1+Q2=-2P. The distinguished factor is
unique, so X has no extra division. Z and X are disjoint exclusions.

T counts signed triples, not unique projected x-tuples. If r(R) is the number
of admissible signed triples summing to R and N is the curve order, then

    E_uniform r(R) = T/(N-1)
    Pr_uniform[r(R)>0] <= min(1, T/(N-1)).

The expected number of projected relations is at most that same mean. This
bound needs no uniform-sum or independence heuristic. An upper bound of one
does not establish high success probability. For independent uniform target
attempts, mean attempts to encounter a decomposable target are at least 1/U
when U is any valid success ceiling; solving cost is an additional question.

## A trace-fiber ceiling

Let k divide n and let the curve be defined over F_(2^k). Write phi for
2^k-power Frobenius and ell=n/k. The group trace is

    tau = 1 + phi + ... + phi^(ell-1): E(F_(2^n)) -> E(F_(2^k)).

It is surjective. Indeed, (phi-1)*tau=phi^ell-1; both difference maps are
separable, so tau is separable too. Its degree, and hence kernel size, is
N_n/N_k. Its kernel lies in E(F_(2^n)), making the image size exactly N_k.
Each trace fiber therefore has h=N_n/N_k points. The finite-target capacities
are h_c=h-1 for the zero trace class and h_c=h otherwise.

For each abscissa x choose one lift P_x. In the integer group ring of the
small curve group, compute the coefficient of t^3 in

    product_x [1 + t*(e_(tau(P_x)) + e_(-tau(P_x)))].

This gives all signed distinct-abscissa triple counts by target trace, using
a degree-three dynamic program. Subtract Z from the zero class and subtract
each target-collision contribution in class tau(-P). Let the resulting counts
be T_c. They sum exactly to T and give the strengthened ceiling

    U_trace = sum_c min(h_c, T_c) / (N_n-1)
            <= min(1, T/(N_n-1)).

A zero T_c certifies that every target in that fiber is impossible in this
chart. A positive T_c does not establish any particular target's solvability.
The script tests k=6 on all ten frozen bases; at n=6 the trace is the identity,
so this ceiling must equal the exhaustively measured support probability.
This is an exact group-trace constraint. It is not an unproved identification
of coordinate trace with point-group trace or a fixed-target orbit quotient.

## The hybrid still visits its conditioned-root branches

For target abscissa r set Z_r=V minus {0,r}. The current candidate generator
visits precisely

    B(r) = sum_(h in V, Tr(h+r)=0)
              (2 - 1_(h=r)) * (|Z_r| - 1_(h in Z_r))

conditioned-root branches on complete enumeration. The absolute field trace
characterizes the two roots of b^2+b=h+r; the correction excludes b=0.
The inner loop excludes z=h. Every surviving branch evaluates a cubic and
performs three additional multiplications before its Artin-Schreier solve.
The frozen Horner routine makes four multiplication calls, giving a lower
bound of 7*B(r) field-multiplication calls, before setup, recovery and checks.
The parity filter executes after this work, so it cannot reduce B(r).

When absolute trace is nonzero on V and r is outside V, with s=|V|,

    B(r) = s*(s-2) + 2*1_(Tr(r)=0).

If trace vanishes on V, there can instead be no branches for half the target
trace classes; the general formula covers this exception. The script checks
the formula by instrumenting the original generator on exhaustive tiny targets.
This is an implementation bound, not a lower bound on all function solvers.

## The pair-table floor and a conditional architecture barrier

The current S3 table visits C(K,2) unordered pairs and stores K*(K-1) entries.
Each pair uses seven field-multiplication calls. A complete target query visits
M signed factors, except that two are skipped when x(R) belongs to the base.
Every residual addition is between distinct abscissas, so it makes at least
two multiplications, excluding the additional inversion work. For b complete
queries with J targets whose abscissas belong to the base,

    multiplications >= 7*C(K,2) + 2*(b*M - 2*J).

This is a deliberately conservative floor in a single primitive, with no
mixing of multiplication, XOR, SAT or curve-addition units. The saved full
batch multiplication vectors are compared against it. It is not a calibrated
total-operation floor, and its ratio is not a speedup or a rho ratio.

For each exact saved eight-target batch, the script also sums the hybrid
floor 7*B(r) over those same eight targets. Comparing that lower bound to
the S3 table's measured total multiplication count can rule out a hybrid
multiplication-count advantage on complete enumeration without pretending
that a timed-out hybrid actually completed. The floor applies to both
hybrid variants. It says nothing corresponding about their first-hit costs.

For uniform affine targets the expected signed-factor visits per complete query
are M - 2*M/(N-1). Since mean projected output is at most T/(N-1), the long-run
query-only multiplication/output ratio is at least

    2*(M - 2*M/(N-1)) / (T/(N-1)),

provided T>0. This is a ratio of expectations for the stated sampling regime;
it is not the expectation of a random finite ratio. Finite batches additionally
pay setup. Supported targets cannot be substituted for uniform ones here.

Suppose this architecture is used to collect at least c*M independent rows
for constant c>0, from independent uniform targets, fully scanning the base
on every query. Even optimistically counting every signed triple as a usable
row gives expected target count Omega(N/M^2). The required query work is then
Omega(N/M), while explicit pair setup is Omega(M^2). Thus

    work >= Omega(M^2 + N/M) >= Omega(N^(2/3)).

The last inequality follows by minimizing at M=Theta(N^(1/3)). It already
ignores rank failures and final linear algebra. This rules out a square-root
exponent for this complete-enumeration architecture under these assumptions.
It does not rule out all Semaev, Nagao, trace-zero, first-hit, target-adaptive,
or subfield index-calculus algorithms. Our measured full groups are composite;
no measured Pollard-rho boundary or complete ECDLP claim follows from this note.

## What would escape these bounds?

The function-first route needs to skip whole coefficient blocks before the
conditioned-root solve, not just reject functions after paying for them. The
pair-table route needs a different time/space/query tradeoff, a justified
target distribution, or fewer unknowns/independent rows. A trace filter could
help only if its construction and application cost are smaller than the
work it removes. Whether the small-group constraints yield useful conditions
directly on function coefficients remains an unproved next hypothesis.

This round changes accounting and the boundary description. It does not
change either solver, claim a performance improvement, or satisfy the broad
regression and full-pipeline promotion gates.
