# Developing a relation encoding without summation polynomials

The implemented starting point is a Nagao-style Riemann–Roch encoding. The
development candidate is **function-coefficient search with an exact factor-base
support constraint**: find a function first, recover its zeroes only after its
residual norm passes the support test. Neither formulation constructs Semaev's
summation polynomial.

This is implementation and toy correctness work, not an ECDLP speedup result.
The target repository is `aburan28/crypto`; the existing Semaev frontend remains
the comparator. The frozen scope is
[`experiments/nagao_relation_contract.json`](experiments/nagao_relation_contract.json).

## Which axis changes?

| Approach | What it changes | Treatment here |
|---|---|---|
| Nagao / Riemann–Roch | Relation representation, using a principal divisor | Implemented incidence and norm encodings |
| Function coefficients + support divisibility | Search variables and how factor-base membership is imposed | Implemented exhaustive coefficient-search prototype |
| Vector-subspace or low-Hamming-weight factor bases | Admissible points | Held fixed within each comparison |
| Fewer summands, including Joux–Vitse variants | Decomposition arity and success probability | Arity fixed at three |
| SAT, Gröbner, XL, triangular methods | Solver | Same SAT solver for the two encodings and chained S3 |
| Elliptic nets | Another algebraic presentation | Related to Semaev; not an independent hardness shortcut [3] |

A vector-space factor base is not the only possible factor base at prime
extension degree: this repository already uses normal-basis Hamming-weight
sets. Prime extension degree removes nontrivial intermediate subfields, not
arbitrary subsets or vector subspaces. Under the usual roughly-uniform-sums
heuristic, a base of about q signed points in a group of size about q^n has
(n−1)-point decomposition probability about 1/(q(n−1)!), with multiplicities
and symmetries accounted for separately. That heuristic is not a universal
measured probability or a guarantee that the cheaper solve wins.

## Exact genus-one construction

Work on the repository's binary curve

    E: y² + xy = x³ + 1 over K = GF(2^n).

Write R=(r,s). Its inverse is −R=(r,r+s). The space L(4O) has basis
1,x,y,x². On the monic pole-order-four chart take

    f(x,y) = x² + a x + c + b y,
    c = r² + a r + b(r+s),       b ≠ 0.

Thus f(−R)=0. If f also vanishes at three affine points P_i whose abscissas
are pairwise distinct and different from r, its four distinct zeroes exhaust
the pole degree. Consequently

    div(f) = (P1)+(P2)+(P3)+(−R)−4(O),
    P1 + P2 + P3 = R.

Conversely, any such decomposition gives this normalized function: the divisor
is principal, its pole order is exactly four, and the x² coefficient is nonzero.
Its y coefficient cannot vanish on four distinct abscissas because a quadratic
in x has at most two distinct roots.

The first frontend keeps x_i,y_i,a,b. It imposes the three curve equations and
the three incidence equations f(P_i)=0. In binary coordinates, squaring is
linear; the arithmetic incidence and curve equations are Boolean quadratic.
The normal basis makes squaring a bit permutation. This degree statement does
not include weight/cardinality constraints and does not bound solving degree.

The second frontend eliminates y locally. With A(X)=X²+aX+c, compute

    N(X) = A(X)² + b X A(X) + b²(X³+1)
         = X⁴ + (b+b²)X³ + (a²+ab)X² + bcX + c²+b².

Match its coefficients to

    (X+r)(X+x1)(X+x2)(X+x3).

Recover y_i=A(x_i)/b, then independently check every point and its signed sum.
This keeps a,b as unknowns; eliminating those coefficients is precisely the
kind of route back to a symmetrized summation polynomial illustrated in [2].
Equivalent solution varieties can still have different circuit sizes and solver
behavior. No equality of runtimes follows from algebraic equivalence.

Both frontends explicitly exclude x_i=0 to match the existing factor base,
and accept valid finite targets including the two-torsion target (0,1).
Repeated abscissas, repeated points, and infinity summands are outside this
first chart. Merely repeating an evaluation-matrix row does not impose a
double zero: a future complete implementation needs local multiplicity/Hasse-jet
conditions. Restricted UNSAT is never reported as unrestricted PDP UNSAT.

## The development candidate: constrain support without root variables

For a d-dimensional binary subspace V of K, define its subspace polynomial

    L_V(T) = product over v in V of (T+v).

It is linearized: only T,T²,T⁴,...,T^(2^d) can occur. Starting at L_{0}(T)=T,
adding an independent basis vector v uses

    L_{W+<v>}(T) = L_W(T)² + L_W(v)L_W(T).

Set H(X)=N(X)/(X+r), which is monic cubic. Impose

    H divides L_V,     H(0) ≠ 0,     H(r) ≠ 0.

Because L_V is a product of distinct linear factors, this is equivalent to H
having three distinct roots in V, none zero or r. Since b≠0, each root has a
unique y from f=0. The norm identity puts that point on E. The incidence at
−R fixes the required sign. Roots are extracted only after admission and the
group sum is checked independently.

The prototype evaluates L_V modulo H using modular squaring and the sparse
linearized coefficients. It never expands a degree-2^d dense polynomial, builds
S4, or introduces x1,x2,x3 as search variables. The initial search simply scans
all a and nonzero b, so it costs q(q−1) coefficient candidates per target.
That is a completeness instrument, not a scalable solver: at n=131 it would
mean approximately 2^262 candidates. In a subspace parameterization three
explicit x variables need only 3d bits; two field coefficients need 2n bits.
Removing root variables therefore does **not** necessarily reduce the number
of unknown bits.

The exact support constraint is an independently derived adaptation of standard
subspace-polynomial algebra. We have not established novelty for its use in
this setting. The research question is whether it yields a cheaper practical
solver than existing Semaev formulations after all costs are counted.

## Boundaries, controls, and falsification

The invariant is the exact set of supported targets at fixed factor base and
arity. If the signed base contains B points, there are at most binomial(B,3)
unordered distinct-point triples. Hence, for a target uniform in the full
group, the success probability is bounded above by
min(1,binomial(B,3)/#E). Conditional on a uniformly sampled affine target, use
#E−1 in the denominator. Our distinct-x restriction can only reduce the count.
This is a counting ceiling on yield, not a runtime lower bound for algebraic
PDP solvers. A changed encoding cannot improve this count.

The implementation gate is zero missing or invalid projected decompositions,
with exhaustive direct group arithmetic as the reference. Every Semaev model
is lifted and verified; any spurious model is blocked and its work retained.
Every RR model is checked against its function certificate and the group law.
Enumeration ends only on solver UNSAT after blocking projected tuples; a
watchdog leaves that comparison incomplete. No sample is dropped.

IR bit-operation counts measure evaluating a constructed arithmetic circuit.
CNF size measures representation size. Neither measures SAT search work. The
coefficient-search prototype reports field-operation categories by phase,
including construction, support checks, extraction and group verification;
these are not silently converted to affine additions. Wall times are diagnostic.

The full attack boundary remains Pollard rho on the same subgroup. The columns
S=all attack operations/sqrt(subgroup order), S/S_rho, and a ratio to an
applicable runtime floor remain **unmeasured** here. Relation-rank collection,
linear algebra and descent have not run. This stage cannot populate a measured
attack row or support an exponent/crossover claim.

## Run the implementation checks

Python 3.10+ and the existing code generator are sufficient for coefficient
search. SAT checks additionally require `python-sat` and `pycryptosat`.

```sh
python -m pip install python-sat pycryptosat
python ecc2k130/codegen/testnagaodecomp.py --solver
python ecc2k130/codegen/testnagaoannihilator.py
python ecc2k130/codegen/nagaocompare.py --out /tmp/nagao-comparison-new.json
python ecc2k130/codegen/nagaoannihilator.py --field-degree 5 --subspace-dimension 3 --out /tmp/nagao-support-new.json
```

The comparison follows the committed contract: all 43 affine targets over
GF(2^5), at weights one and two; 24 uniformly sampled affine targets over
GF(2^9), at weight one. Its outputs contain complete projected solution sets,
input coordinates, source hashes, solver/environment versions, emitted CNF
hashes, and failures. Outputs are written to new paths.

## Next experiments

1. **Compile the support-divisibility condition into the shared IR.** Replace
   exhaustive coefficient enumeration by SAT or a solver that preserves the
   coefficient blocks. Compare incidence, norm-with-roots, and norm-with-support
   on the same V, n, arity, targets, timeout and memory policy. Price construction,
   failed targets, support tests, extraction, verification and repeated use.
   The frozen engineering success target is total verified-relation cost at
   most 0.8 times the best matched Semaev baseline across three increasing sizes.
2. **Exploit coefficient structure before elimination.** Test conditioning one
   coefficient block, incremental row reduction and modular-remainder circuits.
   Include the best symmetrized S4 frontend, not only chained S3, before claiming
   competitive performance. More coefficients or quadratic source equations
   alone do not establish lower solving degree.
3. **Add multiplicity-aware charts and then orbit-aware support.** Compare
   supported-target coverage before and after the additional charts, keeping
   their costs visible. Frobenius transforms R as well as the summands; do not
   quotient individual point orbits while silently holding R fixed. Higher-genus
   divisors are a later extension, with genus-specific bases and reconstruction.

## Observations from the frozen correctness panel

The raw SAT artifact is
[`nagao_relation_comparison_20260913.json`](experiments/nagao_relation_comparison_20260913.json).
Its aggregate `valid: false` means the full panel did not complete: all 24
nine-bit incidence instances reached the ten-second watchdog. There were no
incorrect admitted RR relations. Do not relabel this artifact as a fully passed
campaign. The exact per-variant observations are:

| Variant | Complete target cases / attempted | Watchdogs | Missing or invalid results in completed cases | Full-attack S / rho ratio | Class |
|---|---:|---:|---:|---|---|
| Chained S3 + verified lifts | 110 / 110 | 0 | 0 | Unmeasured / unmeasured | Reference implementation |
| RR incidence | 86 / 110 | 24 | 0 | Unmeasured / unmeasured | Engineering candidate; incomplete panel |
| RR norm | 110 / 110 | 0 | 0 | Unmeasured / unmeasured | Engineering candidate; correctness passed |

Both complete backends recovered the same 947 **projected x-tuples**, summed
over the 110 target cases. The Semaev frontend also produced five non-liftable
models in the n=5, weight=2 panel; they were rejected, blocked and retained in
the raw record. Incidence found six valid projected tuples during the timed-out
cases, but those partial lists are not complete solution sets or UNSAT evidence.

The unoptimized arithmetic IR at n=9 requires 1,587 bit operations per
evaluation for chained S3, 2,847 for incidence and 2,078 for the norm encoding.
These count arithmetic circuit evaluation only. Thus even the norm version
does not win by simply containing fewer gates. Diagnostic solve/enumerate/verify
wall times are recorded, but no SAT operation-count conversion or speedup
headline is inferred from this single toy panel.

The function-first coefficient search uses a different factor-base panel;
its outputs count **signed point triples**, not projected x-tuples:

| n | dim(V) | Signed factor-base points | Affine targets checked | Supported targets | Verified signed triples | Exact agreement with oracle |
|---:|---:|---:|---:|---:|---:|---|
| 5 | 2 | 4 | 43 | 0 | 0 | Yes |
| 5 | 3 | 10 | 43 | 33 | 72 | Yes |

The dimension-two base has only two nonzero curve abscissas, so it cannot
supply three distinct ones; this is a useful negative control, not evidence
against RR search. Each row searched 43×32×31=42,656 coefficient pairs.
Sources: [`n5_d2`](experiments/nagao_annihilator_n5_d2_20260913.json) and
[`n5_d3`](experiments/nagao_annihilator_n5_d3_20260913.json). Their separate
addition/multiplication/squaring ledgers show that exhaustive coefficient search
is not yet an efficient implementation; its purpose is validating the support
certificate before compiling a search solver.

Independent review reconstructed the interpolation equations, examined all
10,640 signed triples with distinct nonzero abscissas over GF(32), and found
9,620 in the declared target domain. Both IR frontends passed all 19,240
evaluations; all 9,360 applicable target-sign mutations were rejected.
Reproduction: `python ecc2k130/codegen/testnagaoreview.py --help`;
[`raw review`](experiments/nagao_independent_review_20260913.json).

The frontend suite passed 10 tests with CryptoMiniSat; the coefficient/support
suite passed five tests. The latter includes full oracle equality over all
affine targets at (n,d)=(3,2),(5,2),(5,3), norm identities, subspace-polynomial
evaluation, repeated-root rejection, and accounting checks.

The contract was committed as `93a09e5` before the campaign; sources and raw
results were snapshotted as `0cf9eaf`. Source hashes in the outputs identify the
exact implementation, including files then untracked. A later maintenance
change preserves partial outputs on exceptions and corrects an ordering
comment; it changes no equations or saved observations. The initial runner
launch failed before any target evaluation because its contract path ascended
one directory too far; the path was fixed before the recorded campaign.

The study note and frozen experiment files were subsequently moved together
under `research/nagao_relations/` at the user's request. Frozen JSON contents
retain their original command paths and source hashes; those describe the
historical run. Current runners load the relocated contract. Run the commands
above from the repository root.

GitHub publication preserves the snapshot trees while assigning new commit
identifiers. The published source/results snapshot is
[`6834d5c`](https://github.com/aburan28/crypto/commit/6834d5c7e44c639b423cb6bdf879933c939528a8).
The [publication record](publication.json) maps the original local IDs to their
published equivalents and retains the original commit objects for reproducing
the identifiers recorded by the historical runs.

The current priority is to compile the exact H-divides-L_V condition and compare
it against the norm backend. Direct incidence remains available as a control
and for alternative solvers; its 24 timeouts apply only to this encoding,
solver and watchdog. No general mathematical avenue is closed.

## Retrieved primary sources

1. Antoine Joux and Vanessa Vitse, *Cover and Decomposition Index Calculus on
   Elliptic Curves Made Practical*, EUROCRYPT 2012, §§2.2–3.1:
   <https://link.springer.com/content/pdf/10.1007/978-3-642-29011-4_3>.
   Describes Nagao's RR/norm approach, the separate splitting check, and reports
   its earlier elliptic formulation as less efficient than Semaev. This supports
   benchmarking the alternative, not assuming it is faster.
2. Vanessa Vitse, *Summation polynomials and symmetries for the ECDLP over
   extension fields*, DLP 2014, slide 10:
   <https://www-fourier.univ-grenoble-alpes.fr/~viva/research/talks/Vitse_talk_DLP14.pdf>.
   Constructs summation polynomials by eliminating RR-function coefficients.
3. Tsunekazu Saito, Shun'ichi Yokoyama, Tetsutaro Kobayashi and Go Yamamoto,
   *Some relations between Semaev's summation polynomials and Stange's elliptic
   nets*, Journal of Math-for-Industry 3 (2011), 89–92:
   <https://api.lib.kyushu-u.ac.jp/opac_download_md/19584/JMI2011A-9.pdf>.

All three sources were retrieved and read during this implementation session.
