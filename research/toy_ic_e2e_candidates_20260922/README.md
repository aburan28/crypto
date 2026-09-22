# ECC2K-130 index-calculus candidates: bounded end-to-end laboratory

This round assembles and measures a complete **toy** index-calculus pipeline.
It accepts two internally generated Koblitz fixtures over GF(2^5) and GF(2^9),
with prime subgroup orders 11 and 127. It is not an ECC2K-130 implementation,
and it supplies no runtime extrapolation to GF(2^131).

The contribution is integration and accounting: the existing Semaev SAT
encodings now run through the same miniature full-DLP pipeline as the Nagao
coefficient solvers. No new decomposition algorithm is claimed. The frozen
contract predates the measurements. See [RESULTS.md](RESULTS.md) for the data.

## Research shortlist

These are candidates for comparative research, not a demonstrated ranking at
ECC2K-130. The literature search did not establish a verified end-to-end
index-calculus result beating the automorphism-aware rho reference there.

| Candidate | Why retain it | Evidence and limitation | This round |
|---|---|---|---|
| Semaev S4 / chained S3, Weil descent, SAT with native XOR reasoning | Directly relevant to binary-field point decomposition; published implementations and reproducible polynomial benchmarks exist. | Trimoska–Ionica–Dequen study SAT and XOR reasoning [1–3]. A polynomial solution still needs rational-point lifting, a verified group relation, independent rows and scalar recovery. | Both encodings connected to CryptoMiniSat and the full toy pipeline. This is **not WDSat**, and it does not reproduce the paper's complete solver configuration. |
| Semaev with Gröbner/F4, symmetry and bounded hybrid guessing | Essential algebraic comparator to SAT; matrix structure and solving degree decide feasibility. | The survey [4] reviews this family. A low first fall degree does not establish a low solving degree; [5] explains why ECDLP descent systems require care. | Not integrated or timed in this matched run. The existing Rust panels use different inputs and are not substituted for a matched F4 result. |
| Nagao / Riemann–Roch coefficient decomposition | Changes the representation of the point-decomposition problem. The existing repository has exact small-field implementations. | Riemann–Roch and summation-polynomial constructions are mathematically related [6]. The repository's quadratic-image, coefficient-pullback and cached-pullback solvers are experimental hybrids, not a proven subexponential algorithm. | All three connected, independently checked and measured. |
| Frobenius-aware factor bases and relation compression | Koblitz structure makes this relevant across decomposition methods. | This is a supporting technique rather than an independent decomposition solver. A chosen support must be compatible with the action. Frobenius-transformed rows are not automatically independent. Rho benefits from automorphisms too [7]. | Reviewed; not implemented in this toy comparison. The current coordinate subspace is not assumed Frobenius-invariant. |
| Quasi-subfield factor bases | A theoretical alternative to ordinary coordinate-subspace supports. | Euler–Petit construct families and prove existence restrictions; their published results do not establish a speedup over previous ECDLP approaches [8]. | Literature candidate only; no concrete admissible ECC2K-130 support or measured advantage established. |

ECC2K-130 is over GF(2^131), and 131 is prime. A subfield GF(2^d) exists
precisely when d divides 131, so there is no intermediate GF(8), GF(16), etc.
Consequently, favorable experiments over small composite-degree extension
fields cannot simply be transferred to this instance. A vector subspace is
not a subfield. This does not rule out every Weil-descent construction; it
rules out treating a nonexistent intermediate field as an available one.

The historical ECC2K-130 project uses rho and publishes an approximately
2^60.9-iteration work estimate [7]. That is a historical algorithmic estimate,
not a measurement of this implementation, a current hardware speed, or proof
that an index-calculus candidate cannot improve. We make no conversion from
the toy timings to that work estimate.

## What runs end to end

The common pipeline checks the tiny group order and generator, generates a
planted target, constructs a factor base, gathers independently verified
relations, performs rank filtering and modular elimination, performs an
individual decomposition, recovers the scalar, and verifies [k]G=Q. Each
candidate receives only the target point and common support, never the
planted scalar or ground-truth decompositions.

Factor points are not assumed to lie in the prime subgroup. Entire relations
are projected by the cofactor, equal or opposite projected columns are
identified, and the scalar matrix is solved modulo the subgroup order.
This scalar field has odd characteristic even though the coordinate field
has characteristic two. The small dense elimination is a correctness
implementation, not a benchmark of large sparse linear algebra.

All adapters use three distinct, nonzero factor abscissae from the same
normal-basis support and exclude the target abscissa. This is a restricted
common domain; repeated summands and exceptional charts are not silently
counted as covered.

The two SAT adapters construct their formulas, invoke CryptoMiniSat with one
thread, independently evaluate returned Boolean models, lift candidate
abscissae to actual curve points, verify the signed group sum, and continue
past rejected algebraic roots. The truth oracle is used only in separate
validation and post-run auditing. In particular, it never supplies a
decomposition or a factor-base logarithm to the timed pipeline.

`pipeline.py` adapts the existing
`research/nagao_relations/coefficient_pullback/e2e.py`. `bindings.py` reuses
the existing solver implementations; its SAT adapter follows
`solver_06/run.py` while removing expected solutions from the computation.
The only accepted parameter pairs are (5,3) and (9,4), and the full pipeline
accepts only the exact two fixed group fixtures. There is no external target,
curve-parameter, or field-size command-line option.

## Accounting and predeclared boundary

For M rational factor abscissae there are at most 8*binom(M,3) signed distinct
triples. Uniform nonidentity targets therefore have decomposition probability
at most

`p_max = min(1, 8*binom(M,3)/(#E-1))`.

With B unknown columns and at most one relation per attempt, the expected
number of collection attempts is at least B/p_max. This is an attempt-count
boundary, not a bound in field multiplications or a universal lower bound on
all index-calculus algorithms. It is weak on these deliberately dense tiny
fixtures. Exact support coverage is independently enumerated and reported.

The reference is the existing **plain** three-partition Floyd rho code on
the same generated target and prime subgroup. It does not implement the
eligible negation/Frobenius quotient, so even beating this reference would
not establish an advantage over the strongest Koblitz rho implementation.

Cold per-fixture time includes setup, exhaustive tiny-curve enumeration,
target generation, factor-base construction, all decomposition attempts,
per-call encoding and preprocessing, lifting, verification, rank filtering,
matrix work, individual recovery, and final verification. Imported module
initialization is reported separately and excluded from per-fixture times.
There is no warm factor-base or solver-state amortization between runs.
Offline exhaustive validation is also reported separately.

Exclusive pipeline phase timers sum to the cold measurement, with a residual
orchestration bucket. Adapter subtimers sit inside the corresponding
collection/descent phase and must not be added to it. Field API additions,
multiplications and squarings are separate from scalar modular operations.
SAT encoding/search has no calibrated operation counter: its field vector
is explicitly **partial**, not a zero-cost search. Common-operation S,
calibrated cost/rho, and calibrated cost/floor remain null for every row.

This is not a performance optimization of the existing solver code, so the
different a=1 WDSat solver-stage regression is not substituted for the a=0
full-DLP experiment. Instead, the equivalent matched toy contract retains
the original support and includes exact historical counter replays. Passing
it does not satisfy the repository's calibrated end-to-end promotion gate.

## Reproduce

Use the existing repository at base commit
`2e2dcbe` with this experiment directory added, Python 3.12 and
`pycryptosat==5.14.7`. No Magma, Sage, GPU, cloud worker or network target
is required. Run sequentially:

```bash
python research/toy_ic_e2e_candidates_20260922/benchmark.py \
  --output research/toy_ic_e2e_candidates_20260922/new_results
```

Use `--smoke` for one seed and one repetition per fixture, or `--validate-only`
for the exhaustive five-bit cross-check. Output paths must be new; existing
evidence is never overwritten. The dependency can be installed in the user's
normal Python environment; this session reused the already installed package.

The main campaign uses the three frozen seeds, seven additional seeds and
three repetitions per variant. Fresh seeds are not necessarily distinct
targets in groups this small; distinctness is reported. Variant order is
deterministically shuffled. Wall-time intervals bootstrap per-seed medians,
not individual repetitions. Results describe these Python/SAT implementations
on this host; they do not rank optimized native implementations universally.

Correctness gates include exhaustive agreement on all 43 five-bit affine
targets for all five relation adapters, independent pair-vs-triple oracle
agreement, rejection of inconsistent relation rows, rejection of other
parameters, post-run auditing of every attempted decomposition including
complete empty answers, and independent verification of every recovered
scalar. Historical frozen scalar and operation records must reproduce.

## Remaining comparison gaps

A matched F4 backend and a matched WDSat backend remain unmeasured. General
large-prime filtering, Frobenius-compatible support compression, scalable
sparse scalar-field linear algebra, and a calibrated common operation unit
are not implemented here. Tiny-field correctness is useful evidence for the
pipeline but cannot establish behavior at the ECC2K-130 scale. No asymptotic
exponent is fitted to two toy field sizes.

## Primary sources

Checked 2026-09-22. Publication dates below describe the work, not search-index
crawl dates. This is a targeted review, not a claim of exhaustive coverage.

1. Trimoska, Ionica, Dequen, *A SAT-Based Approach for Index Calculus on Binary
   Elliptic Curves*, AFRICACRYPT 2020.
   [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7334981/),
   [author publication list](https://mtrimoska.com/).
2. Trimoska, Ionica, Dequen, *Parity (XOR) Reasoning for the Index Calculus
   Attack*, CP 2020. [Author manuscript](https://arxiv.org/abs/2001.11229).
3. Trimoska, *EC-Index-Calculus-Benchmarks*.
   [Source repository](https://github.com/mtrimoska/EC-Index-Calculus-Benchmarks).
   This generates polynomial instances; it is not itself a complete ECDLP timing.
4. Galbraith and Gaudry, *Recent progress on the elliptic curve discrete
   logarithm problem*, DCC 2016; preprint 2015.
   [Author manuscript record](https://eprint.iacr.org/2015/1022).
5. Huang, Kosters, Yang, Yeo, *On the last fall degree of zero-dimensional
   Weil descent systems*, 2015 manuscript.
   [Author manuscript](https://arxiv.org/abs/1505.02532).
6. Vitse, *Summation polynomials and symmetries for the ECDLP over extension
   fields*, DLP 2014, especially slide 10 on Riemann–Roch.
   [Author slides](https://www-fourier.univ-grenoble-alpes.fr/~viva/research/talks/Vitse_talk_DLP14.pdf).
7. *Breaking ECC2K-130*, original project and historical implementation data.
   [Project site](https://ecc-challenge.info/).
8. Euler and Petit, *New results on quasi-subfield polynomials*, accepted
   manuscript 2021. [Author manuscript](https://arxiv.org/abs/1909.11326).
