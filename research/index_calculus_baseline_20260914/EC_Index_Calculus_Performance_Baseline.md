# EC index calculus: initial performance baseline and benchmark contract

Prepared 14 September 2026. Scope: classical ECDLP on ordinary elliptic curves, including binary Koblitz curves, other characteristic-two curves, and generic prime-field curves.

**Use a derived counting floor and a matched, measured reference, then map costs over field size, subgroup order, factor-base size, decomposition length, encoding and solver.** A single field-bit cutoff or a single complexity exponent would hide the main tradeoffs. This report establishes published reference points, a verified local WDSat pilot, an operation-accounting model, and numerical break-even budgets. It does not claim that a complete index-calculus attack has beaten rho, or that all solver/family combinations have now been measured.

The supplied benchmark repository is suitable as one input corpus. It needs a surrounding experiment runner, independent point verification, unbiased target generation, and full relation/matrix accounting before it can serve as the full performance standard.

## Repository integration and admission status

Class: **accounting**. This is an independent imported pilot and planning contract,
not a new solver or an improvement to an existing attack. The repository already
has a [Trimoska corpus audit](../../RESEARCH_TRIMOSKA_BENCHMARKS.md),
[F4 notes](../../RESEARCH_GROEBNER_F4.md),
[Koblitz scaling work](../../RESEARCH_KOBLITZ_SCALING_TARGET.md), and
[matched Nagao controls](../nagao_relations/solver_09/README.md). Statements below
about unmeasured work describe this pilot only. In particular, the existing audit
identifies `n19l6-19-U` as SAT; this six-case pilot does not include that instance.

The subsequent [full-corpus regression batch](regression/RESULTS.md) is the frozen
target for future solver iterations. It retains this initial pilot and adds both
symmetry configurations on all 60 inputs. It also certifies that `n19l6-19-U` is
algebraically SAT but its target does not lift to the curve; it cannot supply a
point relation on this curve over the stated base field.

The primary comparison follows `AGENTS.md`: **S = total common operations / sqrt(N)**,
with all phases included, plus ratios to a derived floor and a measured reference.
Here N is the subgroup order; n below is the binary field degree. Wall time and
solver conflicts are supplementary diagnostics. No binary-field conversion from
conflicts to common operations was measured, and no subgroup order or matched rho
run is in this corpus snapshot. Its full-attack S and ratios are therefore null.
The existing scoreboard verdict and fitted exponents do not change.

| Variant / required reference | Full-DLP S | Cost / matched rho | Cost / applicable floor | Correctness | Class |
|---|---|---|---|---|---|
| Unmodified WDSat, selected S4 corpus | Unmeasured | Unmeasured | Unmeasured | 3 SAT relations and 3 exhaustive UNSAT checks; no DLP | accounting |
| Matched automorphism-aware rho | Not run in this pilot | Unmeasured | Unmeasured | No scalar recovery in this pilot | reference required |

Frozen measurements and admission fields: [summary.json](summary.json).
The counting floor bounds independent uniform-target **attempt counts** in
first-relation mode. It becomes an operation floor only after an applicable
per-attempt operation lower bound is supplied. It is not a theorem ruling out
coordinate-based attacks.

## 1. Compare the right mathematical families

Koblitz curves over binary fields are a subset of characteristic-two curves. Use separate experimental strata so that Frobenius gains are visible.

| Stratum | Definition | Reference attack and index-calculus qualification |
|---|---|---|
| Binary Koblitz | Curve defined over F₂; points over F₂ⁿ | Rho with negation and eligible Frobenius orbits; orbit-aware factor bases deserve a separate IC baseline. |
| Other binary curves | Ordinary curves over F₂ⁿ without the same subfield structure | Rho with negation; Boolean Semaev systems, SAT/XOR and hybrid solvers are applicable. |
| Generic prime field | Ordinary curve over Fₚ with a large prime-order subgroup | Rho with eligible automorphisms; prime-field factor-base membership has to be encoded explicitly. |
| Small-extension positive control | E over F_qᵈ, fixed small d and increasing q | Some index-calculus algorithms have better asymptotic exponents than rho in this regime. Keep this stratum separate. |

Let **N** be the prime subgroup order. It is not necessarily the field cardinality. Record the full curve order and cofactor. Anomalous, pairing-vulnerable, supersingular and other special cases require separate strata; similarly, a composite binary extension degree can enable additional descent constructions.

For fixed small extension degree d, the familiar heuristic index-calculus bound is approximately Õ(q^(2−2/d)), with expensive dependence on d suppressed. In terms of N≈qᵈ, this is Õ(N^(2/d−2/d²)); for d=3 the exponent is 4/9. It cannot be transplanted to q=2 with growing d, or evaluated at d=1 to claim a prime-field speedup. [Trimoska–Ionica–Dequen, introduction](https://eprint.iacr.org/2019/313.pdf).

For the cryptographic-size binary and generic prime-field regimes, these sources do not establish a practical index-calculus replacement for optimized rho. Treat rho as the performance comparator, **not a proven lower bound on every coordinate-aware algorithm**. The uncertainty lies particularly in decomposition-system solving. [Galbraith–Gaudry survey](https://eprint.iacr.org/2015/1022), [subfield-curve analysis](https://eprint.iacr.org/2020/1315).

## 2. The complete cost model

Use distinct symbols for quantities that are often conflated:

| Symbol | Meaning |
|---|---|
| n | Binary extension degree; use log₂p for prime-field size. |
| N | Actual prime subgroup order. |
| m | Number of factor-base summands in a decomposition. |
| M | Number of permitted signed points in the factor-base support. |
| B | Number of unknown factor-base logarithms after sign/orbit compression. |
| ℓ | Dimension of an F₂ subspace defining allowed x-coordinates; M is often on the scale of 2^ℓ but must be counted. |
| p_dec | Probability a random target has an allowed decomposition. |
| u | Fraction of successful attempts producing a fresh usable row, for first-relation mode. |
| R_req | Number of rows actually required to obtain the necessary rank. |
| A | Number of attempted targets, including failures and timeouts. |

The subgroup-support convention matters: intersecting an x-coordinate factor base with the subgroup changes its size. Projecting every point by the cofactor changes its x-coordinates and can destroy the original membership rule. An alternative is to track full-group relations and consistently project the relation itself. State which convention is used.

For a sequential worker, use exclusive timers:

T_cold = T_setup + A · c_attempt + T_filter + T_matrix + T_individual.

c_attempt = E[T_target + T_encode + T_preprocess + T_solve + T_extract + T_lift + T_verify + T_canonicalize].

Point decomposition contains solving; Gröbner/F4 time is nested within solving. **Do not add “relation finding + decomposition + Gröbner” as three independent costs.** If sorting/filtering is charged per attempt, remove it from the global filtering timer.

For first verified relation mode, a useful planning approximation is:

A ≈ R_req / (p_dec · u).

For all-solutions mode, replace p_dec·u by the measured mean number of fresh usable rows per attempt. Multiple decompositions of the same target can provide useful rows, but duplicate rows and rank must be measured. Frobenius images must not be credited as automatically independent relations.

The primary stage metric is calibrated common operations per fresh verified row; the primary whole-method metric is S. Total core-seconds per fresh row and rank gain per core-second are secondary diagnostics. The final gate is successful scalar recovery: verify [k]P=Q. For fixed-base precomputation, report T_precompute and T_individual separately, and name the number K of targets before reporting T_precompute/K+T_individual.

### Counting boundary and success heuristic

For a uniform target in an order-N group and an allowed set of M signed points, the number of unordered m-tuples with repetition is binom(M+m−1,m). Therefore the following support bound is exact:

p_dec ≤ min(1, binom(M+m−1,m)/N).

It can be loose because many tuples have the same sum. When M≫m and tuple sums behave approximately uniformly, a useful **heuristic**, not a theorem, is:

λ ≈ M^m/(m!N), p_dec ≈ 1−exp(−λ).

For λ≪1, this reduces to p_dec≈λ. Estimate actual p_dec from uniform subgroup targets. Planted decompositions cannot estimate it: their distribution favors targets with many decompositions. The same applies to raw random x-coordinates, which need not lift to the curve.

### Typed operation counts

| Stage | Counts to collect | Time accounting |
|---|---|---|
| Factor-base setup | Field tests, root extractions, subgroup checks, stored points and orbit lengths | One-time cost plus memory; distinguish explicit enumeration from implicit membership. |
| Target generation | EC additions, doublings, scalar multiplications, Frobenius evaluations | Measure the actual sampler or walk, including startup and correlations. |
| Encoding | Field operations, created variables, equations, monomials and XOR/CNF clauses | Include Weil descent, basis conversion and all auxiliaries. |
| F4/Gröbner solving | Selected pairs, processed degrees, matrix rows/columns/nonzeros/ranks, coefficient updates | Split symbolic preprocessing, matrix construction and elimination. |
| SAT/XOR solving | Decisions, conflicts, propagations, XOR row operations, restarts | Conflicts are an algorithm-specific work count, not CPU instructions. |
| Extraction and checking | Polynomial factoring, y-lifts, sign combinations, EC additions, subgroup checks | Stop first-solution timer only after a valid point relation is verified. |
| Filtering and matrix | Duplicate rows, row weights, scalar-field operations, matrix-vector products and rank | Include large-prime recombination and fill-in if enabled. |
| Individual logarithm | Attempts, decomposition/verification and scalar recovery | Necessary for a complete attack after fixed-base precomputation. |

**There are two different linear-algebra domains.** Weil-descended binary polynomial systems use F₂. The relation matrix uses the subgroup scalar field F_N, generally a large odd prime field. A fast bit-XOR matrix kernel is not automatically a fast relation-matrix kernel. Coefficients ±1 permit modular additions; Frobenius coefficients or filtered rows can require general modular multiplication.

Record wall time, core-seconds and memory independently. Convert foreign operation units into the chosen common unit using measured conversion factors for that field and implementation. Measured kernel rates can separately translate the vector into a machine-specific time forecast. A rho iteration, field multiplication, SAT conflict, 64-bit XOR and floating-point operation are not interchangeable units.

## 3. Algorithms belong to different layers

Semaev and Nagao describe mathematical formulations. F4 and F5 compute Gröbner bases. WDSat and CryptoMiniSat search Boolean assignments. Compare an **encoding × solver** combination, rather than treating these names as interchangeable algorithms.

| Formulation or engine | Concrete algorithm | Complexity model / caveat |
|---|---|---|
| Direct S₃ | Solve S₃(x₁,x₂,x_R)=0 with both x-values in the allowed support | Two-point decomposition; an important small oracle/control case. |
| Chained S₃ | Replace S_(m+1) with m−1 S₃ equations linked by m−2 intermediate x-values | Binary raw dimension can be mℓ+(m−2)n. Lower input degree costs more variables. Exceptional partial sums need separate handling. |
| Symmetrized S₄ | Solve in elementary symmetric coordinates and enforce their relationship to x₁,x₂,x₃ | Three-point decomposition; recover points and signs after the algebraic solution. |
| Direct higher Semaev | Construct S_(m+1) by resultants and restrict support | Degree in each input is 2^(m−1); symbolic construction itself can dominate. |
| F4 | Batch S-polynomials, symbolically prepare a sparse Macaulay matrix, row-reduce, add new leading terms, repeat | Cost depends on the complete matrix sequence, sparsity and degree growth. |
| F5 / signature methods | Track signatures and syzygies to avoid redundant or zero reductions | Optional comparator; no automatic guarantee of beating a tuned F4 engine. |
| SAT/XOR | Convert to Boolean constraints; branch, propagate and optionally eliminate XOR equations | Measure GE on/off, branching and symmetry independently. |
| Guess-and-solve | Fix h variables, solve all required residual branches, combine verified results | Worst complete cost 2^h·C(v−h) for Boolean guesses; include unsuccessful branches and preprocessing. |
| Exhaustive / meet-in-the-middle | Enumerate group sums and match the target, with exact point verification | Required small-instance oracle; also a real competing decomposition method. |
| Nagao-style function formulation | Search functions whose zero divisor encodes the desired relation | Solve coefficients, support and normalization constraints; extraction and verification remain necessary. |

Magma documents its optimized F4 implementation. It is a useful reference engine, but its licensing is not a mathematical requirement: an explicitly identified open-source Gröbner engine is a valid separate baseline. Do not label a generic `groebner_basis()` call “F4” without checking the backend. [Magma documentation](https://magma.maths.usyd.edu.au/magma/handbook/text/1311).

For WDSat, the published symmetry method has a search bound on the scale of 2^(mℓ)/m! for fixed m, with per-node work requiring separate accounting. This is exponential and is not an end-to-end ECDLP bound. [SAT-based index calculus](https://eprint.iacr.org/2019/313.pdf).

An elementary group oracle enumerates m−1 summands and looks up the last, in O(M^(m−1)) work. Balanced meet-in-the-middle uses tables of partial sums: roughly O(M^ceil(m/2)) work and O(M^floor(m/2)) stored entries for one target. For m=3, a reusable pair table costs O(M²) setup and O(M) per subsequent target. Charge that table as precomputation; otherwise a solver comparison can be badly distorted. Repeated points and sum collisions require retaining the information needed by the selected enumeration mode.

### A concrete function-coefficient lane

The following is an experimental genus-one construction to benchmark alongside Nagao-style methods; it is not a claim of a new published complexity bound. Nagao's decomposition work provides the broader formulation context. [Nagao, decomposition formula](https://eprint.iacr.org/2013/548).

For E: y²+h(x)y=g(x) in characteristic two, write f=A(x)+B(x)y. Its norm under y↦y+h(x) is

Norm(f)=A²+ABh+B²g.

Search a nonzero f in a suitable Riemann–Roch space with pole divisor (m+1)O and f(−R)=0. When its zero divisor is (−R)+P₁+⋯+P_m, the elliptic-curve group law gives P₁+⋯+P_m=R. For the ordinary binary model, h=x and g=x³+ax²+b. In odd characteristic with y²=g(x), the analogous norm is A²−B²g.

After removing the known target factor from the norm, enforce that the residual monic polynomial H has its roots in the x-support V. For distinct residual x-values, this can be expressed as H dividing L_V, where L_V(X)=∏_(v∈V)(X−v). If V is an F₂-subspace, L_V is linearized. Compute X, X², X⁴, … modulo H successively and set L_V mod H=0, preserving the modular-squaring circuit instead of expanding everything.

Necessary controls: reject f=0; cover coefficient-normalization charts; enforce the intended pole/degree conditions; check both curve sheets; handle norm cancellations; and verify the final group relation. Because L_V is squarefree, H|L_V excludes repeated residual x-values. Either add a repeated-root lane or impose the same distinctness restriction on every comparator. This condition also excludes target collisions requiring multiplicities. A lower-degree coefficient system is not sufficient evidence of a faster algorithm: measure total coefficient search, constraint solving, extraction and usable yield.

Over Fₚ there is no proper subfield to supply the same F₂-subspace mechanism. An interval, bit-constrained set or other prime-field factor base needs its own membership encoding, and its degree/circuit cost must be charged.

## 4. Gröbner complexity: the measurable boundary

Record the observed maximum processed degree D_max, per-degree matrix sizes, coefficient operations and peak memory. Keep first fall degree, solving degree, degree of regularity and D_max as distinct quantities. Kosters and Yeo give reasons to doubt heuristics equating a small first fall degree with small regularity. [Notes on summation polynomials](https://arxiv.org/abs/1503.08001).

For v variables, a simple upper universe of monomials through degree D is

C_poly(v,D)=binom(v+D,D), C_bool(v,D)=Σ_(j=0)^min(v,D) binom(v,j).

A rough dense model sometimes writes the elimination contribution as C^ω. This is a scaling diagnostic, not an exact F4 operation count or a certified upper bound for every whole computation. Actual F4 uses selected columns and a sequence of generally rectangular matrices. Measure that sequence.

Here is a concrete size warning calculated directly from the Boolean monomial count:

| v | D | Monomial universe C_bool | Storage for a hypothetical dense square bit matrix |
|---:|---:|---:|---:|
| 30 | 4 | 31,931 | 0.119 GiB |
| 30 | 5 | 174,437 | 3.54 GiB |
| 30 | 6 | 768,212 | 68.7 GiB |

These are not measured Magma memory requirements. They demonstrate why a degree increase of one or two can matter more than a constant-factor kernel improvement. Conversely, adding auxiliary variables to lower input degree can enlarge the monomial universe enough to lose overall.

## 5. Published reference measurements

Times below are published results, not local reproductions or current-machine forecasts. The historical solver comparisons used a 2.40 GHz Intel Xeon E5-2640.

| S₄ parameters (n,ℓ), m=3 | Complete WDSat SAT / UNSAT, seconds | F4 SAT / UNSAT, seconds |
|---|---:|---:|
| (17,6) | 0.220 / 0.605 | 207.220 / 142.119 |
| (23,7) | 3.555 / 7.478 | 3128.844 / 2286.136 |
| (23,8) | 29.584 / 81.767 | exceeded 200 GB memory cap |

WDSat figures are Table 3; F4 figures are Table 2, with different averaging cohorts. They are historical reference points, not a newly paired speedup experiment. [Trimoska–Ionica–Dequen](https://eprint.iacr.org/2019/313.pdf).

Another published comparison for S₃, n=41 and ℓ=20, reports F4 at 16.8/18.7 seconds for SAT/UNSAT and WDSat with extended XOR reasoning and vertex-cover branching at 4.2/13.5 seconds. For a symmetrized S₄ model with 51 variables and 52 equations, WDSat without GE took 0.6/3.8 seconds; enabling GE took 19.0/49.8 seconds. **GE is a benchmarked choice, not a universal improvement.** [Parity reasoning paper, Tables 4–5](https://arxiv.org/pdf/2001.11229).

For Koblitz n=31, m=3, ℓ=5, the subfield study separates setup from Gröbner time:

| Factor-base subspace | Setup, seconds | Gröbner, seconds |
|---|---:|---:|
| Standard | 0.19 | 2.95 |
| Random | 0.19 | 3.15 |
| Frobenius invariant | 0.19 | 3.15 |

These are averages over 200 experiments, not complete DLP timings. Suitable Frobenius-invariant bases can reduce required relations by roughly n and sparse matrix work by roughly n² when long orbits dominate; that does not promise an n-fold faster individual solver call. [Subfield-curve paper, Section 5 and conclusion](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf).

A field size can grow while the allowed support stays tiny. The success probability then collapses even if a planted instance remains solvable. This is why the principal scaling axes must include ℓ, m and M^m/N alongside n.

## 6. Fresh local pilot and repository audit

Pinned corpus commit: `a053971cb155fbab30eaed3e2d11df850b4f3fbb`.
Pinned WDSat commit: `61c6ff3f49445af2729274819edb27b85a89efc1`.

The generator uses E: y²+xy=x³+x²+1, a binary Koblitz curve. It supplies polynomial/ANF/CNF encodings rather than a complete relation-collection and ECDLP pipeline. Its `U` construction samples random candidates, so the filename does not itself certify UNSAT. The Sage path ends with `Variety(I)`, making extraction/enumeration semantics something to inspect before comparing it to first-solution SAT timing. [Pinned repository](https://github.com/mtrimoska/EC-Index-Calculus-Benchmarks/tree/a053971cb155fbab30eaed3e2d11df850b4f3fbb), [generator](https://github.com/mtrimoska/EC-Index-Calculus-Benchmarks/blob/a053971cb155fbab30eaed3e2d11df850b4f3fbb/find_points.sage), [Sage implementation](https://github.com/mtrimoska/EC-Index-Calculus-Benchmarks/blob/a053971cb155fbab30eaed3e2d11df850b4f3fbb/Weil_descent.sage).

I compiled the public WDSat C source with GCC `-O3 -Wall`, without OpenMP, on a shared host reporting AMD EPYC 9V74. The run used the upstream ℓ=6 static configuration, symmetry breaking (`-b`), and GE disabled. The ℓ=5 instance is therefore not tuned to its exact allocation size. Empty lines were removed from the ANF input because the parser assumes every line after the header is an equation. No solver logic was changed. [WDSat configuration](https://github.com/mtrimoska/WDSat/tree/61c6ff3f49445af2729274819edb27b85a89efc1).

| Selected instance | Result | Median process wall time, 3 repetitions | Conflicts |
|---|---|---:|---:|
| n15l5-1-S | verified SAT | 0.00195 s | 22 |
| n15l5-11-U | exhaustively checked UNSAT | 0.03530 s | 5,593 |
| n17l6-1-S | verified SAT | 0.30308 s | 35,231 |
| n17l6-11-U | exhaustively checked UNSAT | 0.38809 s | 43,976 |
| n19l6-1-S | verified SAT | 0.01709 s | 2,793 |
| n19l6-11-U | exhaustively checked UNSAT | 0.41052 s | 44,050 |

Every SAT assignment was checked against all ANF equations, the S₄ expression and independently implemented binary curve addition after lifting the points. For the three UNSAT cases, enumeration covered all x-tuples with replacement up to permutation: 5,984, 45,760 and 45,760 cases, respectively, **97,504 total**. These checks certify the selected support-restricted S₄ instances. They do not prove correctness on arbitrary instances.

The timing includes process startup, input parsing and solver initialization/search. Independent verification was outside that timed region. The post-solve oracle-check time is recorded separately; per-run ANF checks and orchestration are not included in that timer. Input generation was not measured; this pilot consumes pre-generated inputs. Memory peaks and detailed F4-style primitive counts were not measured in this pilot; conflicts are the available solver work count. No GPU was used. These six selected cases do not estimate a random-target success rate or a scaling exponent. In particular, the easy 19-bit SAT case does not show that 19-bit fields are generally easier than 17-bit fields.

The [pilot results](pilot/results/pilot_results.json), compiler output, complete stdout/stderr, verification code and hashes are committed here. Upstream source and inputs are fetched by pinned URL and verified SHA-256 using the [reproduction commands](README.md#reproduce).

## 7. The numerical boundary to work against

For an approximately uniform rho walk in a quotient of average orbit size a, use the ideal birthday estimate

W_rho ≈ sqrt(πN/(2a)) iterations, T_rho=W_rho/rho_rate.

Thus a=1 gives 1.2533√N; negation a=2 gives 0.8862√N. For suitable binary Koblitz subgroups with long Frobenius orbits, a≈2n gives 0.8862√(N/n). Measure canonicalization, short cycles, distinguished-point overhead and the actual iteration rate; the ideal constants are not a performance guarantee. [Subfield-curve discussion of rho](https://eprint.iacr.org/2020/1315).

The maximum admissible mean attempt time for a first-relation strategy is

c_attempt,max = (T_rho − T_setup − T_filter − T_matrix − T_individual) · p_dec · u / R_req.

If the numerator is nonpositive, no decomposition speedup alone can win at those parameters. If p_dec is unmeasured, publish the result as a sensitivity calculation. For parallel runs, calculate worker throughput, total core-seconds and the nonparallel matrix/setup cost explicitly; do not divide the complete attack by worker count.

### Worked 60-bit subgroup scenario — calculated, not measured

Assume N≈2⁶⁰, m=3, M=2²⁰ signed support points, B=2¹⁹ columns after negation compression, and R_req=ceil(1.1B). The 10% row margin is an engineering assumption, not a rank guarantee. The scenario is a size model, not a specified curve.

| Quantity | Calculated value |
|---|---:|
| λ=M³/(6N) | 1/6 |
| Heuristic p_dec | 0.153518 |
| Required rows, under the chosen margin | 576,717 |
| Expected attempts, assuming every success is usable | 3,756,667 |
| Ideal negation-rho iterations | 951,578,915 |
| Rho time at 10⁶ iterations/s | 951.58 s = 15.86 min |
| Maximum mean attempt time if all other costs were zero | **253.30 µs** |
| Same budget in reference-rho time units | 253.30 iteration-equivalents |

At 10⁸ rho iterations/s, the attempt budget becomes 2.533 µs. These rates are hypothetical inputs, not measured kernel throughput, and an iteration-equivalent is a time normalization rather than a claim about the solver's internal operations.

Now account for the matrix. A deliberately simplified sequence of 2B matrix-vector products on a matrix with R_req rows and average row weight 3 gives approximately

2B·R_req·3 = 1.814×10¹² scalar coefficient accumulations.

That is an order-of-work model, not an exact Wiedemann count. At an assumed 10⁹ such accumulations/s it would already take about 1,814 seconds, exceeding the illustrative rho budget before relation finding. Actual sparse methods, coefficients, filtering, dimensions and throughput must be measured. This example shows why the matrix can force a different factor-base size or decomposition length.

Ignoring polynomial factors and constants, M≈N^(1/m) makes success approximately constant for fixed m. If B remains proportional to M, sparse relation-matrix work is on the scale N^(2/m): m=3 gives N^(2/3), m=4 gives N^(1/2), and m≥5 improves that exponent but makes decomposition harder. Large-prime strategies and alternative bases change this balance and need their own cost accounting. This is an algorithm-family tradeoff, not a universal lower bound.

For Koblitz curves, record an actual orbit histogram. If a stable support is preserved while factor-base columns shrink by n, the linear-algebra and relation-count gains can be substantial. The competing rho walk also improves by √n; it is incorrect to compare orbit-aware IC only against unoptimized rho.

On a selected subgroup where Frobenius acts as [λ_F], encode a point τ^j(F) by coefficient λ_F^j modulo N on the representative's column. This is why an orbit supplies fewer unknown logarithms, rather than automatically supplying independent rows. Short orbits, support stability and the Frobenius eigenvalue all belong in the instance record.

## 8. Benchmark execution and acceptance gates

The [machine-readable contract](ec_index_calculus_contract.json) freezes the following initial protocol:

1. **Reproduce and certify.** Begin with the pinned corpus and exact small-field oracles. Keep planted-SAT checks, uniform-target throughput, certified UNSAT, first-relation and all-solutions workloads separate.
2. **Compare matched encodings.** Run direct/chained S₃, symmetrized S₄ and function-coefficient methods against the same curve, support and targets. Compare tuned Magma F4, named open-source Gröbner engines, WDSat, CryptoMiniSat and guess-and-solve variants where the encoding is applicable. Keep a direct group-search/MITM competitor.
3. **Sweep three axes.** Vary field/subgroup size, m∈{2,3,4}, and factor-base dimensions near log₂N/m. Include smaller and larger support offsets. For each size, use three independent curves and initially 30 uniform targets per curve; promote to 100 for confirmation. Koblitz curves have fewer coefficient choices, so use independent field degrees/representations and bases while preserving the family label.
4. **Bound resources.** Start with a 60-second per-target cap and 16 GiB RAM. Confirmation uses 600 seconds and 64 GiB. Larger caps must be declared before the run. These are proposed execution limits, not measured solver capabilities. A full initial Cartesian sweep is unnecessary: expand only the competitive cells.
5. **Keep failures.** TIMEOUT, OOM and ERROR remain separate from UNSAT. Report solved fractions and capped total cost. Do not average only completed successes, and do not fit an uncensored runtime distribution to censored data.
6. **Promote on real work.** Require at least a 20% reduction in paired calibrated common operations per fresh verified relation at fixed support and m, with correctness, yield and memory tradeoffs reported. Require multiple curves, adjacent parameter cells, holdout targets and a 95% paired bootstrap interval excluding no improvement.
7. **Finish the attack on small instances.** Collect sufficient rank, solve the scalar system and verify the recovered logarithm. Compare full-DLP S and total-operation ratios to measured rho on the same subgroup; report cold wall time as a secondary practicality measure. Report precomputation amortization separately.

The important frontiers are: the largest support dimension solving 95% of the chosen workload within its cap; the minimum total cost per fresh relation; matrix/rank growth; and the ratio of total operations to matched rho operations. Runtime ratios are secondary. A method can improve one frontier while worsening another. Preserve the full tradeoff rather than compressing it into “largest field solved.”

## 9. What is established, and what remains unmeasured

| Item | Status in this baseline |
|---|---|
| Corpus scope and target-label audit | Inspected at a pinned commit. |
| WDSat selected-instance timings and conflicts | Measured locally; certificates checked. |
| F4/WDSat historical comparison | Published reference measurements, with parameter and cohort caveats. |
| Koblitz setup versus Gröbner split | Published reference measurements. |
| Success-counting bound, monomial counts and example budgets | Calculated; probability and runtime assumptions labeled. |
| Full stage percentages and primitive counts for tuned F4 | Not measured locally. |
| Matched Nagao/function-coefficient versus Semaev performance | Not measured by this imported pilot; see the existing repository Nagao controls above. |
| Generic-prime and non-Koblitz binary runtime curves | Not measured by this corpus/pilot. |
| Complete ECDLP timings or a rho crossover | Not established by this work. |

The initial practical target is to beat the strongest matched decomposition baseline **after failures, point extraction and verification are charged**, while reducing the predicted full-attack cost. The contract makes that claim testable and prevents a faster polynomial solver from being mistaken for a faster discrete-log algorithm.

To recompute the illustrative boundary, run `python3 ec_index_calculus_budget.py`. Supply the actual subgroup order, signed support size, column count, measured success probability, measured rho rate and all other-stage time before treating the output as a machine-specific budget. The adjacent JSON example is a model, not a benchmark result.
