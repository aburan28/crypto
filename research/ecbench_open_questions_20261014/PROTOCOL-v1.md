# Ten open questions: ecbench campaign, protocol 1

Frozen before execution on 2026-10-14. Source baseline:
`f9011e5375155bd9ef7839abe4eea12718e4b401` of `aburan28/crypto`.
The separate `aburan28/ecbench` repository is empty; this repository contains
the native ecbench harness. The original ten-question literature review is
preserved in `literature/questions.json` and `literature/REPORT.md`.

The requested scope remains all ten open research questions. These experiments
test their executable prerequisites and existing candidates. A completed
diagnostic does not resolve the corresponding mathematical open question.
No dimension, timeout, failed run, or required confidence gate is silently removed.

## Frozen measurement contract

Every measured algorithm runs through native ecbench. The generation, structural
controls, verification and analysis use Rust. Public points are hashed to the
declared subgroup by ecbench; algorithms receive no planted logarithm. One
previously unseen point is each primary workload; four independent workloads
per curve, three paired rounds for prime/isogeny panels, two for binary/orbit/
nonlinear/rank/multi-summand panels, and two workloads/one round for bounded
solver calls. Zero warmup retains cold setup in every execution. Seeds are
20261014 except frozen curve-constructor seed 59297. Inputs are exactly the
JSON files emitted by `examples/ecbench_open_questions.rs freeze`.

Reference: negation rho for prime and ordinary binary curves; signed Frobenius
strong rho for Koblitz curves. A/A repeats use identical target and algorithm
seeds. Setup, oracle construction, relation generation/checks, linear algebra,
recovery and external verification are retained. Reports include cold counted
group-addition equivalents, S = cost/sqrt(r), setup and search phases, native
resource counters, success and exhaustion. Unpriced native operations remain
unpriced. A quotient of incomplete counts is only a counted-work diagnostic,
never a bound on the ratio of full costs. The harness's `floor_s` is a
birthday-style reference constant, not a universal exact Shoup lower bound.

Local runs use `--cpus none --allow-busy`: L0 operation-count diagnostics.
All specs retain required L2; wall results are inadmissible until earned L2+
isolation and the A/A gate pass. No building or testing during a session.
Independent-host replay and GPU measurements require a configured runner.
The same specs can be exported by `ecbench isolab-job` without loosening fidelity.

Each child has the frozen 180-second watchdog; deterministic trial/node/state
caps are in the specs. Timeouts, errors, budget exhaustion, unsupported bases
and invalid results are retained, not replaced or assigned a zero cost.
Audit every deterministic measured execution with `verify --replay-all`.
Record the binary digest, source manifest, Rust version, host capsule, plans,
commands, sealed sessions, audits, tables and resource vectors. No exponent
fit or subexponential claim is authorized by these small fixed-dimension panels.

## Question-specific hypotheses and acceptance criteria

| Q | Frozen experiment | Hypothesis tested | Diagnostic acceptance; unresolved original obligation |
|---|---|---|---|
| 1 | Prime search 16/18/20/22 bits; abscissa 32/64, direct m2/m3 and negation MITM | Sparse-base choice shifts actual complete-pipeline work, not just relation density | Verified recovered scalars, paired counts, all unpriced work exposed; a scalable prime-field base with cheap membership/decomposition and cheaper full cost remains required |
| 2 | Ordinary binary degrees 13/17/19/23; Koblitz degrees 13/19/23/31, summand dimension3, inherited F4/F6 node cap4096/degree4, exact rootfinding control | Solver outcomes need not be predicted by first-fall degree | Preserve infeasible, solved and censored calls; separately record m2 Macaulay ranks, exhaustive roots at dimensions3/4, and actual m3 chain variable/equation counts; asymptotic solving degree as dimensions grow remains open |
| 3 | E0 at degrees19/23/31/37, same compact base4, m3 direct/Frobenius MITM, cap256 | Orbit folding changes cost without changing the point-decomposition problem | Same FB identity and independent witness verification; deriving a compact polynomial witness for relative phases remains open |
| 4 | E1 at degrees11/23, gamma weights1/2/3, cap8 orbits, compact control | Nonlinear stable bases can have different subgroup support and yield | Report raw weight enumeration, projected points, columns and solver outcomes; gamma-weight support is not a normal-basis Hamming ideal, and scalable nonlinear algebraic membership remains required |
| 5 | Exact existing quasi-subfield census at (n,n0)=(5,2),(7,3),(7,4),(11,3), retain3 witnesses | Admissible linearized polynomials sometimes split even at prime degree | Exhaustive coefficient counts; independently enumerate roots of retained witnesses and actual curve point support; a useful general family and cheaper IC pipeline remain open |
| 6 | E1 degree11, exact3/4/6/n−1 summands, subspace3/envelope4, state cap100000 | Higher arity changes coverage and exact-search burden | `n−1` means10, never3; preserve state-budget refusal; symmetric/chained/hybrid asymptotic advantage including lifting remains open |
| 7 | Same compact streams at degrees19/23, target-pin/full-rank/checked-rank/shared-rank; degree11 LP0/1/2 | Raw relations and independent useful rank differ; target-blind setup must be separate | Verified scalars, rank/useful rows/rejections and online phase partition; reusable rank law, descent reliability and general large-prime optimization remain open |
| 8 | CM conductor gap ell263/t241 crater/floor; paired same-model public targets | IC cost can differ among same-order models without proving cross-model transport | Replay model/subgroup/discriminant checks; do not pair points across curves without an isogeny certificate; general Hecke rank and paid map transport remain open |
| 9 | Prime F101 genus2 cover certificate; actual binary genus3 pullbacks at degree7 and registered degree83 E0/E1 | A checked cover can preserve the elliptic subgroup without yielding a cheaper solver | Verify equations, smoothness, actual pullback classes, scalar homomorphism samples, image annihilation and infinity/ramified/off-curve controls; independent exact subgroup-order proof, general norm and faster Jacobian descent remain separate gates |
| 10 | Aggregate every panel, required N83 IC applicability check, exact tiny Boolean QUBO control | Pipeline and encoding costs invalidate inference from isolated solver success | Replayed resource vectors with no lower-bound ratio marketed as speed; N83/N131 full IC, isolated A/A, held-out growing-size crossover, hardware annealing and GPU match remain unresolved |

Degree31 has ord_31(2)=5, unlike ord_53(2)=52, ord_83(2)=82 and
ord_131(2)=130. It is exploratory, never a substitute for the mandatory
degree83 E0 confidence gate (modulus z83+z45+z2+z+1, r=
2417851639230796216685689, cofactor4). The harness currently refuses wide IC;
preserve the refusal and leave that gate blocked. Degree11 multi-summand
results do not transfer to prime degree83 without a new implementation/proof.

The QUBO control uses the repository's Boolean descended Semaev equations,
Rosenberg AND penalties and squared parity residuals. Exhaustively check every
original assignment and minimize the parity slack conditional on a faithful
monomial assignment. Also exhaustively check the AND truth table. This checks
an encoding prerequisite, not an annealer, quantum hardware or a full DLP.

## Stop and decision rules

Stop each individual execution at its frozen deterministic cap or watchdog.
Complete the planned panels and replays, recording infrastructure failures
separately. A diagnostic passes only with its stated correctness controls;
cost comparisons use paired successful executions and explicitly report
censoring. A result that fails the original proof, full-cost, dimension-growth,
or N83 gate leaves the corresponding open question partial/blocked. A failed
finite search does not prove nonexistence. Reports, PDFs and source-linked
graphs go beside the sealed sessions. Ready PR and review/CI gates remain
required; Conductor access is currently unavailable and separately recorded.
