# Ten open questions for elliptic curve index calculus

Evidence snapshot 2026-10-12; final review 2026-10-14.

The ranking prioritizes potential changes to attack complexity and concrete leads in the connected aburan28/cryptanalysis repository and its prime-field companion repositories. It is a research prioritization, not a consensus ranking. Established theoretical problems and repo-derived questions are labeled separately.

Prime-field curves are E(F_p) for prime p. Binary curves are ordinary curves over F_(2^n), with emphasis on prime n. Koblitz means binary subfield curves defined over F_2 and used over F_(2^n), a subset of binary curves. The prime-field j = 0 curve secp256k1 is included in the prime-field scope even though it is sometimes called a Koblitz curve. Throughout, r is the prime DLP subgroup order, n is an extension degree, m is the summand count, B is the actual usable factor-base point count, and C is the effective number of relation-matrix columns.

The reviewed evidence establishes no practical classical index-calculus improvement on the cryptographic prime-field or prime-degree binary/Koblitz instances considered here. Generic square-root lower bounds constrain black-box algorithms; they do not prove that coordinate-based index calculus is impossible. The repo's October 1, 2026 status note explicitly disclaims an ECDLP speedup from its volcano experiments.

Sources: [R1: cryptanalysis status and scope](https://github.com/aburan28/cryptanalysis/blob/main/README.md); [P2: Bernstein and Lange, Safe curves for elliptic-curve cryptography](https://eprint.iacr.org/2024/1265).

Frobenius acceleration itself is already established: Galbraith et al. obtain, under their conditions, about n fewer collected relations and an n^2 saving in final relation algebra. The unresolved issue is whether usable decompositions and the rest of the pipeline make those reductions competitive with Frobenius-accelerated rho. The shifted-base m! mechanism is also prior work, and the linearized two-point splitting idea has prior art.

Source: [P6: Galbraith, Granger, Merz and Petit, On Index Calculus Algorithms for Subfield Curves](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf).

![Research questions mapped to the IC recovery pipeline](research_pipeline.svg)

1. **Can prime-field curves admit a sparse factor base with genuinely cheap decomposition?**

   Prime fields. Broad theoretical open problem.

   For a nonspecial E(F_p) subgroup of order r, can a coordinate-defined factor base support enough efficiently found m-point relations to make complete classical index calculus asymptotically faster than square-root search?

   **Evidence.** The prime-field repo's interval, multiplicative-subgroup, lattice and p-adic experiments expose obstacles in those particular constructions. Its Hecke follow-up also identifies designated relations that additive-closure or Fourier tests would miss. Published prime-field algebraic methods remain limited by decomposition and polynomial solving.

   **What remains open.** Find a target-independent decomposition algorithm, rather than assume a free oracle. A null closure test is not a proof against every non-generic construction. A base of size r^(1/3) and fast relation generation alone do not establish a speedup: final relation solving and target recovery must also fit the budget.

   **What would settle it.** Give a parameterized family, an explicit membership and decomposition algorithm, a natural-query success estimate, and costs for relation rank, matrix solving and target recovery. For an experimental lead, freeze unseen targets and compare the measured exponent against tuned rho.

   Sources: [R2: Prime-field smallness experiments](https://github.com/aburan28/ecdlp-prime-field-smallness/blob/main/README.md); [R3: Hecke relations and the class-number barrier](https://github.com/aburan28/ecdlp-prime-field-smallness/blob/main/LEADS.md); [R16: Representation track and metering](https://github.com/aburan28/ecdlp-cost-challenge/blob/main/REPR_TRACK.md); [P10: Petit, Kosters and Messeng, Algebraic Approaches for the Elliptic Curve Discrete Logarithm Problem over Prime Fields](https://research.birmingham.ac.uk/en/publications/algebraic-approaches-for-the-elliptic-curve-discrete-logarithm-pr/); [P2: Bernstein and Lange, Safe curves for elliptic-curve cryptography](https://eprint.iacr.org/2024/1265).

2. **What is the actual solving degree of binary summation systems as the parameters grow?**

   Binary fields and Koblitz curves. Theoretical open problem with repo evidence.

   Can one prove, or reliably characterize, how the solving degree grows with extension degree n, base dimension l and summand count m in the regimes that could improve ECDLP complexity?

   **Evidence.** Caminata and Gorla distinguish solving degree from first fall degree and related invariants. The repo measures a strong two-summand degree predictor, but labels its degree-131 table as a structural extrapolation and has no corresponding large-parameter m >= 3 result.

   **What remains open.** A small first fall degree, bounded input degree, or a fit at fixed small m does not prove efficient solving when n and m grow together. Satisfiable and unsatisfiable ordinary queries can behave differently.

   **What would settle it.** Establish bounds for a specified formulation and monomial order, including field equations; then test held-out larger systems and retain censored queries. An m = 5 or 6 result needs its own evidence.

   Sources: [R4: Point decomposition degree and factor-base heuristics](https://github.com/aburan28/cryptanalysis/blob/main/experiments/pdp-degree-heuristics/README.md); [P1: Caminata and Gorla, Solving degree, last fall degree, and related invariants](https://arxiv.org/abs/2112.05579).

3. **Can Frobenius orbit decomposition retain a compact algebraic witness?**

   Koblitz curves. Repo-derived research direction.

   Can a low-degree or sparse correspondence represent addition among signed Frobenius orbits while keeping only a small certificate of relative phases?

   **Evidence.** Galbraith et al. already establish Frobenius-based reductions in relation count and final matrix size. The repo's phase-eliminated quotient experiment found that all eight n = 19 target predicates needed degree 9 on its valid-key domain, exactly the degree where every Boolean function on that domain became representable.

   **What remains open.** That negative result applies to the tested coordinates and fully eliminated predicate. Retained relative phases, nonlinear presentations and auxiliary-variable correspondences remain candidates. Orbit folding reduces the number of log unknowns; it does not automatically remove the witness search.

   **What would settle it.** Construct an exact quotient correspondence, independently lift its witnesses to curve points, and measure ordinary-query cost and independent-row yield on the same frozen base.

   Sources: [R5: Orbit-valued quotient summation](https://github.com/aburan28/cryptanalysis/blob/main/experiments/orbit-quotient-summation/README.md); [R6: Frobenius quotient decomposition budgets](https://github.com/aburan28/cryptanalysis/blob/main/experiments/frobenius-quotient-m4/README.md); [P6: Galbraith, Granger, Merz and Petit, On Index Calculus Algorithms for Subfield Curves](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf).

4. **Which Frobenius-stable nonlinear factor bases keep both yield and decomposition cost favorable?**

   Koblitz curves. Repo-derived research direction informed by a recent paper.

   Can nonlinear stable bases combine many usable subgroup points, few orbit columns and cheap membership constraints without moving the cost into auxiliary variables or phase search?

   **Evidence.** The repo combines normal-basis Hamming constraints with the factorized-convolution formulation of La Scala et al. (September 2026) in a complete degree-9 toy IC pipeline. The paper itself concerns syndrome decoding. The repo records no corresponding scaling result at degrees 53, 83 or 131 for that FC-Hamming method.

   **What remains open.** Frobenius invariance and a bounded-degree Hamming encoding are useful ingredients, but their interaction with Semaev equations, actual subgroup density and ordinary-query solver difficulty is unsettled. Hamming weight is not the same as a linear factor-base subspace.

   **What would settle it.** Compare exact point sets and encodings on identical inputs, and measure setup, auxiliary-variable count, failed-query cost, useful rank and one-target recovery. Retain the existing weight encoding as a matched control.

   Sources: [R7: FC-Hamming in a complete toy IC pipeline](https://github.com/aburan28/cryptanalysis/blob/main/experiments/hamming-ic-e2e-20260929/README.md); [R6: Frobenius quotient decomposition budgets](https://github.com/aburan28/cryptanalysis/blob/main/experiments/frobenius-quotient-m4/README.md); [P3: La Scala, Marchesin and Tiwari, Hamming Ideals and Grobner Bases for ISD-like Syndrome Decoding](https://arxiv.org/abs/2609.18866v2); [P6: Galbraith, Granger, Merz and Petit, On Index Calculus Algorithms for Subfield Curves](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf).

5. **Can explicit quasi-subfield constructions become useful at prime extension degrees?**

   Binary fields and Koblitz curves. Published construction problem.

   Are there explicit quasi-subfield polynomial families at the desired prime n whose rational roots yield a usable factor base and a complete attack better than the current alternatives?

   **Evidence.** Huang et al. propose quasi-subfield factor bases. Euler and Petit provide further families and an existence limitation, while explicitly reporting no resulting speedup over prior ECDLP approaches. The repo's explicit-base frontier illustrates why actual usable point counts are essential.

   **What remains open.** The desired regime requires enough rational roots, subgroup-valid points and affordable decomposition. A root-set size or tuple-count upper bound does not establish natural relation yield. Extension-field constructions do not directly give prime-field smallness.

   **What would settle it.** Supply a parameter family meeting all construction hypotheses, prove or measure its usable point count and natural success rate, and account for building the base and solving its relations.

   Sources: [R8: Exact subgroup factor-base feasibility screen](https://github.com/aburan28/cryptanalysis/blob/main/experiments/pdp-scaling/EXPLICIT_BASE_FRONTIER.md); [R4: Point decomposition degree and factor-base heuristics](https://github.com/aburan28/cryptanalysis/blob/main/experiments/pdp-degree-heuristics/README.md); [P8: Euler and Petit, New results on quasi-subfield polynomials](https://arxiv.org/abs/1909.11326v2); [P9: Huang, Kosters, Petit, Yeo and Yun, Quasi-subfield Polynomials and the Elliptic Curve Discrete Logarithm Problem](https://doi.org/10.1515/jmc-2015-0049).

6. **Can multi-summand formulations avoid paying back their degree reduction during lifting?**

   Binary fields and Koblitz curves. Algorithmic open problem.

   Can symmetric coordinates, chained S3 equations or hybrid solvers find verified ordinary m-point decompositions faster than existing methods after splitting and membership checks are charged?

   **Evidence.** The repo shows that symmetric two-summand equations can be linear while the remaining affine solution space is expensive to enumerate. Three-summand behavior depends strongly on the factor-base product spaces. Its linearized-half method improves a decomposition stage without producing a rho-beating complete attack.

   **What remains open.** Solving for elementary symmetric values does not guarantee that the associated polynomial splits into base coordinates with a valid group relation. Planted systems test correctness and satisfiable cost; they do not price ordinary relation collection.

   **What would settle it.** Use independently labeled ordinary and planted corpora at increasing m and l. Count all splitting, sign, rationality, subgroup and failed-query work, and compare direct, chained and symmetric formulations on the same base.

   Sources: [R4: Point decomposition degree and factor-base heuristics](https://github.com/aburan28/cryptanalysis/blob/main/experiments/pdp-degree-heuristics/README.md); [R10: Linearized-half decomposition and its full cost](https://github.com/aburan28/cryptanalysis/blob/main/experiments/linearized-half-decomposition/README.md); [R9: Frobenius-shifted factor-base geometry](https://github.com/aburan28/cryptanalysis/blob/main/experiments/shifted-base-geometry/README.md); [P7: Trimoska, Ionica and Dequen, A SAT-Based Approach for Index Calculus on Binary Elliptic Curves](https://eprint.iacr.org/2019/313); [P1: Caminata and Gorla, Solving degree, last fall degree, and related invariants](https://arxiv.org/abs/2112.05579).

7. **How can a collector buy independent rank and reliable target descent efficiently?**

   Prime fields, binary fields and Koblitz curves. Algorithmic research question.

   Which witness selection, partial-relation recombination and stopping policies minimize cost per new rank contribution and per independently recovered target log?

   **Evidence.** The repo's monitor separates relation hits, attainable rank and descended targets. Its shifted-base audit records increased geometric coverage without a demonstrated faster collector. The published subfield framework also requires a usable relation matrix before recovery.

   **What remains open.** Many hits can be duplicates or dependent, and quotient or cofactor structure can create a known nullspace. More coverage alone need not accelerate target recovery. Final relation algebra is modulo r, even when the PDP solver works over F_2.

   **What would settle it.** Freeze the query stream and rank checkpoints; compare collectors by rank gained per unit work, include partial-relation graph costs, and verify fresh-target descent. State how unavoidable kernels and anchors are handled.

   Sources: [R4: Point decomposition degree and factor-base heuristics](https://github.com/aburan28/cryptanalysis/blob/main/experiments/pdp-degree-heuristics/README.md); [R9: Frobenius-shifted factor-base geometry](https://github.com/aburan28/cryptanalysis/blob/main/experiments/shifted-base-geometry/README.md); [R11: Verified toy IC benchmark and target descent](https://github.com/aburan28/cryptanalysis/blob/main/experiments/ic-bench/README.md); [P6: Galbraith, Granger, Merz and Petit, On Index Calculus Algorithms for Subfield Curves](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf).

8. **Can isogeny and Hecke mechanisms overcome their rank and transport barriers?**

   Prime fields; binary and Koblitz transport variants. Repo conjecture and research direction.

   Can a concrete transport or correspondence produce enough cheap, independent relations to improve the DLP, including any links needed between modular isogeny classes?

   **Evidence.** The prime-field LEADS file verifies toy Hecke relations but reports substantial rank deficiency in its class-walk model. Those ranks use a proxy modulus, so they are not a theorem about every DLP subgroup. The binary volcano study reports nearly identical attempt counts for its source and descendants.

   **What remains open.** For the Hecke mechanism, persistence of the deficit at deeper smooth conductors and a useful class-connecting relation are conjectural. For curve transport, an explicit isogeny preserves the DLP when its degree is coprime to r; any gain must come from a cheaper target representation after route costs.

   **What would settle it.** Recheck rank modulo the actual subgroup prime, prove the appropriate kernel dimension, and charge explicit map and target-transport costs. Compare complete same-target runs rather than infer a benefit from a changed conductor or j-invariant.

   Sources: [R3: Hecke relations and the class-number barrier](https://github.com/aburan28/ecdlp-prime-field-smallness/blob/main/LEADS.md); [R12: Volcano IC over GF(2^19)](https://github.com/aburan28/cryptanalysis/blob/main/experiments/volcano-ic/README.md); [R13: Volcano IC committed numerical summary](https://github.com/aburan28/cryptanalysis/blob/main/experiments/volcano-ic/results/summary.json); [R18: Verified isogeny route registry](https://github.com/aburan28/cryptanalysis/blob/main/experiments/ic-candidate-catalog/isogeny_routes.json); [P2: Bernstein and Lange, Safe curves for elliptic-curve cryptography](https://eprint.iacr.org/2024/1265).

9. **Can a higher-genus cover make the elliptic component cheaper to solve?**

   Prime fields and binary fields, including Koblitz curves. Repo-derived research direction.

   Can a cover or Jacobian construction preserve the desired elliptic subgroup and yield cheaper relation collection and target recovery for that component?

   **Evidence.** The repo has a tiny prime-field quartic norm-projection diagnostic and binary cover transfer infrastructure. It explicitly distinguishes norm relations on E from relations in the full Jacobian, where a Prym component can remain. Diem and Thome's genus-three index calculus supplies a relevant higher-genus benchmark.

   **What remains open.** A map certificate is not an attack advantage. The genus-three heuristic O-tilde(q) bound is faster than rho for a group of size about q^3, but is larger than square-root search on a same-field elliptic subgroup of size about q.

   **What would settle it.** Prove pullback injectivity on the chosen subgroup, provide exact norm or pushforward operations and valid target lifts, identify the required E-component rank, and price the full pipeline against elliptic rho in the original field.

   Sources: [R14: Bielliptic quartic norm-projection diagnostic](https://github.com/aburan28/cryptanalysis/blob/main/experiments/bielliptic-quartic/README.md); [R15: Hyperelliptic cover transfer infrastructure](https://github.com/aburan28/cryptanalysis/blob/main/experiments/hyperelliptic-cover-infrastructure/README.md); [P11: Diem and Thome, Index calculus in class groups of non-hyperelliptic curves of genus three](https://members.loria.fr/EThome/files/non-he-genus3.pdf).

10. **Where is the verified complete-attack crossover with tuned Pollard rho?**

   Prime fields, binary fields and Koblitz curves. Empirical and complexity research question.

   For a precisely specified family, can any current algebraic, SAT, hybrid or quantum-annealing IC method recover an unseen target faster than tuned rho under comparable resources, and does the advantage persist as size increases?

   **Evidence.** The repo's benchmark distinguishes one-target online recovery from reusable setup and cold cost. The 2024 Edwards annealing paper reports a tiny F_1021 example; Lin et al.'s September 2026 paper reports a 16-bit prime-field demonstration. These sizes do not establish an asymptotic or cryptographic-size advantage.

   **What remains open.** The target-dependent PDP, final relation solving when required, descent, verification, memory and failure distribution must be priced coherently. Koblitz rho must receive its signed-Frobenius benefit, and prime-field j = 0 or 1728 curves need their applicable automorphism controls.

   **What would settle it.** Run repeated same-point comparisons over a frozen size ladder and held-out tier. Report online one-target time as primary, reusable preprocessing and cold totals separately, and multi-target workloads with their own target count. For annealing, include embedding, classical preprocessing, readout and repetitions; establish scaling rather than extrapolate a small demonstration.

   Sources: [R11: Verified toy IC benchmark and target descent](https://github.com/aburan28/cryptanalysis/blob/main/experiments/ic-bench/README.md); [R16: Representation track and metering](https://github.com/aburan28/ecdlp-cost-challenge/blob/main/REPR_TRACK.md); [P4: Lin et al., Solving the elliptic curve discrete logarithm problem using quantum tunneling effects on the D-Wave Advantage2 quantum computer](https://www.sciopen.com/article/10.23919/ICN.2026.0017); [P5: Dzierzkowski and Wronski, Solving Elliptic Curve Discrete Logarithm Problem on Twisted Edwards Curves Using Quantum Annealing and Index Calculus Method](https://www.mat.savba.sk/ojs/index.php/TATRA/article/view/945); [P2: Bernstein and Lange, Safe curves for elliptic-curve cryptography](https://eprint.iacr.org/2024/1265).

**A measured negative result that informs question 8.**

![E0 to descendant PDP attempt ratio and reported uncertainty](attempt_ratio.svg)

Committed n = 19 study: E0/descendant mean PDP-attempt ratio = 0.9964, reported 95% interval [0.9765, 1.0166], from 31 curves x 10 scalars with 310/310 verified logs. A ratio of 1 means equal attempts; below 1 favors E0. The interval includes 1. This is an attempt-count result, not a wall-time or IC-versus-rho speedup. The interval is reproduced from the repository summary, not recomputed from raw trials.

Sources: [R12: Volcano IC over GF(2^19)](https://github.com/aburan28/cryptanalysis/blob/main/experiments/volcano-ic/README.md); [R13: Volcano IC committed numerical summary](https://github.com/aburan28/cryptanalysis/blob/main/experiments/volcano-ic/results/summary.json).

**Concrete next work.** For the next bounded repo experiment, prioritize questions 2, 3, 6 and 7: a frozen ordinary-query panel with an exact base, a solver-independent witness check and rank accounting can distinguish a useful algorithmic lead from a faster planted-system solve.

**Evidence scope.** Repository results are reported as committed measurements, controls or conjectures. Published papers were checked at their primary pages or author manuscripts. The recent sources include the 2023 solving-degree paper, the 2024 Safe Curves survey and annealing paper, and September 2026 Hamming-ideal and annealing papers. Older 2020–2021 papers define construction and Frobenius questions that remain relevant. This synthesis adds no new measured attack or proof.

**Curve and route context.** The principal full-width Koblitz reference is EC1N131Ckb1h136f03e58c98: y^2 + xy = x^3 + 1 over F_2[z]/(z^131 + z^13 + z^2 + z + 1), cofactor 4, subgroup order 680564733841876926932320129493409985129. The registered degree-263 descendant is EC1N131Cbinh833014327b07 and the verified route is IW1E263d1hadee4e69fa3d. Neither a route certificate nor these research proposals is a measured full-width IC candidate. The quantitative plot uses the legacy n = 19 volcano study, not a new EC1 or run ID.

Sources: [R17: Exact curve registry](https://github.com/aburan28/cryptanalysis/blob/main/experiments/ic-candidate-catalog/curves.yaml); [R18: Verified isogeny route registry](https://github.com/aburan28/cryptanalysis/blob/main/experiments/ic-candidate-catalog/isogeny_routes.json).

**Graph audit.** Checked the volcano figure inventory (census.svg, ecdlp.svg and tau.png), its canonical analyze.py source and results/summary.json; the factor-base degree/yield figures referenced by pdp-degree-heuristics; and the curve and isogeny route registries. This is a synthesis with no new curve, route or corrected measurement. Existing graph facts remain current. The new pipeline diagram and attempt-ratio plot explain the reviewed evidence and do not alter canonical graph sources.

**Paper dates and scope.**

- [Caminata and Gorla, Solving degree, last fall degree, and related invariants](https://arxiv.org/abs/2112.05579) (2023). published paper; preprint revised June 2022, journal volume 114 January–February 2023.

- [Bernstein and Lange, Safe curves for elliptic-curve cryptography](https://eprint.iacr.org/2024/1265) (2024). survey.

- [La Scala, Marchesin and Tiwari, Hamming Ideals and Grobner Bases for ISD-like Syndrome Decoding](https://arxiv.org/abs/2609.18866v2) (2026). preprint; version 2 dated September 29, 2026; syndrome decoding, not an ECDLP result.

- [Lin et al., Solving the elliptic curve discrete logarithm problem using quantum tunneling effects on the D-Wave Advantage2 quantum computer](https://www.sciopen.com/article/10.23919/ICN.2026.0017) (2026). published September 21, 2026; publisher abstract reports a 16-bit prime-field demonstration.

- [Dzierzkowski and Wronski, Solving Elliptic Curve Discrete Logarithm Problem on Twisted Edwards Curves Using Quantum Annealing and Index Calculus Method](https://www.mat.savba.sk/ojs/index.php/TATRA/article/view/945) (2024). published in volume 87 (2024); tiny prime-field demonstration.

- [Galbraith, Granger, Merz and Petit, On Index Calculus Algorithms for Subfield Curves](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf) (2021). SAC 2020; proceedings published 2021; linked conference version.

- [Trimoska, Ionica and Dequen, A SAT-Based Approach for Index Calculus on Binary Elliptic Curves](https://eprint.iacr.org/2019/313) (2020). AFRICACRYPT 2020; ePrint originally 2019.

- [Euler and Petit, New results on quasi-subfield polynomials](https://arxiv.org/abs/1909.11326v2) (2021). Finite Fields and Their Applications 75, 101881; explicit families and limitations.

- [Huang, Kosters, Petit, Yeo and Yun, Quasi-subfield Polynomials and the Elliptic Curve Discrete Logarithm Problem](https://doi.org/10.1515/jmc-2015-0049) (2020). published paper.

- [Petit, Kosters and Messeng, Algebraic Approaches for the Elliptic Curve Discrete Logarithm Problem over Prime Fields](https://research.birmingham.ac.uk/en/publications/algebraic-approaches-for-the-elliptic-curve-discrete-logarithm-pr/) (2016). PKC 2016; author institution abstract.

- [Diem and Thome, Index calculus in class groups of non-hyperelliptic curves of genus three](https://members.loria.fr/EThome/files/non-he-genus3.pdf) (2008). Journal of Cryptology 2008; linked author preprint March 2007.

Editable figures: [build_report.py](build_report.py). Canonical report content and source index: [questions.json](questions.json). Quantitative plot input: [attempt_ratio_data.json](attempt_ratio_data.json).
