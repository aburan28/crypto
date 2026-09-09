**Factor-base mathematics for binary and Koblitz curves — 2026-09-08**

The most useful construction in this study is a Frobenius-invariant seed union, followed by optional translation closure under rational 2-torsion. The finite experiments below concern sumsets, representation multiplicities, and coefficient-matrix rank. They do not recover unknown scalars. These constructions and their elementary group-theoretic consequences are not claimed as novel.

For the ordinary binary model

\[
E: y^2+xy=x^3+a x^2+b,
\qquad
S_3(x_1,x_2,x_3)=(x_1x_2+x_1x_3+x_2x_3)^2+x_1x_2x_3+b,
\]

one descends the summation equations to Boolean equations after restricting the summand coordinates to the chosen factor base. A polynomial root still needs rational lifts whose group sum is the intended point. The existing symmetric S4 implementation is specialized to b=1. The sparse S3-chain representation handles general b but has a 64-variable representation limit; a wide S4 representation is now used for three-summand nonlinear Koblitz bases. This is a representation choice, not a complexity theorem. [Semaev's original paper](https://eprint.iacr.org/2004/031).

**Invariant linear spaces and invariant sets are different choices.** For odd n, normal-basis coordinates identify the Frobenius module with F2[T]/(T^n-1). Thus invariant linear subspaces correspond to divisors of T^n-1, and their possible dimensions are subset sums of its irreducible-factor degrees. For prime n, the nontrivial degrees are ord_n(2). At n=131 and n=163 those orders are 130 and 162: the only invariant linear dimensions are 0, 1, n-1, n. This conclusion does not forbid small invariant nonlinear sets. Frobenius-invariant factor bases, orbit compression, and their costs are established prior art. [Galbraith–Granger–Merz–Petit](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf).

For a small seed space W, define

\[
U(W)=\bigcup_{j=0}^{n-1}\{w^{2^j}:w\in W\},
\qquad
B(W)=\{P\in E(\mathbb F_{2^n})\setminus\{O\}:x(P)\in U(W)\}.
\]

The set is invariant under Frobenius and negation. Its abscissa count is at most n*2^dim(W), often smaller because of intersections and short orbits. Actual point counts must be computed: rational-lift density varies substantially between seeds. The union is generally nonlinear, so its ambient coordinate dimension must not be reported as the seed dimension. The implementation uses explicit finite-domain prefix constraints for SAT; those constraints cost up to O(n*|U|), and enumerating the base is charged separately. A future selector representation x=sum_i c_i beta_i^(2^j) would preserve the seed/translate parameterization, but no runtime improvement for that representation is established here.

**Reciprocal closure has a precise group interpretation.** On K_a, where b=1 and a is 0 or 1, T=(0,1) has order two. For P=(x,y) with x nonzero, the addition slope is (y+1)/x. Substitution of the curve equation into the binary addition formula gives

\[
x(P+T)=1/x(P).
\]

Consequently

\[
B^\#=(B\cup(B+T))\setminus\{O\}
\]

is obtained by adjoining reciprocals of nonzero abscissae. Frobenius fixes T, so this construction preserves Frobenius invariance. Negation also preserves it. If the cofactor h is even, then [h]T=O and

\[
[h]B^\#=[h]B
\]

for the bases used here, which contain T. Thus the projected signed Frobenius orbit count stays fixed, while B is contained in B# and every existing m-summand representation survives. Coverage can increase, although point count and domain-encoding cost can increase too. Removing O prevents treating an identity summand as an ordinary affine factor-base point.

The same calculation extends to an arbitrary ordinary binary curve with b nonzero: put c=sqrt(b), so T=(0,c), and obtain x(P+T)=c/x(P). Thus torsion closure is available independently of the Koblitz assumption. However, the coordinate-squaring map is an endomorphism of the same curve only when its coefficients are fixed by that Frobenius power. For a curve defined over F_(2^d), use the 2^d-power Frobenius and orbit length dividing n/d. A general binary curve must not be credited with the Koblitz orbit compression merely because its field has characteristic two.

When h=2, this closes every nonzero projected fiber completely. For a nonzero target R in the subgroup, a representation of [2]R over [2]B can be lifted to a representation of R over B#: a lift differs from R by either O or T, and at least one summand has nonzero projection, so its T-translate can correct the difference without becoming O. Hence B# attains the maximum decomposition support possible with that projected image. This is a support statement; it says nothing about the cost of finding the representations.

For larger h, 2-torsion closure need not resolve other cofactor components. More generally B+C, for a Frobenius-stable subgroup C of ker[h], preserves the projected image. Full kernel closure maximizes support for that image but may make the base much larger. Small torsion-subgroup closures are therefore a meaningful family to compare.

**The exact finite measurements.** The degree-19 curve is K_1 over F2[z]/(z^19+z^5+z^2+z+1), with group order 525086=2*262543. Candidates were selected using a 128-target census. Each row below was then evaluated on the same 4,096 distinct nonzero subgroup targets, excluding the original training targets. Targets are explicitly constructed test points, not unknown-scalar instances. Coverage was established by a pair-sum oracle and every selected witness was checked in the group. These are coverage measurements, not SAT timing measurements.

| Seed basis (integer polynomial masks) | Original / reciprocal point counts | Projected signed Frobenius orbits K | Original hits | Reciprocal hits | Reciprocal coverage |
|---|---:|---:|---:|---:|---:|
| 437455, 171788, 13631, 220500 | 191 / 381 | 5 | 3702 / 4096 | 4020 / 4096 | 98.1445% |
| 241511, 473913, 109139 | 229 / 457 | 6 | 3976 / 4096 | 4093 / 4096 | 99.9268% |
| 1, 2, 4, 8 | 229 / 457 | 6 | 3996 / 4096 | 4092 / 4096 | 99.9023% |
| 174944, 378963, 346918, 244674 | 343 / 685 | 9 | 4093 / 4096 | 4096 / 4096 | 100% observed |
| 77539, 461753, 7456 | 77 / 153 | 2 | 877 / 4096 | 910 / 4096 | 22.2168% |

The 381-point base is the best of these five saturated candidates by coverage/K: approximately 0.1963, compared with 0.1665 for the two 457-point bases. This is a screening objective for yield per retained orbit column, not an overall optimum: representation size and solving cost remain separate costs. All five coefficient matrices reached rank K over F_r after cofactor, Frobenius, and negation identifications. Rank was computed directly from known orbit coefficients; no linear system for unknown point scalars was solved. Different sampled equations are not claimed independent merely because they are Frobenius images.

The 229-point random seed originally hit 128/128 training targets but only 3976/4096 held-out targets. The original training result must therefore not be described as complete coverage. Likewise 4096/4096 for the largest saturated base is a sample result, not a proof of universal coverage.

At degree 15, the curve is K_1 over F2[z]/(z^15+z+1), with group order 32494=154*211. Here all 210 nonzero subgroup targets were tested. The divisor base indexed [0,2] has 61 points, two projected signed orbits, and covers 150/210 targets; reciprocal closure has 121 points, the same two projected signed orbits, and covers 210/210. In contrast, divisor [0,3] changes from 31 to 61 points but remains at 150/210 coverage and one projected signed orbit. A Frobenius union seeded by masks 3468,4413 already has 61 points, one projected signed orbit, and covers all 210 targets; its reciprocal closure adds points without improving coverage. Retain the smaller base in that comparison.

**Representation counts can be checked without a uniformity heuristic.** Let H be the cofactor component, and let c_u count the points of B in class u, using the homomorphism P -> [r]P. Define

\[
A=\sum_{u+v+w=0}c_uc_vc_w,\qquad
D=\sum_{2u+v=0}c_uc_v,\qquad
C=\sum_{3u=0}c_u.
\]

Burnside's lemma for the action of S3 on ordered triples gives the exact number of cofactor-admissible unordered triples, with repetitions allowed:

\[
M=(A+3D+2C)/6.
\]

If M_0 of them sum to O and the subgroup order is r, the exact mean number of unordered representations of a nonzero subgroup point is (M-M_0)/(r-1). It does not assume that individual targets receive equal numbers of representations. In all exhaustive degree-15 cases this independently calculated mean equals the observed mean. Histograms, variances, and pair additive energy are retained in the JSONL artifacts.

Writing N(R) for the representation count, Cauchy–Schwarz yields

\[
\frac{\mathbb E[N]^2}{\mathbb E[N^2]}\leq\Pr(N>0)\leq\min(1,\mathbb E[N]).
\]

This explains why maximizing the mean alone is insufficient: concentration can leave many targets uncovered. Pair additive energy is useful context but does not determine the third-sum distribution or its restriction to a cofactor class.

**Four-orbit refinement.** On the degree-19 curve, a generic projected signed Frobenius orbit contains 38 nonzero points. Under a deliberately heuristic Poisson model, K such orbits give approximately lambda=(38K)^3/(6r) unordered three-summand representations per target. Maximizing (1-exp(-lambda))/K leads to exp(lambda)=1+3lambda, whose nonzero solution is 1.9038136944. The predicted optimum is K=3.79497, suggesting a four-orbit base.

All five four-orbit subsets of the 381-point saturated base were then tested. Each retains the zero projected fiber and four of the five nonzero signed Frobenius orbits. The new common set of 4,096 distinct targets excludes both the 128 training targets and the entire earlier 4,096-target holdout. This is an exploratory comparison of all five subsets; the best of them is not claimed a global optimum.

| Omitted projected signed orbit (zero based) | Points | Retained columns | Covered targets | Coverage | Coefficient rank |
|---|---:|---:|---:|---:|---:|
| 0 | 305 | 4 | 3640 / 4096 | 88.8672% | 4 |
| 1 | 305 | 4 | 3541 / 4096 | 86.4502% | 4 |
| 2 | 305 | 4 | 3600 / 4096 | 87.8906% | 4 |
| 3 | 305 | 4 | 3650 / 4096 | 89.1113% | 4 |
| 4 | 305 | 4 | 3585 / 4096 | 87.5244% | 4 |

Orbit order is the lexicographic order of the canonical projected point in each signed Frobenius orbit; the canonical representative minimizes (x,y) over plus/minus Frobenius images. Omit orbit 3, the fourth representative in that order, from the saturation of the seed [437455,171788,13631,220500]. This specifies the selected 305-point subset completely in the stated field model. Its measured coverage/K is 0.2227783, compared with 0.1962891 for the 381-point parent on the previous holdout. The smaller base gives up some coverage while improving the screening objective and reducing point count by 76. Neither metric includes a SAT solving-cost model.

The resulting mathematical choices are the 381-point base when high coverage is the priority and the 305-point subset when coverage per projected column and base size matter more. Confirming the best subset on a further independent sample, and comparing torsion-stable subsets at other field degrees, remain future work.

**SAT and proof boundaries.** The Boolean work completed during this review supplies native XOR equations on the Koblitz path, exact finite-domain membership, and a linear trace constraint. For rational points, delta(P)=Tr(x(P)+a), delta(O)=0, is a homomorphism, giving sum_i Tr(x_i)=Tr(x_R)+(m+1)Tr(a). [Kosters–Yeo](https://arxiv.org/abs/1503.08001) discuss the associated low-degree relation and why low first-fall degree should not be equated with small solving degree. [Trimoska–Ionica–Dequen](https://arxiv.org/abs/2001.11229) study XOR reasoning and branching choices; their related [SAT index-calculus paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7334981/) supplies relevant modelling precedent. None of the sumset results above establishes a SAT complexity bound. The larger nonlinear bases have not been given a completed SAT runtime comparison in this mathematical phase.

A model-blocking defect was also fixed: clauses added after root propagation must discard root-false literals before selecting watches. Exhaustive Boolean model sets, including incremental blocking, agree with direct evaluation. The focused validation log contains 34 passing tests and one pre-existing ignored test. The mathematical artifact checks additionally verify paired target sets, nondecreasing support under reciprocal closure, unchanged projected signed orbit counts, histogram totals, and agreement of exact means in exhaustive cases.

Reproduce the mathematical experiment with:

```sh
cargo run --release --example koblitz_factor_base_yield -- 19 1 128
cargo run --release --example koblitz_factor_base_math -- research/sat_factor_base_review_20260908/yield-n19.jsonl 4096
cargo run --release --example koblitz_factor_base_math -- research/sat_factor_base_review_20260908/yield-n19.jsonl 4096 two-torsion
cargo run --release --example koblitz_factor_base_math -- research/sat_factor_base_review_20260908/yield-n19.jsonl 4096 two-torsion 3
```

Use a new output path for each run. `math-exact-n15.jsonl` and `math-reciprocal-n15.jsonl` contain exhaustive finite results; the degree-19 files contain paired samples. The final manifest records source and artifact hashes. Work occurred across incremental source revisions; final hashes identify the reproducible delivered implementation, not a retroactive binary attestation for every earlier exploratory run.
