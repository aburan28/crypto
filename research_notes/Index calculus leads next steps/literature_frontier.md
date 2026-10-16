# Index Calculus for ECDLP and Relatives: Literature Frontier (as of 2026-10-14)

Scope note: about 17 search and fetch calls. The Galbraith-Gaudry survey PDF (eprint 2015/1022 and the Birmingham mirror) returned HTTP 403, so its claims are taken from abstracts and from secondary quotes in other search results. Items marked **[unverified-recall]** come from background knowledge and were not confirmed against a fetched source in this session. Check them before relying on them. Everything else carries an inline source.

## Q1. Prime fields: Semaev summation polynomials and later work

### Takeaway
No subexponential (proven or heuristic) IC algorithm exists for ECDLP over prime fields F_p. Every prime-field approach replaces the binary/extension "subspace" factor base with a factor base cut out by a rational map, and none beats Pollard rho asymptotically or in practice. Since 2020, prime-field IC activity is close to nil apart from low-impact single-author ePrints (Hu 2024/1923, 2025/015).

### Cited Findings
- Semaev (2004) introduced summation polynomials and their possible use in IC on prime-field curves. Gaudry and Diem then fully developed the framework, combining it with Weil descent for extension fields. — [ePrint 2012/146 (FPPR) intro](https://eprint.iacr.org/2012/146.pdf); [Galbraith-Gaudry 2015/1022](https://eprint.iacr.org/2015/1022)
- **Petit-Kosters-Messeng, PKC 2016**: factor base F = {(x,y) ∈ E(K) : L(x) = 0}, where L is a large-degree rational map that is a composition of small-degree rational maps. Decomposition then becomes a system of low-degree polynomial equations. The authors say it is "only practical for small parameters at the moment" and that the asymptotic analysis is "limited by our understanding of Gröbner basis algorithms". — [PKC16 paper](https://christophe.petit.web.ulb.be/files/16PKC_primeECDLP.pdf); [slides](https://troll.iis.sinica.edu.tw/pkc16/slides/6-1--Algebraic_approaches_for_the_Elliptic_Curve_Discrete_Logarithm_Problem_over_Prime_Fields.pdf); [Birmingham record](https://research.birmingham.ac.uk/en/publications/algebraic-approaches-for-the-elliptic-curve-discrete-logarithm-pr)
- **Amadori-Pintore-Sala (ePrint 2017/609; FFA 51 (2018) 168-182)**: a modified Semaev approach for prime fields that needs only a single Gröbner basis computation. It claims to outperform both Semaev's original method and the PKM16 specialised algorithm. — [ePrint 2017/609](https://eprint.iacr.org/2017/609); [FFA DOI](https://unpaywall.org/10.1016/J.FFA.2018.01.009)
- **McGuire-Mueller (ePrint 2017/1262)**: an IC algorithm using fast summation polynomial evaluation with no Gröbner basis computation, plus a second algorithm that uses neither Gröbner bases nor summation polynomials. It gives complexity estimates and data, and is not peer-reviewed per the listing. — [ePrint 2017/1262](https://eprint.iacr.org/2017/1262)
- **Momonari-Kudo et al. (2018)** improved the PKM16 algorithm, according to Hu's secondary description. I could not retrieve the primary paper. — [ePrint 2024/1923](https://eprint.iacr.org/2024/1923)
- **Hu, ePrint 2024/1923**: reviews PKM16 and Momonari-Kudo, notes that neither PKM nor Amadori et al. gave concrete examples and that the algorithms could not be implemented efficiently, and proposes a "pseudo-index calculus". Its conclusion: no subexponential method is known for prime-field ECDLP. — [ePrint 2024/1923](https://eprint.iacr.org/2024/1923); [search summary](https://eprint.iacr.org/2024/1923.pdf)
- **Hu, ePrint 2025/015** (revised 2025-09): a generic "new IC method". The abstract says EC research has had "little success" and gives no concrete EC complexity. Low evidentiary weight. — [ePrint 2025/015](https://eprint.iacr.org/2025/015)
- **Huang-Petit-Shinohara-Takagi (IWSEC 2013; PJMI 2015)**: symmetrization variant of the Faugère et al. attack in characteristic 2 that reduces time and memory. The prompt lists it under prime fields, but it is binary-field work. — [IWSEC 2013 PDF](https://christophe.petit.web.ulb.be/files/IWSEC2013.pdf); [Birmingham record](https://research.birmingham.ac.uk/en/publications/improvement-of-faug%C3%A8re-et-als-method-to-solve-ecdlp)
- The Christophe Petit mailing-list post (2016) on future ECDLP advances gives author-side framing for PKM16. — [metzdowd 2016-03](https://www.metzdowd.com/pipermail/cryptography/2016-March/028742.html)

### Inferences
- Prime-field IC has three variants in the literature: (a) rational-map factor bases (PKM16, Momonari-Kudo), (b) single-Gröbner-basis reformulation (APS), and (c) Gröbner-free brute evaluation (McGuire-Mueller). All three are exponential in practice. Two levers are unexplored or under-explored: no published work optimises the choice of the composed rational map L (degree pattern, isogeny-induced maps), and none measures cost against rho in a matched unit.
- Since 2020, no strong-venue paper (CRYPTO, EUROCRYPT, ASIACRYPT, PKC, JoC) on prime-field ECDLP IC turned up in my searches. The area looks dormant, not refuted.

### Gaps
- I could not fetch the exact title, venue or contribution of Momonari-Kudo-Takagi (2018). [unverified-recall] It may be "Momonari, Kudo, Takagi, et al., 'Solving the ECDLP using summation polynomials over prime fields' / an improvement of PKM16 via a different decomposition of L", but I could not confirm this.
- Kosters-Yeo on prime fields: what I found (arXiv 1503.08001, "Notes on summation polynomials") concerns first-fall-degree behaviour (see Q7), not a prime-field algorithm.

## Q2. Extension fields F_{q^n} and binary fields F_{2^n}

### Takeaway
For F_{q^n} with n ≥ 2 fixed and q growing, IC (Gaudry/Diem) beats rho, with proven expected time Õ(q^{2-2/n}). For n → ∞ with n/log q → 0, Diem proved subexponential time. For prime-degree binary fields F_{2^n}, every subexponential claim depends on the first-fall-degree heuristic, which is doubted. The consensus (Galbraith-Gaudry) is that these curves are not threatened.

### Cited Findings
- **Diem, Math. Comp. 80 (2011)**: for fixed n ≥ 2, ECDLP over F_{q^n} in expected time Õ(q^{2-2/n}). This is a rigorous version of Gaudry's heuristic. — [AMS Math Comp 2011](https://www.ams.org/mcom/2011-80-273/S0025-5718-2010-02281-1/)
- **Diem, Compositio 147 (2011) 75-104**: for any sequences with n_i → ∞ and n_i/log q_i → 0, ECDLP over F_{q_i^{n_i}} is solvable in subexponential expected time (q_i^{n_i})^{o(1)}. A companion result about F_{q_i} is cut off in the snippet. — [Compositio](https://www.cambridge.org/core/journals/compositio-mathematica/article/on-the-discrete-logarithm-problem-in-elliptic-curves/59B877810708C90F6287972486A5BF0C); [PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/59B877810708C90F6287972486A5BF0C/S0010437X10005075a.pdf/on-the-discrete-logarithm-problem-in-elliptic-curves.pdf)
- **Diem, "On the DLP in elliptic curves II", ANT 7 (2013) 1281-1323** extends the parameter ranges. I did not verify the exact statements. — [Project Euclid](https://www.projecteuclid.org/journals/algebra-and-number-theory/volume-7/issue-6/On-the-discrete-logarithm-problem-in-elliptic-curves-II/10.2140/ant.2013.7.1281.full)
- Small-degree extension fields: ePrint 2010/157, "ECDLP over Small Degree Extension Fields" (Joux-Vitse line). — [ePrint 2010/157](https://eprint.iacr.org/2010/157.pdf)
- **Faugère-Perret-Petit-Renault (EUROCRYPT 2012, ePrint 2012/146)** studies polynomial systems from Weil descent. Heuristically it gives O(2^{c n^{2/3} log n}) bit operations for F_{2^n}, with c < 2, under first-fall-degree-type assumptions. — [ePrint 2012/146](https://eprint.iacr.org/2012/146.pdf); [ia.cr/2012/146](https://ia.cr/2012/146)
- **Semaev, arXiv 1504.01175 (2015)**: a new variant whose heuristic bound for binary fields is 2^{c√(n ln n)} with c ≈ 1.69. Like FPPR, it is heuristic (relies on degree of regularity / first fall degree). — [arXiv 1504.01175](https://arxiv.org/pdf/1504.01175)
- **Galbraith-Gaudry survey (2015/1022; J. Math. Ind. 2016)**: "elliptic curves over characteristic two fields F_{2^n} of prime degree n are not threatened by such methods and are still safe for use" (quoted via search summary). — [ePrint 2015/1022](https://eprint.iacr.org/2015/1022)
- **Galbraith-Granger-Merz-Petit, ePrint 2020/1315 (subfield/Koblitz curves)**: factor bases invariant under the q-power Frobenius reduce the number of polynomial systems by a factor of about 1/n ("the best one could hope for") and speed up linear algebra by n². The polynomial-solving savings are experimental. — [ePrint 2020/1315](https://eprint.iacr.org/2020/1315); [Birmingham record](https://research.birmingham.ac.uk/en/publications/on-index-calculus-algorithms-for-subfield-curves)
- **Huang-Kosters-Petit-Yeo-Yun, "Quasi-subfield polynomials and the ECDLP" (J. Math. Cryptol. 14 (2020) 25-38)** proposes factor bases from quasi-subfield polynomials (roots of x^{q^k} - λ(x)-type polynomials with special splitting). It gives complexity bounds that would be faster if such polynomials exist with suitable parameters. — [JMC 2020 PDF](https://birmingham.elsevierpure.com/files/97895891/Huang_et_al_Quasi_subfield_polynomials_Journal_of_Mathematical_Cryptology_2020.pdf); [DOAJ](https://doaj.org/article/75a4ced306114d568b832c4b20d00fd7); [arXiv 1909.11326](https://arxiv.org/pdf/1909.11326)
- **Euler-Petit, "New results on quasi-subfield polynomials" (FFA 75, 2021)** finds new families and a restriction theorem. Negative result: no speedup over earlier approaches. — [Birmingham record](https://research.birmingham.ac.uk/en/publications/new-results-on-quasi-subfield-polynomials/)
- Linearised trinomials that split completely (FFA 2018 / arXiv 1905.11755) bear on the subspace-style factor base question in binary fields. — [arXiv 1905.11755](https://arxiv.org/pdf/1905.11755)
- Trace-zero variety IC (arXiv 1405.1059) is a related relative. — [arXiv 1405.1059](https://arxiv.org/pdf/1405.1059)

### Inferences
- The asymptotic frontier for composite extension degrees is proven, and Diem's results remain the strongest rigorous ones. The open frontier is entirely prime-degree F_{2^n} and prime fields, where everything depends on unproven Gröbner complexity.
- Quasi-subfield polynomials are a clearly named factor-base family. Existence results so far (Euler-Petit 2021) are negative for speedups, but the existence question is not closed. This lead is concrete and computational.
- Frobenius-invariant factor bases (GGMP 2020) are only experimentally analysed. Extending the invariance idea to other endomorphisms (GLV/GLS, CM, small-degree isogenies) is a natural next step that the literature hints at (see Tsakou-Ionica, Q4).

### Gaps
- [unverified-recall] Gaudry (JSC 2009) heuristic Õ(q^{2-2/n}) for F_{q^n}. Joux-Vitse (EUROCRYPT 2012 / JoC) gave a variant for F_{p^5}/F_{p^6} with (n-1)-decomposition. Petit-Quisquater (ASIACRYPT 2012) gave 2^{ω τ n^{2/3} log n} under first-fall-degree assumptions. Galbraith-Gebregiyorgis (INDOCRYPT 2014) gave summation polynomials invariant under the action of 2-torsion and SAT solvers. Faugère-Gaudry-Huot-Renault (JoC 2014) used symmetries from the 2- and 4-torsion. I did not fetch the sources for any of these.
- GHS/Weil descent (Gaudry-Hess-Smart; Menezes-Teske weak F_{2^n} for composite n): not fetched in this session. [unverified-recall] Only composite-n binary fields are affected, and the attack is now classical.
- I could not verify the exact statements of Diem II (ANT 2013).

## Q3. Lifting and xedni

### Takeaway
Xedni calculus and characteristic-zero lifting are considered dead. Jacobson-Koblitz-Silverman-Stein-Teske (2000) showed the attack is "virtually certain to fail" asymptotically, because of bounded heights / an absolute bound on the coefficients of relations among lifted points. I found no new results revisiting the height obstruction.

### Cited Findings
- Jacobson-Koblitz-Silverman-Stein-Teske, "Analysis of the xedni calculus attack", DCC 20 (2000) 41-64: asymptotically the algorithm is virtually certain to fail, because of an absolute bound on the size of the coefficients of a relation satisfied by the lifted points. Even for small p, experiments show the odds against a suitable lift are prohibitively high. — [xedni.pdf](https://pages.cpsc.ucalgary.ca/~jacobs/PDF/xedni.pdf); [DOI](https://www.doi.org/10.1023/A:1008312401197); [Waterloo CORR 1999-05](https://uwaterloo.ca/combinatorics-and-optimization/sites/default/files/uploads/documents/1999_05.pdf)
- Mahalanobis-Mallick (arXiv 1703.07544, 2017): a Las Vegas ECDLP algorithm, non-generic but not xedni, which says xedni "was found to be no better than exhaustive search". A related ePrint (2018/134) and a 2020 follow-up (arXiv 2005.05039) introduce an "initial minors" conjecture for reducing to linear algebra. — [arXiv 1703.07544](https://arxiv.org/pdf/1703.07544); [ePrint 2018/134](https://eprint.iacr.org/2018/134.pdf); [arXiv 2005.05039](https://arxiv.org/pdf/2005.05039)
- A popular-press blog walkthrough of xedni exists (LeetArxiv Substack). It is non-scholarly. — [Substack](https://leetarxiv.substack.com/p/analysis-of-the-xedni-calculus-attack)

### Inferences
- Lifting is the most firmly closed direction. Any revival would have to defeat the height-growth obstruction, for example via p-adic or non-archimedean lifts rather than lifts to Q. I did not find this proposed in the literature.

### Gaps
- Mahdavi-Shafiee lifting work: not found in my searches.
- [unverified-recall] Silverman-Suzuki and Huang-Kueh-Tan on "lifting and ECDLP" (ANTS 2000) treat lifting obstructions for E(Q) of rank. Not verified this session.

## Q4. Alternative structures: isogenies, higher genus, subgroups and maps, transfers, SAT, symmetries, double large prime

### Takeaway
The most active "alternative factor base" lines are Frobenius- and endomorphism-invariant factor bases (2020-2021) and SAT/XOR-based decomposition (Trimoska-Ionica-Dequen 2020). The higher-genus IC relative (Tsakou-Ionica 2021) uses effective endomorphisms to shrink the factor base. I found no published work on summation ideals via isogenies or on factor bases from isogeny kernels.

### Cited Findings
- **Trimoska-Ionica-Dequen, AFRICACRYPT 2020, "A SAT-based approach for IC on binary elliptic curves"**: WDSat solver with symmetry breaking reduces the work for the (m+1)th summation polynomial by m!, and is up to 300x faster than MAGMA F4, with the factor growing in l and n. Worst case is still exponential in the factor-base dimension l. — [ePrint 2019/313](https://eprint.iacr.org/2019/313); [MATRICS page](https://www.matrics.u-picardie.fr/?p=3157)
- **Trimoska-Ionica-Dequen, "Parity (XOR) reasoning for the IC attack" (arXiv 2001.11229)**: WDSat with XOR reasoning on ANF. The XG-ext method captures a "canceling property" that algebraic methods use but SAT Gaussian elimination misses, though at high cost. Proposes a minimal-vertex-cover preprocessing as the branching rule. — [arXiv 2001.11229](https://arxiv.org/abs/2001.11229v3)
- **Tsakou-Ionica, ePrint 2021/721**: IC on genus-1 and genus-2 Jacobians over F_{q^n} with effective endomorphisms. Shrinking the factor base via endomorphisms gives a net gain despite the extra per-element cost. No explicit complexity formula is given in the abstract. — [ePrint 2021/721](https://eprint.iacr.org/2021/721)
- **Galbraith-Granger-Merz-Petit 2020/1315**: Frobenius-invariant factor bases with Frobenius symmetry breaking (see Q2). — [ePrint 2020/1315](https://eprint.iacr.org/2020/1315)
- **Huang-Petit-Shinohara-Takagi (IWSEC 2013)**: symmetrised variables to cut cost in FPPR systems. — [IWSEC2013](https://christophe.petit.web.ulb.be/files/IWSEC2013.pdf)

### Inferences
- Endomorphism-invariant factor bases have been done for Frobenius (subfield curves) and in genus 2 (Tsakou-Ionica). Doing the same for prime-field curves with efficient endomorphisms (GLV, CM j=0/1728), whose order-3/4/6 automorphisms act on the factor base, appears unexplored in the IC (not rho) setting.
- SAT decomposition has only been benchmarked in binary fields of prime degree. Applying it to prime-field rational-map systems (PKM16-style) looks open.

### Gaps
- [unverified-recall] Faugère-Svartz (ISSAC 2013) give Gröbner bases of ideals invariant under a commutative group (non-modular case), the basis for symmetric summation-polynomial solving. Gaudry-Thomé-Thériault-Diem (Math. Comp. 2007) give a double large prime variation for IC in small genus. The transfer families are MOV/Frey-Rück (pairing to F_{q^k}, small embedding degree) and Smart/Satoh-Araki/Semaev (anomalous curves, #E = p, via p-adic elliptic log). I did not fetch any of these.
- Splitting via genus > 1 covers for prime-field curves, and summation ideals via isogenies: no papers found. Treat as unexplored, not refuted.

## Q5. ML and LLM-guided work relevant to IC (2023-2026)

### Takeaway
ML for Gröbner bases exists (Kera et al., NeurIPS 2024; HATSolver, ePrint 2025/2219 with Lauter). It has been applied to generic or synthetic zero-dimensional systems, not to summation-polynomial or Weil-descent ECDLP systems. I found no paper applying learned orderings or transformers to ECDLP point decomposition.

### Cited Findings
- Kera et al., "Learning to compute Gröbner bases" (NeurIPS 2024): supervised transformers on system-to-Gröbner-basis pairs, using a fast dataset generator for zero-dimensional radical ideals. — [NeurIPS 2024](https://proceedings.neurips.cc/paper_files/paper/2024/hash/3a1de90699eec7d7f42c91d81f94af16-Abstract.html); [Nihon Univ record](https://pure.nihon-u.ac.jp/ja/publications/gr%C3%B6bner-basis-computation-via-learning/)
- Malhou-Perret-Lauter, HATSolver (arXiv 2512.14722; ePrint 2025/2219, Dec 2025): hierarchical attention transformers with curriculum learning solve larger instances than Kera et al. It is categorised under "Attacks and cryptanalysis" but reports no break of a deployed system. — [ePrint 2025/2219](https://eprint.iacr.org/2025/2219); [arXiv 2512.14722](https://arxiv.org/pdf/2512.14722)
- Kambe-Maeda-Vaccon (arXiv 2504.12465, 2025): proves the generality of the Kera et al. dataset generation and extends it. — [arXiv 2504.12465](https://arxiv.org/html/2504.12465v1)

### Inferences
- An open and concrete gap: train or evaluate these models on Semaev S_3/S_4 Weil-descent systems, on PKM16 rational-map systems, or on monomial-order and variable-order selection for F4/F5 on such systems. No literature found.

### Gaps
- I found no LLM-guided ECDLP cryptanalysis papers from 2023-2026. Peifer et al. (2020, RL for Buchberger selection) was mentioned in the HATSolver citations, but I did not fetch it.

## Q6. Recent (2023-2026) venue papers touching ECDLP IC and flagged open problems

### Takeaway
Searches turned up essentially no 2023-2026 CRYPTO, EUROCRYPT or ASIACRYPT papers on ECDLP index calculus. The visible recent items are low-impact ePrints (Hu 2024/1923, 2025/015) and ML-for-Gröbner work (2025/2219). The most recent substantive IC-for-ECDLP papers I found are GGMP 2020/1315, HKPYY 2020, Euler-Petit 2021 and Tsakou-Ionica 2021/721.

### Cited Findings
- Hu 2024/1923 concludes no subexponential method exists for prime-field ECDLP. — [ePrint 2024/1923](https://eprint.iacr.org/2024/1923)
- Hu 2025/015 is a generic IC method with no concrete EC result. — [ePrint 2025/015](https://eprint.iacr.org/2025/015)
- HATSolver 2025/2219 (see Q5). — [ePrint 2025/2219](https://eprint.iacr.org/2025/2219)
- Open problems named in the literature:
  - Existence of good quasi-subfield polynomials. — [HKPYY 2020](https://birmingham.elsevierpure.com/files/97895891/Huang_et_al_Quasi_subfield_polynomials_Journal_of_Mathematical_Cryptology_2020.pdf); [Euler-Petit 2021](https://research.birmingham.ac.uk/en/publications/new-results-on-quasi-subfield-polynomials/)
  - Asymptotic analysis of prime-field algebraic attacks, "limited by our understanding of Gröbner basis algorithms". — [PKM16](https://christophe.petit.web.ulb.be/files/16PKC_primeECDLP.pdf)
  - Full asymptotic analysis of Frobenius-invariant polynomial solving (experimental only). — [2020/1315](https://eprint.iacr.org/2020/1315)
  - SAT worst case exponential in l. — [2019/313](https://eprint.iacr.org/2019/313)

### Inferences
- The field has been quiet since roughly 2021. A careful, unit-matched negative or positive result on any of the named open problems would be novel.

### Gaps
- I could not run an exhaustive per-venue sweep of 2023-2026 proceedings (DBLP or IACR listings) within the tool budget. A 2026 paper I did not see may exist.

## Q7. Barriers: first fall degree controversy, Semaev's conjecture, Diem-style lower bounds

### Takeaway
The binary-field subexponential claims (FPPR 2012, Petit-Quisquater 2012, Semaev 2015) rest on the first fall degree assumption: that the degree of regularity is close to the first fall degree. Huang-Kosters-Yeo (CRYPTO 2015) gave computational evidence and a constructed Weil-descent system where it is unlikely to hold, and proposed the order-independent "last fall degree". Kosters-Yeo and Kousidis-Wiemers refined the bounds. The assumption is regarded as unproven and doubtful, and no proof or definitive refutation for summation systems has appeared since.

### Cited Findings
- Huang-Kosters-Yeo, "Last fall degree, HFE, and Weil descent attacks on ECDLP" (CRYPTO 2015, ePrint 2015/573): computational evidence casting doubt on the first fall degree assumption for characteristic-2 ECDLP, a constructed summation-polynomial Weil-descent system where it likely fails, and the last fall degree as a monomial-order-independent invariant. — [ePrint 2015/573](https://eprint.iacr.org/2015/573)
- Huang-Kosters-Yang-Yeo, "On the last fall degree of zero-dimensional Weil descent systems" (arXiv 1505.02532): a non-heuristic method between full Gröbner basis computation and the first fall degree heuristic. — [arXiv 1505.02532](https://arxiv.org/abs/1505.02532); [ar5iv](https://ar5iv.arxiv.org/html/1505.02532)
- Kosters-Yeo, "Notes on summation polynomials" (arXiv 1503.08001): studies the first fall degree of systems from the third summation polynomial and finds the m = 2 case "somewhat pathological". — [arXiv 1503.08001](https://arxiv.org/pdf/1503.08001); [Kosters ECC 2015 slides](https://lfant.math.u-bordeaux.fr/ecc2015/documents/kosters.pdf)
- Kousidis-Wiemers, "On generalized first fall degree assumptions" (ePrint 2015/358) improves the Petit-Quisquater bound. Kousidis (arXiv 1906.05594) expects the improved degree to be sharp except for the m = 2 case. — [ePrint 2015/358](https://eprint.iacr.org/2015/358); [arXiv 1906.05594](https://arxiv.org/pdf/1906.05594)
- Galbraith-Gaudry (2015/1022) review the complexity claims critically and conclude prime-degree binary curves remain safe. — [ePrint 2015/1022](https://eprint.iacr.org/2015/1022)
- Yokoyama-Yasuda-Takahashi-Kogure analysed Gröbner basis computation in Semaev's naive point decomposition, as cited in search summaries. — [JMC Oct 2020 (DOAJ)](https://doaj.org/article/bb64b7a98fd2462c84a20ffe0d32dd6b)
- Diem's results are upper bounds (Q2). I found no Diem lower-bound theorem for IC on prime fields. — [Compositio 2011](https://www.cambridge.org/core/journals/compositio-mathematica/article/on-the-discrete-logarithm-problem-in-elliptic-curves/59B877810708C90F6287972486A5BF0C)

### Inferences
- The two key barriers are (1) unproven solving-degree behaviour for Weil-descent and summation systems, and (2) for prime fields, no algebraic subspace structure in F_p, so factor bases need rational maps whose decomposition systems grow in degree. Measuring last fall degree or solving degree on small instances (as HKY did) is the accepted way to test claims empirically.

### Gaps
- Semaev's conjecture: [unverified-recall] Semaev (2004) conjectured that the summation polynomial S_m's degree / Gröbner-basis behaviour would allow subexponential decomposition. Its formal status is unproven and widely doubted. I found no explicit source this session stating its current status.
- A 2016-2026 proof or definitive refutation of the first fall degree assumption for summation systems: none found.
