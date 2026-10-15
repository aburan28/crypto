# All ten questions: ecbench execution and remaining obligations

Executed nine ecbench sessions with **937 attempts**, **867 verified scalar recoveries**, **69 exhausted searches** and **1 applicability error**. Nine session audits passed; all 867 recoveries replayed exactly, and the 70 failure outcomes/counts passed a separate native replay. Every original mathematical question remains open.

The native harness is aburan28/crypto; aburan28/ecbench is empty. [Frozen protocol](PROTOCOL-v1.md), [structural amendment](CONTROLS-AMENDMENT-v2.md), [extended-control/replay amendment](CONTROLS-AMENDMENT-v3.md), [source manifest](SOURCE-MANIFEST.json), [commands](execution-manifest.json), [summary](summary.json), [complete outcomes](results.json). Measured binary cc5bba39d33465ddc1b1d31af8ac25c6680d41193fc76012d9a1019b6046f02d, compiled from clean commit c716e462934128c75e99dd458748e08623a2708c. Host ECBENV2h687984d4e73d, virtualized AMD EPYC 9V45, earned L0. Public target law hash_to_subgroup_v1; seed20261014.

All comparisons use one public target per workload with paired rho and A/A controls. Setup is charged in every cold execution. Primary online windows and their exclusive phases are retained in the raw records, but their L0 clocks are inadmissible for a performance claim. Native work remains unpriced. S is counted group-addition equivalents / sqrt(r); its quotient is a counted-work diagnostic, **not a bound on the ratio of complete costs**. The harness floor is a birthday-style reference constant. Successful-only mean S excludes failures; paired ratios are null whenever any execution in that curve/arm is censored. All failures retain their recorded work and resource vector. No scaling exponent was fitted.

## Results by original question

1. **Can prime-field curves admit a sparse factor base with genuinely cheap decomposition?** — partial.

   Evidence: q01-prime. Remaining: A scalable prime-field base with executable cheap membership/decomposition and a fully priced held-out crossover.

   Literature/repo basis: [R2](https://github.com/aburan28/ecdlp-prime-field-smallness/blob/main/README.md); [R3](https://github.com/aburan28/ecdlp-prime-field-smallness/blob/main/LEADS.md); [R16](https://github.com/aburan28/ecdlp-cost-challenge/blob/main/REPR_TRACK.md); [P10](https://research.birmingham.ac.uk/en/publications/algebraic-approaches-for-the-elliptic-curve-discrete-logarithm-pr/); [P2](https://eprint.iacr.org/2024/1265)

2. **What is the actual solving degree of binary summation systems as the parameters grow?** — partial.

   Evidence: q02-binary, q02-solvers. Remaining: Growing-dimension solving-degree bounds; correct chained-system shape metadata and held-out satisfiable/infeasible queries.

   Literature/repo basis: [R4](https://github.com/aburan28/cryptanalysis/blob/main/experiments/pdp-degree-heuristics/README.md); [P1](https://arxiv.org/abs/2112.05579)

3. **Can Frobenius orbit decomposition retain a compact algebraic witness?** — partial.

   Evidence: q03-frobenius. Remaining: Compact algebraic witnesses for relative Frobenius phases, including membership, splitting and lift costs.

   Literature/repo basis: [R5](https://github.com/aburan28/cryptanalysis/blob/main/experiments/orbit-quotient-summation/README.md); [R6](https://github.com/aburan28/cryptanalysis/blob/main/experiments/frobenius-quotient-m4/README.md); [P6](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf)

4. **Which Frobenius-stable nonlinear factor bases keep both yield and decomposition cost favorable?** — partial.

   Evidence: q04-nonlinear. Remaining: Scalable nonlinear algebraic membership/solving; gamma weight is not a normal-basis Hamming ideal.

   Literature/repo basis: [R7](https://github.com/aburan28/cryptanalysis/blob/main/experiments/hamming-ic-e2e-20260929/README.md); [R6](https://github.com/aburan28/cryptanalysis/blob/main/experiments/frobenius-quotient-m4/README.md); [P3](https://arxiv.org/abs/2609.18866v2); [P6](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf)

5. **Can explicit quasi-subfield constructions become useful at prime extension degrees?** — partial.

   Evidence: controls-v3.json / quasi_subfield. Remaining: A useful prime-degree parameter family, subgroup yield law and complete measured IC integration.

   Literature/repo basis: [R8](https://github.com/aburan28/cryptanalysis/blob/main/experiments/pdp-scaling/EXPLICIT_BASE_FRONTIER.md); [R4](https://github.com/aburan28/cryptanalysis/blob/main/experiments/pdp-degree-heuristics/README.md); [P8](https://arxiv.org/abs/1909.11326v2); [P9](https://doi.org/10.1515/jmc-2015-0049)

6. **Can multi-summand formulations avoid paying back their degree reduction during lifting?** — partial.

   Evidence: q06-many-summand. Remaining: Growing arity with full splitting/lifting/matrix costs; literal ten-summand toy results do not transfer to N83.

   Literature/repo basis: [R4](https://github.com/aburan28/cryptanalysis/blob/main/experiments/pdp-degree-heuristics/README.md); [R10](https://github.com/aburan28/cryptanalysis/blob/main/experiments/linearized-half-decomposition/README.md); [R9](https://github.com/aburan28/cryptanalysis/blob/main/experiments/shifted-base-geometry/README.md); [P7](https://eprint.iacr.org/2019/313); [P1](https://arxiv.org/abs/2112.05579)

7. **How can a collector buy independent rank and reliable target descent efficiently?** — partial.

   Evidence: q07-rank, q06-many-summand. Remaining: Rank and descent reliability laws, optimized large-prime policies and fully priced reusable setup/online windows.

   Literature/repo basis: [R4](https://github.com/aburan28/cryptanalysis/blob/main/experiments/pdp-degree-heuristics/README.md); [R9](https://github.com/aburan28/cryptanalysis/blob/main/experiments/shifted-base-geometry/README.md); [R11](https://github.com/aburan28/cryptanalysis/blob/main/experiments/ic-bench/README.md); [P6](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf)

8. **Can isogeny and Hecke mechanisms overcome their rank and transport barriers?** — partial.

   Evidence: q08-isogeny. Remaining: Explicit paid isogeny transport and a Hecke-rank theorem or certified counterexample in its stated model.

   Literature/repo basis: [R3](https://github.com/aburan28/ecdlp-prime-field-smallness/blob/main/LEADS.md); [R12](https://github.com/aburan28/cryptanalysis/blob/main/experiments/volcano-ic/README.md); [R13](https://github.com/aburan28/cryptanalysis/blob/main/experiments/volcano-ic/results/summary.json); [R18](https://github.com/aburan28/cryptanalysis/blob/main/experiments/ic-candidate-catalog/isogeny_routes.json); [P2](https://eprint.iacr.org/2024/1265)

9. **Can a higher-genus cover make the elliptic component cheaper to solve?** — partial.

   Evidence: controls-v3.json / covers; prime-cover-replay.json. Remaining: Independent exact source-order proofs, a general Jacobian norm and a cheaper higher-genus relation/descent solver.

   Literature/repo basis: [R14](https://github.com/aburan28/cryptanalysis/blob/main/experiments/bielliptic-quartic/README.md); [R15](https://github.com/aburan28/cryptanalysis/blob/main/experiments/hyperelliptic-cover-infrastructure/README.md); [P11](https://members.loria.fr/EThome/files/non-he-genus3.pdf)

10. **Where is the verified complete-attack crossover with tuned Pollard rho?** — partial.

   Evidence: q10-n83-gate. Remaining: Wide IC at N83, N131 transfer justification, isolated A/A, independent host, complete pricing, held-out scaling and matched GPU/annealing hardware.

   Literature/repo basis: [R11](https://github.com/aburan28/cryptanalysis/blob/main/experiments/ic-bench/README.md); [R16](https://github.com/aburan28/ecdlp-cost-challenge/blob/main/REPR_TRACK.md); [P4](https://www.sciopen.com/article/10.23919/ICN.2026.0017); [P5](https://www.mat.savba.sk/ojs/index.php/TATRA/article/view/945); [P2](https://eprint.iacr.org/2024/1265)

## Structural observations

The complete quasi-subfield coefficient searches examined32,16384,2097152 and2048 vectors at (n,n0)=(5,2),(7,3),(7,4),(11,3), finding0,127,127 and0 splitting polynomials. Retained witnesses passed independent root enumeration; their actual point/subgroup support is recorded. This finite census is not a general family or complete IC attack.

All six gamma-weight base constructions passed on-curve, subgroup-annihilation, Frobenius-closure and coefficient checks. The cap is eight projected orbits; raw weight domains and projection counts are retained. Hamming ideals for syndrome decoding do not automatically apply to this base.

The degree controls measure m2 Macaulay matrices through degree5, at dimensions3/4, and exhaust their Boolean assignments. They also record actual m3 chain sizes. Ecbench PDP3 shape currently reports only summand coordinates; intermediate field variables and the second equation block must be included before interpreting solver shape/growth. First fall degree is not solving degree. Pipeline trial-cap exhaustion does not by itself imply an individual algebraic query hit its node budget.

The degree7 and exact degree83 E0/E1 genus3 pullbacks passed generator/image annihilation, sampled scalar homomorphism, typed push-pull, infinity, ramified and off-curve controls. An additional degree7 public-point Jacobian BSGS control recovered the scalar and replayed it on the elliptic curve. This is functional generic recovery, with no measured speed claim. The F101 genus2 cover certificate replay passed; two-infinity Jacobian arithmetic/pullback/norm remain unsupported. E1 degree83 uses catalogued r=8569786107849059 and cofactor1128547018; it does not replace the E0 N83 subgroup. The first E1 cofactor2 constructor failure is preserved.

The QUBO encodes the repository's descended Semaev system. All8192 original assignments passed faithful monomial/slack completion checks; a small independent fixture exhausts every auxiliary assignment. This does not enumerate every auxiliary assignment of the full encoding, prove an x-only zero-energy solution lifts to a point relation, or run hardware annealing.

## Session ledger and counted costs

| Session | Attempts | Verified | Exact verified replays | Audit |
|---|---:|---:|---:|---|
| [q01-prime](sessions/q01-prime/session.json) | 288 | 288 | 288 | true |
| [q02-binary](sessions/q02-binary/session.json) | 128 | 109 | 109 | true |
| [q02-solvers](sessions/q02-solvers/session.json) | 40 | 16 | 16 | true |
| [q03-frobenius](sessions/q03-frobenius/session.json) | 128 | 118 | 118 | true |
| [q04-nonlinear](sessions/q04-nonlinear/session.json) | 96 | 89 | 89 | true |
| [q06-many-summand](sessions/q06-many-summand/session.json) | 64 | 55 | 55 | true |
| [q07-rank](sessions/q07-rank/session.json) | 96 | 96 | 96 | true |
| [q08-isogeny](sessions/q08-isogeny/session.json) | 96 | 96 | 96 | true |
| [q10-n83-gate](sessions/q10-n83-gate/session.json) | 1 | 0 | 0 | true |

| Curve | Arm | Verified / attempts | Mean cold S | Counted IC / rho | Status |
|---|---|---:|---:|---:|---|
| icv1-fp16-t295-8d3c3165 | abs32-subtract2 | 12 / 12 | 146.3892 | 75.5481 | {"verified":12} |
| icv1-fp16-t295-8d3c3165 | abs64-mitm2 | 12 / 12 | 22.7599 | 11.7458 | {"verified":12} |
| icv1-fp16-t295-8d3c3165 | abs64-subtract2 | 12 / 12 | 137.3309 | 70.8733 | {"verified":12} |
| icv1-fp16-t295-8d3c3165 | abs64-subtract3 | 12 / 12 | 23.4979 | 12.1267 | {"verified":12} |
| icv1-fp16-t295-8d3c3165 | rho | 12 / 12 | 1.9377 | 1.0000 | {"verified":12} |
| icv1-fp16-t295-8d3c3165 | rho-aa | 12 / 12 | 1.9377 | 1.0000 | {"verified":12} |
| icv1-fp18-tm337-d28d5e09 | abs32-subtract2 | 12 / 12 | 282.4994 | 177.1919 | {"verified":12} |
| icv1-fp18-tm337-d28d5e09 | abs64-mitm2 | 12 / 12 | 13.4152 | 8.4144 | {"verified":12} |
| icv1-fp18-tm337-d28d5e09 | abs64-subtract2 | 12 / 12 | 278.5684 | 174.7262 | {"verified":12} |
| icv1-fp18-tm337-d28d5e09 | abs64-subtract3 | 12 / 12 | 17.5888 | 11.0322 | {"verified":12} |
| icv1-fp18-tm337-d28d5e09 | rho | 12 / 12 | 1.5943 | 1.0000 | {"verified":12} |
| icv1-fp18-tm337-d28d5e09 | rho-aa | 12 / 12 | 1.5943 | 1.0000 | {"verified":12} |
| icv1-fp20-t727-cd198a38 | abs32-subtract2 | 12 / 12 | 524.8381 | 536.5484 | {"verified":12} |
| icv1-fp20-t727-cd198a38 | abs64-mitm2 | 12 / 12 | 11.2446 | 11.4955 | {"verified":12} |
| icv1-fp20-t727-cd198a38 | abs64-subtract2 | 12 / 12 | 510.0127 | 521.3923 | {"verified":12} |
| icv1-fp20-t727-cd198a38 | abs64-subtract3 | 12 / 12 | 20.6335 | 21.0938 | {"verified":12} |
| icv1-fp20-t727-cd198a38 | rho | 12 / 12 | 0.9782 | 1.0000 | {"verified":12} |
| icv1-fp20-t727-cd198a38 | rho-aa | 12 / 12 | 0.9782 | 1.0000 | {"verified":12} |
| icv1-fp22-tm1385-475dcb5f | abs32-subtract2 | 12 / 12 | 1081.0897 | 716.3848 | {"verified":12} |
| icv1-fp22-tm1385-475dcb5f | abs64-mitm2 | 12 / 12 | 10.1679 | 6.7377 | {"verified":12} |
| icv1-fp22-tm1385-475dcb5f | abs64-subtract2 | 12 / 12 | 819.1889 | 542.8360 | {"verified":12} |
| icv1-fp22-tm1385-475dcb5f | abs64-subtract3 | 12 / 12 | 35.3830 | 23.4466 | {"verified":12} |
| icv1-fp22-tm1385-475dcb5f | rho | 12 / 12 | 1.5091 | 1.0000 | {"verified":12} |
| icv1-fp22-tm1385-475dcb5f | rho-aa | 12 / 12 | 1.5091 | 1.0000 | {"verified":12} |
| icv1-f2m13-tm9-6863a838 | rho | 8 / 8 | 3.5026 | 1.0000 | {"verified":8} |
| icv1-f2m13-tm9-6863a838 | rho-aa | 8 / 8 | 3.5026 | 1.0000 | {"verified":8} |
| icv1-f2m13-tm9-6863a838 | sub3-mitm3 | 8 / 8 | 127.2247 | 36.3234 | {"verified":8} |
| icv1-f2m13-tm9-6863a838 | sub4-mitm3 | 8 / 8 | 46.7808 | 13.3562 | {"verified":8} |
| icv1-f2m17-t531-2048fb7c | rho | 8 / 8 | 2.0296 | 1.0000 | {"verified":8} |
| icv1-f2m17-t531-2048fb7c | rho-aa | 8 / 8 | 2.0296 | 1.0000 | {"verified":8} |
| icv1-f2m17-t531-2048fb7c | sub3-mitm3 | 8 / 8 | 311.3590 | 153.4071 | {"verified":8} |
| icv1-f2m17-t531-2048fb7c | sub4-mitm3 | 8 / 8 | 159.0557 | 78.3670 | {"verified":8} |
| icv1-f2m19-tm649-a9133fe6 | rho | 8 / 8 | 1.5420 | 1.0000 | {"verified":8} |
| icv1-f2m19-tm649-a9133fe6 | rho-aa | 8 / 8 | 1.5420 | 1.0000 | {"verified":8} |
| icv1-f2m19-tm649-a9133fe6 | sub3-mitm3 | 5 / 8 | 236.8360 | — | {"exhausted":3,"verified":5} |
| icv1-f2m19-tm649-a9133fe6 | sub4-mitm3 | 8 / 8 | 201.9391 | 130.9581 | {"verified":8} |
| icv1-f2m23-tm4129-45500fce | rho | 8 / 8 | 1.3823 | 1.0000 | {"verified":8} |
| icv1-f2m23-tm4129-45500fce | rho-aa | 8 / 8 | 1.3823 | 1.0000 | {"verified":8} |
| icv1-f2m23-tm4129-45500fce | sub3-mitm3 | 0 / 8 | — | — | {"exhausted":8} |
| icv1-f2m23-tm4129-45500fce | sub4-mitm3 | 0 / 8 | — | — | {"exhausted":8} |
| icv1-f2m13-t181-515ee569 | f6-ic | 0 / 2 | — | — | {"exhausted":2} |
| icv1-f2m13-t181-515ee569 | inherited-f4 | 0 / 2 | — | — | {"exhausted":2} |
| icv1-f2m13-t181-515ee569 | rho | 2 / 2 | 14.6759 | 1.0000 | {"verified":2} |
| icv1-f2m13-t181-515ee569 | rho-aa | 2 / 2 | 14.6759 | 1.0000 | {"verified":2} |
| icv1-f2m13-t181-515ee569 | rootfind-control | 0 / 2 | — | — | {"exhausted":2} |
| icv1-f2m19-t797-b6cf2467 | f6-ic | 0 / 2 | — | — | {"exhausted":2} |
| icv1-f2m19-t797-b6cf2467 | inherited-f4 | 0 / 2 | — | — | {"exhausted":2} |
| icv1-f2m19-t797-b6cf2467 | rho | 2 / 2 | 3.0459 | 1.0000 | {"verified":2} |
| icv1-f2m19-t797-b6cf2467 | rho-aa | 2 / 2 | 3.0459 | 1.0000 | {"verified":2} |
| icv1-f2m19-t797-b6cf2467 | rootfind-control | 0 / 2 | — | — | {"exhausted":2} |
| icv1-f2m23-t5197-69e76b73 | f6-ic | 0 / 2 | — | — | {"exhausted":2} |
| icv1-f2m23-t5197-69e76b73 | inherited-f4 | 0 / 2 | — | — | {"exhausted":2} |
| icv1-f2m23-t5197-69e76b73 | rho | 2 / 2 | 1.1604 | 1.0000 | {"verified":2} |
| icv1-f2m23-t5197-69e76b73 | rho-aa | 2 / 2 | 1.1604 | 1.0000 | {"verified":2} |
| icv1-f2m23-t5197-69e76b73 | rootfind-control | 0 / 2 | — | — | {"exhausted":2} |
| icv1-f2m31-tm90707-c95f16f5 | f6-ic | 0 / 2 | — | — | {"exhausted":2} |
| icv1-f2m31-tm90707-c95f16f5 | inherited-f4 | 0 / 2 | — | — | {"exhausted":2} |
| icv1-f2m31-tm90707-c95f16f5 | rho | 2 / 2 | 1.2694 | 1.0000 | {"verified":2} |
| icv1-f2m31-tm90707-c95f16f5 | rho-aa | 2 / 2 | 1.2694 | 1.0000 | {"verified":2} |
| icv1-f2m31-tm90707-c95f16f5 | rootfind-control | 0 / 2 | — | — | {"exhausted":2} |
| icv1-f2m19-t797-b6cf2467 | direct | 8 / 8 | 19.2284 | 6.2226 | {"verified":8} |
| icv1-f2m19-t797-b6cf2467 | orbit-folded | 8 / 8 | 4.1025 | 1.3276 | {"verified":8} |
| icv1-f2m19-t797-b6cf2467 | rho | 8 / 8 | 3.0901 | 1.0000 | {"verified":8} |
| icv1-f2m19-t797-b6cf2467 | rho-aa | 8 / 8 | 3.0901 | 1.0000 | {"verified":8} |
| icv1-f2m23-t5197-69e76b73 | direct | 8 / 8 | 7.8934 | 6.9527 | {"verified":8} |
| icv1-f2m23-t5197-69e76b73 | orbit-folded | 8 / 8 | 2.3315 | 2.0536 | {"verified":8} |
| icv1-f2m23-t5197-69e76b73 | rho | 8 / 8 | 1.1353 | 1.0000 | {"verified":8} |
| icv1-f2m23-t5197-69e76b73 | rho-aa | 8 / 8 | 1.1353 | 1.0000 | {"verified":8} |
| icv1-f2m31-tm90707-c95f16f5 | direct | 8 / 8 | 14.6069 | 12.3489 | {"verified":8} |
| icv1-f2m31-tm90707-c95f16f5 | orbit-folded | 8 / 8 | 2.2417 | 1.8952 | {"verified":8} |
| icv1-f2m31-tm90707-c95f16f5 | rho | 8 / 8 | 1.1828 | 1.0000 | {"verified":8} |
| icv1-f2m31-tm90707-c95f16f5 | rho-aa | 8 / 8 | 1.1828 | 1.0000 | {"verified":8} |
| icv1-f2m37-tm534059-32aad96b | direct | 3 / 8 | 5.2241 | — | {"exhausted":5,"verified":3} |
| icv1-f2m37-tm534059-32aad96b | orbit-folded | 3 / 8 | 4.0678 | — | {"exhausted":5,"verified":3} |
| icv1-f2m37-tm534059-32aad96b | rho | 8 / 8 | 0.2574 | 1.0000 | {"verified":8} |
| icv1-f2m37-tm534059-32aad96b | rho-aa | 8 / 8 | 0.2574 | 1.0000 | {"verified":8} |
| icv1-f2m11-t67-05f5aa36 | compact-control | 8 / 8 | 48.0597 | 2.0277 | {"verified":8} |
| icv1-f2m11-t67-05f5aa36 | gamma-weight-1 | 8 / 8 | 17.3176 | 0.7306 | {"verified":8} |
| icv1-f2m11-t67-05f5aa36 | gamma-weight-2 | 8 / 8 | 46.1181 | 1.9457 | {"verified":8} |
| icv1-f2m11-t67-05f5aa36 | gamma-weight-3 | 8 / 8 | 44.2287 | 1.8660 | {"verified":8} |
| icv1-f2m11-t67-05f5aa36 | rho | 8 / 8 | 23.7020 | 1.0000 | {"verified":8} |
| icv1-f2m11-t67-05f5aa36 | rho-aa | 8 / 8 | 23.7020 | 1.0000 | {"verified":8} |
| icv1-f2m23-tm5197-1f85e9e1 | compact-control | 8 / 8 | 2.2833 | 2.6757 | {"verified":8} |
| icv1-f2m23-tm5197-1f85e9e1 | gamma-weight-1 | 1 / 8 | 3.1227 | — | {"exhausted":7,"verified":1} |
| icv1-f2m23-tm5197-1f85e9e1 | gamma-weight-2 | 8 / 8 | 2.0658 | 2.4209 | {"verified":8} |
| icv1-f2m23-tm5197-1f85e9e1 | gamma-weight-3 | 8 / 8 | 2.0951 | 2.4552 | {"verified":8} |
| icv1-f2m23-tm5197-1f85e9e1 | rho | 8 / 8 | 0.8533 | 1.0000 | {"verified":8} |
| icv1-f2m23-tm5197-1f85e9e1 | rho-aa | 8 / 8 | 0.8533 | 1.0000 | {"verified":8} |
| icv1-f2m11-t67-05f5aa36 | m3-lp0 | 0 / 8 | — | — | {"exhausted":8} |
| icv1-f2m11-t67-05f5aa36 | m3-lp1 | 7 / 8 | 126.1067 | — | {"exhausted":1,"verified":7} |
| icv1-f2m11-t67-05f5aa36 | m3-lp2 | 8 / 8 | 55.3047 | 2.3333 | {"verified":8} |
| icv1-f2m11-t67-05f5aa36 | m4-lp0 | 8 / 8 | 374.7798 | 15.8121 | {"verified":8} |
| icv1-f2m11-t67-05f5aa36 | m6-lp0 | 8 / 8 | 1161.9226 | 49.0220 | {"verified":8} |
| icv1-f2m11-t67-05f5aa36 | mn-1-lp0 | 8 / 8 | 8797.0309 | 371.1508 | {"verified":8} |
| icv1-f2m11-t67-05f5aa36 | rho | 8 / 8 | 23.7020 | 1.0000 | {"verified":8} |
| icv1-f2m11-t67-05f5aa36 | rho-aa | 8 / 8 | 23.7020 | 1.0000 | {"verified":8} |
| icv1-f2m19-t797-b6cf2467 | checked-rank | 8 / 8 | 4.3208 | 1.3983 | {"verified":8} |
| icv1-f2m19-t797-b6cf2467 | full-rank | 8 / 8 | 4.1135 | 1.3312 | {"verified":8} |
| icv1-f2m19-t797-b6cf2467 | rho | 8 / 8 | 3.0901 | 1.0000 | {"verified":8} |
| icv1-f2m19-t797-b6cf2467 | rho-aa | 8 / 8 | 3.0901 | 1.0000 | {"verified":8} |
| icv1-f2m19-t797-b6cf2467 | shared-rank | 8 / 8 | 2.8644 | 0.9270 | {"verified":8} |
| icv1-f2m19-t797-b6cf2467 | target-pin | 8 / 8 | 4.1025 | 1.3276 | {"verified":8} |
| icv1-f2m23-t5197-69e76b73 | checked-rank | 8 / 8 | 2.3994 | 2.1134 | {"verified":8} |
| icv1-f2m23-t5197-69e76b73 | full-rank | 8 / 8 | 2.3331 | 2.0550 | {"verified":8} |
| icv1-f2m23-t5197-69e76b73 | rho | 8 / 8 | 1.1353 | 1.0000 | {"verified":8} |
| icv1-f2m23-t5197-69e76b73 | rho-aa | 8 / 8 | 1.1353 | 1.0000 | {"verified":8} |
| icv1-f2m23-t5197-69e76b73 | shared-rank | 8 / 8 | 1.9003 | 1.6738 | {"verified":8} |
| icv1-f2m23-t5197-69e76b73 | target-pin | 8 / 8 | 2.3315 | 2.0536 | {"verified":8} |
| icv1-fp19-t241-6012dab8 | abscissa-mitm | 12 / 12 | 11.4171 | 9.6626 | {"verified":12} |
| icv1-fp19-t241-6012dab8 | abscissa-subtract | 12 / 12 | 317.7792 | 268.9462 | {"verified":12} |
| icv1-fp19-t241-6012dab8 | rho | 12 / 12 | 1.1816 | 1.0000 | {"verified":12} |
| icv1-fp19-t241-6012dab8 | rho-aa | 12 / 12 | 1.1816 | 1.0000 | {"verified":12} |
| icv1-fp19-t241-b992f548 | abscissa-mitm | 12 / 12 | 11.5799 | 10.1450 | {"verified":12} |
| icv1-fp19-t241-b992f548 | abscissa-subtract | 12 / 12 | 335.5223 | 293.9474 | {"verified":12} |
| icv1-fp19-t241-b992f548 | rho | 12 / 12 | 1.1414 | 1.0000 | {"verified":12} |
| icv1-fp19-t241-b992f548 | rho-aa | 12 / 12 | 1.1414 | 1.0000 | {"verified":12} |
| icv1-f2m83-tm6151469093347-debefd74 | ic-n83 | 0 / 1 | — | — | {"error":1} |

## Evidence limits and delivery

All ten diagnostic tracks were attempted; all ten original resolution criteria remain partial. No IC speedup, new asymptotic result, exact universal lower-bound crossing or new mathematical theorem is established. N83 IC was actually attempted and returned the missing-wide-implementation error. Degree31 has ord_31(2)=5 and cannot replace degree83/131 structure. Independent-host/L2/GPU/annealer gates remain blocked by missing configured runner access; [access problem](../../problems/ecbench-ten-question-runner-access.md), [open resolution gates](../../problems/ecbench-ten-question-resolution-gates.md). Local unreviewed evidence is not a merged or independently accepted result.

Build the native ecbench_open_questions, ecbench_open_questions_run, ecbench_open_questions_replay_failures and ecbench_open_questions_publish examples plus the ecbench binary. Freeze a **new directory**, run the native dispatcher, audit/supplement failures, analyze, and render the recorded output. Never reuse or overwrite a sealed session. Exact commands are in execution-manifest.json and COMMANDS.md.
