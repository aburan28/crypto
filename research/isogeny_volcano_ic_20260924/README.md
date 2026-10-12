# Isogeny-volcano index-calculus campaign

This campaign turns the isogeny/Gröbner ideas into a large, bounded experiment portfolio in the **research** repository. It does not add production attack code and does not claim that volcano descent makes ECDLP easier.

## Naming convention

A curve instance gets a stable, machine-readable identity independent of its nickname:

```
ICV1:<field>:<trace>:<order>:<j>:<end>:<level>:<path>:<model>
```

where:

- `field` = `fp-<p>` or `f2m-<m>-<modhash8>`;
- `trace` = signed Frobenius trace in decimal;
- `order` = `#E(F_q)` in decimal;
- `j` = canonical field-element encoding of the j-invariant;
- `end` = endomorphism-order discriminant when certified, otherwise `unk`;
- `level` = certified volcano level relative to the chosen ell-volcano, otherwise `unk`;
- `path` = root-relative sequence such as `r`, `r.d3.d3.h5`; each edge records direction and prime;
- `model` = first 12 hex chars of SHA-256 over the canonical field+curve model JSON.

The human slug is shorter:

```
icv1-f2m53-t<signed-trace>-ell3-L2-p<path>-<modelhash>
```

Do **not** use only `j`: twists can share j. Do **not** use only group order: an isogeny class shares it. The canonical JSON stores the full coefficients, field modulus/basis, subgroup, maps, conductor evidence, and parent edge certificate. Unknown mathematical facts stay `null`, never guessed.

Candidate algorithms use the existing tournament idea but a stricter ID:

```
ICCAN1/<curve-id>/<fb>/<relation>/<solver>/<symmetry>/<revision>
```

Example: `ICCAN1/icv1-f2m53-.../orbit8/chained-s3/f4/frob-r3/r01`.

This is compatible with `docs/ic/boundary_targets.json` schema v2: every candidate emits the required factor-base, decomposition, relation-yield, rank, end-to-end and vs-rho fields for the stage it claims. Stage-only candidates leave full-DLP ratios null.

## Experiment portfolio

The machine-readable registry is `experiments.json`. It deliberately contains many orthogonal cells. Workers should start with the first wave and only promote a mechanism when correctness gates pass.

### Wave A — topology and attribution
1. Enumerate certified ell=2,3,5,7 neighborhoods on odd-prime ordinary toy classes.
2. Enumerate separable ell!=2 neighborhoods on binary ordinary toy classes.
3. Certify horizontal/ascending/descending edges from endomorphism-order evidence.
4. Compare vertical pairs against horizontal same-level pairs.
5. Compare every neighbor against eight same-curve coordinate/basis changes.
6. Compare rebuilt codomain systems against explicit isogeny pullbacks.
7. Measure map numerator/denominator degrees and exceptional charts.
8. Search for graph-local cost valleys without using holdouts for selection.

### Wave B — factor bases
9. Transport the exact point base through the isogeny; require target-wise decomposition-count conservation.
10. Fresh equal-cardinality vector-subspace bases.
11. Frobenius-orbit union bases.
12. Uniform random point bases.
13. Low-Hamming-weight coordinate supports where meaningful.
14. Torsion-stable bases.
15. Native bases selected for sparse membership equations.
16. Native bases selected for low resultant growth.
17. Quasi-subfield-like supports as a separately labelled experimental family.
18. Factor-base cardinality sweep at fixed curve and m.
19. B^m/r density stratification.
20. Rank-increasing-row yield, not raw relation count.

### Wave C — relation formulations
21. Direct S3/S4/S5 where constructible.
22. Chained S3 for m=4 and m=5.
23. Symmetrized S4 versus unsymmetrized control.
24. Function-first/Nagao coefficient solver.
25. Fixed-span quotient-space rank testing.
26. Lazy support transformation.
27. Isogeny-transported support constraints.
28. Graph-ideal encoding versus rebuilt codomain equations.
29. Resultant-tree shape sweep.
30. Elimination-order sweep with frozen variable semantics.

### Wave D — symmetry
31. Permutation canonicalization.
32. Small-order torsion symmetrization.
33. Frobenius orbit canonicalization.
34. Fixed-target stabilizer versus moving-target orbit accounting.
35. Orbit-aware Gröbner variable ordering.
36. Orbit-block Macaulay rank profiles.
37. Circulant/normal-basis constraint representation.
38. Symmetry-aware sparse linear algebra.
39. Quotient solution lifting cost.
40. Short-orbit/repeated-point correction.

### Wave E — Gröbner/F4/F5
41. F4 versus F5/signature backend on identical serialized systems.
42. Degree-filtration full-column rank profile.
43. First observed fall versus completed solving-degree telemetry.
44. Macaulay sparsity/nonzero growth.
45. Matrix ordering and pivot strategy.
46. CPU dense/sparse crossover.
47. NVIDIA GPU elimination adapter when available.
48. AMD GPU adapter when available.
49. Apple Silicon adapter when available.
50. Preprocessing reuse across targets, cold and K=10/100 amortized.

### Wave F — Koblitz/subfield
51. Koblitz root versus isogenous model retaining a certified Frobenius action.
52. Neighbor where the cheap action is lost or maps to a conjugate curve.
53. Orbit-aware relation reuse.
54. Subfield k=2 and k=4 curves with relative odd extension.
55. Frobenius-invariant base versus transported base.
56. Actual orbit-length stratification.
57. Tau/Frobenius representation cost accounting.
58. Same experiment at n=9,11,13 then larger only after gates.
59. Binary model coordinate-basis controls.
60. Net benefit after automorphism-aware rho normalization.

### Wave G — cover/descent
61. F_(2^6)/F_(2^2) relative-degree-3 cover feasibility.
62. F_(2^12)/F_(2^4) counterpart.
63. Original versus isogenous source curve.
64. Explicit correspondence construction cost.
65. Genus and map-degree accounting.
66. Relation-generation cost on the descended object where implemented.

### Wave H — replication/accounting
67. Whole-trace-class holdouts.
68. Six held-out classes overall before replicated label.
69. Five seeds and three timing repetitions.
70. Full failed-target accounting.
71. Candidate-search/isogeny-construction cost.
72. Memory Pareto frontier.
73. Frozen-suite comparison under `research/index_calculus_baseline_20260914`.
74. Boundary ledger schema-v2 report.
75. Matched automorphism-aware rho reference.
76. Complete toy-DLP recovery and scalar verification.
77. Four-size whole-method exponent fit only after complete costs exist.
78. Negative-result archive for missing edges/bases/backends.
79. Mechanism witness extraction for every interesting degree drop.
80. Independent replay of promoted candidates.

## First worker wave

Launch independent workers on:

- **W1 topology:** fixtures, canonical IDs, certified vertical/horizontal edges.
- **W2 transport:** exact point-base transport and decomposition-count oracle.
- **W3 algebra meter:** F4/F5 trace semantics and full-column Macaulay profiles.
- **W4 factor bases:** native equal-cardinality structured/random/orbit bases.
- **W5 function-first:** Nagao/support/quotient variants against Semaev.
- **W6 symmetry:** torsion/Frobenius/permutation encodings and reconstruction.
- **W7 Koblitz:** Frobenius-retention/loss experiments on small binary fields.
- **W8 accounting:** schema-v2 boundary reports, frozen-suite adapter, matched rho.

Workers must preserve failures and emit PR-ready evidence. No worker may call a stage diagnostic a full speedup.

## Promotion gate

The existing repository contract remains authoritative. A candidate needs >=20% lower paired common operations per fresh verified relation at fixed support and m to pass the engineering gate, with correctness preserved, failed attempts charged, memory reported, independent curves/adjacent cells and holdouts. A full attack claim additionally requires all phases, the same-series rho comparison, and final scalar verification.

