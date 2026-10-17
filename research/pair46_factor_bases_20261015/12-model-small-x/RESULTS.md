# 12-model-small-x: small-coordinate factor bases in other curve models (pair 46 campaign)

2026-10-15. Protocol `PROTOCOL.md` (sha256 aedc3b41...a3f), frozen before data and not edited. Instrument: Rust binary `src/main.rs` (sha256 1dccee3b...484), one core (taskset -c 3), 596 s wall, run 2026-10-15T16:45:28Z to 16:55:24Z (`run.log`). Raw records and summaries: `results.json`. A smoke run at b=16 (same code and seeds) was made during debugging; its output was discarded and the full run regenerated it.

## Verdict

**Exact obstruction for pair 46. On toy curves the result is null: no model is distinguishable from a random factor base.**
For pair 46 itself (crater and floor curves and their twists), Montgomery, twisted Edwards and Hessian models do not exist over F_p. On 240 toy curves where the models exist, small-coordinate factor bases in the Montgomery u, Edwards y, Edwards x, Hessian x and Hessian x+y coordinates have the same m=2 decomposition coverage as Weierstrass small-x and as random sets of the same size. Every one of the 21 preregistered tests is "indistinguishable": its Bonferroni CI lies inside the ±1% margin. At b=22 the largest |mean| is 0.04%. Both positive controls are detected at every size.

## Exact statements

- **A1 (exact computation, pair 46).** For both curves, N = 5947697996161643 is odd and N mod 3 = 2. For the twist order 5v^2 = 5947698065141045 (= p+1+t, checked), the order is odd and 5v^2 mod 3 = 2. gcd(x^p - x, x^3+ax+b) = 1 for the crater and the floor curve, so neither has a rational 2-torsion point. N·P = O was checked on 3 random points of each curve.
- **A2 (theorem, standard).** A Montgomery or twisted Edwards model over F_p needs 4 | #E(F_p). Jacobi-quartic and Huff models need rational 2-torsion. Hessian, twisted Hessian and generalized Hessian models need a rational point of order 3. A1 and A2 together mean that none of these models exists for the pair-46 curves or their twists over F_p. The 2-torsion appears only over F_{p^3}, because the cubic is irreducible, and that is a Weil-descent setting outside this item. So the pair-46 instance of this item cannot run. Part B tests the general question on toy curves where the models exist.
- **A3 (identity).** For any set F of K points, the mean number of ordered representations R = P1+P2 over all R in E is exactly K^2/N. A choice of model can therefore change only the spread of the representations: their coverage and collisions.
- **A4 (exact, explains the only visible structure).** Negation-closed sets (W-x, M-u, E-y, H-s, R-sym) spend K/2 of their ~K^2/2 pairs on P + (-P) = O. Their coverage is therefore about 1/K below that of an unstructured set: 0.39%, 0.10% and 0.05% at K = 256, 1024 and 2048. The measured "vs R-asym" gaps match this: -0.23 to -0.42%, -0.06 to -0.13% and -0.02 to -0.07%. The E-x fibre {x_E = c} is {P, T-P} with T = (0,-1) of order 2. Since σ(P1)+σ(P2) = -(P1+P2), F+F is negation-closed and the pairs P+σ(P) = T are wasted. That is why E-x sits 1/K below R-asym (-0.35%, -0.11%, -0.03%) even though it is not negation-closed. A Hessian x-fibre is a collinear triple with sum O, and it creates no extra collisions; H-x vs R-asym is ≤ 0.05% in absolute value.
- **Map verification (exact, every curve).** Montgomery: B v^2 = u^3+Au^2+u was checked on every point. Twisted Edwards: a x^2+y^2 = 1+d x^2y^2 was checked on every affine-chart point. Hessian: all affine points were enumerated natively. The fitted cubic relation was verified on every point, the map is injective, #H = #W, collinear x_H-triples sum to O, and (x,y)→(y,x) maps to -P. No curve failed (`failures: []`).

## Comparative results (relative coverage C_arm/C_ctrl − 1, mean over curves, Bonferroni bootstrap CI, 20000 resamples)

| pop | b (K, curves) | arm vs control | mean | Bonferroni CI | verdict |
|---|---|---|---|---|---|
| M | 16 (256, 60) | W-x vs R-sym | +0.03% | [-0.24, +0.30]% | indistinguishable |
| M | 16 | M-u vs R-sym | -0.13% | [-0.36, +0.11]% | indistinguishable |
| M | 16 | E-y vs R-sym | -0.10% | [-0.41, +0.20]% | indistinguishable |
| M | 16 | E-x vs R-asym | -0.35% | [-0.57, -0.15]% | indistinguishable (A4) |
| M | 16 | R-sym2 vs R-sym (null) | -0.14% | [-0.43, +0.16]% | indistinguishable |
| H | 16 (256, 60) | W-x vs R-sym | -0.14% | [-0.41, +0.13]% | indistinguishable |
| H | 16 | H-x vs R-asym | +0.05% | [-0.08, +0.17]% | indistinguishable |
| H | 16 | H-s vs R-sym | -0.17% | [-0.50, +0.14]% | indistinguishable |
| M | 20 (1024, 40) | W-x / M-u / E-y | +0.01 / -0.01 / +0.02% | all within [-0.10, +0.13]% | indistinguishable |
| M | 20 | E-x vs R-asym | -0.11% | [-0.18, -0.05]% | indistinguishable (A4) |
| H | 20 (1024, 40) | W-x / H-x / H-s | -0.03 / -0.01 / -0.04% | all within [-0.14, +0.08]% | indistinguishable |
| M | 22 (2048, 20) | W-x / M-u / E-y / E-x | +0.04 / +0.01 / -0.00 / -0.03% | all within [-0.07, +0.10]% | indistinguishable |
| H | 22 (2048, 20) | W-x / H-x / H-s | -0.03 / +0.01 / +0.02% | all within [-0.09, +0.07]% | indistinguishable |

Mean coverage (measured): C ≈ 0.394 on every arm at every size, against the heuristic 1 - e^{-1/2} = 0.393. The null control R-sym2 is indistinguishable at all 6 (pop, size) cells.

Positive controls (detector power, vs R-sym): PC1 (DL interval ±kG) gives -98.0%, -99.50% and -99.75% at b = 16, 20 and 22. PC2 (half DL interval) gives -19.6% to -20.4%, with Bonferroni CIs of width ≤ 1.5%. Both are "different" in every cell. The instrument can therefore see a structured factor base whose effect is ~20× the equivalence margin, and ~200× the half-width of the b=22 CIs.

Semaev cross-check: after mapping each negation-closed arm back to Weierstrass X, the S3(x1, x2, X(R)) root test agreed with the enumeration bitmap on 312000 of 312000 targets, with 0 degenerate cases.

## Effect size statement (empirical)
At b=22 (K=2048, 20 curves per population), every model arm's coverage equals its matched random control to within |Δ| ≤ 0.10% (Bonferroni CI bounds). The only systematic offset across arms is the exact 1/K negation-pair effect (A4), and it is identical for Weierstrass small-x and the random negation-closed control. So the change of model does nothing measurable to m=2 decomposition.

## What was cut
Nothing from the preregistered design: all 3 sizes and all 120 + 120 curves ran within budget. Not in the design: m ≥ 3 decompositions, Jacobi-quartic and Huff models, and models over extension fields (for pair 46 those exist only over F_{p^3}, see A2).

## Reproduce
```
cd research/pair46_factor_bases_20261015/12-model-small-x
cargo build --release --offline && ./target/release/model_small_x 16,20,22 2100   # writes results.json (~10 min, 1 core)
python3 plot.py && python3 make_report.py                                        # figure + report.pdf (layout only)
```
