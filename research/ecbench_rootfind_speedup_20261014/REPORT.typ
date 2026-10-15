#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set par(justify: true)
= Cheaper targets for the ecbench rootfind engine

2026-10-14 · study `research/ecbench_rootfind_speedup_20261014` · Conductor T-1 · *class: engineering* · host 4 vCPU Xeon 2.10 GHz, shared, isolation L0

== Summary

About half of the counted cost of the best complete rootfind arm (`pdp3-koblitz:m=2,engine=rootfind-negfold`) was the relation loop's target generator. The r-adding walk sets up 16 jumps and up to 16 restart offsets at two scalar multiplications each. `targets=stride` presents $R_0 + [k]S$. `targets=stride-renew` also renews the stride ($S <- S + T$) after every relation and steps past repeats.

Results:
- *n = 13:* stride-renew halves S against the walk, $0.531 [0.506, 0.557]$, and puts complete index calculus below rho-strong at d = 4, $0.770 [0.720, 0.824]$.
- *n = 19:* it still cuts the walk, $0.835 [0.741, 0.943]$, but index calculus stays $4.7 times$ above rho-strong.
- *rho-frob* (S ≈ 1.9 at n = 13) is 6× below every index-calculus arm.

No exponent moved; no security claim.

== Sessions (every record verified, every replay matched; A/A exactly 1.000)

#table(columns: 5, [id], [session], [curve], [arms], [records],
  [A], [`ECBS1h66ed5f6ccf76`], [Koblitz a=0, n=13, r=2003], [rho-strong, rho-frob, d5 walk/stride, d4 stride, A/A], [768],
  [B], [`ECBS1ha187c9428699`], [n=19, r=130873], [rho-strong, d6 walk/stride], [144],
  [C], [`ECBS1h97bbbf34584f`], [n=13], [rho-strong, rho-frob, d5 walk/stride/renew, d4 renew, A/A], [896],
  [D], [`ECBS1hb8d9ed954ea3`], [n=19], [rho-strong, d6 walk/stride/renew], [256])

== Hypotheses (paired candidate ÷ baseline, `ecbench compare`, cluster 95 % interval)

#table(columns: 4, [], [comparison], [ratio], [verdict],
  [H1 (A)], [d5 stride ÷ d5 walk], [0.584 [0.548, 0.629]], [pass],
  [H2 (A)], [d4 stride ÷ rho-strong], [0.745 [0.618, 0.845]], [pass (reference is a lower bound)],
  [H3 (B)], [d6 stride ÷ d6 walk], [1.146 [0.960, 1.363]], [contains 1; above the predicted band],
  [H4 (C)], [d5 renew ÷ d5 walk], [0.531 [0.506, 0.557]], [pass],
  [H5 (D)], [d6 renew ÷ d6 walk], [0.835 [0.741, 0.943]], [pass (both arms lower bounds)],
  [H6 (C)], [d4 renew ÷ rho-strong], [0.770 [0.720, 0.824]], [pass (reference is a lower bound)])

Descriptive comparisons:
- *renew ÷ stride:* 0.957 [0.910, 1.004] (C); 0.626 [0.545, 0.713] (D).
- *index calculus ÷ rho-strong at n = 19:* renew 4.68 [4.19, 5.20]; walk 5.61 [5.12, 6.08].

#figure(image("figures/s-by-arm.svg", width: 100%), caption: [Mean S (GAE/√r) with the recorded 95 % interval for every arm of sessions A–D. Circles are fully priced; triangles are lower bounds.])

#figure(image("figures/phases.svg", width: 92%), caption: [Mean counted work per complete solve by phase, sessions C and D (warm-ups excluded).])

== Mechanism

Diagnosis: in the inherited n = 13 claim, about 500 of the 703 relation-phase GAE were the walk's 46 scalar multiplications, spread over about 140 one-addition targets.

The plain stride saves that set-up, but it needs about 1.7× the trials at the same hit rate. Under one fixed $S$, two pairs of relations the same number of steps apart, with coinciding factor-base differences, give a dependency with $sum c = sum c k = 0$. That dependency cancels the target and pins nothing, so rows exceed rank:
- A: 14.9 rows at rank 12.1;
- B: 29.6 rows at rank 23.4.

Random targets end with rows = rank, which singles out the stride. Renewing the stride per relation (one addition each) removes these dependencies, and rows = rank again in C and D. Trials return to the walk's level (C 128.5 vs 123.6; D 3,988 vs 3,732), at about 6 scalar multiplications against the walk's 42.7 (C) and 63.4 (D).

== Scope and limits

- *rho-strong is not the best rho here.* At n = 13 its fixed jump-table set-up is about 526 of its roughly 640 GAE, and rho-frob reaches S ≈ 1.9 on the same targets. H2 and H6 do not show index calculus below rho in general, and B and D show it is not.
- *Lower bounds.* rho rows leave canonicalisations unpriced (rho-strong also table operations). n = 19 index-calculus rows leave AS solves, lookups and row operations unpriced (no calibration for that curve).
- *Exponent.* m = 2 decomposition is at best $N^(2\/3)$ against rho's $N^(1\/2)$. The target source moves a constant.
- *Selection.* d = 4 and stride-renew were chosen on exploratory seeds; every claim session uses fresh seeds.
- *Wall time.* L0, descriptive only.

Run reports: `reports/run-20261014T1700Z/` (A, B) and `reports/run-20261014T1728Z/` (C, D). Protocol: `PROTOCOL.md`, `AMENDMENT_1.md`.
