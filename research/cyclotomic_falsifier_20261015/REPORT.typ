#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9pt)
#set heading(numbering: "1.")
#align(center)[#text(size: 14pt, weight: "bold")[The cyclotomic γ-weight base's falsifier] \ Encyclopedia §8 entry D4, lead A2 — 2026-10-15 \ `research/cyclotomic_falsifier_20261015`]

*Verdict (pre-registered rule, PROTOCOL.md frozen at `9972ea0e8`).* The *yield clause* of A2 ("fail if its yield per column is lower") *survives, yield-neutral*: on all 15 $n = 23$ tests the cyclotomic base $C_w$ is inside the distribution of 1999 random Frobenius-stable bases with the same columns (effects −0.49% … +0.65%, smallest p = 0.363, critical 0.00333), while a planted positive control is detected on every cell (p = 0.001; −35% … −64% coverage, +246% … +293% energy). D4 therefore *survives* only in the sense that its yield is not lower; it is not higher either. Any value of the lead rests on the untested *solver clause*, whose only record ($n = 11$) points the wrong way. Cost (ecbench, operations only, L0): cyclotomic vs compact-orbit-scan at equal columns is costlier in 2 of 4 pairs (+2.2%, +12.2%, factor-base build phase) and unresolved in 2.

= Question
A2 (`docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md` L187–221): $C_w = {sum_(k in S) gamma^k : |S| = w}$, $gamma$ a primitive $(2n+1)$-th root of unity, is a Frobenius-stable base with an algebraic membership test. Falsifier: at $n = 23$, "fail if the cyclotomic system's solving degree is not lower, or its yield per column is lower". The encyclopedia recorded the falsifier as not established; the 2026-10-12 census sampled 32 targets with no null population and its $n = 23$ solver cells timed out.

= Method
Every base is moved to discrete-log space (full table of $[k]G$) and measured exactly over all $r$ targets: $"cov"_2 = |B+B|$, $"cov"_3 = |B+B+B|$, and the additive energy. The mean number of ordered $m$-decompositions per target is $|B|^m\/r$ for every base (identity), so at matched columns yield can differ only through sum collisions. Cyclotomic and compact bases are the `ic_framework` plug-ins' own output. Null: 1999 bases from uniform abscissae with identical lifting, cofactor projection and Frobenius/negation fold. Positive control: the λ-progression base. Cells at the transition size $K_m = "round"((m! r ln 2)^(1\/m)\/2n)$ capped at the base size. One-column bases are all unit multiples of ${plus.minus lambda^j}$ and hence identical (amendment A1): $C_1$ cannot differ from any 1-column base.

= Yield (measured, exact; N = 1999; critical p = 0.05/15)
#table(columns: 8, inset: 3.5pt,
[*rung*], [*w*], [*K*], [*metric*], [*$C_w$ effect*], [*p*], [*null 2.5–97.5%*], [*pos. control*],
[K1 n23], [2], [30], [cov2], [−0.022%], [0.951], [−0.76, +0.73%], [−40.1%],
[], [2], [30], [energy], [+0.176%], [0.704], [−0.94, +1.00%], [+279%],
[], [2], [6], [cov3], [−0.351%], [0.728], [−2.60, +1.44%], [−64.0%],
[], [3], [52], [cov2], [+0.060%], [0.842], [−0.58, +0.57%], [−40.8%],
[], [3], [52], [energy], [−0.210%], [0.553], [−0.71, +0.75%], [+293%],
[], [3], [6], [cov3], [−0.389%], [0.700], [−2.60, +1.44%], [−64.0%],
[K0 n23], [2], [17], [cov2], [−0.488%], [0.392], [−1.12, +1.01%], [−39.3%],
[], [2], [17], [energy], [+0.653%], [0.363], [−1.36, +1.47%], [+246%],
[], [2], [4], [cov3], [−0.387%], [0.828], [−4.70, +2.10%], [−58.5%],
[], [3], [37], [cov2], [−0.238%], [0.582], [−0.84, +0.86%], [−35.5%],
[], [3], [37], [energy], [+0.252%], [0.644], [−1.05, +1.08%], [+250%],
[], [3], [4], [cov3], [+0.018%], [0.996], [−4.70, +2.10%], [−58.5%],
)
$w = 1$: one column at $K_1$ (identical to every 1-column base, p = 1); no subgroup point at $K_0$ (no base). Compact-orbit-scan and the ONB comparator are also inside the null on every cell (README §3). At $n = 11$ ($r = 991$) the positive control is not detected at $K = 2$, so those cells are uninformative.

#figure(image("figures/yield-vs-null.svg", width: 100%), caption: [Exact yield relative to the null median, $n = 23$ cells; thin line = null min–max, bar = 2.5–97.5%; black edge = significant. Source `yield/yield-*.json`.])

= Cost (measured, ecbench, operations only, L0; IC S are lower bounds)
#table(columns: 4, inset: 3.5pt,
[*curve*], [*cyclotomic / compact, equal columns*], [*ratio of S [95%]*], [*verdict*],
[K1 n23], [`w=1` / `columns=1`], [0.884 [0.544, 1.413]], [unresolved],
[K1 n23], [`w=2` (30) / `columns=30`], [1.022 [1.021, 1.023]], [cyclotomic costlier],
[K0 n23], [`w=2` (17) / `columns=17`], [1.122 [1.120, 1.125]], [cyclotomic costlier],
[K0 n23], [`w=3,columns=3` / `columns=3`], [0.977 [0.818, 1.161]], [unresolved],
)
Sessions `ECBS1h6e7055fb8555`, `ECBS1h6c2295d00f35`: 16 public targets × 2 rounds, 240/240 verified, replay-all 160/160 identical each. The costlier pairs differ in the factor-base phase (2247 vs 1661 and 2186 vs 1079 GAE in one record each) with identical oracle setup — build engineering, not yield (derived). Every IC arm is 2.7–30× rho-frob.

#figure(image("figures/cost-s-by-arm.svg", width: 100%), caption: [Mean S per arm with recorded 95% intervals.])

= Limits
The solver clause was not re-run; its only record (`a2_summary_n11.md`: F4 max degree 4 vs 3, 754.6 vs 7.3 ms, Macaulay 0/16 vs 3/16 resolved at $n = 11$, $w = 1$, $m = 2$) goes against the lead but is not at A2's $n = 23$. $n = 83$: no honest cheap cell (single-word plug-in; $r approx 2^81$). Metrics cover $m <= 3$ in the projected-subgroup model at transition sizes; not raw bases, $m >= 4$, or other sizes. No security or asymptotic claim. Run report: `reports/run-20261015T0728Z/report.pdf`.
