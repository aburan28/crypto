#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)
#show raw: set text(size: 8.5pt)
#show figure.caption: it => it.body

= Factor base → degree of regularity: yield and correlations on the ECC2K-130 family

2026-10-14. Study directory `research/factor_base_dor_campaign_20261013/`; preregistered hypotheses in `PREREGISTRATION.md`; run report in `reports/run-20261014T0600Z/RUN_REPORT.md`. Class: stage diagnostic of the decomposition oracle (`D*` and work-per-relation proxies); no end-to-end speedup is claimed. Curve family: Koblitz $y^2 + x y = x^3 + 1$ over $F_(2^n)$, toy degrees $n = 9..15$ plus $n = 131$ structural checks; the challenge member is $n = 131$ (ECC2K-130). Engine: `sage/factor_base_dor.sage`, selftest clean.

== The question

Does the factor base — subspace dimension $n'$, algebraic structure, and Shantz–Teske slack $delta$ — change the degree of regularity $D^*$ (Macaulay refutation degree of the Weil-descended Semaev system on non-decomposable targets) and the decomposable yield enough to point at a drastic index-calculus speedup at $n = 131$?

== Findings (measured)

*F1 — multiplicative closure is the only structure that lowers $D^*$.* At matched $(n, n')$ the subfield base (closed under products) refutes at $D^* = 2$ exactly, every seed, every cell (n=9 and n=15, both dims), against 2.35–3.00 for standard/random; in the symmetric-$S_4$ m=3 encoding at $n = 15, l = 3$ the subfield base refutes at $D^* = 3$ with 9 variables against $D^* = 5$ with 15 variables for standard and for the Galbraith–Gebregiyorgis geometric special subspace — two degrees, i.e. roughly two orders of magnitude of Macaulay work at fixed nvars. Normal-basis windows and symmetry breaking (`nb_window`, `nb_break`) do *not* lower $D^*$ (they track random within ±0.2 at every matched cell).

*F2 — $D^*$ is set by the base dimension, not the field degree.* Pearson $r(D^*, n') = 0.87$ pooled over the 45 non-subfield m=2 aggregates; $r(D^*, n) = -0.04$. $D^*$ rises about one degree per $+2$ of $n'$ across n=11 and n=13 at every structure.

*F3 — the partial-closure gradient is visible below the subfield threshold.* The standard base $theta dot "span"{1, z, dots, z^(n'-1)}$ (products stay mostly inside $"span"{1..z^(2n')}$, partially inside the base) refutes measurably lower than a random base of the same dimension at every matched cell: 2.23 vs 2.79 (n=13, $n'=4$), 2.57 vs 2.88 (n=11, $n'=4$), 2.86 vs 3.00 (n=13, $n'=5$), 2.36 vs 2.99 (n=15, $n'=5$). Consistent across all three seeds.

*F4 — work per relation has an interior optimum in $n'$.* With $"W\/rel"$ = "Macaulay rows x cols at " $D^*$ / "decomposable fraction", the optimum sits at small dimensions and the largest cells are 20–300x worse: n=13 standard $n'=3$: $1.6 times 10^4$ against $n'=6$: $4.5 times 10^6$ (286x); n=11 standard $n'=4$: $1.0 times 10^5$ against $n'=5$: $4.0 times 10^5$; n=9 subfield $n'=3$: $7.8 times 10^3$ against random $n'=3$: $2.7 times 10^4$. Yield grows $4x$ per $+1$ of $n'$ (counting floor $|V|^2\/2^n$), but $D^*$ grows one degree per $+2$, and the Macaulay size at degree $D$ grows combinatorially with $2n'$ — the product is minimized around $n' approx n\/3..n\/4$ at toy scale.

*F5 — Shantz–Teske $delta$ buys 6–11x per relation at toy scale.* Sweeping $delta$ (base shrunk so $m n' = n - delta$) on the standard m=2 base: work per relation falls $7.9 times 10^5 -> 1.2 times 10^5$ (6.5x, n=11) and $5.4 times 10^6 -> 4.9 times 10^5$ (11x, n=13) from $delta = -1$ to $delta = +3$, while $D^*$ falls 4.00 → 2.65 and 4.00 → 2.87. Yield falls by $2^(-delta)$ as the model predicts, and the net still favours the shrink at these sizes.

*F6 — the symmetric $S_4$ encoding does not beat the chain on degree, only closure does.* Geometric == standard at both cells ($D^* = 4$ at n=9 l=2; $D^* = 5$ at n=15 l=3); the subfield base is the only one at $D^* = 3$. The special subspace buys variables (9 vs 15), not degree — confirming the Galbraith–Gebregiyorgis caveat at the next size up.

*F7 — $n = 131$ applicability.* The Frobenius-window fold was verified at the challenge degree itself: a normal element exists, $pi^("dim" V_i) = V_(i+1)$ holds, the windows are disjoint (`raw/frobwin_n131.log`). Weight-class invariant bases against windows at n=11/13 cut trials by 10x and 38x respectively (stacking tables) — an unknowns/relation-count lever, not a $D^*$ lever. But $x^131 - 1 = (x+1) dot ("irreducible of degree 130")$: *there is no proper subfield at n=131, so the one structure that measurably crushes $D^*$ (F1) is unavailable on the challenge curve itself.*

#figure(image("figures/dstar_vs_dim.svg", width: 100%), caption: [Figure 1. Mean refutation degree $D^*$ against factor-base dimension $n'$ (3 seeds, 24 targets/cell, m=2). Left: n=11; right: n=13. Structures: standard, random, nb_window, nb_break.])

#figure(image("figures/wpr_vs_dim.svg", width: 100%), caption: [Figure 2. Work-per-relation proxy against $n'$ (log scale). The optimum is interior; the largest bases are 20–300x worse per relation than the optimum.])

#figure(image("figures/delta_sweep.svg", width: 62%), caption: [Figure 3. Shantz–Teske $delta$ sweep, standard base, m=2. Per-relation work falls 6.5x (n=11) and 11x (n=13) from $delta = -1$ to $delta = +3$.])

== What this says about a drastic ECC2K-130 speedup

Every measured degree lever (F1, F3, F6) runs through multiplicative closure of the base, and none of the closure-free structures (windows, geometric special subspaces, weight classes) lower $D^*$. At $n = 131$ no proper subfield exists, so on the challenge curve the Gröbner-oracle track's per-system cost can only be trimmed by *size* levers (F4, F5: smaller effective bases, interior optimum near $n' approx n\/3$, $delta approx n\/3$) worth 10–300x per relation at toy scale, and by *arity* levers (the m=5/6-sum and 6-sum yield models of `RESEARCH_ECC2K130_IC_FEASIBILITY.md` §4), while the Frobenius-invariant weight-class bases (F7) cut the number of systems by ~$n$ with the windows verified at 131. The counting floor (yield $approx |V|^m\/2^n$) fixes how small a base can be before relations vanish; at $n = 131$ that floor is what forces $K >= 2 times 10^7$ columns in the 4-sum ladder, so the size levers trade against it and must be re-priced at m=83 under the repository's confidence gate before any transfer claim.

== Classification and falsification record

H1 confirmed (closure lowers $D^*$; windows do not). H2 confirmed ($r = 0.87$ with $n'$, interior optimum of $"W\/rel"$). H3 confirmed (6.5–11x at $delta = +3$; the $2^(-delta)$ yield cost is visible in the recorded fractions). H4 confirmed (geometric $=$ standard on $D^*$). H5 confirmed ($r(D^*, n) = -0.04$). Chained m=3 re-measurement: not completed (engine budget, preserved failures); prior $D^* = 5$ figure stands as prior evidence. All numbers are toy-scale stage diagnostics; nothing here is an end-to-end or sub-ρ claim, and the m=83 gate (AGENTS.md §8a) is unmet by this study.
