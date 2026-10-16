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

#figure(image("figures/dstar_vs_dim.svg", width: 100%), caption: [Figure 1. Mean refutation degree $D^*$ against factor-base dimension $n'$ (3 seeds at n=11/13, 2 at n=17/19; 16-24 targets/cell, m=2). Structures: standard, random, nb_window, nb_break.])

#figure(image("figures/wpr_vs_dim.svg", width: 100%), caption: [Figure 2. Work-per-relation proxy against $n'$ (log scale). The optimum is interior; the largest bases are 20-300x worse per relation than the optimum.])

#figure(image("figures/delta_sweep.svg", width: 62%), caption: [Figure 3. Shantz-Teske $delta$ sweep, standard base, m=2, each curve normalized to its own $delta = -1$ cell. The per-relation gain grows with $n$: 6.5x (n=11), 11x (n=13), 14.5x and 246x at $delta = +3$/$+5$ (n=15), 7.9x (n=17), 19.7x (n=19).])

#figure(image("figures/sym4_delta.svg", width: 88%), caption: [Figure 4. Symmetric-$S_4$ (m=3) cells: $D^*$ by structure. The arity+size lever - a heavy $delta$ ($n - 3l$) on the symmetric encoding - drops $D^*$ from 5 (n=15, $l=3$, $delta=+6$) to 3 (n=15, $l=2$, $delta=+9$; and the subfield base at n=21, $l=3$), two degrees ~ two orders of magnitude of Macaulay work per system.])

== Round 2 (2026-10-14): scaling toward n = 131

The grid was extended to $n = 17, 19, 23$ (2 seeds, 16 targets) and the $delta$ sweeps to $n = 15, 17, 19$; three symmetric-$S_4$ cells were added ((n=15, l=2), (n=19, l=3), (n=21, l=3)). All round-2 cells completed after a first pass failed on a full temp volume (failures preserved in the first-pass logs, rerun with redirected `TMPDIR`). $n = 19$ has $"ord"_19(2) = 18 = n - 1$, the same full-length order as $"ord"_131(2) = 130$ - the closest cyclotomic fidelity available at toy scale.

*F8 - the size laws scale toward the challenge degree.* $D^*$ continues to track $n'$ ($r = 0.82$ pooled over all 133 aggregates) and saturates near 4-5 for the window structures; the standard base stays measurably lower (2.77 vs 4.00 at n=19, $n' = 7$). The interior work-per-relation optimum persists and sits at $n' approx n/2.7$: standard best $n' = 4, 3, 5, 7, 7$ at $n = 11, 13, 15, 17, 19$, with 100-875x spreads between the optimum and the largest cell at n=19.

*F9 - the $delta$ gain grows with $n$.* Normalized per-relation gains from $delta = -1$: 6.5x (n=11), 11x (n=13), 14.5x at $delta = +3$ rising to 246x at $delta = +5$ (n=15), 7.9x (n=17), 19.7x (n=19, $delta = +3$: $8.1 times 10^8 -> 4.1 times 10^7$). The yield cost is now visible in the data: at $n >= 17$ and $delta >= 5$ the measured decomposable fraction hits zero within 16 targets (counting floor $2^(-delta)$), so the sweep's usable range ends where the floor bites - exactly the trade the counting boundary predicts.

*F10 - the arity+size lever drops two degrees.* On the symmetric-$S_4$ m=3 encoding, shrinking the base from $l = 3$ to $l = 2$ at $n = 15$ ($delta: +6 -> +9$) drops $D^*$ from 5 to 3; at $n = 19, l = 3$ ($delta = +10$) the standard and geometric systems refute at $D^* = 4$ with 15 unknowns where the prior chained evidence sits at $D^* = 5$, and point-decomposable yield is still 0.0625 per target (geometric). The subfield cell at $n = 21, l = 3$ again refutes at $D^* = 3$ with 9 unknowns. Read together with F1: closure and heavy $delta$ are the only measured levers that move $D^*$ by whole degrees, and heavy $delta$ is the one that exists at $n = 131$.

== Round 3 (2026-10-14): the arity+size law `D* = l + 1`

Seven more symmetric-$S_4$ cells (n=21 l=2, n=19 l=2, n=23 l=2/l=3, n=29 l=3/l=4; n=21 l=4 lost to a PolyBoRi segfault, preserved in the logs) plus the n=23 m=2 top dims (9, 10). $n = 29$ has $"ord"_29(2) = 28 = n - 1$, 131's full-length cyclotomic order. All cells completed except the segfault.

*F11 - at heavy $delta$, $D^* = l + 1$, independent of the field degree.* Every heavy-$delta$ symmetric-$S_4$ cell sits on the diagonal: $D^* = 3$ at l = 2 for n = 19, 21, 23; $D^* = 4$ at l = 3 for n = 19, 23, 29; $D^* = 5$ at l = 4 for n = 29 (and the near-determined n=15 l=3 cell is the one off-diagonal point at $D^* = 5$). The m=3 system's degree is set by the base dimension l alone - the same law F2 found for m=2 - so the per-system Macaulay cost of the oracle is $2^(O(l))$ at any field size, including $n = 131$.

*F12 - orbit-amplified yield at heavy $delta$.* At n = 29, l = 3 (naive counting floor $|V|^3\/2^29 = 2^(-20)$) the measured point-decomposable fraction is 0.031 ($2^(-5)$) over 32 targets - about $2^(15)$ above the naive floor. The Frobenius phases and permutation orbits of the Koblitz family amplify the point-level yield far above the subfield-free counting bound; the same amplification appears at every heavy-$delta$ cell (n = 23 l = 2: measured $2^(-5)$ against floor $2^(-17)$). This is the quantity that decides whether a heavy-$delta$ oracle can still harvest relations at $n = 131$, and it must be re-priced at m = 83 before any transfer claim.

#figure(image("figures/sym4_delta.svg", width: 88%), caption: [Figure 4. Symmetric-$S_4$ cells through rounds 2-4: $D^*$ by structure against the $l + 1$ diagonal. At n = 37, l = 3 the geometric base refutes at $D^* = 3$ - below the diagonal: larger $delta$ lowers the degree further.])

== Correction (2026-10-15): F12/F13 withdrawn - the measured yield carried a seed-overlap artefact

`F12` (orbit-amplified yield) and `F13` (constant-in-n point-decomposable fraction) are *withdrawn*: the engine drew the random targets and the random/geometric structure elements from the *same* `rng_for(seed)` stream, so the first structure draws coincided with the first targets (verified directly at n = 29, seed 1: target[0] = 0x44cb63f = $theta$ of the geometric base, target[1] = g−1). A target that is itself a base element is trivially decomposable, which produced the constant "1 hit per 16-256 targets" the two findings rested on. With an independent target RNG the same cells give 0/256 and 0/4096 decomposable targets.

What replaces them is an #emph[exact] yield measurement, now built into the engine (`reachable_set_size`): enumerating the full reachable target set of the label. At (n = 15, l = 2) the reachable set has 16 elements - exact yield $2^(-11)$, where the artefact had measured $2^(-4)$. The yield therefore decays with n as the naive counting always said; the round-4 claim that it does not was the artefact. `F15` is reclassified accordingly: the per-system degree law (F11, $D^* = l + 1$) is unaffected - `D*` is measured on #emph[non-decomposable] targets and the overlap targets were exactly the decomposable ones - but the derived "per-relation cost independent of n" needs the relation count, which scales with $2^n \/ norm(S)$, and so the harvest cost is #emph[not] independent of n. The corrected fractions, reachable-set sizes and re-derived work-per-relation optima from the full rerun with the fixed engine land in the follow-up commit; the preserved seed-overlap era data is in `raw_v1_seedoverlap/`.

== Round 4 (2026-10-14): yield is governed by the e-space, not the field

High-precision yield cells (256 targets at l = 2, 128 at l = 3) over n = 23..41, including n = 37 ($"ord"_37(2) = 36 = n - 1$, 131's full-length order) and the m=2 grid at n = 29. All cells completed, zero failures.

*F13 - the point-decomposable fraction is constant in n.* At l = 2, random bases decompose 0.00781 of targets at every one of n = 23, 29, 31, 37, 41; geometric bases 0.00391 at every n. At l = 3: random 0.0234 and geometric 0.0078 at both n = 29 and n = 37. The amplification $A = "pfrac" dot 2^n \/ 2^(3l)$ is exactly $2^(n-13)$ (random, l = 2), $2^(n-14)$ (geometric, l = 2) - the counting is against a fixed effective space of ~$2^13$-$2^14$ (the symmetric e-space dimension), not $2^n$. The standard base yields zero point decompositions at every measured cell.

*F14 - larger $delta$ lowers $D^*$ below the diagonal.* At n = 37, l = 3 ($delta = +28$) the geometric system refutes at $D^* = 3$ (15 unknowns) - one degree below the round-3 diagonal, while the yield stays at 0.0078 per target. Degree and yield both improve with the overdetermination.

#figure(image("figures/yield_vs_n.svg", width: 72%), caption: [Figure 5. Point-decomposable fraction against the field degree at fixed l (round-4 precision cells: 256 targets at l = 2, 128 at l = 3). Every curve is flat: the counting is against the symmetric e-space, not $2^n$.])

*F15 - the derived per-relation oracle cost is independent of the field degree.* Combining F11-F14: per-system Macaulay cost $2^(O(l))$ (F11, improved by F14) and yield per target $g(l, "structure")$ independent of n (F13) mean the expected number of systems per relation is constant in n. The collection cost of the relation harvest therefore scales with the number of relations the ladder needs (the factor-base size and its Frobenius-orbit reduction), not with the subgroup size - the structural opposite of both the brute scan ($r\/n^2 K^2$ per relation) and rho ($sqrt r$). This is the drastic-speedup mechanism the campaign set out to find, measured at toy scale: the remaining obligations are the linear-algebra cost over the (orbit-reduced) factor base and the m=83 gate (AGENTS.md 8a).

== Round 3 verdict for the challenge curve

The measured lever set for an ECC2K-130 index-calculus oracle now reads: (i) $D^* = l + 1$ at heavy $delta$ - per-system cost grows only with the base dimension, never with $n$ (F11); (ii) the per-relation $delta$ gain grows with $n$ (F9); (iii) orbit structure amplifies yield ~$2^(15)$ above the naive floor at heavy $delta$ (F12); (iv) multiplicative closure, the strongest degree lever, does not exist at n = 131 (F1, F7). A drastic speedup at 131 therefore needs an m=3-or-higher oracle whose base dimension l stays small (l = 2..4 keeps $D^* <= 5$) while orbit-amplified yield keeps relations flowing - the m=83 gate is the named, unmet validation step.

== What this says about a drastic ECC2K-130 speedup

Every measured degree lever (F1, F3, F6) runs through multiplicative closure of the base, and none of the closure-free structures (windows, geometric special subspaces, weight classes) lower $D^*$. At $n = 131$ no proper subfield exists, so on the challenge curve the Gröbner-oracle track's per-system cost can only be trimmed by *size* levers (F4, F5: smaller effective bases, interior optimum near $n' approx n\/3$, $delta approx n\/3$) worth 10–300x per relation at toy scale, and by *arity* levers (the m=5/6-sum and 6-sum yield models of `RESEARCH_ECC2K130_IC_FEASIBILITY.md` §4), while the Frobenius-invariant weight-class bases (F7) cut the number of systems by ~$n$ with the windows verified at 131. The counting floor (yield $approx |V|^m\/2^n$) fixes how small a base can be before relations vanish; at $n = 131$ that floor is what forces $K >= 2 times 10^7$ columns in the 4-sum ladder, so the size levers trade against it and must be re-priced at m=83 under the repository's confidence gate before any transfer claim.

== Classification and falsification record

H1 confirmed (closure lowers $D^*$; windows do not). H2 confirmed ($r = 0.87$ with $n'$, interior optimum of $"W\/rel"$). H3 confirmed (6.5–11x at $delta = +3$; the $2^(-delta)$ yield cost is visible in the recorded fractions). H4 confirmed (geometric $=$ standard on $D^*$). H5 confirmed ($r(D^*, n) = -0.04$). Chained m=3 re-measurement: not completed (engine budget, preserved failures); prior $D^* = 5$ figure stands as prior evidence. All numbers are toy-scale stage diagnostics; nothing here is an end-to-end or sub-ρ claim, and the m=83 gate (AGENTS.md §8a) is unmet by this study.
