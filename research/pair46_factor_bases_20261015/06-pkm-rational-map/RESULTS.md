# 06-pkm-rational-map: results

Pair-46 factor-base campaign, item 06. Run 2026-10-15, one core (`taskset -c 3`), 1579 s wall,
pure Python. Protocol: `PROTOCOL.md`. Original frozen hash `8d8ea86d…81b5`; Amendment 1
(adding the lcm-degree metric) was made before the recorded run, and the amended hash
`baa68336…e7c5` is stored in `results.json`. A code smoke test at n1 <= 3 ran before the
amendment. Its output was discarded.

## Verdict

**Negative (empirical, toy sizes; heuristic for pair 46).** This run reproduces the defining
property of the Petit-Kosters-Messeng chained factor bases: the solving degree stays low.
The lcm step degree was 4-5 on all three chained arms (positive control: matches PKM
Table 4). Even so, the cost of one relation system grows like p^1.9-2.1 in this solver. PKM's
own published Magma times grow like p^1.51. Matched rho grows like p^0.52. So the projected
index-calculus relation stage has exponent about 2.4-2.6 here, or about 2.0 on PKM's own data,
against 0.5 for rho. It is 10^3 to 10^7 times slower than rho over the measured range, and the
gap widens with p. The chained bases do not change the relation yield (an exact result). They
were also slower than the same set given by its flat equation in this instrument.

## Source

C. Petit, M. Kosters, A. Messeng, "Algebraic approaches for the Elliptic Curve Discrete Logarithm
Problem over prime fields", PKC 2016, LNCS 9615, pp. 3-18, doi:10.1007/978-3-662-49387-8_1.
The author PDF was read in this session. It defines F = {(x,y) : L(x) = 0} with
L = L_n0 o ... o L_1 of small-degree maps (System (4)). Sec. 3.3 covers p-1 smooth subgroups
and Sec. 3.4 covers auxiliary-curve isogeny kernels, with Velu chains and their Lemma 1. Sec. 4
has Magma experiments with m = 2, |V| = 2^n1 and p >= 2^(2n1). The paper itself states that the
algorithms "are only practical for small parameters" and that "all attacks are outperformed by
generic discrete logarithm algorithms for practically relevant parameters". It makes no
complexity claim; it leaves the complexity open.

## Exact statements (checked by computation)

- **X1 (exact).** For the same set V, the chained and flat systems have the same solutions.
  Relation yields were identical target by target: A1 = A2 and A4 = A5 at every level, for
  example 12/12, 5/5 and 6/6 valid relations. A rational-map factor base changes only the
  algebra of membership, not which relations exist.
- **X2 (exact).** Every solution lies in V x V ⊂ F_p^2. For each GB, dim F_p[x]/I is at least
  the number of enumerated solutions. All GBs passed this check (`all_ok`). At n1 <= 3 all
  168 GBs matched sympy's reduced Groebner basis exactly.
- **X3 (exact, pair 46).** p - 1 = 2·3^4·7·241·631·34489703 and
  p + 1 = 2^4·11·21061·1604564779.
  - Instance 1 (Sec. 3.3), m = 2: the only subgroup of size at least p^(1/2) ≈ 2^26.2 has order
    2·3^4·7·241·631 ≈ 2^27.4. Its chain needs maps of degree 241 and 631, not small degree.
    With degree <= 7 the subgroup has order 1134 ≈ 2^10.1, which would need m ≈ 5.
  - Dickson/p+1 variant: the smooth part with B = 21061 is about 2^21.8.
  - Instance 2 (Sec. 3.4): an auxiliary curve with a 2^26-smooth subgroup is a precomputation
    of roughly p^(1/m) curve trials. It is not an obstruction.
  - The Frobenius prime v = 34489703 divides p - 1. This is a known pair-46 fact, and it is
    useless here because a degree-v map is not small.

## Measured (empirical), m = 2, 12 targets per level

Mean CPU seconds T for one relation system, using one pure-Python Buchberger
(grevlex, sugar, Gebauer-Moeller) on every arm. Dl = mean max lcm degree of a
nonzero-reducing pair (Amendment 1 metric). D = mean sugar degree (preregistered metric).

| n1 | p | A1 PKM p-1 chain T / Dl / D | A3 PKM isogeny chain T / Dl | A4 Dickson chain T / Dl | A2 flat x^d=1 T / Dl | A6 small-x T | A7 random T | rho DLP time (ops) |
|--|--|--|--|--|--|--|--|--|
| 2 | 17 | 0.0037 / 4 / 8.1 | 0.0071 / 5 | 0.0035 / 4 | 0.0017 / 6 | 0.0017 | 0.0020 | 1.2e-5 (12.8) |
| 3 | 73 | 0.074 / 4 / 9.0 | 0.100 / 5 | 0.072 / 4 | 0.012 / 10 | 0.0088 | 0.0105 | 3.3e-5 (25.5) |
| 4 | 257 | 1.20 / 4 / 10.0 | 1.62 / 5 | 1.23 / 4 | 0.045 / 18 | 0.043 | 0.044 | 4.3e-5 (50.5) |
| 5 | 1153 | 17.4 / 5 / 12.0 | 20.1 / 5 | 17.5 / 5 | 0.29 / 34 | 0.31 | 0.30 | 1.2e-4 (117) |
| 6 | 4289 | stopped | stopped | stopped | 2.08 / 66 | 1.96 | 1.94 | 3.9e-4 (226) |
| 7 | 17921 | stopped | stopped | stopped | 14.3 / 130 | 14.1 | 14.7 | 4.7e-4 (446) |

The Dickson arms (A4/A5) use their own prime p' = 19, 71, 271, 1087, 4159 and 17791. The
chained arms stopped after n1 = 5 under the preregistered rule (mean T > 15 s). The flat arms
have Dl = 2d + 2 exactly at every level, regardless of V.

Growth exponents (OLS slope of log2 mean T on log2 p, 95% bootstrap CI over targets):

| arm | slope | 95% CI | IC exponent (slope + 1/2) |
|--|--|--|--|
| A1 PKM p-1 chain (n1 2-5) | 2.02 | [2.00, 2.04] | 2.52 |
| A3 PKM isogeny chain (n1 2-5) | 1.91 | [1.88, 1.95] | 2.41 |
| A4 Dickson chain (n1 2-5) | 2.11 | [2.08, 2.13] | 2.61 |
| A2 subgroup flat (n1 2-7) | 1.29 | [1.28, 1.31] | 1.79 |
| A5 Dickson flat | 1.30 | [1.29, 1.32] | 1.80 |
| A6 small-x | 1.31 | [1.30, 1.32] | 1.81 |
| A7 random-x (null) | 1.28 | [1.26, 1.30] | 1.78 |
| PKM Table 4, published Magma, n1 7-11 (not our data) | 1.51 | n/a | ≈ 2.0 |
| rho group operations (ours) | 0.52 | n/a | 0.52 |

**IC/rho.** The projected IC relation stage ((|F|+10) / mean valid relations × mean T, with
linear algebra ignored, which favours IC) divided by the mean rho DLP time was:

- best chained arm at n1 = 5: 4.0e6, 95% CI [2.4e6, 7.6e6] (A4 Dickson chain);
- A1 at n1 = 5: 7.0e6, CI [4.1e6, 1.7e7];
- best flat arm at n1 = 7: 5.8e6, CI [3.3e6, 1.6e7] (A7).

At n1 = 2 every arm was between 1.7e3 and 4e4, except A7, which found 0 relations (ratio infinite). In the instrument's absolute timings this is a
Python constant times an exponent gap. Only the exponent gap carries information.

**Relation yield (M1).** This is valid relations divided by the heuristic 2|F|^2/N, pooled over
12 targets with a 95% Garwood interval. The 95% intervals of all seven arms A1-A7 contain 1 at
every level, except A7 at n1 = 5 (ratio 1.95, CI [1.01, 3.41]). A7 is the random null arm,
so with 42 intervals one such excursion is expected. No factor-base family shows a yield
effect.

## Controls

- **Degree positive control, preregistered sugar metric: FAILED as stated.** The sugar degree
  D_corr of A1 was 8-12 against PKM Table 4's 4-4-4-4-4 for n1 = 1..5. The sugar degree
  measures something different from Magma's F4 step degree.
- **Degree positive control, Amendment 1 lcm metric: PASSED.** A1 measured Dl = 4, 4, 4, 5
  against PKM's D^corr = 4, 4, 4, 4 for n1 = 2..5, within ±1 at every level. A3 measured 5 at
  every level, against PKM Table 5's isogeny instance at 4-5.
- **Relation-detector positive control (A8, planted known-log base).** The ratio was 1.0, 2.3,
  4.4, 8.4, 13.3 and 39 for n1 = 2..7. The preregistered threshold of 10x was met only at
  n1 = 6 and 7, so the yield detector has demonstrated power only at those two levels.
- **Negative control (A7, random-x, flat).** Its yield ratio was about 1, and its degrees and
  times were indistinguishable from small-x and the other flat arms. So in this solver the
  flat cost depends on d = |V| only.

## What was cut, and limits

- Chained arms at n1 = 6 and 7 were not run (stop rule: mean T > 15 s at n1 = 5).
- Pair 46 itself (p ≈ 2^52.4) was not run, because one system is out of budget by many orders
  of magnitude. Pair-46 statements are X3 (exact) plus extrapolation (heuristic).
  Extrapolating PKM's own published slope of 1.51 from T(2^22) ≈ 5e3 s gives
  T(2^52.4) ≈ 2^46 × 5e3 s per relation system, against about 2^26.4 group operations for rho.
- There was no F4 or msolve implementation (neither msolve nor Magma is installed). Our
  Buchberger is slower than PKM's Magma F4 (p^2.0 against p^1.5 on the chained systems), and in
  our solver the flat systems were faster than the chained ones. That ordering may not hold
  under F4. Neither solver gives an exponent below 1, and m = 2 needs one below 0 to beat rho.
- The rho baseline is the Floyd r-adding walk, charged at 3 group operations per step:
  446 measured operations against sqrt(πN/2) = 169 at n1 = 7. This slows rho and favours IC.

## Files and reproduce

- `run_pkm.py` (driver), `gb.py` (Groebner instrument), `results.json` (raw per-system records),
  `summary.json`, `analyze.py`, `fig_pkm.svg/png`, `fig_yield.svg/png`, `run.log`,
  `report.pdf` (`make_pdf.py`).
- To reproduce, from this directory:
  `taskset -c 3 python3 run_pkm.py 7 > run.log && python3 analyze.py && python3 make_pdf.py`

## Figures

The figures (`fig_pkm.svg` and `fig_yield.svg`) are embedded in `report.pdf`.
