
### 7.A1-R. Experiment A1 at m = 31: measured (2026-10-10)

Cell `icv1-f2m31-tm90707-c95f16f5` (K₀/F_{2³¹}, r = 1 439 393, h = 1492),
`ic bench`, seed 123212651130, three matched repeats, eight counted rho runs
(rho S = 1.055). Bases: the ledger's Frobenius-invariant x-frame base on the
divisor 0;1;2 subspace (dimension 11), the +T-closed u-frame base on the same
V (the quotient-isogeny base of §4.4, i.e. F_u = V⁻¹(F₂)), and the trace-zero
base on divisor 1;2 (dimension 10, every point in [2]E). Each base was run
with the ledger's signed-Frobenius-orbit columns and with a new fold added
for this experiment, `projected=1`, which merges orbits whose cofactor
projections coincide (so P and P + T share one unknown; this is Prop. 4.1
applied to the column map). Oracle: `mitm-frobenius`, exact and complete,
so the yield and the column count are measured without a solver constant.
Full table: `research/isogeny_fibers_20261010/a1_m31/RESULTS.md`.

| arm | pts | cols | hit rate | GAE / indep. rel | S | S / S_x (same fold) |
|:--|--:|--:|--:|--:|--:|--:|
| x, m = 2 | 2049 | 35 → 33 proj | 0.0025 | 1047 | 34.4 → 31.0 | 1 |
| u (+T-closed), m = 2 | 2357 | 39 → **19** proj | 0.0021 / 0.0025 | 945 / 949 | 43.8 → **22.7** | 1.28 → **0.73** |
| x, m = 3 | 2049 | 35 → 33 | 0.444 | 3503 | 122 → 119 | 1 |
| u (+T-closed), m = 3 | 2357 | 39 → **19** | 0.207 / 0.219 | 10340 / 10387 | 235 → 140 | 1.92 → 1.18 |
| trace-zero (dim 10), m = 3 | 1179 | 20 → 19 | 0.232 | 4740 | 84 → 85 | **0.69** |

Reading:

1. **The projected fold is a free factor on the symmetrised base.** It halves
   the unknowns (39 → 19) at identical hit rate and identical cost per
   relation, and S falls by 0.52 (m = 2) and 0.60 (m = 3). On the x-frame base
   it changes nothing (35 → 33). The prior round's u-arm columns were counted
   before projection, so its K₀/2³¹ verdict was carrying a factor two it did
   not need to; with the fold the combinatorial u-arm *beats* the x-arm at
   m = 2 (0.73) where it lost before (1.28).
2. **The predicted ×k on relation collection does not appear.** §4.4 predicted
   unchanged yield per target and k = 2 fewer relations needed. The relations
   needed do halve, but the hit rate per target halves with them at m = 3
   (0.207 vs 0.444) and the cost per trial is 15 % higher (2357 vs 2049
   lookups), so the primary metric is **2.95× worse** per independent relation
   and S is 1.18× worse. The counting is exact: representations over a
   +T-closed base come in 2^{m−1} copies, so the distinct decomposable targets
   scale as |F_u|^m / 2^{m−1}; at m = 3 that is 2357³/4, against 2049³ for the
   x-base — ratio 0.38, measured 0.47.
3. **At equal unknowns the quotient base is dominated by the trace-zero base.**
   u-proj (19 columns, 2357 points) and trace-zero (19–20 columns, 1179 points)
   have the same hit rate (0.22 vs 0.23) — (2|F|)³/4 = 2·|F|³, exactly the
   class-cancellation gain the trace-zero base gets for free — but the
   trace-zero base pays half the lookups per trial, so it is 1.65× cheaper in
   S and 2.2× cheaper per relation. This is the §3.5 statement in numbers: the
   torsion structure the quotient isogeny exposes is already available as the
   subgroup [2]E, and a base chosen inside the subgroup costs less than a base
   closed under the kernel.
4. **Against rho** every arm is 22–223× worse at this size, unchanged in kind
   from the prior round.

Verdict for H1 on the combinatorial oracle at m = 31: the mechanism is real
(unknowns halve), its collection gain is cancelled by the representation
multiplicity, and the trace-zero base is the better use of the same
structure. The 1.5× kill criterion of §7.A1 is not met (0.34× per relation
at m = 3; 1.1× at m = 2). The F4 arms (`symmetrised` oracle with the
projected fold, m = 2 and m = 3) are reported below when they finish.
