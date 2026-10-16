#set page(paper: "a4", margin: (x: 2cm, y: 2cm), numbering: "1")
#set text(size: 10pt)
#set par(justify: true)
#show heading.where(level: 1): set text(size: 14pt)
#show heading.where(level: 2): set text(size: 11.5pt)

= ISO-1 continuation (2026-10-16): certified labels, the 8th-power mass, and the residual zero classes

Research evidence report. Study: `research/iso1_weak_classes_20261007`. Markdown version with links: `residual_model_20261016/REPORT.md`. Branch `aburan28/quirky-feynman-prppv2`, PR aburan28/crypto\#1954. Host: 4-core Linux container, shared; wall times are not isolated measurements.

== Summary

+ *Certified census labels through p = 79.* A new counter returns a group order only when it is the unique Hasse-interval candidate compatible with the tested points on the curve and its quadratic twist. Reruns at p = 7, 13, 37, 41, 43, 47, 53, 59 are byte-identical to the retained CSVs, so their zero labels are exact. p = 61, 67, 71, 73 and 79 are new.
+ *Exact 8th-power mass (conjecture).* The number $N_8(t)$ of 8th-power Legendre parameters in a class is a closed 2-adic factor times an odd-conductor Hurwitz number. It matches exhaustive censuses at p = 7, 11, 13, 17 and the registered holdouts p = 19 and p = 23 with 0 mismatches.
+ *Random thinning.* Beyond that structure, the norm-one torus samples each class at random. The parameter-free model predicts the residual zero count at the pristine holdouts p = 53, 59 and at the forward predictions p = 61, 67, 71, 73, 79, with binomial dispersion $1.00 plus.minus 0.01$.

Beyond the proved obstruction $v_2(f_pi) >= 2$, whether a class is weak is a draw with probability $1 - (1 - 8/(p^2-1))^(N_8(t)/12)$. The evidence does not support any compact deterministic criterion beyond that.

== 1. Certified point counts

$Q = p^6$, so the Hasse interval is exactly $[Q+1-2p^3, Q+1+2p^3]$. Every curve has full rational 2-torsion, so $4 | \#E$, and the twist has order $2Q+2-\#E$. A full baby-step/giant-step scan of $4P$ lists every candidate. Further points on the curve and the twist filter the list, and the order is returned only when one candidate remains. Source: `census_runtime/src/certified.rs`. Unit test: 12 random p = 7 curves against brute-force counts.

#table(
  columns: (auto, auto, auto, auto, 1fr),
  align: (right, right, right, right, left),
  table.header([*p*], [*point counts*], [*extra points*], [*wall s*], [*result*]),
  [7], [213], [5], [0.01], [identical to exact square-table census],
  [13], [2,423], [8], [0.2], [identical to retained CSV],
  [37], [156,523], [16], [32], [identical],
  [41], [235,901], [20], [56], [identical],
  [43], [285,363], [27], [110], [identical],
  [47], [407,193], [35], [118], [identical],
  [53], [658,243], [27], [347], [identical],
  [59], [1,010,651], [23], [941], [identical to retained p59 run],
  [61], [1,154,751], [41], [1,265], [new: 0 / 111,630 depth-1 weak; 546 / 111,630 high-depth zero],
  [67], [1,680,383], [30], [1,448], [new: 0 / 148,138 depth-1 weak; 640 / 148,136 high-depth zero],
  [71], [2,118,901], [34], [resumed], [new: 0 / 176,436 depth-1 weak; 634 / 176,434 high-depth zero],
  [73], [2,367,853], [27], [resumed], [new: 0 / 191,844 depth-1 weak; 622 / 191,844 high-depth zero],
  [79], [3,247,401], [35], [resumed], [new: 0 / 243,400 depth-1 weak; 746 / 243,398 high-depth zero],
)

Receipts: `certified_20261016/` (comparison receipt, stderr receipts, CSV hashes, queue script and status).

== 2. Exhaustive λ census and the 8th-power mass

`iso1_lambda_census` enumerates every $lambda$ in $F_Q$ or in a 2-power subgroup. For each it records the certified trace of $y^2 = x(x-1)(x-lambda)$, whether $lambda$ is a $2^k$-th power, and whether it lies in the norm-one torus $T$ (the weak curves).

- At p = 11 it agrees with the orbit census on all 666 rows. In each class, membership in $T$ thins $Omega = (F_Q^*)^(2^e)$ at random: dispersion 1.012, with 31 zero classes against 32.4 expected.
- 32nd powers at p = 17 are not 2-adically local; 8th powers are, at every prime tested. This matches the fact, quoted here and not reproved, that the cover $lambda = mu^N$ is a congruence modular curve only for $N | 8$.

With $Delta = t^2 - 4Q = f_pi^2 D_K$, $k = v_2(f_pi) >= 2$, $s = (D_K / 2)$ and $H_"odd"$ the Hurwitz sum with the 2-part of the conductor removed:

$ N_8(t) = H_"odd"(t) dot cases(
  (2-s) 2^(k-2) & "if" t equiv.not plus.minus (Q+1) med (mod 32),
  5(2-s) 2^(k-2) - 6(1-s) & "if" t equiv plus.minus (Q+1) med (mod 32),
) $

#table(
  columns: (1fr, auto, auto, auto),
  table.header([*test*], [*primes*], [*rows*], [*mismatches*]),
  [derivation (exhaustive 8th-power census)], [7, 11, 13, 17], [4,078], [0],
  [registered exactness holdout], [19], [3,249], [0],
  [registered exactness holdout], [23], [5,819], [0],
)

The formula is read off the tables and is not proved.

== 3. Residual zero classes: registered model V2

V1, which used class size only, was registered and failed on its training primes (dispersion about 5). V2 uses only small-prime information: weak orbits follow $"Binomial"(N_8/12, 8/(p^2-1))$, with nothing fitted. It is evaluated natively by `iso1_residual_model`.

#figure(
  image("figures/zero_rows.svg", width: 100%),
  caption: [Observed high-depth zero classes ($t > 0$) against registered V2 predictions with 95% Poisson-binomial intervals. Filled markers: certified censuses (grey: seen during V1 diagnosis; blue: pristine holdouts; orange: forward predictions, scored). Hollow: forward predictions whose census is running. Data: `figures/zero_rows.csv`; sources `model_v2_results.txt`, `predictions_v2.txt`.],
)

#table(
  columns: (auto, 1fr, auto, auto, auto, auto, auto),
  align: (right, left, right, right, left, right, right),
  table.header([*p*], [*role*], [*rows*], [*zeros*], [*predicted [95%]*], [*disp.*], [*AUC*]),
  [37], [seen], [12,321], [145], [168.8 [146.5, 191.1]], [0.981], [0.959],
  [41], [seen], [16,810], [198], [186.2 [162.8, 209.6]], [0.998], [0.969],
  [43], [seen], [19,414], [188], [195.9 [171.9, 219.9]], [1.003], [0.972],
  [47], [seen], [25,403], [212], [212.6 [187.6, 237.6]], [0.997], [0.978],
  [53], [pristine holdout], [36,517], [247], [240.4 [213.8, 266.9]], [1.005], [0.982],
  [59], [pristine holdout], [50,474], [270], [265.9 [238.0, 293.8]], [1.001], [0.986],
  [61], [forward], [55,815], [273], [275.3 [247.0, 303.6]], [0.992], [0.984],
  [67], [forward], [74,068], [320], [301.1 [271.5, 330.8]], [0.996], [0.987],
  [71], [forward, e = 4], [88,217], [317], [317.9 [287.5, 348.3]], [1.000], [0.989],
  [73], [forward, e = 4], [95,922], [311], [327.1 [296.3, 357.9]], [1.011], [0.989],
  [79], [forward, e = 5], [121,699], [373], [353.5 [321.5, 385.5]], [0.998], [0.993],
)

p = 37 lies 1.5 below its interval; every other prime is inside. $h_K$ alone ranks zero rows with AUC 0.71–0.78.

== 4. Scope and limits

- *Proved* (`THEOREM.pdf`): no ordinary class with $v_2(f_pi) = 1$ is weak, for every odd $p$.
- *Exact* for $p <= 79$: the certified class labels.
- *Conjecture*, exact on six primes: the $N_8$ formula.
- *Measured*: the binomial thinning model, at two pristine holdouts and five forward predictions. It describes the data statistically and decides no individual class. It does not show that no exotic arithmetic criterion exists.
- Class numbers come from PARI/GP 2.15.4 `quadclassunit`, which is GRH-conditional for large $|D_K|$. This affects the predictions, not the census labels.
- Pending: the censuses for p = 83–199. p = 71, 73 and 79 ran across several container sessions with checkpointing, so it has no single wall time. p = 199 alone needs about 4 days on 4 threads.
