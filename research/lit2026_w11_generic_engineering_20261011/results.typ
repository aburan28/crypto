#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)
#set heading(numbering: "1.")
#show table: set text(size: 7.5pt)
#let R = "reports/run-20261011T2049Z/"

= lit2026 W11: generic-method engineering — results

Workstream W11 of the lit2026 speedup program, 2026-10-11; stacked on
aburan28/crypto\#1681. Run `run-20261011T2049Z`; host Apple M4 Pro (14 cores, shared,
load 72–155), *isolation L0*. Arms: `original` = `4b6e09b78`; `previous`/`candidate` =
this branch with its switch off/on. Deterministic counters are the gate; wall ratios are
ratio of medians with [min, max] of per-round ratios. No exponent, deployed-curve or
IC-versus-rho claim.

#table(columns: 3,
  [*item (paper)*], [*built*], [*verdict*],
  [(a) ePrint 2026/2102], [`bsgs_fast` kernels `paired`, `coz` (`BSGS_STEP`)], [speedup: paired 1.15–1.67× wall vs original, 40–56 bits; 16.75 → 3.71 M-eq per point],
  [(b) ePrint 2026/704], [x-only binary evaluators; one-inversion pair step (`BTW_XEVAL`)], [step 1.50× fewer `F_(q^8)` mults, 1.29× wall; class walk end to end: no change],
  [(c) ePrint 2025/1652], [cubical ladder for the k = 2 Tate pairing], [16–19% fewer `F_p` M-eq than Jacobian Miller (28–52-bit r); identical values],
  [(d) ePrint 2026/1823], [DP trail recording in `koblitz_strong_rho`], [diagnostic: competing-risk geometric law; fruitless abandonment truncates trails at dp = 12],
  [(e) arXiv 2508.14011], [Algorithm 3 in Rust + `prime_fast` rho], [diagnostic: 1.077 ± 0.020 × √(πn/4) over 720 solves; paper's 8/12/16-bit cards use an ascending search],
)

Papers: abstracts of 2026/2102, 2026/704, 2025/1652, 2026/1823 (PDFs returned HTTP 403
from this host; formulas here are our own derivations), full PDF of 2508.14011.

= (a) BSGS: co-Z / paired giant steps, windowed batch inversion

`paired`: affine centre C, offsets jA (j ≤ K), C ± jA share 1/(x(jA) − x(C)); non-centre
points x-only; one batch inversion over W(K+1) denominators. `coz`: Jacobian centre,
co-Z conjugate pairs, Z's inverted per window (our reading of 2026/2102).

#table(columns: 7,
  [curve], [kernel], [K], [M/pt], [S/pt], [I/pt], [M-eq/pt],
  [secp48k1], [chain], [–], [4.625], [1.000], [0.1250], [16.75],
  [secp48k1], [paired], [32], [2.539], [1.000], [0.00193], [3.71],
  [secp48k1], [coz], [32], [5.117], [2.054], [0.00193], [7.34],
  [secp56k1], [chain], [–], [4.625], [1.000], [0.1250], [18.88],
  [secp56k1], [paired], [32], [2.539], [1.000], [0.00193], [3.74],
)

#image(R + "bsgs_opcount.svg", width: 100%)

#table(columns: 8,
  [curve], [orig. s], [prev. s], [paired s], [coz s], [orig/paired], [prev/paired], [orig/prev],
  [secp40k1], [4.67], [4.28], [3.20], [3.49], [1.459 [1.054, 1.917]], [1.337 [1.223, 2.037]], [1.091 [0.655, 1.091]],
  [secp44k1], [5.63], [6.12], [4.27], [4.22], [1.316 [1.015, 2.002]], [1.430 [1.262, 1.737]], [0.920 [0.804, 1.399]],
  [secp48k1], [8.63], [9.55], [6.64], [8.16], [1.301 [1.299, 1.603]], [1.439 [1.203, 1.849]], [0.904 [0.847, 1.157]],
  [secp52k1], [15.55], [19.40], [13.58], [14.77], [1.145 [0.996, 1.976]], [1.429 [1.204, 2.466]], [0.801 [0.726, 1.002]],
  [secp56k1], [45.84], [48.99], [27.44], [35.32], [1.670 [1.349, 1.670]], [1.785 [1.190, 1.849]], [0.936 [0.860, 1.327]],
)

#image(R + "bsgs_wall.svg", width: 85%)

Answers and table sizes identical across arms. The `previous` arm (same binary, chain
kernel, now generic over a field trait) ran 0.80–0.94× of `original` at four sizes, with
intervals crossing 1; recorded, not explained.

= (b) x-only binary isogeny evaluation

Counted `F_(2^163)` mul / inv per step (two points); all evaluators equal `velu_x_map`,
the Vélu point map's abscissa lies on the `velu_codomain` curve.

#table(columns: 7,
  [ℓ (s)], [batch inv.], [kernel poly], [fraction], [fraction, 1 inv.], [product], [projective],
  [3 (1)], [2 / 2], [6 / 2], [4 / 2], [7 / 1], [8 / 2], [14 / 0],
  [5 (2)], [8 / 2], [8 / 2], [8 / 2], [11 / 1], [14 / 2], [20 / 0],
  [7 (3)], [14 / 2], [14 / 2], [12 / 2], [15 / 1], [20 / 2], [26 / 0],
)

#image(R + "binary_xonly.svg", width: 100%)

The product form x·Π((x xᵢ + √b)/(x + xᵢ))² equals (Π xᵢ)²·φ(x) (checked on every
kernel). Torsion-walk step at (19, 1), ℓ = 457 in `F_(q^8)`: 627 460 → 419 070 mults,
916 → 458 inversions, 1564 → 1211 µs per kernel. End-to-end class walk (12 rounds,
identical snapshots): walk seconds 1.280 / 1.164 / 1.198 (original / previous /
candidate), original/candidate 1.069 [0.714, 1.831]: *no change*.

= (c) k = 2 pairing: cubical ladder vs Miller

#table(columns: 7,
  [r bits], [affine Miller M+S / I], [Jacobian Miller M+S], [cubical M+S / I], [unsquare M+S], [(cub.+unsq.)/Jac. M+S], [same, M-eq],
  [16], [438 / 52], [506], [363 / 1], [64], [0.844], [0.909],
  [32], [996 / 120], [1148], [737 / 1], [147], [0.770], [0.830],
  [40], [1215 / 146], [1401], [912 / 1], [181], [0.780], [0.840],
  [52], [1611 / 194], [1857], [1188 / 1], [241], [0.770], [0.828],
)

#image(R + "biext_k2.svg", width: 75%)

The cubical ladder's reduced value is the square of the pairing on every instance;
`unsquare` (exponent (r+1)/2, counted) restores it. Values identical to Miller and to
`eds_mov::modified_tate`.

= (d) Distinguished-point statistics

#table(columns: 9,
  [curve], [dp], [DPs], [hits/θ], [trail/geo], [trail/c-r], [hazard], [KS geo / c-r / 95%], [steps/ideal],
  [K0/F2^53], [4], [1 945 590], [1.066], [0.993], [0.997], [2.97e-4], [0.0023 / 0.0006 / 0.0010], [0.82],
  [K0/F2^53], [8], [165 190], [1.001], [0.921], [0.991], [2.96e-4], [0.031 / 0.0042 / 0.0033], [1.18],
  [K0/F2^53], [12], [13 518], [0.982], [0.433], [0.950], [2.91e-4], [0.301 / 0.022 / 0.012], [1.58],
  [K0/F2^61], [4], [2 465 006], [1.066], [0.994], [0.998], [2.55e-4], [0.0019 / 0.0005 / 0.0009], [0.80],
  [K0/F2^61], [8], [183 269], [1.003], [0.932], [0.993], [2.56e-4], [0.026 / 0.0025 / 0.0032], [1.01],
  [K0/F2^61], [12], [15 276], [0.986], [0.482], [0.983], [2.54e-4], [0.264 / 0.0088 / 0.011], [1.37],
)

#image(R + "dp_stats.svg", width: 100%)

= (e) The secpkk1 challenge ladder

Cards 6/24/32/48 = Algorithm 3 as printed (descending); cards 8/12/16 = ascending search
from 2^(k−1). 720 verified solves at 30–60 bits: iterations/√(πn/4) mean 1.077
(s.d. 0.546, s.e. 0.020), median 0.998; DP start-up overhead estimate 0.12–0.18.

#image(R + "challenge_ladder.svg", width: 80%)

= Limits

L0 on a heavily shared host; (a) is not an ecbench method-level comparison (ecbench's
BSGS arms do not call `bsgs_fast`), and the rho walk does not share the BSGS step; (b)
scanned sect163-field curves, s ≤ 3, paper formulas unavailable; (c) u64 fields only,
squared pairing corrected by `unsquare`; (d) two curves, one walk configuration; (e)
8–16 seeds per level above 48 bits.
