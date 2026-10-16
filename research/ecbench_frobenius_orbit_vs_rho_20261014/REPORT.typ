#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show link: underline

= Frobenius-orbit index calculus on Koblitz curves, fully priced, against the best rho

2026-10-14 · Conductor T-1 · Class: *engineering* · Study `research/ecbench_frobenius_orbit_vs_rho_20261014` · Protocol frozen in commit `e1708eef8` before claim data.

*Verdict: constant-only.* With every native priced, Frobenius-orbit IC ($m = 2$, `compact-orbit-scan` + `mitm-frobenius-counted`, `stride-renew` targets) is not below rho-frob on any of eight Koblitz curves (log#sub[2] $r$ = 11–39). The paired ratio grows with $r$: slope $beta = 0.165 [0.117, 0.213]$ per $ln r$ with $ln n$ in the fit (model: $1/6$). Orbits are 2–15× cheaper than the subspace base (constant), and beat rho-strong only where rho-strong's fixed set-up dominates (log#sub[2] $r <= 17$). No security claim.

== What was built
- *Full pricing.* New opt-in methods `ic.pipeline_full`, `rho.signed_frobenius_full`, `rho.signed_frobenius_strong_full` run unchanged solvers, then price every native counter from `src/cryptanalysis/ecbench/full_prices.json` (ratios to one affine addition measured by `ecbench price-units`, median of 15 passes, this host). Unknown counters stay unpriced. Existing methods and `docs/ic/calibration.json` unchanged; committed sessions replay as before.
- *`ecbench trend`*: weighted least squares of $ln$(paired ratio) on $ln r$ and $ln n$ from `compare --save` files.

== Panel and verification
Eight sessions, 16 fresh public targets × 3 rounds (+1 warm-up), L0, 4 vCPU Xeon 2.80 GHz shared. 2,816/2,816 records verified; 2,112/2,112 measured records replayed identically; A/A exactly 1.000 on every curve; no record a lower bound.

#table(columns: 8, align: (left, right, right, right, left, left, left, left),
  [*curve*], [*log#sub[2] r*], [*h*], [*c*], [*orbit ÷ rho-frob*], [*orbit ÷ rho-strong*], [*renew ÷ walk*], [*orbit ÷ negfold*],
  [a0 n13], [10.97], [4], [1], [2.956 [2.552, 3.242]], [0.411 [0.397, 0.423]], [0.403 [0.395, 0.412]], [0.511 [0.465, 0.562]],
  [a1 n17], [16.00], [2], [1], [1.719 [0.810, 2.980]], [0.631 [0.558, 0.711]], [0.503 [0.448, 0.562]], [0.234 [0.203, 0.270]],
  [a0 n19], [17.00], [4], [2], [2.805 [2.295, 3.236]], [0.739 [0.699, 0.782]], [0.556 [0.525, 0.591]], [0.148 [0.131, 0.170]],
  [a0 n23], [21.00], [4], [3], [3.847 [3.035, 4.819]], [1.379 [1.249, 1.521]], [0.753 [0.678, 0.834]], [0.065 [0.059, 0.071]],
  [a1 n29], [15.37], [12646], [1], [1.566 [0.992, 2.347]], [0.574 [0.556, 0.593]], [0.480 [0.464, 0.499]], [–],
  [a0 n31], [20.46], [1492], [2], [3.664 [2.907, 4.476]], [1.116 [1.022, 1.226]], [0.668 [0.607, 0.733]], [–],
  [a0 n37], [27.78], [596], [8], [5.048 [4.335, 6.086]], [4.676 [4.080, 5.310]], [1.029 [0.903, 1.176]], [–],
  [a0 n41], [39.00], [4], [90], [18.16 [15.67, 21.55]], [15.34 [13.16, 18.28]], [1.009 [0.938, 1.082]], [–],
)
Paired B ÷ A ratios of S from `ecbench compare`, cluster-bootstrap 95 % intervals. $h$ = cofactor, $c$ = orbit columns.

#figure(image("figures/ratio-trend.svg", width: 100%), caption: [Paired ratios against log#sub[2] $r$ with the `free_r_n` and fixed-model fitted values. Sessions as listed in README; `trend/trend-orbit-vs-rho-frob.json`.])

== Trend (H6)
#table(columns: 3, [*fit*], [*coefficients (95 %)*], [*χ²/dof*],
  [$c + beta ln r$], [$beta = 0.093 [0.062, 0.123]$], [9.9],
  [$c + beta ln r + gamma ln n$], [$beta = 0.165 [0.117, 0.213]$, $gamma = -1.18 [-1.90, -0.46]$], [3.9],
  [model $c + ln r\/6 - ln n\/2$], [$e^c = 1.84 [1.41, 2.40]$], [15.7])
$beta > 0$ in both free fits: constant-only by the frozen rule. The fixed model fails the χ²/dof ≤ 2 rule; the free $gamma$ rests on the two large-cofactor small-$r$ curves where $c$ is at its integer floor and rho-frob is heavy-tailed, so it is descriptive. Extrapolating `free_r_n` to $r approx 2^n$ (derived, not measured): the ratio is smallest near $n approx 10$ at about 3.3 and never below 1.

#figure(image("figures/s-by-arm.svg", width: 100%), caption: [Mean S of every arm on every curve with recorded 95 % intervals (`ecbench table`).])

== Pricing
The orbit base's formerly `*_uncharged` construction work costs 9.4 %, 6.9 %, 3.8 % and 0.16 % of the IC total at n = 13, 19, 23, 41: the old lower bound hid little. All natives are 36 % (n = 13) to 16 % (n = 41) of the IC total; hash probes are priced on a $2^20$-entry map, about 11× the pinned ratio, which leans against IC. With every native priced at zero on both sides (derived, ratio of means), IC ÷ rho-frob is 1.98, 1.28, 1.84, 2.77, 0.96, 2.60, 4.09, 16.3 in panel order, so the verdict does not hinge on prices.

#figure(image("figures/orbit-phase-split.svg", width: 100%), caption: [Phase shares of the orbit-renew arm (means of recorded phase GAE, `tables/derived-phase-means.jsonl`).])

== Limits
Toy and medium sizes on one shared host (L0; wall descriptive); prices are one host's microbenchmark ratios; $c$ and $d$ chosen on exploratory seeds (disclosed); negfold only at n ≤ 23 with cofactor ≤ 4; rho exponents fitted on these sizes are pre-asymptotic (rho-frob α = 0.238 [0.135, 0.303]); two a1 sessions record `git_dirty` from a whitespace-only `rustfmt` during the run (binary built from `e1708eef8`, unchanged).
