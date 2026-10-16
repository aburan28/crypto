#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set par(justify: true)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 14pt, weight: "bold")[`k = 4` index calculus on $E(FF_(p^4))$, measured]\
  #v(2pt)
  relations $prop p^(1.00)$, Wiedemann $prop p^(2.04)$, rho $prop p^(1.93)$;
  crossover from the fits $2^(100)$ ($h = 4$), sensitivity $2^(86)$–$2^(125)$\
  #v(2pt)
  `research/gaudry_quartic_d4_exponent_20261015` · 2026-10-15 · class: measurement
]

*Summary.* Gaudry/Diem index calculus with FGHR's $D_4$ oracle (PR \#1863) and rho
with the negation map were run end to end in ecbench on `fp4_torsion` curves for
12 primes $p = 29 dots 1009$: 288 records, all verified ($[d]G = Q$), `S` a count.
Fitted over the rungs, the relation phase grows like $p^(1.00)$ (95 %: 0.74–1.26),
the linear algebra (Wiedemann, new here) like $p^(2.04)$ (1.85–2.22), and rho like
$p^(1.93)$ (1.77–2.08), as the analysis says ($q^1$, $q^2$, $q^2$). Index calculus
still loses to rho at every measured size ($9.6 dot 10^5 times$ at $p = 29$,
$6.7 dot 10^4 times$ at $p = 1009$). From the measured fits, it meets rho at a group
of $2^(100.4)$ ($h = 4$; bootstrap $2^(98.0)$–$2^(103.6)$); with the relations slope
at its interval ends, $2^(85.9)$–$2^(124.9)$. The §11.16 model gave
$2^(105.4)$–$2^(106.5)$. *This extrapolates 15 bits in $p$, 60 bits in group
size, beyond the top rung.* At that point the relation phase is still 97 % of the
cost. Double large primes cut the linear-algebra slope to 1.57 (theory 1.5) but
cost 1.1–1.8× more on paired walks, and move the fitted crossover later
($2^(124)$).

= Method

Protocol `PROTOCOL.md`, registered at `08722a2b7` before any session; amendments:
(1) queue reordered and stop extended to 16 h (runtime, before any IC record);
(2) audit scope by size; (3) relaunch of 20 jobs after the container's CPU set
shrank. Code `0bea75038`. Arms (all `fghr-y-line:zero=1`, oracle
`gaudry-d4:macaulay=26,weighted=1,chains=1,sparsest=1,replay=1,charge_unsolved=1`):
`ic-d4-w` (`linalg=wiedemann`), `ic-d4-lp2`/`ic-d4-lp1` (`large_primes=2`/`1`,
small base the first $ceil(N^theta)$ columns, $theta = 3/4$, $7/9$), `rho-neg`.
One spec per (arm, $p$, target); the IC arms on a target share the walk. 59 of 72
jobs ran; 13 (third plain targets at five small rungs, second LP2 targets) did not
start before the stop.

*New code.* `linalg=wiedemann` solves for the target through a left-kernel vector
of the relation matrix, $O(m dot "nnz" + m^2)$ multiplications mod $r$, priced as
`row_ops`; the default dense elimination ran on the same rows as an unpriced
shadow and *gave the same logarithm in all 48 IC runs*. The large-prime wrapper
works on any summand oracle; the existing large-prime code is binary-field only.

*Fits.* OLS of $ln("rung mean")$ on $ln p$; interval = union of the Student-$t$
95 % interval and a bootstrap (10,000, targets resampled within rungs). The
crossover solves $sum_k A_k p^(beta_k) = overline(S)_"rho" sqrt(p^4 slash h)$.

#figure(image("figures/phases_vs_p.svg", width: 92%),
  caption: [Cost by phase and rho against $p$ (rung means; light dots: single
  targets). The dense elimination is derived from its `row_ops` at the same
  price and is not part of `S`.])

= Results

#table(columns: (auto, auto, auto, auto), stroke: 0.3pt, inset: 4pt,
  [*quantity*], [*$beta$*], [*95 % interval*], [*expected*],
  [relations phase (= total)], [*1.00*], [0.74–1.26 (bootstrap 0.96–1.04)], [$q^1$],
  [trials], [0.90], [0.64–1.15], [],
  [linear algebra, Wiedemann], [*2.04*], [1.85–2.22], [$q^2$],
  [dense elimination, same rows (derived)], [2.34], [2.11–2.57], [$q^3$ worst case],
  [rho (219 distinct runs)], [*1.93*], [1.77–2.08], [$q^2$],
)

The relations interval is wide because trials per relation vary between curves
(60 at $p = 29$, 76 at 61, 310–340 at 37 and 53, 81 at 1009), which the
within-rung bootstrap does not see. One $D_4$ solve costs $6.1 dot 10^5$ GAE at
$p = 29$ and $8.3 dot 10^5$ at $p = 1009$; 106 of 202,677 systems (0.05 %) were
degenerate and are charged. $overline(S)_"rho" = 1.08$ (1.02–1.15), drifting from
1.6 at $p = 29$ to 0.9–1.0 at the top.

#table(columns: (auto, auto, auto, auto), stroke: 0.3pt, inset: 4pt,
  [*crossover (group size)*], [*$h = 2$*], [*$h = 4$*], [*$h = 16$*],
  [measured fits, point], [$2^(98.4)$], [*$2^(100.4)$* ($p^* = 2^(25.1)$)], [$2^(104.6)$],
  [bootstrap within rungs], [$2^(96.0)$–$2^(101.4)$], [$2^(98.0)$–$2^(103.6)$], [$2^(102.0)$–$2^(107.9)$],
  [relations slope at 0.74 / 1.26], [$2^(84.3)$–$2^(122.1)$], [*$2^(85.9)$–$2^(124.9)$*], [$2^(89.3)$–$2^(130.6)$],
  [§11.16 model (\#1863)], [$2^(103.0)$–$2^(103.7)$], [$2^(105.4)$–$2^(106.5)$], [—],
)

With rho's own fitted power law in $p$: $2^(108.7)$ ($2^(103.5)$–$2^(116.1)$). The
top rung is a $2^(39.9)$ group, so the crossover is *60 bits of group size (15.1
bits of $p$) beyond the measurements*. The honest interval contains the model's
figure. With the dense elimination priced instead, the cost meets rho at
$2^(101.8)$ ($h = 4$) and rises above it again at $2^(126)$; at $h = 16$ it never
meets it: Wiedemann's $q^2$ is what keeps the window open.

#figure(image("figures/crossover.svg", width: 92%),
  caption: [Extrapolation of the fitted phases to the crossover with rho. Grey:
  the measured range. Blue band: bootstrap; hatched: relations slope at its
  interval ends; orange: the §11.16 model at $h = 4$.])

= Large primes

#table(columns: (auto, auto, auto, auto, auto), stroke: 0.3pt, inset: 4pt,
  [*arm*], [*runs (rungs)*], [*lin. alg. $beta$*], [*relations $beta$*], [*crossover $h = 4$*],
  [none], [27 (12)], [2.04 (1.85–2.22)], [1.00 (0.74–1.26)], [$2^(100.4)$],
  [two, $theta = 3/4$], [9 (8, $p <= 293$)], [*1.57* (1.12–2.01)], [1.23 (0.55–1.91)], [$2^(123.9)$],
  [one, $theta = 7/9$], [12 (7, $p <= 149$)], [1.62 (1.11–2.14)], [1.49 (0.43–2.55)], [$2^(178)$ (153–214)],
)

Paired $S$ ratios (same target and walk), two large primes: 1.00, 1.00, 1.12,
1.36, 1.32/1.09, 1.22, 1.80, 1.44 for $p = 29 dots 293$; one: 0.89–1.12 at
$p <= 53$, 1.78–3.45 at $p >= 61$. The matrix shrinks as the analysis says
(Wiedemann 4.7× cheaper at $p = 293$) but it is $10^(-6)$ of the cost here and
3 % at the crossover; the extra trials land on the relations phase. Large primes
would pay only where the linear algebra dominates, past about $2^(120)$.

#figure(image("figures/large_primes.svg", width: 92%),
  caption: [Linear algebra and total cost with and without large primes.])

#figure(image("figures/slopes.svg", width: 80%),
  caption: [All fitted exponents with their 95 % intervals.])

= Predictions

P1 held: 288/288 verified, 268/268 replays reproduced. P2 held (1.00 ∈ [0.9, 1.3]). P3 held (2.04).
*P4 not held* (dense elimination 2.34 < 2.4). P5 held (1.93; $S_"rho"$ slope
−0.15). P6: crossover $2^(100.4) in [2^(88), 2^(115)]$ held; *"no dense crossover
below $2^(200)$" not held* ($2^(101.8)$, recrossing at $2^(126)$). P7 held (1.57;
1.23 > 1.00; every paired ratio at $p >= 101$ above 1). P8 held.

= Audit

`ecbench verify` on all 60 sessions (receipts `audits/<job>.audit.json`; scope
fixed by amendment 2): *audit OK on all 60, 288 records, 0 problems.* Full replay
(`--replay-all`) on 50 sessions, every session with $p <= 197$, every rho-only
session and the spot checks `main-p293-t0` and `main-p401-t0`: *268 replays
reproduced exactly*. Integrity only (stored $[d]G = Q$ checks, every planned run
recorded, hashes) on the other 10: `main-p1009-t0/t1`, `main-p701-t0/t1`,
`main-p509-t0/t1`, `main-p401-t1`, `main-p293-t1/t2`, `lp2-p293-t0`.

= Limits

The exponents are fitted over $p = 29 dots 1009$ (groups $2^(19)$–$2^(40)$); the
crossover is an extrapolation of those fits, not a measurement, and its width is
set by curve-to-curve scatter in the relation rate. The large-prime fits rest on
one or two targets per rung over $p <= 293$. $S$ is a count on one host class;
wall time is not admitted. No `docs/bounds` domain exists for $FF_(p^4)$ curves,
so there is no challenge verdict. Study files: `README.md`, `raw/summary.json`,
`sessions/`, run report `reports/run-20261015T0457Z/`.
