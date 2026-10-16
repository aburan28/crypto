#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set par(justify: true)
#show heading: set text(size: 11pt)

= PKM's isogeny-tower factor base in ecbench, with its auxiliary-curve search charged

2026-10-15 · `research/pkm_isogeny_tower_ecbench_20261015` · protocol committed in `43986439e` before any measured run · class: *accounting*.

== Result

`prime-isogeny:t=T[,seed=S]` (`ic_framework::prime_tower`) builds PKM's isogeny tower $V = x(R' + ⟨Q⟩)$ for $Q$ of order $2^t$ on an auxiliary curve $E' \/ FF_p$ found from a public seed, and charges the search (residue-class baby-step giant-step filter on every curve, Mestre's order on survivors, then $Q$, $⟨Q⟩$, $R'$, $V$, the Vélu chain) to the factor-base phase in group-addition equivalents. All 1,080 runs (720 measured) verified; 8 audits `OK`, 720/720 replays identical.

*Measured.* IC on the isogeny base costs *1.00–1.33×* IC on a plain base of equal size (paired $S$); the excess is the search's share of $S$, 0.8–25 %, set by how many curves the seed tried (3–744 at $t = 8..11$), not by curve size. On the three $v_2(p-1) = 1$ curves (P-192 shape, no Kummer tower) it is 13–25 %. Every IC arm is 17–137× rho with the negation map. Panel 3: curves tried $approx C dot 2^t$ with pooled $C = 1.05$ (224 builds, $t = 8..14$), spread of 2–3 orders of magnitude between seeds.

*Derived.* At P-192 size ($log_2 p = 192$, $t = 48$) the search would cost $approx 2^74$ GAE ($2^82$ at $t = 64$) against rho's $2^95.8$; the `u64` code cannot run it.

#figure(image("figures/ratios.svg", width: 92%), caption: [Paired ratios of $S$ (95 % cluster bootstrap over 8 targets, 16 pairs) and the search's share of the isogeny arm's $S$. Sessions in `sessions/`, L0, counts only.])

== Paired ratios (ecbench compare)

#table(columns: 7, align: (left, right, right, right, right, right, right), stroke: 0.4pt,
  [*curve*], [*isogeny ÷ plain*], [*isogeny ÷ Kummer*], [*aux share*], [*isogeny ÷ rho*], [*rootfind ÷ mitm*], [*A/A*],
  [p20 `icv1-fp20-tm728-d215a2a1`], [1.035 [1.005, 1.064]], [1.031 [1.006, 1.056]], [5.1 %], [38.5], [1.196], [1.000],
  [p22 `icv1-fp22-tm2947-e15fed69`], [1.140 [1.134, 1.147]], [1.196 [1.188, 1.204]], [12.4 %], [44.4], [1.115], [1.000],
  [p24 `icv1-fp24-t2705-27cbb01b`], [1.002 [0.979, 1.028]], [0.976 [0.951, 1.002]], [1.3 %], [29.6], [1.078], [1.000],
  [p26 `icv1-fp26-t2838-8e779ca8`], [1.078 [1.070, 1.085]], [1.115 [1.097, 1.129]], [7.4 %], [113.5], [1.044], [1.000],
  [p28 `icv1-fp28-t4155-50f492b1`], [1.008 [1.004, 1.011]], [1.050 [1.047, 1.054]], [0.8 %], [136.5], [1.025], [1.000],
  [q20 `icv1-fp20-t199-d1684bf5`], [1.204 [1.190, 1.218]], [—], [17.4 %], [21.0], [1.168], [1.000],
  [q22 `icv1-fp22-t285-9fc7e3e5`], [1.331 [1.321, 1.343]], [—], [25.0 %], [52.2], [1.090], [1.000],
  [q24 `icv1-fp24-t4451-d5368b2d`], [1.172 [1.142, 1.206]], [—], [13.5 %], [38.3], [1.058], [1.000],
)
Every $S$ is a lower bound (prime-curve Legendre symbols, square roots, lookups and row operations counted, unpriced) for all IC arms alike. No wall time is reported (L0, shared host).

#figure(image("figures/s_by_arm.svg", width: 88%), caption: [Mean $S$ per arm, 16 measured runs each, every run verified.])

#figure(image("figures/panel3.svg", width: 92%), caption: [Panel 3: `ecbench fb` builds of `prime-isogeny:t=T,seed=S`, 16 seeds per $t$, two 24-bit curves.])

== Decision against the protocol

- P1 held (all verified, all replays identical). P2 held (tests; run `fb_id` = dump `fb_id` on all 8 curves).
- P3 failed on single draws (p20, p24, p28: tried/$2^t$ = 0.09, 0.006, 0.05), held for the 16-seed mean at every $t$ (0.61–1.27).
- P4 held on 7 of 8 curves; p24 not resolved (1.3 % share).
- P5 held by the letter (5.1 % → 0.8 %), not as a trend: the share follows the seed's draw.
- P6 held on 4 of 5 (p24: 0.976, Kummer base 4 abscissae larger). P7 held.
- Controls: A/A 1.000; rootfind ÷ mitm reproduces the 2026-10-13 values (1.212 … 1.026).

== Projection (derived)

GAE per curve $approx 1.5(b - t) + t + 2 dot 2^((b\/2 + 2 - t)\/2)$ (matches panel 3: 46 vs 50 at $t = 8$, 30 vs 29 at $t = 14$). With $C dot 2^t$ curves: P-192, $t = 48$: $2^48.1$ curves × $2^26$ = $approx 2^74$ GAE, baby table $2^25$; $t = 64$: $approx 2^82$. Rho: $sqrt(pi r\/4) approx 2^95.8$.

== Not done

Wide integers ($p < 2^62$, $t <= 20$); SEA, CM construction of $E'$, 2-Sylow by halving, $X_1(8)$ families; Vélu evaluations uncharged ($2t$ inversions); seed variance not in the per-curve intervals (panel 3 is the distribution); two rounds instead of three, L0 instead of L2; scoreboard and encyclopedia not updated. Conductor control plane unreachable.
