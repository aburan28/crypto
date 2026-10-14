#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)
#show figure.caption: it => it.body
#show table: set text(size: 8.5pt)

= The $k = 4$ Macaulay solve, 6× faster in wall time, the same arithmetic

2026-10-13. Commit `029309579`, `src/cryptanalysis/gaudry_quartic.rs`; study `research/gaudry_quartic_speed_20261013/`. *Class: engineering* — the counted cost $C_4$ is unchanged, so the crossover with rho does not move.

== Change

Characteristic polynomial and eigenvectors moved from `u128` remainders, one entry at a time, to `u16` rows with the vectorised `axpy`, rows in parallel (the Hessenberg step's elementary similarities commute, so its row operations precede its column operations with the same result). Normal forms resolved in parallel blocks of 128 columns. `axpy` compiled for AVX2/AVX-512 and chosen at run time. Same operations, counts and outputs: identical on the committed $d = 3 dots 5$ panel, pinned by a unit test, equal in `total_muls` on all 72 solves here, and $1.2269 dot 10^(12)$ on the $S_5$ solve, §11.17's $C_4$.

== Result

Interleaved before/after, generic planted systems at $p = 269$, 12 solves per arm and degree, median wall time in ms.

#table(columns: 7, align: right,
  [$d$], [dim], [elimination], [normal forms], [char. poly.], [eigenvectors], [whole solve],
  [4], [256], [155 → 128], [91 → 15], [70 → 17], [33 → 14], [354 → 184 (*1.9×*)],
  [5], [625], [1,074 → 690], [1,456 → 179], [1,192 → 93], [608 → 73], [4,368 → 1,067 (*4.1×*)],
  [6], [1,296], [7,116 → 3,682], [13,961 → 1,412], [12,342 → 507], [3,364 → 169], [37,433 → 5,862 (*6.4×*)],
)

#figure(image("figures/phase_wall.svg", width: 100%), caption: [Median wall time of one full solve by phase. $S_5$ ($d = 8$, dimension 4,096, three systems) after only: 164–203 s a solve; before not timed here (about 30 min a solve in earlier runs on this host).])

The gain grows with size; the elimination is now about 64 % of the $S_5$ solve and the next lever.

Environment: Claude Code cloud container, Intel Xeon \@ 2.10 GHz, 4 vCPU with AVX-512, 15 GiB RAM, rustc 1.97.0; shared host, otherwise idle; four rayon threads. Conductor control plane unreachable: no task card.
