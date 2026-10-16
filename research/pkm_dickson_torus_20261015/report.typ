#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 8pt)

#align(center)[
  #text(size: 14pt, weight: "bold")[The Dickson torus tower, and a deeper-2-adic Kummer control]\
  research/pkm_dickson_torus_20261015 · 2026-10-15 · lane dickson-torus\
  Protocol committed before any measured cell: `9776493c4` (PROTOCOL.md sha256 `43e0071a…fb39fb`)
]

*What this is.* A stage diagnostic of the solver axis (AGENTS.md §2). It measures the F4 solving degree
$D$ of PKM tower systems with `f4_fp_tower`. It prices no pipeline, reports no $S$ and claims no speed.
Labels: _measured_ is a committed row; _derived_ is a calculation or argument; _recorded_ is a value from an
earlier round (PKM note §10–§16). The source of truth is README.md in the study directory.

= Answers

#table(columns: (auto, 1fr, auto), align: left,
  [*question*], [*answer*], [*status*],
  [Q-A: does the torus Dickson tower ($2^t | p+1$) give the split $D$?],
  [*H_same.* At both torus primes, torus $D$ = split $D$ in all 16 cells: $m=2$, $N=2..20$; $m=3$, $N=6,9,12,15$ (one target per prime at 15); $m=4$, $N=8,12$. Width, step count and (at the 31-bit primes) multiply-adds also equal the split Dickson's, so it is the same F4 computation on different numbers.], [measured],
  [Q-A∞: degenerate at infinity?], [*Yes.* The equations are the split Dickson ones with another constant $c$, so the top parts are $y_(i,j)^2$. §3.4 gives a BFS degree of regularity of $N+1$.], [derived],
  [Q-B: does deeper 2-adicity change Kummer $D$?], [*H2_same at $v_2 = 30$* ($p_2 = 3 dot 2^30+1$). $D$ and width are equal to $p_0$ ($v_2=18$) and $p_1$ ($v_2=27$) at every $(m,N)$, and so are the multiply-adds of $p_1$ on this build. *$v_2 >= 40$ was not run.* Both engines require $p < 2^32$ (32-bit coefficients), and no prime below $2^32$ has $v_2(p-1) > 30$.], [measured / blocked],
  [Q-C: $m=3, N=21$ and $m=4, N=20$?], [About 130–150 GB and $2.5$–$7.4 dot 10^15$ multiply-adds; about 240–380 GB and $0.85$–$2 dot 10^16$. They need a 256 GB or 512 GB host, or a disk-backed $B'$. Not run.], [derived],
)

= What ran

Shared 4-core, 15 GB container; each run single-threaded, `nice -n 10`, 2.5 GB address-space cap. Two random
targets per cell; seed `0x504B4D54`; flags as round 2 (`--engine tower --budget 7200 --max-nnz 2e9`).
Primes: torus $T_24 = 7 dot 2^21 - 1$ (14,680,063) and $T_31 = 61 dot 2^25 - 1$ (2,046,820,351); split Dickson at
$p_1 = 15 dot 2^27+1$; deep Kummer at $p_2 = 3 dot 2^30+1$; positive control (isogeny) at $p_0 = 3 dot 2^18+1$.
The torus builder checks that $zeta$ has exact order $2^t$ in $T$, that $u in T without F_p$ and $u^2 in.not mu_(2^t)$,
that $V subset F_p$ (zero $i$-part of $w + w^(-1)$), $|V| = 2^t$, that $y^2-4$ is a non-residue on $V$, that no level
is ramified, and that $c != plus.minus 2$.

- *Positive control (R0): passed.* Isogeny $D = 5$ at $N = 4..10$, where the split towers give 4.
- *Exhaustive check:* `pkm_tower_check verify`, 127 finished tower rows, 127 agree.
- *Run status:* no system timed out or ran out of memory. The session's 2-hour task limit stopped the run script
  during R14 before it wrote a row. R14 was restarted whole and finished. R15 (the second targets at $m=3$,
  $N=15$) was not run, because the protocol deadline had passed.
- *Build provenance:* R0–R13 used a binary built in the lanes' shared target directory. A rebuild in this lane's
  own target directory re-measured a subset of every (kind, $m$, control) series, 144 rows. Every field except
  wall time was identical (`compare-rows`, 0 differences).

= Task A: torus against split, measured ($D$ / width / steps)

#table(columns: 6, align: center,
  [$m$, $N$], [torus $T_31$], [torus $T_24$], [split Dickson $p_1$], [split recorded (Kummer)], [torus null $T_31$],
  [2, 2], [3/4/2], [3/4/2], [3/4/2], [3], [3/4/2],
  [2, 4], [4/16/4], [4/16/4], [4/16/4], [4], [4/16/4],
  [2, 6], [4/56/6], [4/56/6], [4/56/6], [4], [4/56/6],
  [2, 8], [4/152/8], [4/152/8], [4/152/8], [4], [4/152/8],
  [2, 10], [4/336/10], [4/336/10], [4/336/10], [4], [4/336/10],
  [2, 12], [5/1223/13], [5/1223/13], [5/1223/13], [5 (1231)], [5/1223/13],
  [2, 14], [5/2303/15], [5/2303/15], [5/2303/15], [5 (2304)], [5/2303/15],
  [2, 16], [5/4786/17], [5/4786/17], [5/4786/17], [5 (4789)], [5/4786/17],
  [2, 18], [5/8639/20], [5/8639/20], [5/8639/20], [5 (8644)], [5/8639/20],
  [2, 20], [6/33426/24], [6/33426/24], [6/33426/24], [6 (33446)], [6/33426/24],
  [3, 6], [5/164/10], [5/164/10], [5/164/10], [5], [5/164/10],
  [3, 9], [6/1184/16], [6/1184/16], [6/1184/16], [6 (1186)], [6/1184/16],
  [3, 12], [6/5018/18], [6/5018/18], [6/5018/18], [6 (5019)], [6/5066/16],
  [3, 15], [7/31378/23 (1)], [7/31378/23 (1)], [7/31378/23 (1)], [7 (31397)], [not run],
  [4, 8], [6/1080/21], [6/1080/21], [6/1080/21], [6 (1049)], [6/952/20],
  [4, 12], [7/9451/29], [7/9451/29], [7/9451/29], [7 (9452)], [7/8944/33],
)
The $T_24$ null equals the $T_31$ null in every cell. (1) marks one target. Under rule A1 the torus $m=2$ line
(3 4 4 4 4 5 5 5 5 6) reads inconclusive, as Kummer's does. The shape-matched null gives the torus $D$ everywhere
(H1b: no curve-specific effect).

= Task B: Kummer at $v_2(p-1) = 30$, measured

$p_2$, $m=2$, $N = 2..18$: $D$ = 3 4 4 4 4 5 5 5 5, widths 4, 16, 56, 152, 337, 1231, 2304, 4789, 8644. These equal the
recorded Kummer values at $p_0$ and $p_1$. $m=3$, $N=6,9,12$: $D$ = 5, 6, 6, widths 164, 1186, 5019, also equal. The null
equals the tower in every cell. _Derived:_ the Kummer system's coefficients ($g$, $a$, $b$, $x_R$, and the tower
constants 1) do not involve $v_2(p-1)$ beyond the condition $2^t | p-1$, so 2-adicity can act only through the
prime's size.

= Task C (derived, not run)

#table(columns: 6, align: center,
  [cell], [ratio basis], [width], [peak memory], [multiply-adds], [time (recorded host)],
  [$m=3$, $N=21$], [last step], [≈ 360k], [≈ 146 GB], [$2.5 dot 10^15$], [≈ 86 h],
  [$m=3$, $N=21$], [geometric mean], [≈ 475k], [≈ 133 GB], [$7.4 dot 10^15$], [≈ 138 h],
  [$m=4$, $N=20$], [last step], [≈ 625k], [≈ 384 GB], [$8.5 dot 10^15$], [≈ 327 h],
  [$m=4$, $N=20$], [geometric mean], [≈ 658k], [≈ 243 GB], [$2.0 dot 10^16$], [≈ 168 h],
)
The inputs are the recorded Kummer rows at $p_1$ (requirements README §3.2). The projections assume $D$ does not
rise again. If it does, they are too low. The round-6 early exit removes up to 39% of the multiply-adds and does
not change memory. The dense reducer block $B'$ holds about 93% of peak memory, so a disk-backed $B'$ is the engine
lever. Without one, $m = 3$, $N = 21$ needs about a 256 GB host and $m = 4$, $N = 20$ about a 512 GB host.

#figure(image("d_vs_n.svg", width: 100%),
  caption: [$D$ against $N$: torus ($T_31$, $T_24$), split Dickson ($p_1$, this build), torus null, deep-2-adic Kummer
  ($p_2$) and the isogeny positive control. The grey line is the recorded split value. Points are offset horizontally;
  two random targets per cell, one per prime at $m = 3$, $N = 15$. Runs R0–R14.])

#figure(image("width_vs_n.svg", width: 70%),
  caption: [Widest F4 matrix up to the solution (columns, log scale), same runs.])

= Scope

All primes are below $2^32$ and $t <= 10$, so nothing here is a statement about P-192, P-224 or P-256. Equal $D$
means the torus tower inherits every open question of the split towers (lead 1). It is not evidence that $D$ is
bounded. The torus $m = 2$ line rises to 6 at $N = 20$, exactly as Kummer's does.
