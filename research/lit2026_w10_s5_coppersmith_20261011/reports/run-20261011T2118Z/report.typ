#set page(paper: "a4", margin: 1.6cm, numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 8pt)

= Run report: lit2026 W10, run-20261011T2118Z

*What ran.* `examples/lit2026_w10_coppersmith_bound.rs` at commit `a21b5de50` (branch
`lit2026/w10-s5-coppersmith`, stacked on the program commit `3c239129d`), built with
`CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/lit2026-w10 cargo build --release --example lit2026_w10_coppersmith_bound`,
run as `./lit2026_w10_coppersmith_bound --part all --max-ell 3 --json run.json` (wrapper `run_w10.sh`).
Start 2026-10-11 21:18 UTC, end 22:01 UTC; 2564 s wall, 301 s user, peak RSS 1.45 GB
(`time_and_exit.txt`). Host: shared Apple M4 Pro, load average 108 at start, 177 at end (`host_load.txt`).
Exit status 0. No run failed or timed out.

*Isolation.* None claimed. Every reported quantity except the `ms` columns is a deterministic
count or an exact rational (closure optimum, monomial counts, LLL vanishing counts with fixed seeds).
Wall times are L0, single sample, informational only.

*Raw records.* `run.md` (all tables), `run.json` (schema `lit2026-w10-coppersmith/v1`, 181 records).
Graphs are drawn from `run.json` by `../../plot_w10.py` (layout only).

== Parts and what each measured

+ *Known bounds* (part `known`). Univariate degree $d in {2,3,5}$, $ell <= 16$: exact closure
  optimum equals $ell/(d ell + d - 1)$ in all 30 cells and tends to $1/d$. Boneh–Durfee
  $f = 1 + x(A+y) mod e$, $X = e^delta$, $Y = e^(1/2)$: the exact closure optimum equals the
  Boneh–Durfee $x$+$y$-shift lattice with the best $t$ at every $ell in {1,2,3,4,6,8,12,16,24,32}$
  (0.27949 at $ell = 32$); the closed form at $ell = 2000$, $tau = 0.215$ gives 0.28466 against
  the literature value $(7 - 2 sqrt(7))/6 = 0.28475$. Jochemsz–May basic lattice on the box
  $[0,d]^m$: 1/18, 1/10, 0.11905 at $ell = 1, 2, 3$ for $S_3$ (the repo's probe-B column) and
  limits $1/6, 1/24, 1/80$.
+ *Supports* (part `support`, 61-bit prime, random $a, b$, random target $z$). $S_3, S_4, S_5$
  with $x_(m+1) = z$ have 9/9, 125/125, 6561/6561 monomials of the box $[0, 2^(m-1)]^m$, all box
  vertices, total degree 4, 12, 32; symmetric; the interpolant equals the direct resultant at 8/8
  off-grid points, Semaev vanishes at 8/8 random relations and the interpolant at 4/4 planted
  decompositions of $z$. Elementary-symmetric rewriting: 6, 35, 495 monomials, weighted degree 4, 12, 32.
+ *Box families* (part `families`): JM basic and best-extension box for $ell <= 128$.
+ *Closure optimum* (part `closure`): exact maximum-weight-closure optimum over all monomial sets
  in $[0,L]^m$ for each leading vertex class.
+ *Finite lattices* (part `lattice`): 32-bit curve `a3-fast-32bit-p256class`
  ($p = 4294967291$, $a = -3$, $b = 31$, prime order $n = 4294959973$), repo `lll_reduce`
  ($<= 30$ dims) / `lll_reduce_hp`, Lovász 0.99.

== Closure optimum (selected rows; all rows in `run.md` §4)

#table(columns: 8, align: (left, right, right, right, right, right, left, right),
  [poly], [$ell$], [vertex classes], [$L$], [$delta_ell$ exact], [$delta_ell$], [closure = best box $<= L$], [dim],
  [$S_3$], [1], [0,1,2], [4], [9/100], [0.090000], [yes], [25],
  [$S_3$], [12], [0,1,2], [48], [202/1305], [0.154789], [yes (L′=29)], [900],
  [$S_4$], [1], [0..3], [18], [1331/70875], [0.018780], [yes (L′=14)], [3375],
  [$S_4$], [5], [0..3], [41], [25090/761349], [0.032955], [yes (L′=37)], [54872],
  [$S_5$], [1], [0,4], [44], [923521/175821516], [0.005253], [yes (L′=38)], [2313441],
  [$S_5$], [1], [0..4], [24], [83521/18750000], [0.004454], [yes (bounded)], [390625],
  [$S_5$], [2], [0,4], [32], [237073/37949472], [0.006247], [yes (bounded)], [1185921],
  [$S_5$], [2], [1..4], [24], [45041/9375000], [0.004804], [yes (bounded)], [390625],
  [$S_5$], [3], [4], [32], [480707/75898944], [0.006334], [yes (bounded)], [1185921],
)
In every one of the 70 closure cells the optimum equals the best box $[0,L']^m$, $L' <= L$, and
every vertex class gives the same value. The charging check $Sigma k <= Sigma_(u in M) J(u)$ held
in all cells (with equality). "Bounded" means the best box is the bounding box itself.

Symmetric encoding $G(e_1, dots, e_m)$ with weights $i delta$ (§4b): $S_3$ 0.100000, 0.125000,
0.136111 at $ell = 1,2,3$; $S_4$ 0.018176, 0.024903 at $ell = 1, 2$; $S_5$ 0.004061 at $ell = 1$
(bounding box $[0, 2 dot "max deg"]$; best vertex the constant term).

== Finite lattices

#table(columns: 7, align: (left, right, right, right, right, right, right),
  [poly], [$ell$], [dim], [$delta$ (B)], [$>= m$ rows vanish], [HG first $m$ rows], [mean vanishing],
  [$S_3$], [1], [9], [0.10 (9)], [6/6], [3/6], [2.2],
  [$S_3$], [1], [9], [0.14 (22)], [3/8], [0/8], [1.5],
  [$S_3$], [1], [9], [0.16 (34)], [0/7], [0/7], [0.4],
  [$S_3$], [2], [25], [0.12 (14)], [7/7], [1/7], [12.1],
  [$S_3$], [2], [25], [0.16 (34)], [4/8], [0/8], [1.4],
  [$S_3$], [2], [25], [0.20 (84)], [0/8], [0/8], [0.1],
  [$S_4$], [1], [125], [0.035 (2) … 0.125 (15)], [0/36], [0/36], [0.0–0.7],
)
$S_4$ at $ell = 1$ recovered no planted root at any tested bound, including the smallest
($B = 2$); its asymptotic-in-$p$ reach at $ell = 1$ is $delta_1 = 0.00133$, i.e. $X < 1.03$ at
32 bits. Trials per cell: 3–8 for $S_3$ (a trial is dropped, not counted, when the planted sum is the point at infinity — frequent at $B <= 3$ — or LLL returns an error), 6 for $S_4$.

#figure(image("delta_vs_arity.svg", width: 100%), caption: [δ against arity: upper bound, JM limit, closure optimum, needed δ, repo S3 probes. Source `run.json`.])
#figure(image("delta_vs_multiplicity.svg", width: 100%), caption: [δ_ℓ against multiplicity ℓ. Source `run.json`.])
#figure(image("finite_lattice.svg", width: 85%), caption: [Finite lattices on the 32-bit curve. Source `run.json`.])
