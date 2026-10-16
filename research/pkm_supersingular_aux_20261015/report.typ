#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 10pt)
#set heading(numbering: "1.")
#show table: set text(size: 8.5pt)

#align(center)[
  #text(size: 15pt, weight: "bold")[A supersingular auxiliary curve for PKM's isogeny tower]
  #v(2pt)
  research/pkm_supersingular_aux_20261015 · 2026-10-15 · branch `claude/pkm-supersingular-aux` \
  Protocol `PROTOCOL.md` frozen in `8d719fddf` (sha256 `db109df9…c1ad65`) before any measured cell
]

= Question and verdicts

Over the NIST P-192 prime ($p equiv 3 mod 4$, $v_2(p+1) = 64$) the supersingular curve
$y^2 = x^3 + x$ has $p + 1$ points and a cyclic 2-Sylow subgroup, so a point of order $2^48$
costs one draw, against about 11.2 core-years for the ordinary curve the backlog card
`XP-20261015-p192-isogeny-aux-search-t48` prices. Does PKM's isogeny-tower factor base work as
well on it? Three pre-registered questions:

#table(
  columns: (auto, 1fr, auto),
  [*question*], [*verdict*], [*class*],
  [Q1: is $V_"ss" = x(R' + chevron.l Q chevron.r)$ the Dickson torus set up to a Möbius map?],
  [*No; by the rule's letter, partial.* At $t = 3 dots 5$ over 103 primes, 0 of 6,507 sets of $y^2 = x^3 + x$, 0 of 6,768 $j = 0$ sets and 0 of 24,580 ordinary sets match a torus set; 2 of 12,969 sets of $y^2 = x^3 - x$ do, at one cell ($p = 383$, $t = 3$). At $t = 2$ every family matches about 75 % (one modulus).],
  [exact],
  [Q2: F4 solving degree $D$, supersingular vs ordinary vs torus],
  [*H_same fails; supersingular = ordinary.* All four supersingular curves give the ordinary isogeny
   tower's $D$ in every cell ($m = 2$, $N <= 20$; $m = 3, 4$, $N <= 12$; three primes); the torus differs
   (4 against 5 at $m = 2$, $N = 4 dots 10$; 5 against 6 at $m = 3, N = 6$).],
  [measured],
  [Q3: relation phase, charged (ecbench)],
  [*The search cost goes.* Aux share 0.9–1.9 % of $S$ (build only) against 10–27 % for the searched
   ordinary curve; yield unchanged. Under the AGENTS.md §5a band (unpriced work at 0, 0.1, 1 GAE per
   operation): ss < ordinary on s22, s24 at every weight (0.82–0.92×), ordinary < ss on s20 (larger
   base); ss vs a plain base of its size *unresolved* on s20, s24 and dearer by 1–7 % on s22.],
  [measured],
)

*For the card (derived).* Nothing measured distinguishes a supersingular $E'$ from an ordinary one
except the cost of finding it, so the ordinary search is not needed to feed the relation phase. The
gate of any P-192 use, lead 1 (bounded $D$ as $N$ grows), is untouched; nothing here concerns P-192's
security.

= Q1: the structure of $V$ (exact)

For every prime $p equiv 3 mod 4$ in $[100, 3000]$ and $t = 2 dots 5$ with $t + 1 <= v_2(p+1)$, every
tower set of every family (all Dickson torus cosets; every $x(R' + chevron.l Q chevron.r)$ on
$y^2 = x^3 plus.minus x$, and on $y^2 = x^3 + 1$, $y^2 = x^3 + d^3$ when $p equiv 2 mod 3$; three
ordinary curves; the positive control $y^2 = x^3 + 16 x tilde.equiv y^2 = x^3 + x$) is reduced to a
complete $"PGL"_2(FF_p)$ invariant and compared.

#figure(image("figures/q1_equivalence.svg", width: 100%),
  caption: [Fraction of each family's sets equivalent to some torus set (positive control: to some
  $y^2 = x^3 + x$ set), pooled over primes; labels are counts. Source `q1/pooled.tsv`.])

#table(
  columns: 7,
  [$t$], [cells], [$x^3 + x$], [$x^3 - x$], [$j = 0$ (both)], [ordinary], [pos. control],
  [2], [103], [13965 / 18617], [27883 / 37181], [14407 / 19276], [54738 / 72689], [18617 / 18617],
  [3], [53], [0 / 4993], [*2* / 9959], [0 / 5186], [0 / 18763], [4993 / 4993],
  [4], [27], [0 / 1256], [0 / 2500], [0 / 1300], [0 / 4944], [1256 / 1256],
  [5], [12], [0 / 258], [0 / 510], [0 / 282], [0 / 873], [258 / 258],
)
Sets equivalent to a torus set / all sets (positive control: to a $y^2 = x^3 + x$ set). The canonical
form agreed with brute force over $"PGL"_2(FF_p)$ on all 906 checked pairs ($p <= 131$). The one
exception at $t >= 3$ is a single shared class at $p = 383$; the pre-registered rule therefore reads
*partial*, but no family-wide map exists on the data, and none for the P-192 curve $y^2 = x^3 + x$.

*Why.* The torus set is a fibre of the polynomial $D_(2^t)$, totally ramified over $infinity$; the
supersingular set is a fibre of a composed Vélu x-map, whose ramification indices are at most 2.
The maps are not Möbius-conjugate; single fibres can still match by coincidence, which happens only
at $t = 2$ (one cross-ratio modulus) and for every pair of families alike.

= Q2: F4 solving degree (measured)

`f4_fp_tower` (reduced presentation, full-rank exit, degree bound $n + d + 6$), two random targets
per cell, one thread; 336 systems, none stopped; 296 tower rows ($m <= 3$) agree with an exhaustive
recount over $V^m$. Primes $T_24 = 7 dot 2^21 - 1$, $T_31 = 61 dot 2^25 - 1$, $S_25 = 9 dot 2^21 - 1$,
and the positive control at $p_0 = 3 dot 2^18 + 1$.

#figure(image("figures/q2_d_vs_n.svg", width: 100%),
  caption: [$D$ against $N$ at $m = 2$; both targets equal in every cell. Source `q2/cells.tsv`.])

#table(
  columns: 7,
  [$N$], [ss-j1728], [ss twists, j=0], [ordinary], [torus], [null on ss], [$p_0$: ordinary / split],
  [2], [3/4/2], [3/4/2], [3/4/2], [3/4/2], [3/4/2], [3 / 3],
  [4], [5/16/4], [5/16/4], [5/16/4], [*4* /16/4], [5/16/4], [5 / 4],
  [6], [5/57/7], [5/57/7], [5/57/7], [*4* /56/6], [5/57/7], [5 / 4],
  [8], [5/150/9], [5/150/9], [5/150/9], [*4* /152/8], [5/150/9], [5 / 4],
  [10], [5/346/12], [5/346/12], [5/346/12], [*4* /336/10], [5/346/11], [5 / 4],
  [12], [5/1182/16], [5/1181/16], [5/1181/16], [5/1223/13], [5/1182/15], [],
  [14], [5/2629/19], [5/2629/19], [5/2629/19], [5/2303/15], [5/2628/18], [],
  [16], [5/5034/22], [5/5034/22], [5/5034/22], [5/4786/17], [5/5034/21], [],
  [18], [5/9220/24], [5/9221/24], [5/9221/24], [5/8639/20], [5/9184/23], [],
  [20], [6/35157/28], [—], [6/35158/28], [6/33426/24], [—], [],
)
Entries $D$ / width / steps ($m = 2$), identical at every prime where the arm ran. At $m = 3$:
ss 6, 6, 6 and ordinary 6, 6, 6 against torus 5, 6, 6 ($N = 6, 9, 12$); at $m = 4$: all 6, 7
($N = 8, 12$). The torus rows reproduce round 8 (at $N = 20$ its multiply-add count to the unit).
*Why:* every Vélu tower has top forms $y_j y_(j+1) - y_j^2$; Dickson towers have $-y_j^2$.

= Q3: relation phase, charged (measured)

Three P-192-shape curves ($v_2(p - 1) = 1$, $v_2(p+1) = 10$, 20–24 bits), one ecbench session each,
7 arms × 8 workloads × 2 measured rounds; 504/504 runs verified, every replay identical, L0
(operation counts only).

#figure(image("figures/q3_s_and_ratios.svg", width: 100%),
  caption: [Left: mean $S$ per arm at $w = 0$ (16 runs), hatched = auxiliary search + build. Middle: paired ratios at $w = 0$
  with ecbench's 95 % cluster-bootstrap intervals. Right: the §5a band (addendum). Sources `q3/tables/`, `q3/band/`.])

#table(
  columns: 7,
  [curve], [ss ÷ plain(ss)], [ord ÷ plain(ord)], [ss ÷ ord], [plain(ss) ÷ plain(ord)], [aux share ss / ord], [rootfind ÷ mitm],
  [s20], [1.022 [1.010, 1.035]], [1.115 [1.098, 1.132]], [1.108 [1.092, 1.122]], [1.209], [1.9 % / 10.4 %], [1.256],
  [s22], [1.013 [1.006, 1.019]], [1.322 [1.311, 1.332]], [0.824 [0.820, 0.827]], [1.075], [1.4 % / 24.5 %], [1.138],
  [s24], [1.020 [0.992, 1.048]], [1.388 [1.354, 1.421]], [0.736 [0.715, 0.764]], [1.002], [0.9 % / 27.1 %], [1.074],
)

Rules at weight 0, as pre-registered: P1 held; H3a held; H3b held on s22, s24 and failed on s20; H3c held
on s22, s24 and failed on s20 (base 141 against 123 abscissae); H3d held (yield within 1.8 %); controls
passed.

*Addendum 2026-10-16, AGENTS.md §5a.* Every arm's unpriced counters carry counts (IC: row operations,
lookups, Legendre symbols, square roots, target-guard probes, Vélu evaluations; rho: canonicalisations;
`q3/band/unpriced.tsv`). `ecbench band` at $w = 0, 0.1, 1$ GAE per operation (ratios of mean priced $S$):

#table(
  columns: 6,
  [curve], [pair], [$w = 0$], [$w = 0.1$], [$w = 1$], [order across the band],
  [s20], [ss ÷ plain(ss)], [1.022], [1.015], [0.981], [*unresolved*],
  [s22], [ss ÷ plain(ss)], [1.013], [1.023], [1.074], [plain < ss],
  [s24], [ss ÷ plain(ss)], [1.020], [1.016], [0.999], [*unresolved*],
  [s20], [ss ÷ ordinary], [1.108], [1.102], [1.072], [ordinary < ss],
  [s22], [ss ÷ ordinary], [0.824], [0.838], [0.916], [ss < ordinary],
  [s24], [ss ÷ ordinary], [0.736], [0.762], [0.882], [ss < ordinary],
  [all], [rootfind ÷ mitm], [1.07–1.26], [1.07–1.25], [1.05–1.25], [mitm < rootfind],
  [all], [ss ÷ rho (lower bound)], [31.6–77.3], [33.0–77.5], [39.5–78.2], [rho < ss],
)
So under §5a: H3a, H3c (held s22/s24, failed s20), H3d and both controls stand at every weight; H3b is
unresolved on s20 and s24 and not held on s22. Every IC arm is 36–85× rho at $w = 1$ (29–94× at $w = 0$),
lower-bound comparisons on both sides.

= Scope

Small fields only ($p < 2^31$ for F4, $p < 2^24$ for ecbench); $N <= 20$; one seed per arm in Q3.
Q1 is exact on the primes run, not a proof for all $p$. The registry gained three curves; its derived
artifacts were regenerated on this branch in `b90b36670`. Two container restarts are recorded in `progress.txt`.
