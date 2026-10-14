#set page(paper: "a4", margin: 1.7cm, numbering: "1")
#set text(size: 9.5pt)
#show table: set text(size: 8.3pt)

#align(center)[
  #text(size: 14pt, weight: "bold")[lit2026 W10: Coppersmith reach for the prime-field $S_5$ small-root oracle] \
  2026-10-11 · branch `lit2026/w10-s5-coppersmith` · engine `examples/lit2026_w10_coppersmith_bound.rs` @ `a21b5de50` · run `run-20261011T2118Z`
]

*Question.* Index calculus over $FF_p$ first beats rho at arity $m = 4$. That needs roots of
$f = S_5(x_1, dots, x_4, z) equiv 0 mod p$ with $|x_i| < p^delta$, $delta approx 1/5$ (repo
`RESEARCH_PRIME_FIELD_BREAKTHROUGH_PROGRAM.md` §2.4). Do the 2024–26 automated-Coppersmith tools
(Ryan ePrint 2024/1577 "cuso"; Ding et al. 2026/1027, an LP for Meers–Nowakowski asymptotics;
Ding et al. 2026/423) reach that? Only the abstracts were read: the 2026/1027 PDF returned HTTP
403 (bot protection, not bypassed). The optimiser below is our own exact formulation.

*Problem.* $f$ has degree $d = 2^(m-1)$ per variable. Measured support with $z$ substituted
(61-bit field): $S_3$ 9/9, $S_4$ 125/125, $S_5$ 6561/6561 monomials of the box $[0,d]^m$, all
of them present. Total degree $m d$ = 4, 12, 32. Checked against genuine relations and planted
decompositions.

*Family $cal(F)$.* Lattices spanned by $p^(ell-k(beta)) x^(beta - k(beta) a) f^(k(beta))$,
$beta in M$, with every support inside $M$; $a$ is a vertex of $"Newt"(f)$. As $p -> oo$ the
lattice succeeds iff $sum_M k(beta) > sum_M lambda dot beta$. This family contains Jochemsz–May
basic and extended shifts and every monomial-set choice for one polynomial over $(p,f)^ell$.

*Exact optimiser.* For fixed $(a, ell)$ and a bounding box, the choice of $(M, k)$ is a
maximum-weight closure (totally unimodular LP). We solve it by min-cut inside Dinkelbach
iteration and get an exact rational.

*Ceiling.* Charge the $j$-th unit of $k(beta)$ to $beta + j(c - a)$. This gives
$sum k <= sum_M (lambda dot u)/(lambda dot c)$, so success needs $lambda dot c < 1$ for the top
vertex $c$, i.e. $delta < 1/(m d)$ at every $ell$. The bound holds for the $e$-encoding as well,
whose weighted degree is the same $m d$.

#table(columns: 4, align: (left, right, right, right),
  [], [$S_3$ ($m=2$)], [$S_4$ ($m=3$)], [$S_5$ ($m=4$)],
  [$delta^*$ in $cal(F)$ (box limit $2\/(m(m+1)d)$)], [*1/6*], [*1/24*], [*1/80 = 0.0125*],
  [ceiling of $cal(F)$, every $ell$: $1\/(m d)$], [1/4], [1/12], [1/32 = 0.03125],
  [exact closure optimum (largest $ell$)], [0.154789 ($ell=12$)], [0.032955 ($ell=5$)], [0.006334 ($ell=3$)],
  [closure, $e$-encoding], [0.136111 ($ell=3$)], [0.024903 ($ell=2$)], [0.004061 ($ell=1$)],
  [$delta$ needed (repo §2.4)], [1/2], [1/4 (tie)], [*1/5*],
  [repo / this run, finite lattices], [HG ≤ 0.10, heuristic to 0.14–0.16], [ℓ=1, dim 125: 0/36 at B = 2…15], [—],
)

*Validation.* Univariate: the closure optimum equals $ell\/(d ell + d - 1)$ in 30/30 cells and
tends to $1\/d$. Boneh–Durfee: the closure optimum equals the $x$+$y$-shift lattice with the best
$t$ at every $ell <= 32$, and that lattice reaches 0.28466 at $ell = 2000$ against the literature
value 0.28475. The 0.2929 bound lies outside $cal(F)$ (it needs a sublattice). JM basic gives
1/18, 1/10, 0.119 for $S_3$, the repo's column. In 70/70 closure cells
($S_3$ $ell <= 12$, $S_4$ $ell <= 5$, $S_5$ $ell <= 3$, every vertex class) the optimum over all
monomial sets is a box $[0,L']^m$.

*Verdict.* $delta^*(S_5) = 1/80 < 1/5$, and the ceiling of $cal(F)$ is $1/32 < 1/5$. *The
Coppersmith route to the $m = 4$ oracle is closed for full-rank triangular single-polynomial
lattices over $(p,f)^ell$*: JM basic and extended shifts, all monomial sets, all leading
vertices, in both the $x$-encoding and the $e$-encoding. *Not covered:* non-full-rank
sublattices and relation-based linearisation (the Boneh–Durfee 0.292 / Herrmann–May type);
cuso's Gröbner-optimal combinations that cancel out-of-$M$ monomials; multi-polynomial systems
(curve equations, shared unknowns); unequal or non-interval bounds. *Correction to the repo:*
the "measured 1/80" in `RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md` is the derived JM-basic
limit. Only $S_3$ was ever measured.

#figure(image("reports/run-20261011T2118Z/delta_vs_arity.svg", width: 92%), caption: [δ against arity. Source `reports/run-20261011T2118Z/run.json`.])
#figure(image("reports/run-20261011T2118Z/delta_vs_multiplicity.svg", width: 100%), caption: [Convergence in multiplicity ℓ. Source `run.json`.])
#figure(image("reports/run-20261011T2118Z/finite_lattice.svg", width: 80%), caption: [Finite lattices, 32-bit curve `a3-fast-32bit-p256class`. Source `run.json`.])
