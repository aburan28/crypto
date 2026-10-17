#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)
#set par(leading: 0.64em)

= Torsion-coordinate DLP over quadratic extensions

*Eight-arm same-target comparison with a public scaling-isomorphism control · 2026-10-14*

The frozen groups are prime-order subgroups of $E(F_(p^2))$, with $E:y^2=x^3+1$ and orders 19, 61, 211, 811. Translation by rational three-torsion acts on $chi=(y+1)/(y-1)$. The direct degree-three image is $E':y^2=x^3-27$. In every exact field tested, $u^2=-3$ has a base-field solution, so $psi_u(x,y)=(u^2 x,u^3 y)$ is also a public map to $E'$. The selected values of $u$ were 27, 57, 28 and 264, each checked on 64 point pairs plus all frozen inputs.

The table uses the *whole cold current-thread CPU interval*, divided by measured source-addition time and $sqrt(r)$, denoted $S_("CPU/GAE")$. It charges setup, queries, relation recovery, map work and scalar verification. Every value is a three-round median on a shared Darwin arm64 host (L0); the χ/rho and scaling/rho columns are medians of paired full-cost ratios.

#table(
  columns: (auto, auto, auto, auto, auto, auto, auto),
  table.header([Order $r$], [Projected], [$chi$ index], [Degree-three], [Scaling], [Rho], [$chi$/rho]),
  [19], [202.85], [256.79], [239.44], [231.65], [27.30], [14.53],
  [61], [138.26], [165.85], [156.96], [252.68], [15.28], [11.09],
  [211], [157.59], [190.15], [156.02], [154.41], [9.70], [21.46],
  [811], [198.17], [277.34], [200.41], [204.56], [6.91], [50.27],
)

The four projected arms had *identical relation count, query count, rank and recovered scalar in all twelve paired size-by-round cells*. The 96 run rows have 87 scalar recoveries checked by a second native $F_(p^2)$ and curve-arithmetic implementation; nine raw-pair rows reached their cap. The $chi$ *orbit index* increased median CPU over the plain projected table at every size. This result says nothing about the still-unimplemented $chi$ polynomial decomposition oracle. The three-round ranges, including a wide scaling-arm spread at order 61, prevent map-specific wall-performance attribution.

#pagebreak()

== Full arm ranges

#image("GRAPH_V4.svg", width: 100%)

Dots are medians; whiskers span the three measured rounds. Hollow raw-pair points are capped work with no scalar answer. The dashed generic counted-operation floor is a reference marker; the CPU-equivalent axis does not make it a universal lower bound in CPU time.

#pagebreak()

== Accounting, proof and open work

The cold interval includes line lifting, square roots, cofactor projection, factor-base construction, pair-table allocation, failed and successful queries, exact source relation checks, modular elimination, target recovery and internal verification. The degree-three and scaling arms charge their map work; scaling constructs $u$ by an ascending search in measured setup. Six-class rho charges eigenvalue discovery. Calibration and output serialization are outside method intervals. A fixed two-class generic counted-operation floor is $sqrt(pi/4)=0.886226925$ per $sqrt(r)$; converting CPU to local addition equivalents does not create a universal lower bound. Per-arm peak RSS, portable primitive conversion, L2 isolation and independent-host timing remain open.

The exact raw torsion-shift census found two extra subgroup support points at $p=787$, but no extra cofactor-projected pair sum and no representation of that frozen target. The projected support sizes were 18/18, 58/60, 196/210 and 742/810 nonidentity subgroup points. These are finite observations on four frozen bases. Four exact coefficient witnesses are proven in a pinned Lean 4 build with no axioms. The universal sixth-power factorization is also proven in Lean using its standard foundations and checked by a separate native identity certificate. The prior conditional quotient-support theorem and polynomial χ certificates are linked from the report. A faithful curve-level Lean isomorphism and χ theorem is still open.

For this fixed source and destination, the scaling isomorphism exists over every $F_(p^2)$ with odd $p != 3$: each nonzero $a in F_p$ has order dividing $p-1$, which divides $(p^2-1)/2$, so $a$ is a square in $F_(p^2)$. Taking $a=-3$ gives the required $u$. The next non-isomorphic search must therefore use a different curve/isogeny pair. The finite-field existence theorem is written in the Markdown report but is not yet formalized in Lean.

The embedding degrees $"ord"_r(p)$ were 9, 10, 105 and 810, or 9, 5, 105 and 405 relative to $p^2$. The root-of-unity output fields have those first four degrees over $F_p$, while common fields containing the frozen $F_(p^2)$ points and the roots have degrees 18, 10, 210 and 810. A second pairing input may require more. A pairing transfer and multiplicative extension-field DLP were screened by those degrees only: they were neither constructed nor priced. The table instead compares complete elliptic-curve DLP methods on the same extension-field groups.

The corrected seven-arm V3 run and its 75-scalar replay remain in `REPORT_V3.md` and `REPORT_V3.pdf`. The final eight-arm bytes are in `RUN_V4.jsonl`; `REPLAY_RUN_V4.json` contains the independent 87-scalar check, `PAIRED_AUDIT_V4.jsonl` the twelve exact outcome comparisons, `SUMMARY_V4.jsonl` all 32 table rows, and `ARTIFACTS.sha256` the full manifest. The pre-merge release-library suite passed 4,009 tests with zero failures and 113 ignored; the merged PR tree passed 4,083 with zero failures and 117 ignored. Earlier failed environment receipts are retained. The Markdown `REPORT.md` gives the complete cost table, failures, proof scope and next gate.
