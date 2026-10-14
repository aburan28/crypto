#set page(paper: "a4", margin: 1.7cm, numbering: "1")
#set text(size: 9.2pt)
#set heading(numbering: "1.")
#show table: set text(size: 7.8pt)

#align(center)[
  #text(size: 16pt, weight: "bold")[The ECBench ECDLP performance frontier] \
  #text(size: 10pt)[Factor bases, symmetry, F4/F5/F6, and root finding] \
  #text(size: 9pt)[Evidence synthesis and current-binary replay audit · 2026-10-13 · `run-20261013T2025Z`]
]

#v(5pt)
#box(fill: rgb("eef5f7"), inset: 10pt, radius: 3pt, width: 100%)[
  *Boundary.* Generic rho remains the cold one-target, setup-inclusive
  frontier. Compact-orbit index calculus owns the admitted retained-setup
  online point at n=61: rho/IC 13.9x, 95% interval [9.8, 19.1], 64 unseen
  targets on an isolated x86_64 host. The newer N83 dense-root point is
  18.528x online but exploratory; its complete IC process is 339.068 s
  against 1.173 s rho for that point. Solver-stage wins do not change the
  cold frontier.
]

= Scope and evidence class

This document reconciles committed measurements at repository commit
`759978250c1ce53d7b87bd8a1c2304781fcc4fc4`. It is not a new cross-method
session and creates no bound or challenge verdict. Machine-readable values,
source paths, derivations, and classifications are in `frontier-data.json`.

#table(
  columns: (1.2fr, 1.45fr, 2.4fr),
  align: (left, left, left),
  table.header([requested line], [status], [evidence boundary]),
  [factor bases], [complete for catalog; 7 measured arms], [10 sealed bases; four native sessions; operation smokes at L0],
  [symmetry], [partial], [signed Frobenius is native; torsion-symmetrised u-frame remains a legacy non-ECBench study],
  [F4/F5/F6], [committed cells covered], [K13/R13/K19/R19; matrix m19 and F6 m31 timeouts retained],
  [native root finder], [covered], [pair/subtraction ECBench engines; current root-finder rows replay exactly],
  [new N83 root finder], [not ECBench-integrated], [native independent study; 8 GiB table; shared-host online observation],
)

= Formal boundary

A performance point is keyed by curve/subgroup, exact targets, cold/online/stage
interval, counted unit and unpriced work, amortisation target count, memory and
resource envelope, and correctness/replay/isolation class. Ratios are formed
only when that key matches apart from method identity. A point is Pareto-optimal
only if no admitted point is no worse on every declared axis and better on one.

This separates four quantities that cannot be collapsed: setup-excluded online
wall and cold end-to-end cost; word XORs/probes/field operations/GAE; lower
bounds and complete counts; exploratory L0 wall and admitted isolated wall.

For q targets and reusable preparation B:

#align(center)[$ T_("IC")(q) = B_("IC") + q O_("IC"), quad T_rho(q) = B_rho + q O_rho $]

At N83, the reported reusable components sum to 339.004 s. With
`O_IC=0.063324125 s` and `O_rho=1.17329325 s`, the arithmetic crossover is
305.418, or 306 targets. This is a derived what-if, not a measured batch
frontier: it assumes fully reusable setup, constant online cost, no capacity
change, and zero rho setup.

#pagebreak()

= Canonical cold frontier

#table(
  columns: (1.2fr, 1.6fr, 1.45fr, 1.4fr),
  align: (left, left, left, left),
  table.header([axis], [incumbent], [measurement], [coverage]),
  [Koblitz medium ops], [`rho.signed_frobenius_strong` / `ECBND1h05a2ea918f3c`], [1.089x floor [0.927, 1.261]; alpha 0.451 [0.399, 0.512]], [144 runs, 6 sizes],
  [Koblitz medium memory], [`rho.signed_frobenius` / `ECBND1h46fbc9140ac3`], [1.188x floor [0.954, 1.464]; 0.000706 sqrt(r) memory], [72 runs, 3 sizes; constant bound],
  [retained-setup online], [compact-orbit IC n61 K500], [rho/IC 13.9x [9.8, 19.1]], [64 unseen points, isolated x86_64],
  [newer online observation], [dense root N83 K1798], [rho/IC 18.528377], [one point, shared host, exploratory],
)

#figure(
  image("cold_ops_gap.svg", width: 100%),
  caption: [Cold whole-pipeline operation gaps to each study's matched generic
  rho, retaining its reported estimator. Curves, references, and estimators
  differ, so parity—not bar-to-bar raw cost—is the shared comparison. Hatched
  bars retain unpriced work and are lower bounds.],
)

Every selected IC bar is above parity. The smallest selected gap is matrix
F4/F5 at 6.97x rho-negation on K13; compact orbit is 9.985x its matched
signed-Frobenius rho on the n17 factor-base cell. No candidate displaces rho.

#pagebreak()

= Factor bases and symmetry

#table(
  columns: (1.15fr, 1.55fr, 0.85fr, 0.75fr, 0.75fr, 1.6fr),
  align: (left, left, right, right, right, left),
  table.header([cell], [candidate], [points / cols], [mean S], [S/rho], [reading]),
  [K17 a=1], [compact orbit], [272 / 8], [9.114], [9.985], [best in matched cell; lower bound],
  [K17 a=1], [full orbit], [443 / 14], [17.081], [18.714], [complete count],
  [K17 a=1], [binary subspace], [10 / 5], [43.103], [47.224], [complete count],
  [K17 a=1], [standard subspace], [11 / 6], [43.150], [47.275], [complete count],
  [K23 a=0], [compact orbit], [460 / 10], [2.996], [8.360], [separate curve],
  [K23 a=1], [cyclotomic w=1], [46 / 1], [3.427], [9.933], [separate curve],
  [prime 20], [abscissa], [64 / 32], [11.329], [12.300], [better than GLV here; lower bound],
  [prime 20], [GLV orbit], [32 / 16], [22.390], [24.309], [lower bound],
)

The native catalog contains ten sealed bases and reports 10/10 remote
verification; seven IC arms have performance rows. Compact orbit is the local
Koblitz winner, and ordinary abscissa is the local prime winner. Neither beats
rho.

The legacy torsion-symmetrised u-frame has whole-pipeline `S_sym/S_x` 0.84 at
n17 and 0.77/0.76 on the two n23 curves, then 4.01 at n31. Its combinatorial
controls are 0.89, 0.90, 1.09, and 1.28. The small-size gain is not persistent.
The current ECBench registry has no torsion-symmetrised factor-base plug-in.

Signed Frobenius is a separate symmetry: it is native in rho and owns the
canonical Koblitz operations frontier. Factor-base orbit folding must not be
credited with the walk's automorphism discount.

= F4, F5, F6-IC

#table(
  columns: (1.35fr, 2fr, 2.1fr),
  align: (left, left, left),
  table.header([comparison], [observation], [boundary verdict]),
  [matrix F5 / matrix F4], [word-XOR ratio 0.6465 K13, 0.6414 R13; same trees], [solver-stage engineering win],
  [matrix F4 / inherited F4], [198.6x K13, 222.1x R13], [inherited F4 is the practical algebraic baseline],
  [F6-IC / inherited F4], [solver-XOR ratio about 0.82 at m13 and 0.77-0.78 at m19; relation wall 1.7-1.9x faster], [solver-stage win],
  [charged F6-IC S], [64.6-68.8x inherited at m13; 1,303-1,353x at m19], [relabelling/accounting loss],
  [matrix m19], [four 1,800 s timeouts], [censored],
  [F6 ladder m31], [both arms hit 3,600 s on target 1], [protocol decision open],
)

Every IC arm is above matched rho on every target in the cross-family cells.
F5 and F6 improve a subroutine, not the cold ECDLP frontier.

#pagebreak()

= Root finding

Three root-finder lines are distinct.

1. *Native ECBench pair/subtraction engines.* K13 pair-table root finding has
   complete S 28.383 (m2) and 38.822 (m3): 15.170x and 20.749x matched
   signed-Frobenius rho. It matches MITM operations exactly. Descriptive L0
   wall favours it over inherited F4, but a complete count cannot be divided
   by an F4 lower bound to claim an operation speedup.
2. *PDP root-finding study.* The table/enumeration oracles beat F4 within the
   decomposition stage, but every measured slope equals its generic
   pair-table or enumeration floor. Classification: engineering, not an
   exponent change.
3. *N83 dense five-byte index.* K1798 lowers mean probes from 300,021.90625
   to 169,069.8125. The ratio is 0.563525; the paired reduction is
   130,952.09375 [52,057.703125, 209,990.828125]. It uses an 8 GiB root table
   and is not registered in `ecbench methods`.

#figure(
  image("online_boundary.svg", width: 92%),
  caption: [The same methods cross parity only after reusable preparation is
  excluded. Both series use rho cost / IC cost. The n61 online interval and
  cold ratio are reported measurements. N83 is a one-point shared-host
  observation; its lower bar divides rho online by the complete IC process and
  is a setup diagnostic, not a matched cold ratio.],
)

#pagebreak()

= Current-binary replay audit

Build: Rust 1.94.1, current commit, auditor SHA-256
`53654ef1c0e93ed02025acd8188514d2674c4da2f57f25f61db40a510a422950`,
hardware class `ECBHW1h7373da65e433`.

#table(
  columns: (1.7fr, 1.55fr, 1.05fr, 2fr),
  align: (left, left, left, left),
  table.header([session], [audit], [receipt prefix], [reading]),
  [factor bases K17], [8/8 identical], [`9306d33e…`], [current-binary pass],
  [F4/F5/F6 K13], [7/7 identical], [`0304c496…`], [representative pass covering all six arms; not replay-all],
  [rootfinder K13 sample], [8/8 identical], [`4814136e…`], [sampler selected rho only],
  [rootfinder K13 replay-all], [120/168 identical; 48 failures], [`6110727a…`], [all rootfinder/rho rows pass; every old inherited-F4 row differs in solver phase counts],
)

The replay-all failure is preserved. It changes no scalar, factor base, trial
count, charged GAE, or whole-pipeline S. Seq 58, for example, retains
`total_gae=741.091311`, but current x86_64 counts 301,961 solver word XORs
where the archived Apple record counts 300,531. The committed same-host audit
made with the measuring binary reports 168/168 identical. Until this
cross-host counter drift is explained or the method identity is versioned,
the root-finder study's F4 stage-count comparison is historical evidence. The
root-finder rows themselves replay identically.

= Exact frontier-moving obligations

1. Register torsion-symmetrised bases and N83 dense root finding under new,
   immutable ECBench method IDs, with cold and five-phase online records.
2. Resolve inherited-F4 word-XOR cross-host replay drift before fitting a bound.
3. Run a native challenge at four or more sizes against the standing rho bound.
4. Repeat N83 on a Tier 1/0 isolated host, charging every cold phase and byte.
5. Pin a conversion for F6's batched geometry before interpreting wall/GAE
   disagreement.

*Non-claims.* Public synthetic known-answer inputs only; no deployed-curve or
key-recovery claim; no ECC2K-130 transfer; no new exponent; no canonical
scoreboard promotion. Existing canonical graphs were checked and left
unchanged. The required Conductor executable and workflow file were absent, so
work stayed inside this report, benchmark copies, and benchmark index.
