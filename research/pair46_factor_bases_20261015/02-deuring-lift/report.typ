#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set par(justify: true)

= 02-deuring-lift: canonical CM(−19) lift of the pair-46 crater

Pair-46 factor-base campaign, item 02 · 2026-10-15 · protocol `PROTOCOL.md` (frozen before data) · data `results.json`, `summary.json`, `run.log`.

*Verdict: exact obstruction (prediction "fails by rank obstruction" confirmed).*

== Exact statements

+ The Deuring lift of the crater is Cremona *361a1*, $y^2+y=x^3-38x+90$, with $j=-96^3$ and conductor 361. Its $a_p$ is $t = 34489701$, so $\#E(bb(F)_p)=N$, and the explicit isomorphism to the crater (with $u^2$ a square mod $p$) is verified. No twist is needed.
+ $E(bb(Q)) = bb(Z) G$ with $G=(0,9)$. The rank is 1 unconditionally (ellrank bounds [1,1]; analytic rank 1), torsion is trivial, and $hat(h)(G)=0.181292$.
+ $G mod p$ has order $N$, so reduction is onto. The lift of a target $T=k G$ is $k G$, with height $k^2 hat(h)(G)$, and its expected height is $approx 2.1 dot 10^30$ nats. Lifting is therefore the DLP itself, and the only relations that involve the target are the trivial ones once $k$ is known.
+ Twists $E_d$ ($1<|d| <= 200$, one per ${d,-19d}$ class, 115 curves): 51 have rank 0, 54 rank 1, 8 rank 2 and 2 rank 3. Generators on distinct twists are $bb(Z)$-independent (they lie in distinct Galois eigenspaces), so a twist factor base has *0* $bb(Q)$-relations. Every relation among them exists only mod $p$, as a point of a determinant-$N$ kernel lattice.
+ Over $bb(Q)(sqrt(-19))$ the only extra relation is $w G equiv lambda G$, the known GLV relation, which the glv-orbit item already covers. This one is a derivation and was not computed here.

== Empirical (24 toy primes at 16/20/24 bits, $R=6$)

#table(columns: 4, align: left,
  [*metric*], [*A2 canonical*], [*NC random*], [*PC planted*],
  [$lambda_1 \/ "GH"$ median], [1.393], [1.250], [0.290 (max 0.461)],
  [A2 − NC (95% CI)], table.cell(colspan: 3)[+0.143 [−0.180, +0.262], Mann–Whitney p = 0.50],
  [coverage diff], table.cell(colspan: 3)[+0.0008 [−0.052, +0.056]],
  [lifted Q-relations], [0 / 24 primes], [n/a], [48 (planted 2 per prime)],
  [Xedni generic $r=2..5$], [0 / 96 trials], [n/a], [24 / 24 ($P,-P$ planted)],
)

#figure(image("fig_relations.svg", width: 100%), caption: [Left: shortest mod-$p$ relation among $R=6$ factor-base points, normalized by the Gaussian heuristic. The canonical-lift generators (A2) cannot be distinguished from random points (NC), while the planted dependency (PC) is caught at every prime. Right: characteristic-0 relations (|F| − rank of the canonical-height Gram matrix) against the number a factor base needs. 24 toy primes; run `run_deuring_lift.py`, seed 20261015.])

#figure(image("fig_twist_ranks.svg", width: 90%), caption: [Mordell–Weil ranks of the quadratic twists of 361a1 for $1<|d|<=200$ (PARI ellrank, every bound closed). The largest rank is 3.])

== Deviations and scope

- *Precision fix.* cypari2 defaults to 15 digits of precision. At that precision the positive control (the planted $G,2G,3G$ Gram matrix) failed the frozen $10^(-20)$ tolerance, so I raised realprecision to 80 for every arm. The tolerance itself was not changed.
- *Nothing cut.* Nothing in the preregistered design was cut, and the whole run took 2.4 s on 1 core.
- *Limits.* One generator per twist, $|d|<=200$, and no DLP at pair-46 size. The pair-46 statements are exact.

*Reproduce:* `python3 run_deuring_lift.py results.json && python3 analyze.py && python3 -c "import typst; typst.compile('report.typ', output='report.pdf')"`
