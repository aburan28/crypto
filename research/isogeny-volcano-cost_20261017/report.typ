#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)
#set table(stroke: 0.4pt)
#show image: it => block(width: 100%, it)
= Isogenous curves and ECDLP cost: vertical, horizontal and coordinate
effects
<isogenous-curves-and-ecdlp-cost-vertical-horizontal-and-coordinate-effects>
#strong[Question.] On the curves of one ℓ-isogeny volcano, with one
ECDLP instance carried to all of them, does solving cost differ:

- #strong[vertically], between levels (crater against floor);
- #strong[horizontally], between curves of one level (one endomorphism
  ring, hence horizontally isogenous);
- #strong[by coordinates], between isomorphic models of one curve?

The methods are Pollard rho, BSGS and index calculus. The protocol was
frozen before any data (#link("PROTOCOL.md")[`PROTOCOL.md`], v1). The
design is in #link("PLAN.md")[`PLAN.md`]. The success criterion
(CLAUDE.md) asks for the gap to rho, reported as the distance remaining.

== Answer
<answer>
#strong[Measured, PROTOCOL v1.1 (§ "Amendment"); 6 classes, 200 curves ×
4 models: 43,200 measured runs (46,912 with the two control classes,
each after one warm-up), every one verified; 34 sessions audited,
408/408 replays exact.]

#align(center)[#table(
  columns: 5,
  align: (col, row) => (left,left,left,left,left,).at(col),
  inset: 6pt,
  [method], [vertical: floor / crater, 90% CI], [verdict (fixed
  wording)], [horizontal: τ between curves of a level (null 95th %ile),
  p], [coordinates: τ between isomorphic models (null 95th %ile), p],
  [BSGS (`bsgs.negation`)],
  [1.000 \[1.000, 1.000\]],
  [identical cost on every curve and model for every paired target
  (exact, 12/12 targets, 800 workload curves)],
  [0],
  [0],
  [rho (`rho.negation`)],
  [1.022 \[0.983, 1.062\], p \= 0.28],
  [no level effect detected; #strong[not equivalent within ±5%] (the
  interval reaches +6.2%)],
  [1.7% (2.6%), Q p \= 0.28, max|z| p \= 0.73],
  [2.6% (4.6%), p \= 0.70],
  [rho, A/A seeds],
  [1.001 \[0.965, 1.039\], p \= 0.94],
  [the noise floor of the same design],
  [0% (2.7%), p \= 0.87],
  [2.5% (4.6%), p \= 0.75],
  [IC relation yield (`ic.pipeline`, 32 smallest abscissae, MITM)],
  [0.980 \[0.957, 1.004\], p \= 0.20],
  [#strong[inconclusive]: the planted +5% floor control was not detected
  (§ Controls), though the interval lies inside ±5%],
  [2.2% (1.4%), p \= 0.0005],
  [#strong[4.7% (2.5%), p \= 0.0005]],
  [IC cost S (reported, not tested)],
  [1.039 \[1.008, 1.071\], p \= 0.061],
  [],
  [2.9% (2.1%), p \= 0.002],
  [5.6% (3.7%), p \= 0.0005],
)
]

What it says:

+ #strong[BSGS:] the cost is a function of the instance alone. Carrying
  the instance by isogeny, and changing coordinates, leaves it exactly
  unchanged on 800 workload curves (measured; it also follows from BSGS
  being generic).
+ #strong[Rho:]
  - No vertical, horizontal or coordinate effect was detected. Every
    statistic sits inside its null band and the A/A arm reads the same.
  - The vertical bound is ±6.2%, not the ±5% the protocol aimed for. The
    design’s predicted ±2% did not hold: with 1–11 crater curves per
    class, the floor/crater ratio is noisier than pooled solve counts
    suggest.
+ #strong[Index calculus:]
  - #strong[Horizontal (measured):] curves of one level differ in
    relation yield by τ ≈ 2.2%. That is above the null, so it is
    significant.
  - #strong[Coordinates (measured):] isomorphic models of one curve
    differ #strong[more], τ ≈ 4.7%. By the protocol’s rule this is a
    #strong[coordinate effect], not an isogeny effect. The factor base
    "the 32 points of smallest abscissa" is a different set of points in
    each model, and the yield follows that set, not the curve.
  - #strong[Vertical:] the pooled test cannot see a 5% effect (failed
    analysis control), so it is inconclusive. Exploratory per-class
    contrasts (§ Exploratory) suggest an ℓ \= 2-only effect of about
    −6.5%.
+ #strong[Gap to rho (distance remaining):] on every one of the 800
  workload curves, rho is below IC at every weight of the §5a band (0,
  0.1 and 1 GAE per unpriced native operation; `band/`).
  - Per curve, IC costs S ≈ 27–96 at w \= 0 against rho’s 0.8–1.6.
  - So IC is about 30–60× above rho at 24 bits.
  - No isogenous curve or model closes any part of that gap beyond the
    few-percent effects above.

Scope: 24-bit `p`, depth-3 volcanoes, ℓ ∈ {2, 3}, the three methods at
these parameters. No security statement. Nothing above 24 bits, and no
other factor base, is covered.

== Controls
<controls>
#align(center)[#table(
  columns: 4,
  align: (col, row) => (left,left,left,left,).at(col),
  inset: 6pt,
  [control], [measured], [predicted], [verdict],
  [planted IC: factor base 33 vs 32 points, yield ratio matched over
  3,200 workloads],
  [1.0613 \[1.0583, 1.0644\]],
  [(33/32)² ≈ 1.063],
  [detected; instrument sees +6%],
  [j \= 0 crater: `rho.cubic` / `rho-neg-crater`, 64 workloads],
  [0.697 \[0.589, 0.849\]],
  [√(1/3) ≈ 0.58],
  [detected (excludes 1)],
  [j \= 1728 crater: `rho.quartic` / `rho-neg-crater`],
  [0.707 \[0.604, 0.815\]],
  [√(1/2) ≈ 0.71],
  [detected (reported only)],
  [BSGS exactness, all classes],
  [1 distinct cost per paired target],
  [1],
  [exact],
  [A/A rho (independent measurement seed), matched pair],
  [0.9915 \[0.9795, 1.0032\]],
  [1],
  [within noise],
  [analysis: +5% planted on every floor curve, rho],
  [V 1.022 → 1.073, p \= 0.0006],
  [detect],
  [#strong[detected]],
  [analysis: +5% planted on every floor curve, IC yield],
  [V 0.980 → 1.029, p \= 0.073],
  [detect at 0.0125],
  [#strong[not detected] → IC vertical inconclusive],
  [analysis: +10% planted on one curve per class, rho],
  [max|z| p 0.73 → 0.73],
  [detect],
  [#strong[not detected] → single-curve rho outliers below ~±15% are not
  excluded],
  [analysis: +10% planted on one curve per class, IC yield],
  [max|z| p 0.0105 → 0.0005],
  [detect],
  [detected],
)
]

The planted-control runs are in `analysis-controls/`, for both the v1
and v1.1 statistics.

== Amendment v1.1 (after seeing data, stated as such)
<amendment-v1.1-after-seeing-data-stated-as-such>
BSGS cost is identical on every curve, but the frozen v1 statistic gave
its floor/crater ratio a bootstrap interval of \[0.88, 1.22\].

- #strong[Why:] v1 pools all classes into one floor mean and one crater
  mean. The classes contribute different shares of floor and crater
  curves. So when a bootstrap replicate redraws a class’s paired targets
  and shifts that class’s mean, the ratio moves even though no curve
  differs from another.
- #strong[Fix:] v1.1 forms the ratio inside each class and combines the
  six ratios as an equal-weight mean of logs
  (`volcano-stats --vertical   stratified`). The permutation tests are
  valid under either statistic.

Both versions are reported:

#align(center)[#table(
  columns: 3,
  align: (col, row) => (left,left,left,).at(col),
  inset: 6pt,
  [arm · metric], [v1 pooled: V, 90% CI, p], [v1.1 stratified: V, 90%
  CI, p],
  [BSGS S],
  [1.000 \[0.881, 1.224\], 1.0],
  [1.000 \[1.000, 1.000\], 1.0],
  [rho S],
  [1.024 \[0.988, 1.062\], 0.071],
  [1.022 \[0.983, 1.062\], 0.275],
  [rho S, A/A],
  [0.978 \[0.947, 1.010\], 0.069],
  [1.001 \[0.965, 1.039\], 0.938],
  [IC yield],
  [#strong[0.972 \[0.951, 0.993\], 0.0012]],
  [0.980 \[0.957, 1.004\], 0.202],
  [IC S],
  [1.029 \[0.998, 1.060\], 0.024],
  [1.039 \[1.008, 1.071\], 0.061],
)
]

Under v1 the IC-yield vertical test is nominally significant (p \=
0.0012 \< 0.0125). The same planted analysis control fails under v1 too
(+5% detected only at p \= 0.081). Under v1.1, which removes the
class-mix artefact the BSGS check exposed, it is not significant. By the
frozen v1 rule the v1 result would read "level changes IC yield by
−2.8%". This report does not claim that, because the statistic behind it
is shown wrong by an exact control. The v1 numbers stay here and in
`analysis.json`.

== Exploratory (not tested; six unadjusted contrasts)
<exploratory-not-tested-six-unadjusted-contrasts>
Per-class floor/crater, v1.1 statistic (`analysis-exploratory/`):

#align(center)[#table(
  columns: 4,
  align: (col, row) => (left,left,left,left,).at(col),
  inset: 6pt,
  [class], [E(F\_p) 2-part: floor vs above], [rho S], [IC relation
  yield],
  [`l2-ramified`],
  [cyclic vs Z/2 × Z/2],
  [0.996 \[0.889, 1.124\]],
  [#strong[0.930 \[0.863, 0.994\]], p \= 0.033],
  [`l2-split`],
  [cyclic vs Z/2 × Z/2],
  [1.043 \[0.960, 1.139\]],
  [#strong[0.937 \[0.888, 0.987\]], p \= 0.001],
  [`l2-inert` (1 crater curve)],
  [cyclic vs Z/2 × Z/2],
  [1.074 \[0.931, 1.225\]],
  [0.988 \[0.906, 1.086\]],
  [`l3-ramified`],
  [all cyclic],
  [1.012 \[0.897, 1.130\]],
  [0.964 \[0.913, 1.018\]],
  [`l3-split`],
  [all cyclic],
  [1.019 \[0.940, 1.110\]],
  [1.003 \[0.958, 1.052\]],
  [`l3-inert`],
  [all cyclic],
  [0.991 \[0.860, 1.140\]],
  [1.066 \[0.983, 1.165\]],
)
]

#strong[Hypothesis for PROTOCOL v2 (not a finding).]

- #strong[Pattern:] in ℓ \= 2 volcanoes, floor curves (cyclic 2-part,
  one rational 2-torsion point) give about 6.5% fewer IC relations per
  trial than crater curves (full rational 2-torsion). ℓ \= 3 volcanoes,
  whose groups are cyclic at every level, show nothing.
- #strong[Mechanism to test:] the small-abscissa factor base is drawn
  from the whole group, and how its points split across the 2-power
  torsion cosets differs with the group structure.
- #strong[What would confirm it:] fresh ℓ \= 2 classes with several
  crater curves each, a pre-registered one-sided test, and the IC cost S
  on the same curves.
- #strong[What would make it an improvement:] a confirmed yield gain
  also has to show in S and survive the §5a band before it counts as an
  index-calculus improvement (CLAUDE.md success criterion).

== Population and instrument
<population-and-instrument>
Six ordinary classes plus two positive-control classes, all 24-bit,
depth 3. The table (p, root model, curves per level, group structures)
is in #link("PROTOCOL.md")[`PROTOCOL.md`] §2.

- #strong[Curves and transport:] `volcano-rho --ecbench` exported every
  curve of each volcano, up to 12 per level, with every crater curve
  kept. Each curve came in 4 isomorphic models, with the
  isogeny-transported generator and `planted_shared` targets, so every
  curve and model solved the same k.
- #strong[Runs:] `ecbench` (L0, counted units only); each session’s
  command and spec are in `classes/` and `sessions/`.
- #strong[Analysis:] `volcano-stats` (native), with bootstrap 2,000 and
  permutations 10,000 (2,000 for the heterogeneity nulls). Files:
  `analysis-v1.1.json` (headline), `analysis.json` (v1 as frozen) and
  `analysis-j0.json` / `analysis-j1728.json` (controls).
- #strong[Figures:] `plot.py` → `figures/`: one volcano heat map per
  class with rho cost and IC yield panels, `vertical_forest` and
  `heterogeneity`.

#figure([#image("figures/vertical_forest.svg")],
  caption: [
    vertical
  ]
)

#figure([#image("figures/heterogeneity.svg")],
  caption: [
    heterogeneity
  ]
)

#figure([#image("figures/volcano-l2-split.svg")],
  caption: [
    l2-split
  ]
)

#figure([#image("figures/volcano-l2-ramified.svg")],
  caption: [
    l2-ramified
  ]
)

#figure([#image("figures/volcano-l3-split.svg")],
  caption: [
    l3-split
  ]
)

== Evidence and limits
<evidence-and-limits>
- #strong[Sessions:] 34 complete sessions, 46,912 measured runs (each
  after one warm-up), every run `verified`. `ecbench verify --replay 12`
  passed on every session: 408/408 replays identical. Audits are in
  `audits/`.
  - `l3-split/spec` was interrupted at 7,784/12,032 runs by the host’s
    2-hour background limit, not by its result. It stays as evidence;
    the analysis uses its complete re-run `l3-split/spec-r2`.
  - The two BSGS sessions the limit prevented from starting ran
    afterwards.
- #strong[Storage:] records and plans are gzipped (`*.jsonl.gz`,
  `plan.json.gz`; `MANIFEST.sha256` holds the uncompressed hashes).
  Decompress before `ecbench verify`.
- #strong[Isolation:] the host is a shared 4-vCPU VM at L0. No wall time
  is cited. Counted units are host-independent. Independent-host replay
  is not done, so this is not a claim under the IC claim rules, and it
  makes no IC-beats-rho claim.
- #strong[Power:]
  - The rho vertical bound is ±6.2%.
  - IC vertical is inconclusive.
  - Rho single-curve outliers below about ±15% are not excluded:
    per-curve rho has 64 runs, about ±11% at 95%.
  - More crater curves (more classes, more models) is the lever for all
    three.
- #strong[Not run (protocol v2):] the 20- and 28-bit size arms, and rho
  at 32 and 40 bits (PLAN §4).
- #strong[Conductor:] not installed in this container, so no task was
  reserved.
