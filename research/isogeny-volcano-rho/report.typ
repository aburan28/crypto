#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)
#set table(stroke: 0.4pt)
#show image: it => block(width: 100%, it)
= Pollard rho across the strata of an isogeny volcano
<pollard-rho-across-the-strata-of-an-isogeny-volcano>
#strong[Question.] Take one ECDLP instance, carry it by explicit
isogenies to every level of an ℓ-isogeny volcano, and run Pollard rho on
each curve: does the level (crater → floor), or the particular curve,
change how fast rho solves it? Report min, max and mean.

#strong[Update:] curves D (j \= 0 crater, `run-20261017T0140Z`) and E (j
\= 1728 crater, `run-20261017T0305Z`) show the one place the curve
matters. Rho folded by the crater’s automorphisms needs 1.89 ± 0.09× (j
\= 0, √3 predicted) and 1.32 ± 0.06× (j \= 1728, √2 predicted) fewer
operations. In wall time that is 1.33× faster on j \= 0 and 2.3× slower
on j \= 1728, because `ca` has a tuned walk only for j \= 0. A separate
control (`run-20261017T0415Z-instance`, 5,000 solves) shows that the
generator and the secret k do not change rho’s cost either: redrawing k,
P, both, or using k ≤ 2^16 all match the walk-seed spread (every |z| \<
0.8, KS p \> 0.6).

#strong[Answer (measured on three curves: A, `run-20261016T2302Z`, ℓ \=
3, ramified crater; B, `run-20261016T2315Z`, ℓ \= 3, split crater; C,
`run-20261017T0049Z`, ℓ \= 2, a whole depth-3 2-volcano).] No measurable
effect on any of them. Curve A in detail: All 29 isogenous curves and 3
isomorphic controls took the same work to within sampling noise. Pooled
over the 2,900 isogenous solves: min 8,477, mean 659,263 ± 12,422 (95 %
CI), max 2,099,707 group operations, against √(πn/2) \= 654,547. Every
level’s 95 % CI includes or touches that value, and per-step cost is
flat (641–661 ns/op on every curve). That is what theory predicts:
isogenous curves over F\_p have the same group order, and generic rho
depends only on the order of ⟨P⟩. The spread #emph[between] curves
(per-curve means 572,846 – 733,141) matches the spread of 100-trial
means of a quantity whose single-solve sd is about 341,000 (sd of a mean
≈ 34,000).

== Method
<method>
`volcano-rho` (`src/bin/volcano_rho.rs`, native Rust, all statistics
computed there; `plot.py` only lays out recorded values):

+ #strong[Curve search.] Random curves `y² = x³ + ax + b` over random
  40-bit primes `p ≡ 3 (mod 4)`, point count by baby-step giant-step,
  until `D = t² − 4p = f²·D_K` has `v₃(f) = 3` (a 3-volcano with three
  levels below the crater) and `#E = h·n` with `n` prime, `h ≤ 64`.
+ #strong[Walk.] Breadth-first over 3-isogenies. Each edge’s kernel is a
  rational root `x₀` of the 3-division polynomial
  `ψ₃ = 3x⁴ + 6ax² + 12bx − a²`, and the map is Vélu’s, applied
  explicitly to points. Vertices are keyed by j-invariant. A vertex’s
  level is `depth − dist(v, floor)`, the floor being the vertices with
  exactly one rational 3-kernel; the tool asserts that no edge skips a
  level.
+ #strong[One problem everywhere.] `P` of order `n` and `Q = kP` on the
  root curve are pushed along the BFS tree to every curve. Each edge is
  checked: images lie on the codomain, `φ(P+Q) = φ(P) + φ(Q)`, and
  `n·φ(P) = O` with `φ(P) ≠ O` (3 ∤ n, so every map is injective on
  ⟨P⟩). Every rho solve on every curve returned the same `k`.
+ #strong[Rho.] Single walker, 32-adding walk (Teske), distinguished
  points with θ \= 2⁻⁸, no negation map; restarts on a DP-free run
  longer than 40/θ. Trial seed `s` uses the same walk constants on every
  curve (the partition hashes `x`, so walks differ). Measured unit:
  group additions to the solve; also wall time per solve.
+ #strong[Control.] Three random isomorphic re-scalings of the root,
  `(u⁴a, u⁶b)`, `(x, y) ↦ (u²x, u³y)`: the same curve in new
  coordinates, showing the spread that a pure change of representation
  produces.
+ All 3,200 (arm, seed) jobs shuffled and run on 4 threads.

== Curve A: the instance (`run-20261016T2302Z`, ramified crater)
<curve-a-the-instance-run-20261016t2302z-ramified-crater>
#align(center)[#table(
  columns: 2,
  align: (col, row) => (left,left,).at(col),
  inset: 6pt,
  [p],
  [1090993606999 (40 bits)],
  [root E],
  [a \= 690282649366, b \= 1050606678233],
  [\#E],
  [1090993806052 \= 4 · n, n \= 272748451513 (prime, 38 bits)],
  [trace t],
  [−199052],
  [D \= t² − 4p],
  [−4324352729292 \= 54² · (−1482974187); v₃(f) \= 3],
  [crater],
  [D\_K ≡ 0 (mod 3): ramified, 2 crater vertices joined by one edge],
  [k],
  [49631396114 (Q \= kP)],
  [volcano],
  [80 curves \= 2 + 6 + 18 + 54 (crater, L1, L2, floor), 79 edges, all
  verified],
)
]

== Curve A results (iterations \= group operations; time on a shared
4-vCPU VM)
<curve-a-results-iterations-group-operations-time-on-a-shared-4-vcpu-vm>
Pooled by level:

#align(center)[#table(
  columns: 8,
  align: (col, row) => (left,right,right,right,right,right,right,right,).at(col),
  inset: 6pt,
  [level], [curves], [solves], [min], [mean ± 95 % CI], [median], [max],
  [mean ms ± CI],
  [0 crater],
  [2],
  [200],
  [17,415],
  [700,483 ± 47,849],
  [662,463],
  [1,621,252],
  [456.9 ± 32.0],
  [1],
  [6],
  [600],
  [48,186],
  [681,766 ± 26,525],
  [656,746],
  [1,861,764],
  [440.7 ± 17.6],
  [2],
  [10],
  [1,000],
  [9,649],
  [637,158 ± 20,410],
  [598,424],
  [1,792,889],
  [414.4 ± 13.6],
  [3 floor],
  [11],
  [1,100],
  [8,477],
  [659,589 ± 20,993],
  [620,458],
  [2,099,707],
  [426.0 ± 13.9],
  [all isogenous],
  [29],
  [2,900],
  [8,477],
  [659,263 ± 12,422],
  [622,386],
  [2,099,707],
  [427.1 ± 8.2],
  [isomorphic control],
  [3],
  [300],
  [23,698],
  [674,042 ± 40,816],
  [628,245],
  [1,820,972],
  [434.0 ± 27.0],
)
]

Reference: √(πn/2) \= 654,547 (+ about 2⁸ for the distinguished-point
tail).

Per curve (100 solves each), iterations:

#align(center)[#table(
  columns: 11,
  align: (col, row) => (left,right,right,right,right,center,left,right,right,right,right,).at(col),
  inset: 6pt,
  [curve], [level], [min], [mean], [max], [], [curve], [level], [min],
  [mean], [max],
  [v5],
  [0],
  [113,746],
  [732,044],
  [1,621,252],
  [],
  [v18],
  [2],
  [44,056],
  [649,444],
  [1,724,073],
  [v10],
  [0],
  [17,415],
  [668,921],
  [1,565,882],
  [],
  [v19],
  [2],
  [106,910],
  [621,995],
  [1,638,981],
  [v4],
  [1],
  [108,640],
  [695,082],
  [1,655,814],
  [],
  [v20],
  [2],
  [117,650],
  [667,890],
  [1,600,331],
  [v8],
  [1],
  [115,038],
  [665,553],
  [1,802,980],
  [],
  [v21],
  [2],
  [72,048],
  [692,313],
  [1,527,838],
  [v9],
  [1],
  [48,186],
  [618,080],
  [1,325,879],
  [],
  [v22],
  [2],
  [71,551],
  [620,647],
  [1,792,889],
  [v23],
  [1],
  [134,281],
  [716,324],
  [1,422,523],
  [],
  [v44],
  [2],
  [54,056],
  [674,911],
  [1,686,754],
  [v24],
  [1],
  [123,525],
  [733,141],
  [1,861,764],
  [],
  [v0 (root)],
  [3],
  [55,979],
  [658,723],
  [1,634,032],
  [v25],
  [1],
  [55,435],
  [662,413],
  [1,730,218],
  [],
  [v2],
  [3],
  [54,631],
  [659,817],
  [1,854,099],
  [v1],
  [2],
  [69,112],
  [572,846],
  [1,308,207],
  [],
  [v3],
  [3],
  [107,168],
  [629,575],
  [1,601,367],
  [v6],
  [2],
  [38,512],
  [628,345],
  [1,584,239],
  [],
  [v11],
  [3],
  [94,240],
  [678,053],
  [2,094,241],
  [v7],
  [2],
  [9,649],
  [616,240],
  [1,283,603],
  [],
  [v12],
  [3],
  [64,211],
  [634,612],
  [1,509,154],
  [v17],
  [2],
  [100,833],
  [626,950],
  [1,696,010],
  [],
  [v13],
  [3],
  [8,477],
  [704,738],
  [1,838,374],
  [iso0],
  [ctrl],
  [23,698],
  [672,888],
  [1,820,972],
  [],
  [v14],
  [3],
  [17,758],
  [676,654],
  [1,460,561],
  [iso1],
  [ctrl],
  [44,112],
  [701,414],
  [1,670,056],
  [],
  [v15],
  [3],
  [49,801],
  [663,488],
  [1,678,869],
  [iso2],
  [ctrl],
  [109,711],
  [647,822],
  [1,591,405],
  [],
  [v16],
  [3],
  [68,503],
  [659,946],
  [2,099,707],
  [],
  [],
  [],
  [],
  [],
  [],
  [v26],
  [3],
  [51,815],
  [655,015],
  [2,017,181],
  [],
  [],
  [],
  [],
  [],
  [],
  [v27],
  [3],
  [46,289],
  [634,861],
  [1,689,685],
)
]

Full per-curve rows (with times, ns/op, medians, p10/p90) are in
`runs/run-20261016T2302Z/results.json`; every solve is a row of
`trials.csv`.

#image("runs/run-20261016T2302Z/volcano.svg")
#image("runs/run-20261016T2302Z/per_curve.svg")
#image("runs/run-20261016T2302Z/per_level.svg")

== Curve B: a second volcano with a split crater (`run-20261016T2315Z`)
<curve-b-a-second-volcano-with-a-split-crater-run-20261016t2315z>
Same tool, `--crater split --seed 20261017`. This volcano’s crater is a
cycle, so the walk reaches a long run of crater vertices, each with its
own depth-3 tree. The cycle is truncated at the walk radius.

#align(center)[#table(
  columns: 2,
  align: (col, row) => (left,left,).at(col),
  inset: 6pt,
  [p],
  [714095314843 (40 bits)],
  [root E],
  [a \= 354455170800, b \= 66680642576],
  [\#E],
  [714094928137 \= n, prime (cofactor 1)],
  [trace t],
  [386707],
  [D \= t² − 4p],
  [−2706838955523 \= 27² · (−3713084987); v₃(f) \= 3],
  [crater],
  [D\_K ≡ 1 (mod 3): split, the crater is a cycle],
  [k],
  [661235192333],
  [walk],
  [323 curves, 322 edges verified; inside the measured radius 11 crater,
  18, 42 and 90 floor curves by level],
  [rho],
  [41 curves + 3 controls × 100 solves \= 4,400, all correct; 754 s],
)
]

#align(center)[#table(
  columns: 8,
  align: (col, row) => (left,right,right,right,right,right,right,right,).at(col),
  inset: 6pt,
  [level], [curves], [solves], [min], [mean ± 95 % CI], [median], [max],
  [mean ms ± CI],
  [0 crater],
  [10],
  [1,000],
  [42,088],
  [1,091,654 ± 35,041],
  [1,032,925],
  [3,188,276],
  [685.7 ± 22.6],
  [1],
  [10],
  [1,000],
  [42,018],
  [1,090,391 ± 34,841],
  [1,046,357],
  [3,327,413],
  [681.2 ± 22.0],
  [2],
  [10],
  [1,000],
  [44,771],
  [1,081,516 ± 35,498],
  [1,012,223],
  [3,437,779],
  [678.1 ± 22.7],
  [3 floor],
  [11],
  [1,100],
  [36,930],
  [1,099,508 ± 33,756],
  [1,050,592],
  [3,110,345],
  [686.8 ± 21.1],
  [all isogenous],
  [41],
  [4,100],
  [36,930],
  [1,090,981 ± 17,380],
  [1,034,403],
  [3,437,779],
  [683.0 ± 11.0],
  [isomorphic control],
  [3],
  [300],
  [51,690],
  [1,097,651 ± 60,787],
  [1,094,657],
  [2,961,196],
  [688.9 ± 38.9],
)
]

√(πn/2) \= 1,059,102. Per-curve means span 951,361 (v22) to 1,210,896
(v2); the per-op cost is 622–635 ns on every curve. Per-curve rows are
in `runs/run-20261016T2315Z/results.json`.

#image("runs/run-20261016T2315Z/volcano.svg")
#image("runs/run-20261016T2315Z/per_curve.svg")
#image("runs/run-20261016T2315Z/per_level.svg")

== Curve C: a whole 2-isogeny volcano (`run-20261017T0049Z`)
<curve-c-a-whole-2-isogeny-volcano-run-20261017t0049z>
Same tool with 2-isogenies:
`--ell 2 --bits 48 --depth 3 --max-cofactor 4096 --seed 20261018`. A
2-volcano deeper than one level needs `p ≡ 1 (mod 8)`; with
`p ≡ 3 (mod 4)`, `v₂(t² − 4p) ≤ 3`. ℓ \= 2 searches draw such primes.
Depth 3 also puts at least 2⁶ into `#E`, hence the larger cofactor
bound. The walk reaches the whole volcano, a binary tree under a
two-vertex crater, and every curve in it is measured.

#align(center)[#table(
  columns: 2,
  align: (col, row) => (left,left,).at(col),
  inset: 6pt,
  [p],
  [146833024135393 (48 bits, ≡ 1 mod 8)],
  [root E],
  [a \= 14391545271884, b \= 141956518099848],
  [\#E],
  [146833030777920 \= 960 · n, n \= 152951073727 (prime, 37.2 bits)],
  [trace t],
  [−6642526],
  [D \= t² − 4p],
  [−543208944880896 \= 56² · (−173217138036); v₂(f) \= 3],
  [crater],
  [ramified (D\_K even), 2 crater vertices joined by one edge],
  [k],
  [50935020165],
  [walk],
  [30 curves \= 2 + 4 + 8 + 16 (the whole volcano), 29 edges verified],
  [rho],
  [30 curves + 3 controls × 100 solves \= 3,300, all correct; 310.7 s],
  [cross-check],
  [`ca volcano --ell 2` (aburan28/cryptanalysis) counts the same \#E and
  finds the same 30 curves, 2/4/8/16],
)
]

#align(center)[#table(
  columns: 8,
  align: (col, row) => (left,right,right,right,right,right,right,right,).at(col),
  inset: 6pt,
  [level], [curves], [solves], [min], [mean ± 95 % CI], [median], [max],
  [mean ms ± CI],
  [0 crater],
  [2],
  [200],
  [36,433],
  [519,550 ± 39,163],
  [479,542],
  [1,313,611],
  [394.4 ± 31.3],
  [1],
  [4],
  [400],
  [18,197],
  [513,250 ± 26,834],
  [487,539],
  [1,399,367],
  [392.0 ± 20.7],
  [2],
  [8],
  [800],
  [17,768],
  [486,931 ± 18,140],
  [465,983],
  [1,429,985],
  [362.9 ± 14.1],
  [3 floor],
  [16],
  [1,600],
  [26,502],
  [496,402 ± 12,454],
  [466,178],
  [1,586,151],
  [372.3 ± 9.7],
  [all isogenous],
  [30],
  [3,000],
  [17,768],
  [497,666 ± 9,336],
  [467,154],
  [1,586,151],
  [373.9 ± 7.3],
  [isomorphic control],
  [3],
  [300],
  [21,764],
  [478,319 ± 28,790],
  [447,625],
  [1,308,949],
  [366.9 ± 24.1],
)
]

√(πn/2) \= 490,158 (×1.016 for the r \= 32 walk: 498,000). Per-curve
means span 438,599 (v17) to 561,536 (v8). Per-op cost is 736–800 ns on
every curve. Every level’s CI contains √(πn/2) or its r-adding value,
and every pair of level CIs overlaps.

#image("runs/run-20261017T0049Z/volcano.svg")
#image("runs/run-20261017T0049Z/per_curve.svg")
#image("runs/run-20261017T0049Z/per_level.svg")

== Why min and max differ so much
<why-min-and-max-differ-so-much>
Pollard rho is randomized: it walks pseudo-randomly through ⟨P⟩ until
the walk meets a point it has visited before. The step at which that
happens is a birthday-problem collision time:
`P(no collision after k steps) ≈ exp(−k²/2n)`. That gives a mean of
`√(πn/2) ≈ 1.25√n` and a standard deviation of `√((4 − π)n/2) ≈ 0.66√n`,
52 % of the mean. The spread between min and max is that distribution’s
width, not anything about a curve.

#strong[Observed vs the birthday model] (derived from the recorded
`trials.csv`, all isogenous solves):

#align(center)[#table(
  columns: 7,
  align: (col, row) => (left,right,right,right,right,right,right,).at(col),
  inset: 6pt,
  [curve], [solves], [mean (model)], [sd (model)], [1st pct (model)],
  [median (model)], [99th pct (model)],
  [A],
  [2,900],
  [659,263 (654,547)],
  [341,296 (342,147)],
  [71,551 (74,043)],
  [622,025 (614,906)],
  [1,584,239 (1,584,962)],
  [B],
  [4,100],
  [1,090,981 (1,059,102)],
  [567,779 (553,617)],
  [124,423 (119,807)],
  [1,034,344 (994,960)],
  [2,609,159 (2,564,577)],
  [C],
  [3,000],
  [497,666 (490,158)],
  [260,902 (256,217)],
  [59,046 (55,447)],
  [467,109 (460,473)],
  [1,190,999 (1,186,900)],
)
]

- #strong[One curve, 100 solves:] the model’s typical fastest solve is
  `√(2n·ln(101/100))` and its typical slowest `√(2n·ln 101)`, about
  0.14√n and 3.0√n. That is a ~20× spread on every curve, the controls
  included.
- #strong[Pooled extremes:] pools of 1,000+ solves reach further into
  both tails, which is where single minima like 8,477 come from.
- #strong[Per-curve means:] each mean averages 100 such draws, so it
  carries about `0.66√n / 10` of noise: model 34,215 / 55,362 / 25,622
  for A / B / C. The measured spread of the per-curve means is 36,686 /
  62,351 / 26,723, so the curves differ by sampling noise alone.
- #strong[Why the curve doesn’t matter:] an isogenous curve changes the
  coordinates, so the x-hash that drives the walk takes different paths.
  The group has the same order, so the collision-time distribution is
  the same.
- #strong[What a real effect would need:] a curve effect would have to
  move the mean by more than this noise. Detecting an x % difference
  takes roughly `(2 · 0.52 · 100/x)²` solves per arm, about 1,100 for 3
  %.

== Same walk on every curve: the pull-back control
(`run-20261017T0111Z`)
<same-walk-on-every-curve-the-pull-back-control-run-20261017t0111z>
The per-curve differences come from the walk, not the curve. Rho picks
each step, decides distinguished points and keys its table by a hash of
the current point’s x-coordinate. An isogeny or a rescaling changes
those coordinates, so each curve runs a different pseudo-random walk on
the same group.

#strong[Prediction.] If every curve read the root curve’s coordinates
instead, the walks would coincide and every curve would take exactly the
same number of steps for each seed.

#strong[How the walk reads the root’s coordinates (`--walk pullback`).]
- #strong[The canonical point.] Each curve m isogeny steps from the root
reads `[ℓ^(M−m)] · φ̂₁ ∘ … ∘ φ̂_m (R)` on the root curve, which equals
`±ℓ^M φ⁻¹(R)`. Here M is the largest distance among the measured curves,
7 for curve C. - #strong[The dual edge.] Each dual φ̂ is the rational
kernel whose codomain has the parent’s j-invariant, composed with the
isomorphism onto the parent’s model. It is accepted only if
`φ̂(φ(P)) = ±[ℓ]P`. - #strong[Controls.] These pull back by undoing their
rescaling.

#strong[Result] on curve C: 30 curves + 3 controls × 20 seeds.

#align(center)[#table(
  columns: 2,
  align: (col, row) => (left,left,).at(col),
  inset: 6pt,
  [], [result],
  [seeds that took the same number of steps on all 33 arms],
  [#strong[20 / 20]],
  [all 660 solves],
  [correct],
  [steps per seed (sorted)],
  [67,490 · 165,318 · 202,642 · 212,430 · 226,841 · 274,404 · 369,316 ·
  497,854 · 524,385 · 593,906 · 621,975 · 624,229 · 722,575 · 737,697 ·
  745,537 · 755,347 · 767,003 · 839,832 · 903,373 · 1,021,850],
  [mean over seeds],
  [543,700 (√(πn/2) \= 490,158; one draw per seed, sd of the mean ≈
  57,000)],
  [cost per step],
  [6,425 ns (root, chain length 0, ℓ^7 lift) rising to 7,562 ns (chain
  length 7); native walk 757 ns],
)
]

The native walk on the same curves gives per-seed counts that are
uncorrelated across curves (mean correlation +0.002 on A, −0.025 on C).
The pulled-back walk gives the same count everywhere. The step count is
a property of the walk and the instance; the curve only supplies the
coordinates the walk hashes. Pulling back makes every step about 9×
dearer and saves no steps: whichever walk a curve runs, its collision
time follows the same birthday distribution.

#figure([#image("runs/run-20261017T0111Z/per_curve.svg")],
  caption: [
    pull-back per curve
  ]
)

== Where the curve does matter: the j \= 0 crater (`run-20261017T0140Z`)
<where-the-curve-does-matter-the-j-0-crater-run-20261017t0140z>
On a prime field, the one property that can differ between isogenous
curves and change rho’s step count is the automorphism group. A curve
with j \= 0 has six automorphisms; every other curve has two (±1). Rho
can walk on classes of six points instead of two, which predicts √(6/2)
\= √3 ≈ 1.73× fewer steps than the negation-map walk. Such curves exist
only in classes with D\_K \= −3, where j \= 0 sits on the crater. For a
D\_K \= −3 class with a 2-volcano, (−3/2) \= −1, so the crater is the
single j \= 0 vertex, and its three 2-isogenies are permuted by the
order-3 automorphism and land on one curve.

#strong[Curve D.] It was found by
`volcano-rho --j0 --ell 2 --bits 48 --depth 5 --max-cofactor 100000 --seed 20261015`,
which searches y² \= x³ + b with p ≡ 1 (mod 24). It was measured by
`ca volcano … --automorphisms` from aburan28/cryptanalysis at `e4f9e25`;
exact commands are in `runs/run-20261017T0140Z/command.txt`. Every curve
runs `ca`’s plain negation-map rho. The j \= 0 crater curve also runs
`ca_curve_solve`’s automorphism-folded rho on the same instance with the
same seeds.

#align(center)[#table(
  columns: 2,
  align: (col, row) => (left,left,).at(col),
  inset: 6pt,
  [p],
  [173805494607697 (48 bits, ≡ 1 mod 24)],
  [root E],
  [y² \= x³ + 26953713445357, j \= 0],
  [\#E],
  [173805517224804 \= 10524 · n, n \= 16515157471 (prime, 33.9 bits)],
  [trace t],
  [−22617106; D\_K \= −3, 2-volcano depth 5, inert crater],
  [walk],
  [32 curves \= 1 + 1 + 2 + 4 + 8 + 16, 31 edges verified (the whole
  volcano)],
  [rho],
  [32 curves + 3 controls (plain) + the crater curve (folded), 400 seeds
  each: 14,400 solves, all correct; 124.0 s],
)
]

#align(center)[#table(
  columns: 6,
  align: (col, row) => (left,right,right,right,right,right,).at(col),
  inset: 6pt,
  [arm], [curves], [solves], [min], [mean ± 95 % CI], [max],
  [L0, the j \= 0 crater, plain],
  [1],
  [400],
  [32,500],
  [138,688 ± 5,881],
  [323,112],
  [L1],
  [1],
  [400],
  [29,091],
  [140,115 ± 5,964],
  [345,705],
  [L2],
  [2],
  [800],
  [27,442],
  [141,521 ± 4,151],
  [351,803],
  [L3],
  [4],
  [1,600],
  [24,572],
  [141,185 ± 2,972],
  [360,062],
  [L4],
  [8],
  [3,200],
  [24,150],
  [140,314 ± 2,092],
  [370,427],
  [L5, the floor],
  [16],
  [6,400],
  [26,545],
  [139,916 ± 1,498],
  [416,873],
  [all isogenous, plain],
  [32],
  [12,800],
  [24,150],
  [140,242 ± 1,053],
  [416,873],
  [isomorphic controls, plain],
  [3],
  [1,200],
  [23,098],
  [142,544 ± 3,447],
  [390,953],
  [#strong[j \= 0 crater, automorphism-folded]],
  [1],
  [400],
  [10,573],
  [#strong[74,024 ± 3,560]],
  [231,754],
)
]

- #strong[Operations (measured):] the folded walk needs #strong[1.89 ±
  0.09×] fewer group operations than the plain walks, 140,242 against
  74,024. On the crater curve alone, same instance and seeds, it’s
  138,688 against 74,024 (1.87×). Every plain level is flat at about
  140k, as on curves A to C.
- #strong[Against the √3 prediction:] the interval excludes 1.73. Both
  arms sit above their ideal values: plain at 1.23× √(πn/4) \= 113,890,
  folded at 1.13× √(πn/12) \= 65,755. The extra is `ca`’s walk overhead
  (its r \= 1024 negation walk with fruitless-cycle handling), which
  falls more on the plain walk. The algorithmic factor is √3; this
  implementation’s measured factor is 1.89.
- #strong[Wall time (shared host, indicative):] 6.52 ms against 8.65 ms
  per solve, #strong[1.33×]. A folded step costs 88 ns against 62 ns,
  because it picks a canonical representative among six points.
- #strong[What this means for the volcano (derived):] an instance on any
  of the 31 lower curves can be carried up to the j \= 0 crater by at
  most five 2-isogeny evaluations, which costs nothing next to 10⁵ group
  operations, and solved there about 1.9× cheaper in operations (1.3× in
  time). This is the one place in a prime-field volcano where the curves
  are not interchangeable. It exists only for D\_K \= −3 (or −4, j \=
  1728: 2× fewer, √2 beyond negation) classes, and a random curve’s
  class has |D\_K| ≈ p, so there is nothing to climb to.
- #strong[A recorded-output note:] the `ca` build that ran wrote the
  folded arm’s mean into the crater vertex of `volcano.diagram.json`,
  giving a duplicate key. `plot_ca.py` takes vertex values from
  `results.json`’s plain arms, and `ca` is fixed (`31b72d2`).

#image("runs/run-20261017T0140Z/volcano.svg")
#image("runs/run-20261017T0140Z/per_level.svg")
#image("runs/run-20261017T0140Z/crater_cdf.svg")

== The j \= 1728 crater (`run-20261017T0305Z`)
<the-j-1728-crater-run-20261017t0305z>
The other prime-field case is j \= 1728, y² \= x³ + ax over p ≡ 1 (mod
4), whose four automorphisms give classes of four points. That predicts
√(4/2) \= √2 ≈ 1.41× fewer steps than the negation-map walk. Its class
has D\_K \= −4, which is ramified at 2: the crater is the single j \=
1728 vertex, and its horizontal 2-isogeny is the degree-2 endomorphism
`1 + i`, a loop. The two other 2-kernels are swapped by the automorphism
and land on one curve below.

#strong[Curve E.] It was found by
`volcano-rho --j1728 --ell 2 --bits 48 --depth 5 --max-cofactor 100000 --seed 20261025`
and measured by `ca volcano … --automorphisms` (aburan28/cryptanalysis
`5841fd7`); exact commands are in `runs/run-20261017T0305Z/command.txt`.

#align(center)[#table(
  columns: 2,
  align: (col, row) => (left,left,).at(col),
  inset: 6pt,
  [p],
  [213656134960913 (48 bits, ≡ 1 mod 8)],
  [root E],
  [y² \= x³ + 122917133887808x, j \= 1728],
  [\#E],
  [213656154587108 \= 116 · n, n \= 1841863401613 (prime, 40.7 bits)],
  [trace t],
  [−19626194; D\_K \= −4, 2-volcano depth 5, ramified crater],
  [walk],
  [32 curves \= 1 + 1 + 2 + 4 + 8 + 16, 31 edges verified (the whole
  volcano)],
  [rho],
  [32 curves + 3 controls (plain) + the crater curve (folded), 400 seeds
  each: 14,400 solves, all correct; 1,147.7 s],
)
]

#align(center)[#table(
  columns: 7,
  align: (col, row) => (left,right,right,right,right,right,right,).at(col),
  inset: 6pt,
  [arm], [curves], [solves], [min], [mean ± 95 % CI], [max], [mean ms],
  [L0, the j \= 1728 crater, plain],
  [1],
  [400],
  [156,527],
  [1,360,140 ± 66,271],
  [3,954,060],
  [],
  [L1],
  [1],
  [400],
  [137,160],
  [1,333,500 ± 65,462],
  [3,891,950],
  [],
  [L2],
  [2],
  [800],
  [127,700],
  [1,353,380 ± 44,854],
  [4,141,490],
  [],
  [L3],
  [4],
  [1,600],
  [135,195],
  [1,332,150 ± 30,476],
  [3,972,090],
  [],
  [L4],
  [8],
  [3,200],
  [122,902],
  [1,347,930 ± 22,393],
  [4,807,580],
  [],
  [L5, the floor],
  [16],
  [6,400],
  [113,777],
  [1,334,080 ± 15,370],
  [3,946,280],
  [],
  [all isogenous, plain],
  [32],
  [12,800],
  [113,777],
  [1,339,300 ± 11,010],
  [4,807,580],
  [77.0],
  [isomorphic controls, plain],
  [3],
  [1,200],
  [162,329],
  [1,320,980 ± 35,656],
  [4,098,680],
  [76.1],
  [#strong[j \= 1728 crater, automorphism-folded]],
  [1],
  [400],
  [161,794],
  [#strong[1,018,110 ± 47,342]],
  [2,604,700],
  [176.9],
)
]

#strong[Curves D and E side by side]

#align(center)[#table(
  columns: 3,
  align: (col, row) => (left,right,right,).at(col),
  inset: 6pt,
  [], [j \= 0 (curve D)], [j \= 1728 (curve E)],
  [automorphism order],
  [6],
  [4],
  [predicted operation speed-up],
  [√3 \= 1.73],
  [√2 \= 1.41],
  [measured operation speed-up],
  [#strong[1.89 ± 0.09]],
  [#strong[1.32 ± 0.06]],
  [plain rho / its ideal √(πn/4)],
  [1.23],
  [1.11],
  [folded rho / its ideal √(πn/(2m))],
  [1.13],
  [1.20],
  [ns per operation, plain / folded],
  [62 / 88],
  [57 / 174],
  [wall-time speed-up],
  [1.33×],
  [#strong[0.44× (2.3× slower)]],
  [`ca` code path for the folded arm],
  [`glv_covariant_rho_solve` (tuned j \= 0 walk)],
  [`glv_rho_solve` (generic)],
)
]

- #strong[Measured:] plain rho is again flat on every level, and the
  controls match. The folded crater arm needs 1.32 ± 0.06× fewer
  operations, and its interval excludes √2. It is 2.3× slower in wall
  time.
- #strong[Why the two craters differ (derived from reading `ca`’s
  source, not measured separately):] the comparison is between
  implementations as much as automorphism groups. `ca`’s default build
  routes j \= 0 to a dedicated covariant walk
  (`CA_J0_RHO_COVARIANT_WALK=ON`), which labels the six-point class
  without canonicalising the point. j \= 1728 runs the generic
  `glv_rho_solve`, which on every step:
  - applies the automorphism three times and hashes four points to pick
    the class representative, consistent with the 3× cost per operation;
  - walks with at most r \= 32 multipliers (plain rho uses 1,024) and
    escapes only 2-cycles by look-ahead. Longer fruitless cycles wait
    for a 24·2^dp-step abandon and restart, which is consistent with
    landing 20 % above the folded ideal.

  So curve E measures this implementation, not the √2 limit. A j \= 1728
  counterpart of the covariant walk is the obvious follow-up in `ca`.
- #strong[Implication for the volcano (derived):] the j \= 1728 crater
  lowers the operation count as theory says it should, but with this
  solver, moving there costs wall time. On curve D the same move gained
  1.33× in wall time.

#image("runs/run-20261017T0305Z/volcano.svg")
#image("runs/run-20261017T0305Z/per_level.svg")
#image("runs/run-20261017T0305Z/crater_cdf.svg")

== Generator and k: does the instance itself matter?
(`run-20261017T0415Z-instance`)
<generator-and-k-does-the-instance-itself-matter-run-20261017t0415z-instance>
Every volcano run above solved one fixed (P, k) per curve and varied
only the walk’s seed. This run isolates the instance. On the root of
curve C it runs five arms of 1,000 solves each, all with the same binary
(`volcano-rho --instance-study`, native 32-adding DP walk, no negation
map), and redraws a different part of the problem per trial. The
protocol and the decision rule were frozen before the run
(#link("runs/run-20261017T0415Z-instance/PROTOCOL.md")[`PROTOCOL.md`]).
An arm counts as different from `walk` only if |z| \> 2.58 or KS p \<
0.0125.

#align(center)[#table(
  columns: 10,
  align: (col, row) => (left,left,right,right,right,right,right,right,right,right,).at(col),
  inset: 6pt,
  [arm], [varies per trial], [min], [mean ± 95 % CI], [median], [max],
  [sd/mean], [vs `walk`: Δ mean ± 95 % CI], [z], [KS p],
  [`walk`],
  [walk seed (P, k fixed)],
  [17,176],
  [487,658 ± 16,110],
  [454,066],
  [1,567,146],
  [0.53],
  [],
  [],
  [],
  [`k`],
  [k uniform in \[1, n)],
  [22,331],
  [480,254 ± 16,127],
  [451,165],
  [1,468,003],
  [0.54],
  [−7,404 ± 22,795],
  [−0.64],
  [0.75],
  [`generator`],
  [P, a fresh point of order n],
  [30,819],
  [492,522 ± 16,136],
  [460,232],
  [1,530,409],
  [0.53],
  [+4,864 ± 22,802],
  [+0.42],
  [0.68],
  [`all`],
  [walk seed, k and P],
  [23,729],
  [496,427 ± 16,165],
  [461,552],
  [1,561,929],
  [0.53],
  [+8,769 ± 22,822],
  [+0.75],
  [0.64],
  [`small_k`],
  [k uniform in \[1, 2^16\]],
  [10,099],
  [487,206 ± 15,362],
  [459,656],
  [1,528,969],
  [0.51],
  [−452 ± 22,260],
  [−0.04],
  [0.68],
)
]

Reference: √(πn/2) / √(1 − 1/32) \= 498,001 iterations.

- #strong[Measured:] no arm differs from `walk` in mean or in shape.
  Every |z| is below 0.8 and every KS p is above 0.6. Each interval on
  the mean difference excludes an effect larger than about ±4.7 % of the
  mean. All 5,000 solves returned their own k, and replaying the first
  five trials of every arm matches 25/25 iteration counts.
- #strong[Measured:] a small secret is not cheaper. k ≤ 2^16, about
  2^−21 of the group, costs the same as a uniform k, and cost shows no
  trend against k
  (#link("runs/run-20261017T0415Z-instance/instance_vs_k.png")[`instance_vs_k`]).
- #strong[Derived:] the generator doesn’t matter because any other
  generator is P′ \= mP: solving against P′ is the same problem up to a
  known factor m. The walk mixes P and Q \= kP with random coefficients,
  so the size of k never enters; only the collision time of a
  pseudo-random walk does. That differs from kangaroo/lambda, whose cost
  grows with the interval known to hold k.
- #strong[Derived:] with the walk seed fixed, a new k or P still changes
  every walk point, so each trial is a fresh draw from the same
  collision-time distribution. The four arms with redraws cannot be told
  apart from the walk-seed arm.

#figure([#image("runs/run-20261017T0415Z-instance/instance_arms.png")],
  caption: [
    every solve per arm
  ]
)

#figure([#image("runs/run-20261017T0415Z-instance/instance_cdf.png")],
  caption: [
    cumulative distributions
  ]
)

#figure([#image("runs/run-20261017T0415Z-instance/instance_vs_k.png")],
  caption: [
    iterations against k
  ]
)

== Reading it
<reading-it>
- #strong[Observed:] level means differ by at most 63,000 iterations (L0
  vs L2), and the 95 % CIs of L0 and L2 overlap. The CIs of L0, L2 and
  L3 contain √(πn/2); L1’s lower bound (655,241) sits 700 above it, with
  four level intervals and no correction for multiple comparisons.
  Single solves range over ~250× (8,477 to 2,099,707) on every level and
  on the control, so min and max are the extremes of a wide distribution
  and say little about a curve. The mean is the meaningful statistic.
- #strong[Curve B] (observed): the level means lie within 1.7 % of each
  other (1,081,516 to 1,099,508), and every level’s CI overlaps every
  other’s. All of them, and the control, sit about 3 % above √(πn/2).
  The pooled mean, 1,090,981 ± 17,380, excludes it. Curve A’s pooled
  mean sits 0.7 % above (inside its CI). A known heuristic (derived, not
  measured here) is that an r-adding walk needs about 1/√(1 − 1/r) times
  a truly random walk’s steps: ×1.016 for r \= 32, or 1,076,000 on B and
  665,000 on A. Both pooled CIs contain those values. The offset is a
  property of the walk, shared by every curve and the control, not of
  any stratum.
- #strong[Interpretation (derived, not measured):] an ℓ-isogeny with ℓ ∤
  n is a group isomorphism ⟨P⟩ → ⟨φ(P)⟩, so generic rho sees the same
  group of order n on every stratum; what changes between levels is the
  endomorphism ring, which generic rho does not use. The per-curve
  differences are what re-drawing the walk’s random partition gives, as
  the isomorphic control (same curve, new coordinates) shows. Walking
  the volcano is not a way to make this generic DLP cheaper.
- #strong[Gap to matched rho:] this #emph[is] rho on every arm; no
  index-calculus arm is involved, so there is no IC baseline to improve
  and no §5a band to rank across.

== Limits
<limits>
- One instance per curve, three curves:
  - A: 40-bit p, ℓ \= 3, depth 3, ramified crater with 2 vertices, so
    its L0 has only 200 solves;
  - B: 40-bit p, ℓ \= 3, depth 3, split crater, truncated at the walk
    radius;
  - C: 48-bit p, ℓ \= 2, depth 3, ramified crater, the whole volcano.
- Not an ecbench session: the counted unit is this tool’s group-addition
  count; no isolation (shared 4-vCPU VM, 4 threads, jobs interleaved),
  so wall times are indicative only. Per-op time includes affine
  inversion by `i128` extended Euclid.
- Curves measured: A 29 of 80 (the whole volcano), B 41 of the 323
  walked (root plus up to 10 per level in BFS order), C all 30.
- No negation map and no endomorphism-accelerated rho; a curve at the
  crater with a larger endomorphism ring would only help an algorithm
  that uses that ring (GLV-style speedups exist only for tiny |D\_K|,
  and this D\_K \= −1482974187 is not tiny).

== Run it on your own curve
<run-it-on-your-own-curve>
```sh
cargo build --release --bin volcano-rho
target/release/volcano-rho --curve P,A,B [--order N] [--ell 3] [--seed S] [--trials 100] --out OUT
target/release/volcano-rho --curve P,A,B --engine ca --ca /path/to/ca ...   # cryptanalysis's C rho
python research/isogeny-volcano-rho/plot.py OUT    # optional: the standard graphs
```

The cryptanalysis repository has the same measurement built into its own
binary, with no Rust involved:
`ca volcano --p P --a A --b B [--trials T] [--diagram FILE]`
(`ca_volcano.h`). It reproduces this study’s curve A volcano exactly (80
curves, 2 / 6 / 18 / 54).

`y² = x³ + Ax + B` over `F_P`, any prime `P < 2^62`. Without `--order`
the tool counts points (baby-step giant-step). It checks the order
against random points, takes `n` the largest prime factor, draws a
generator (cofactor × a seeded random point) and a secret `k` from
`--seed`, and reads the volcano depth `v_ℓ(f)` and the crater shape from
`D = t² − 4p`. Then it walks, carries `(P, kP)` to every curve and
measures. A curve whose conductor ℓ does not divide has depth 0, so the
walk covers crater curves only. Outputs go to `--out`: `results.json`,
`trials.csv`, `instances.jsonl` (every curve’s instance in `ca` format,
ready for `ca solve`) and the standard graphs.

- #strong[Engines.] `native` is this binary’s single-walker 32-adding
  walk with distinguished points. The DP rate defaults to
  `min(8, log₂n/2 − 6)` bits, so small groups keep every point. With
  every point distinguished, the x-keyed table also catches `±X`
  matches, so expect about √(πn/4). `--engine ca` runs
  `ca solve --alg rho` from aburan28/cryptanalysis on each curve, with
  its defaults (negation map on) plus any `--ca-arg`. Iterations are
  then `ca`’s reported `ops` and times its reported `seconds`.
- #strong[Cost.] About √n group operations per solve: n ≈ 2^40 takes a
  second or so per solve on the native engine. Choose `--trials` with
  that in mind.
- #strong[Determinism.] The same command and seed give the same curves,
  instance and iteration counts; wall times vary.
- #strong[Instance study.] `--instance-study` skips the walk and
  measures the root curve alone, `--trials` solves per arm. The five
  arms redraw the walk seed, k, the generator, all three, or k from \[1,
  2^16\]. It writes `instance_study.json` and `instance_trials.csv`;
  `plot_instance.py OUT` draws them.

== Reproduce
<reproduce>
```sh
cargo build --release --bin volcano-rho
./target/release/volcano-rho --bits 40 --ell 3 --depth 3 --seed 20261016 \
    --trials 100 --per-level 10 --controls 3 --out OUT
./target/release/volcano-rho --bits 40 --ell 3 --depth 3 --crater split --seed 20261017 \
    --trials 100 --per-level 10 --controls 3 --out OUT_B
./target/release/volcano-rho --bits 48 --ell 2 --depth 3 --max-cofactor 4096 --seed 20261018 \
    --trials 100 --per-level 16 --controls 3 --out OUT_C
./target/release/volcano-rho --bits 48 --ell 2 --depth 3 --max-cofactor 4096 --seed 20261018 \
    --trials 20 --per-level 16 --controls 3 --walk pullback --out OUT_PULLBACK
python plot.py OUT     # matplotlib, plotting only
```

The curve, volcano, instance and every solve’s iteration count are a
deterministic function of `--seed`; wall times are not. Run from the
then-uncommitted `src/bin/volcano_rho.rs` on top of `9b09d308e`
(committed in this PR with formatting and lint-only changes), Intel Xeon
\@ 2.80 GHz, 4 threads, 342.5 s for the rho phase.

Curve B ran on the committed binary (with `--crater`), on the same host,
in 754.0 s for the rho phase.

Curve C ran on the committed binary (with the `p ≡ 1 mod 8` draw for ℓ
\= 2) in 310.7 s for the rho phase. Replaying its first five seeds
reproduces 165/165 iteration counts, and curve A still replays 160/160
on this binary.

#strong[Replay.] The committed binary (formatting and lint fixes only)
re-run with `--trials 5` reproduces 160/160 iteration counts of the
first five seeds of every arm exactly.

Report PDF: #link("report.pdf")[`report.pdf`].
