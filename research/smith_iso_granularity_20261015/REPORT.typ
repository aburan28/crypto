// Results report. Every number is read from data/*.json; nothing here
// computes a result (rounding is presentation only).
#set page(paper: "a4", margin: (x: 1.6cm, y: 1.5cm), numbering: "1")
#set text(size: 9pt)
#set par(justify: true, leading: 0.55em)
#set heading(numbering: "1.")
#show heading: set text(size: 10.5pt)
#let s7 = json("/data/summary_p7.json")
#let s11 = json("/data/summary_p11.json")
#let vd = json("/data/verdict.json")
#let f(x, n: 4) = str(calc.round(x, digits: n))
#let tbl(..args) = text(size: 7.8pt, table(inset: 3pt, stroke: 0.4pt + luma(170), ..args))
#align(center)[#text(size: 13pt, weight: "bold")[Smith's class-level census, counted by isomorphism class] \
`research/smith_iso_granularity_20261015` · 2026-10-15 · protocol `PROTOCOL.md` frozen at `05f3d9203` before any data]

= Question
The 2026-10-14 census (`research/smith_prop24_native_20261014`) found that a tractable $L$-class usually holds a
model on which Smith's construction fully succeeds, but that the class-level share $C$ lies *below* an independence
null that scattered successes over census *models*. The open item in
`problems/smith-kernel-isogeny-class-walk.md` asked whether the deficit is a structural invariant of the class or the
sparsity of small classes. This study tests the simplest sparsity explanation: `full_success` is a property of a curve
up to $FF_p$-isomorphism, and one curve accounts for many census models, so the model-level null over-counts the
independent coins per class.

= Verdicts (frozen wording)
#tbl(columns: (11%, 15%, 30%, 30%, 14%),
  [*prime*], [*observed $C$*], [*iso-level null: median [2.5%, 97.5%] (999 perms)*],
  [*stratified by orbit size too (secondary)*], [*verdict*],
  ..(s7, s11).map(s => ([$p = #s.p$], [#f(s.C)],
    [#f(s.null_primary.median) [#f(s.null_primary.lo), #f(s.null_primary.hi)]; $P("null" <= C) = #f(s.null_primary.frac_le_observed, n: 3)$],
    [#f(s.null_secondary.median) [#f(s.null_secondary.lo), #f(s.null_secondary.hi)]],
    [*#s.verdict*])).flatten(),
  [overall], [], [], [], [*#vd.overall*],
)

#figure(image("/figures/iso_granularity.svg", width: 100%),
  caption: [(a) Observed $C$ against the model-level null of 2026-10-14 (square; its 95% range is a single value) and
  the iso-level null of this study (circle, 95% range); dotted lines are the direct share $D$ and the ceiling $U$.
  (b) Share of tractable $L$-classes with no `full_success` isomorphism class, by the number $k$ of tractable
  isomorphism classes they hold; bars observed, ticks the expectation under independent isomorphism classes with the
  per-stratum success rate. Source `figures/plot.py`, data `data/summary_p*.json`.])

= Measured
#tbl(columns: 3,
  [*quantity (exact)*], [*p = 7*], [*p = 11*],
  [census models / isomorphism classes / $L$-classes], [#s7.models / #s7.iso_classes / #s7.l_classes],
  [#s11.models / #s11.iso_classes / #s11.l_classes],
  [tractable isomorphism classes / with `full_success`], [#s7.tractable_iso_classes / #s7.full_iso_classes],
  [#s11.tractable_iso_classes / #s11.full_iso_classes],
  [median tractable isomorphism classes per tractable $L$-class], [#s7.median_tractable_iso_classes_per_tractable_l_class],
  [#s11.median_tractable_iso_classes_per_tractable_l_class],
  [$D$ / $C$ / $U$ (reproduce 2026-10-14)], [#f(s7.D) / #f(s7.C) / #f(s7.U)], [#f(s11.D) / #f(s11.C) / #f(s11.U)],
  [tractable $L$-classes with no `full_success`: observed / expected (independent iso classes)],
  [#s7.tractable_l_classes_without_full.observed / #f(s7.tractable_l_classes_without_full.expected_independent_iso, n: 1)],
  [#s11.tractable_l_classes_without_full.observed / #f(s11.tractable_l_classes_without_full.expected_independent_iso, n: 1)],
)

#include "findings.typ"

= Validation (controls of the protocol; all must pass)
#tbl(columns: (auto, 1fr),
  [*control*], [*result*],
  [C1: per-$L$-class model, tractable and `full_success` counts equal `census_p*_classes.csv`],
  [$p = 7$: #s7.controls.c1_mismatched_rows of #s7.controls.c1_census_rows rows differ; $p = 11$:
   #s11.controls.c1_mismatched_rows of #s11.controls.c1_census_rows differ.],
  [C2: `full_success`, stratum and $L$-key constant on every isomorphism class (every model checked)],
  [inconsistent models: #s7.controls.c2_inconsistent_models ($p = 7$), #s11.controls.c2_inconsistent_models ($p = 11$);
   orbit images outside the squarefree census: #s7.controls.c2_orbit_hits_nonsquarefree, #s11.controls.c2_orbit_hits_nonsquarefree.],
  [C3: orbits partition the census], [unvisited models: #s7.controls.c3_unvisited_models, #s11.controls.c3_unvisited_models.],
  [instrument], [#if s7.instrument_valid and s11.instrument_valid [*valid*] else [*invalid*]],
)
Orbits are enumerated exactly: every non-Weierstrass point of $PP^1(FF_p)$ is sent to infinity, the leading
coefficient is normalised by a square (its class decides the twist), the octic is depressed and $x -> u x$ is
applied for every $u in FF_p^*$; this covers all of $"PGL"_2(FF_p)$ modulo translations. C2 is the check that
`full_success` is a curve invariant on every one of the census models, not only on sampled ones.

= Scope
Exact enumeration at $p = 7, 11$, degree-8 models only. Nothing here builds, prices or bounds an isogeny walk, or
says anything at cryptographic sizes or in characteristic 2. A rejected granularity hypothesis means the deficit is
not explained by counting curves instead of models; it does not identify the structural cause.

= Reproduce
```
cargo build --release --example smith_iso_granularity
B=target/release/examples/smith_iso_granularity; D=research/smith_iso_granularity_20261015/data
$B --p 7 --out $D; $B --p 11 --out $D        # 17 s and 15 min on 2 vCPU
$B summarize --dir $D 7 11
python research/smith_iso_granularity_20261015/figures/plot.py   # layout only
```
Source: B. Smith, J. Cryptology 22 (2009), arXiv:0806.2995. Prior: `research/smith_prop24_native_20261014/`.
