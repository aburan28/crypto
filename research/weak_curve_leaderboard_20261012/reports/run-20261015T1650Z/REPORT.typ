#set page(paper: "a4", margin: 18mm)
#set text(font: "Libertinus Serif", size: 10pt)
= Weighted relation support for geometry-guided curve queries

A finite-seed nonuniform query law can place at most $binom(B+m-1,m) M$ seeds on the unordered positive $m$-sum support of a fixed base of $B$ usable classes, if each accepted target has at most $M$ preimage seeds. Lean checks the full-support and accepted-set bounds.

#image("WEIGHTED_SUPPORT.svg", width: 100%)

== Proof and probability interpretation

Let $w(a)$ count seeds mapping to target class $a$. The fixed-base support $S$ has at most $C=binom(B+m-1,m)$ classes. On an accepted set $E$, $sum_(a in E ∩ S) w(a) <= abs(E ∩ S) M <= C M$. Divide by the attempted seed count for a per-attempt bound, or by accepted seed mass for a conditional bound. The theorem requires an exact maximum-preimage certificate and no uniform-target premise.

== Curve-screening consequence

For the P-192 Type I 16-class, three-sum base, $C=816$. Its earlier character law admits 224 support targets and has a separately certified uniform accepted-target hit bound below $2^(-182)$. The new result does not assign a preimage cap to a different sampler; that cap, verified relation yield, useful rank, and target descent are concrete gates for each guided collector.

The source-linked report is `research/weak_curve_leaderboard_20261012/reports/run-20261015T1650Z/REPORT.md`. Lean module: `WeightedSupport.lean` (Lean 4.34.1, pinned Mathlib). The SVG is the editable diagram source.
