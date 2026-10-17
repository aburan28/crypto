#set page(paper: "a4", margin: 16mm)
#set text(size: 10pt)
#set heading(numbering: "1.")

= Prebuilt factor bases at degree 83

*Review date:* 2026-10-16  \
*Primary curve:* `icv1-f2m83-tm6151469093347-debefd74`  \
*Question:* construction and reuse of symmetric bases for two, three and four summands.

#image("symmetry-flow.svg", width: 100%)

The archived native sweep already constructed and replayed 57 signed-Frobenius
bases, including a 54-object panel with 45 distinct point sets. Their
content-addressed objects were stored and fetched from S3. Each retained full
orbit contributes 166 labeled subgroup points but one relation column.
The newly frozen local protocol measures construction and checked import
separately; its run receipts are in `RESULT.md`.

== Exact structural screen

For a base of $B$ distinct points and a uniform nonidentity subgroup target,
the full-smooth hit probability is at most
$min(1, binom(B + m - 1, m) / (r - 1))$.
This is an upper bound, not observed yield.

#table(
  columns: (auto, auto, auto, auto, auto),
  inset: 5pt,
  align: (right, right, right, right, right),
  [*Orbit columns*], [*Points*], [*$m=2$*], [*$m=3$*], [*$m=4$*],
  [64], [10,624], [2.334e-17], [8.268e-14], [2.197e-10],
  [256], [42,496], [3.735e-16], [5.290e-12], [5.621e-8],
  [600], [99,600], [2.051e-15], [6.811e-11], [1.696e-6],
)

The dimension-12 projected source closure is another verified stored base:
332,166 points in 2,001 signed-Frobenius columns. Its exact four-summand
count/order ratio is 0.000209791390. The current materialized root-index
design requires 22,945,941,707 states at K=16,627, the first primary-curve
orbit count where the four-summand counting ceiling reaches one. These
counts do not price a solver or establish actual support.

== Interpretation and next gate

Prebuilding and verified reuse are demonstrated at the factor-base stage.
On the primary 81-bit subgroup, the archived two-summand probes found no
relation in their bounded fixed-point panel; the three-summand attempt was
capped, and there is no matched natural four-summand relation/rank panel for
these bases. The separate degree-83 `a=1` four-summand cold measurements use
a 53-bit subgroup and cannot be transferred to this curve.

The next discriminating comparison uses the same public primary-curve query
corpus for the K=2,048 explicit-orbit object and the 2,001-column structured
closure. It requires a full-width native four-summand oracle, group-lift and
rank replay, then matched cold accounting with A/A and rho controls. The
checkable obligation is recorded in
`problems/n83-prebuilt-base-quality-m2-m4.md`.

== Sources

- #link("https://github.com/aburan28/crypto/blob/main/research/koblitz_n83_factor_base_sweep_20261008/RESULT.md")[N83 stored-base sweep and receipts]
- #link("https://github.com/aburan28/crypto/blob/main/research/koblitz_n83_factor_base_sweep_20261008/SUPPORT_MOMENTS.md")[Exact multiset support screen]
- #link("https://github.com/aburan28/crypto/blob/main/research/koblitz_n83_factor_base_sweep_20261008/S5_CAPACITY.md")[Four-summand index capacity]
- #link("https://github.com/aburan28/crypto/blob/main/research/f6_n83_mixed_base_20261004/ORBIT_RESULT.md")[Dimension-12 closure census]
- #link("https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf")[Galbraith, Granger, Merz and Petit: Frobenius-invariant bases]
