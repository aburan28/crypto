# Protocol I-3: model-dependent screens over a whole isogeny class, and whether large-degree neighbours look different

Frozen 2026-10-08, before any instrument is built.  **Boundary; engineering
if a screen proves cheap enough to run over a class.**  `S`, end-to-end
cost and speedup are **unset**.

## Derivation (stated before measuring)

Three things a non-generic method reads off a curve's *equation* are not
isogeny invariants:

1. **The GHS magic number** `m(b)` (binary, extension degree `N = n·l`):
   `dim_{F_2} span{σ^i(√b)}`.  Galbraith–Hess–Smart (2002) extended GHS
   precisely by walking the class to a neighbour with small `m`; Menezes–
   Teske quantified how many classes over `F_{2^{155}}`, `F_{2^{161}}`
   contain such a curve.  `audit_curve` in `ec_trapdoor.rs` computes `m`
   for every factorisation of `N`.
2. **Cover existence and genus** for `E/F_{q^n}` (`curve_cover_check`,
   `jv_cover`): the Joux–Vitse family `y² = x(x−α)(x−σ(α))` is a model
   condition, which is why #1556 walked depth-1 neighbours looking for it.
3. **Subfield-definedness** of `j` or of the coefficients, and `|Aut|`.

Hypothesis H3a (positive control): over composite-degree binary fields
the distribution of `min_{(n,l)} m(b)` over a class is non-degenerate and
its minimum reproduces the literature's weak-curve rates.
Hypothesis H3b (prime-degree control): over `F_{2^n}` with `n` prime,
`m(b) = n` for every ordinary curve, so the screen is constant.
Hypothesis H3c (the large-`ℓ` question): the screen values of a class's
large-`ℓ` neighbours are drawn from the same distribution as those of its
small-`ℓ`-reachable members, because a horizontal large-`ℓ` neighbour *is*
a small-`ℓ`-reachable member (class-group generation).  A difference
would mean either a reach the class-group argument misses (impossible for
horizontal edges under GRH) or that the sample of large-`ℓ` edges is
biased by how they were found (the rational-torsion window selects by
eigenvalue order, which correlates with nothing in the equation).

## Instrument (Rust, follow-on PR)

- `isogeny_walk walk` (prime field) and a `binary_isogeny` class walk
  (binary; the `Φ_ℓ mod 2` root path exists) to enumerate the class at
  toy size to exhaustion where the class number allows, else to a fixed
  node budget with the frontier recorded.
- A `screen` subcommand writing one row per node: slug, `j`, edge path,
  `min m(b)` with its factorisation, cover status (`curve_cover_check`
  certificate or `unsupported`), `j ∈ F_{q^d}` for proper `d | n`,
  coefficients in a proper subfield, `|Aut|`.  Every unknown stays `null`.
- Large-`ℓ` neighbours tagged by how they were found (rational-torsion,
  `Φ_ℓ`, class-group route, vertical).

## Frozen inputs

| item | value |
|:--|:--|
| binary composite | three classes each over `F_{2^{12}}`, `F_{2^{15}}`, `F_{2^{21}}`, registered in the follow-on PR |
| binary prime | `icv1-f2m13-t181-515ee569`, `icv1-f2m19-t797-b6cf2467`, `icv1-f2m23-t5197-69e76b73` |
| prime field | the I-1 classes at `p ≈ 2^{20}` and `2^{24}` |
| node budget | full class where `h ≤ 20,000`, else 20,000 nodes |
| large `ℓ` | the I-2 sweep, tagged by discovery route |
| seed | 20261013 |

## Predictions (pass/fail)

- **Q1 (positive control).**  Over each composite-degree binary class the
  screen `min m(b)` takes at least two distinct values, and the fraction of
  classes whose minimum gives genus `≤ 2^{m−1}` with `m ≤ 4` is within a
  factor 2 of the Menezes–Teske rate for that `(n, l)` shape.
- **Q2 (negative control).**  Over every prime-degree binary class the
  screen is constant: `m(b) = n` on every node.
- **Q3 (large `ℓ`).**  For every class and every screen feature, a
  two-sample Kolmogorov–Smirnov test between large-`ℓ` neighbours and
  small-`ℓ`-reachable nodes does not reject at `p < 0.05` after
  Benjamini–Hochberg over the features.
- **Q4 (prime field).**  Over prime-field classes every screen except
  `|Aut|` is constant (`null` cover, no subfield structure), and `|Aut| >
  2` occurs only on the `D ∈ {−3, −4}` classes.

## Decision rule (registered)

- Q1–Q4 pass: **boundary**.  Screens reproduce the known weak-curve
  mechanism where it exists, return constants where theory says they
  must, and large-`ℓ` neighbours add no reach.  The screen table becomes
  step 2 of the methodology.  If the screen's cost per node is under
  `10⁴` field operations, it is also **engineering**: cheap enough to run
  over a class before any DLP work.
- Q1 fails: the instrument disagrees with the literature; a bug report
  against `audit_curve` or the walk, not a result.
- Q3 fails on a feature: the discovery route is biased or the class-group
  argument is being misapplied (a vertical edge mislabelled horizontal).
  Trace the cells; if the difference survives with correctly labelled
  edges on a fresh class, it is a **reproducible unexplained anomaly** and
  gets its own protocol.  It cannot be classed higher from here.
- Q4 fails: same treatment as Q3; a prime-field model screen that varies
  is exactly what the directory is for and is held to every control.

## Stop condition and inadmissible moves

Bounded: the classes above, one walk each, one screen pass.

Inadmissible: treating a screen value as a DLP cost; counting a node
twice through two routes; mixing `F_{2^n}` representations within a class
without the conversion recorded; reporting the class minimum of a screen
without the route cost to reach it (I-2).
