# Protocol: which isogeny classes over F_{p⁶} hold a Joux–Vitse weak curve, and where in the 2-volcano the weak curves sit

Frozen 2026-10-07, before the census is built or run.  **Stage diagnostic,
toy sizes.**  `S`, end-to-end cost and speedup are **unset**.  Nothing here
concerns a prime-field or deployed curve: the weak class exists only over
an extension field of degree divisible by 3, and the route it serves is
[JV12]'s, reproduced in the cover ledger on branch
`research/jv-cover-end-to-end-20261007`
(`research/notes/index-calculus/RESEARCH_COVER_DECOMPOSITION_LEDGER.md`,
§§13, 17, 18).  Status: **RUN** 2026-10-07; results and the reading against every prediction are in `README.md`.  `p = 19` stopped at the wall limit (partial).

## What is already known, and the gap

- The weak class (ledger §13): a curve with full rational 2-torsion over
  `F_{q³}`, `q = p²`, is weak exactly when one of its three 2-torsion
  cross-ratios has `N_{F_{q³}/F_q} = 1`; density `3/q` among
  full-2-torsion curves; an isomorphism invariant, and twist-invariant,
  since a quadratic twist keeps the abscissae of the 2-torsion.
- The walk (ledger §17.5): random 2,3-isogeny components reach 5 to 180
  curves of a class of about `q^{3/2}`; 60 to 90% of walks end in a
  component with no weak curve; success falls with `p`.  "The cause is
  where the weak curves are."
- The sampled census (ledger §17.5.3, §18.5): weak curves fall in far
  fewer classes than random curves, some classes hold many, and the
  fraction of random full-2-torsion curves whose class holds a weak curve
  is `0.53`–`0.61` at `p = 7`–`17`.  No residue of `t` mod 4 or 8 separates
  weak classes.  A sampled census cannot count the classes with none.
- §18.3 registers the characterization as exploratory and post hoc, with
  candidates `t` modulo small primes, class size, the square-free part of
  `t² − 4q³`, its 3-adic valuation.

A random model is already refuted by the data: a class of `q^{3/2}` curves
at density `3/q` would hold about `3√q` weak curves and almost never none.
The weak curves are therefore placed by structure.  Weakness is a
2-torsion property, so the structure to test first is 2-adic: the
2-volcano of the class and the 2-part of the conductor of `Z[π]`.

## Derivation (stated before measuring)

Write a weak curve as `y² = x(x − α)(x − σα)` with `α ∈ F_{q³} ∖ F_q`.
Its 2-torsion points are `T₀ = (0, 0)`, `T₁ = (α, 0)`, `T₂ = (σα, 0)`.  A
2-torsion point `(e₁, 0)` on `y² = (x − e₁)(x − e₂)(x − e₃)` is in
`2E(F_{q³})` exactly when `e₁ − e₂` and `e₁ − e₃` are squares.  Since
`−1` is a square in `F_{q³}` (`q ≡ 1 mod 4`) and `σ` preserves squares:

- `T₀` is halvable iff `α` is a square;
- `T₁` is halvable iff `α` and `α − σα` are squares;
- `T₂` is halvable iff `T₁` is.

So a weak curve has either no halvable 2-torsion point (`α` a non-square),
exactly one (`α` a square, `α − σα` not), or all three.  A generic
full-2-torsion curve has independent conditions on its three points.
Halvability governs which of the three rational 2-isogenies descend the
2-volcano, so the weak curves should be distributed over the volcano's
levels differently from the class as a whole, and classes whose volcanoes
have the wrong shape at 2 should hold none.  The 2-adic shape of a class
is read off `D = t² − 4q³ = f²·D_K`: the depth is `v₂(f)`, and the crater
is a cycle, a point or a segment according to `D_K mod 8`.

## Instrument (Rust, standalone, `census.rs`)

`F_{p⁶} = F_{p²}[s]/(s³ − c)` with `F_{p²} = F_p[i]/(i² − n)`, `n` a
quadratic non-residue, `c` a non-cube in `F_{p²}`; `σ` is `s ↦ ωs` with
`ω = c^{(p²−1)/3}`.  For each `p`:

1. Enumerate every `λ ∈ F_{p⁶} ∖ {0, 1}`; `j(λ)`; keep one `λ` per `j`.
   Every full-2-torsion curve over `F_{p⁶}` is `F_{p⁶}`-isomorphic to a
   Legendre curve or its quadratic twist, so the `j` set is complete.
2. Weak flag per `j`: `N(c) = 1` for some `c ∈ {λ, 1 − λ, λ/(λ − 1)}`.
3. The three rational 2-isogenies per `j` by Vélu on each 2-torsion
   point: codomain `j′` and whether it keeps full 2-torsion (the product
   of the other two roots is a square).  A codomain without full
   2-torsion is a floor node.
4. Trace `t` of the Legendre curve by baby-step giant-step on the order of
   two random points in the Hasse interval, once per connected component
   of the 2-isogeny graph (the trace is an isogeny invariant), with a
   second independent BSGS on another member of every component as a
   check.  The twist has trace `−t`, the same `j`, the same weak flag and
   the same graph, so classes are recorded by `|t|`.
5. Height of each node: the shortest path to a floor node in the
   2-isogeny graph.
6. Per class `|t|`: the number of full-2-torsion `j`, the number of weak
   `j`, `D`, `D_K`, `f`, `v₂(f)`, `D_K mod 8`, `|t| mod 16`, and the
   histogram of all and weak nodes by height.

Sizes: `p ∈ {5, 7, 11, 13, 17}`, exhaustive; `p = 19` if the run fits in
the session.  Seed for the random points: `20261007`.  Outputs:
`results/census_p{p}.jsonl` (one line per class), `results/summary.json`.

## Predictions (pass/fail)

- **W1 (density).**  Weak `j` over full-2-torsion `j` lies in
  `[2.5/q, 3.5/q]` at every `p ≥ 7`.
- **W2 (reach).**  The fraction of full-2-torsion `j` whose class holds a
  weak `j` lies in `[0.45, 0.70]` at every `p ≥ 7`.  Reproduces the
  sampled census.
- **W3 (2-adic law, the registered hypothesis).**  At every `p ≥ 11`, the
  indicator "class holds a weak `j`" is a function of the triple
  `(v₂(f), D_K mod 8, |t| mod 16)`: no two classes with the same triple
  disagree.  *Falsified by one disagreeing pair.*
- **W4 (placement).**  Within weak classes, weak nodes are not uniform
  over height: at every `p ≥ 11` the weak fraction at height 1 is below
  `0.5×` the class-pooled weak fraction, or the weak fraction at the
  maximum height is above `1.5×` it.  *Falsified if neither holds.*
- **W5 (consistency).**  Every component's second BSGS agrees with the
  first; the maximum height in a class equals `v₂(f)` whenever the class
  has a floor node.

## Decision rule (registered)

- W3 passes: the walk is replaced by a decision.  Read `t` from the
  public order, compute the triple, and either refuse the class or
  navigate to the height W4 singles out by halvability tests, then walk
  horizontally.  Expected steps fall from `q/3` to the inverse of the
  weak fraction at that height, measured here.  Class: **advance for the
  walk stage only**, still an accounting change for the route, since the
  walk is below 1% of its cost at `p ≥ 251` (§18.4).
- W3 fails and W4 passes: navigation by height still shortens the walk by
  the measured ratio; the class law stays open.  Class: **engineering**.
- Both fail: the 2-adic structure does not place the weak curves, and
  §18.3's candidates are next.  Class: **boundary**.

## Stop condition and inadmissible moves

Bounded by the sizes listed.  A size that does not finish in 30 minutes
of wall time is stopped and reported as partial.

Inadmissible: fitting the triple after seeing the data (the triple is
fixed above); dropping classes with few curves; counting the twist as a
second class; reading any of this as a statement about a prime-field
curve or about the route's `S / rho`.

[JV12] A. Joux, V. Vitse, *Cover and decomposition index calculus on
elliptic curves made practical*, EUROCRYPT 2012.

## Addendum A — registered 2026-10-07 after reading p = 5, 7, 11 and before reading p = 13, 17, 19

The runs at `p ≤ 11` falsify W3 as written (the triple predicts but is not
a function) and show one exact pattern and one graded one.  Both are
registered here, before the held-out sizes are read, as predictions on
`p ∈ {13, 17, 19}`:

- **A1 (depth-1 exclusion).**  Every class with `v₂(f) = 1` holds zero
  weak curves.  *Falsified by one weak `j` in any such class.*
- **A2 (depth ≥ 2 admission).**  Among classes with `v₂(f) ≥ 2`, the
  fraction holding at least one weak `j` is at least `0.75` at each size.
- **A3 (density by height, pooled over all classes).**  The weak
  fraction times `q` lies in `[1.5, 2.5]` at height 1, in `[3, 5.5]` at
  height 2, and in `[10, 17]` at height 3.
- **A4 (open, no prediction).**  Whether heights 4 and above carry weak
  curves: `p = 7` has them (25 to 60% weak), `p = 11` has none in 2,353
  nodes.  Reported, not predicted.

Only the `census_p{13,17,19}.jsonl` files decide A1 to A3; the
`tabulate` mode of `census.rs` prints the tables they are read from.

## Addendum B — registered 2026-10-07 after reading p = 13 and before reading p = 17, 19

`p = 13` passes A1, A2 and A3 and, like `p = 11`, has no weak curve above
height 3, while `p = 7` has weak curves at heights 4 to 6.  Of the sizes
run, `p = 7` is the only one with `p ≡ ±1 (mod 8)`, which is exactly the
condition for `F_{p²}` to contain the 8th roots of unity (`i` is a square
in `F_{p²}` iff `p² ≡ 1 mod 16`), and the halvability conditions behind
the heights are square tests.  Registered on the two remaining sizes:

- **B1.**  At `p = 17` (`17 ≡ 1 mod 8`, `v₂(p² − 1) = 5`) weak curves occur
  at height 4 or above.
- **B2.**  At `p = 19` (`19 ≡ 3 mod 8`, `v₂(p² − 1) = 3`) no weak curve
  occurs above height 3.
- A1 to A3 as in Addendum A, at both sizes.
