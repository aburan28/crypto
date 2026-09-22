# The descent crossover: overdetermination and yield are one parameter

**Class:** derivation. No run, no variant, no scoreboard row (§ "What this
owes the scoreboard").
**Background:**
[`research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md`](../ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md)
(the yield law and the free-oracle floor this note reads against),
[`research/notes/ecc2k130/RESEARCH_ECC2K130_RR_SOLVER_PANEL.md`](../ecc2k130/RESEARCH_ECC2K130_RR_SOLVER_PANEL.md)
§9 (the linear `NO`-certificate and its death at `d = 7`),
[`RESEARCH_DREG_MEASUREMENT.md`](RESEARCH_DREG_MEASUREMENT.md)
(solving degree against first fall degree; the Kosters–Yeo decoupling),
[`research/notes/ecc2k130/RESEARCH_ECC2K130_ROUTES.md`](../ecc2k130/RESEARCH_ECC2K130_ROUTES.md)
Route 4 (reach `m = 4`, which this note scopes).

**The question.**  The Weil descent of a Semaev system is cheap to refute
when it is massively overdetermined, and every cheap-refutation result this
repository has — the degree-6 refutation with multiplier set `{1}`, the
linear `NO`-certificate, the flat `m = 2` first-fall degrees — was measured
in that regime.  The RR panel already says why the certificate cannot
survive: *"the collapse is an artefact of `3d ≪ 131`, and it necessarily
disappears exactly when `d` grows enough for the problem to become
interesting."*  That is stated for `m = 3`, as a remark.  **Is it a general
law, and at what rate?**

**Bottom line.  It is an exact identity, it holds for every arity, and the
rate is one bit of yield per surplus equation.**  The overdetermination of
the descended system and the expected number of decompositions are not two
quantities that happen to move together — they are the same parameter
written two ways, so no choice of `m` or `l` can buy a cheap regime and a
non-empty one at once.  What that does *not* say is in §6, and one thing it
corrects is in §5.

## 1. The boundary, stated before anything is derived

The identity below is a bookkeeping statement about two counts.  It bounds
no algorithm, so it has no floor and no reference of its own; the boundaries
it is read against are the ones already in force
(`RESEARCH_ECC2K130_DECOMPOSITION.md` §0):

```text
#E(F_2^131) = 4·r,  r prime, 2^129.0000
rho reference:      2^60.8090      S = 0.077430
free-oracle floor:  2^56.40 at m = 4
```

**Falsification target.**  This note is wrong if a solver exists whose
refutation cost on these systems, at fixed `n`, is flat in the surplus
`n − m·l` — or, weaker and more interesting, if the deciding degree does not
grow as the surplus falls to zero.  Either would break the scoping rule in
§4 and make the cheap-regime measurements quotable about the attack's
regime after all.  §7 says how to test it with runs this repository can
afford.

Inadmissible: reading §4 as a lower bound on any algorithm, reading the
table in §3 as a cost, or treating the `λ` column as measured at `n = 131`
(it is derived; the yield law behind it is measured at `n ≤ 19` and its
conclusion spot-checked at full size by the §4 witnesses of the background
note).

## 2. The identity

Fix a curve `E/F_2^n`, a factor base `F = {P : x(P) ∈ V}` with `V` an
`F_2`-subspace of dimension `l`, and an arity `m`.  Two counts:

**The descended system.**  Semaev's `S_{m+1}(x_1, …, x_m, x(R)) = 0` is one
equation over `F_2^n`.  Writing each summand abscissa in a basis of `V`,
`x_i = Σ_j v_ij b_j` with `v_ij ∈ F_2`, and expanding the value in a basis
of `F_2^n` over `F_2`, gives

```text
    n Boolean equations   in   m·l Boolean unknowns.
```

That is the system `src/cryptanalysis/ic_descent_degrees.rs` builds — its
table prints `n` and `vars = m·n'` in adjacent columns, so the surplus below
is already in its output, unnamed.

**The yield.**  The `m`-subsets of `F` number `C(|F|, m)` and their sums
spread over the group, so a fixed target is hit

```text
    λ = C(|F|, m) / #E  ≈  2^{ml} / (m! · 2^n)
```

times on average (`RESEARCH_ECC2K130_DECOMPOSITION.md` §2, measured on twelve
toy cells to within a factor 1.09).

Take logs of the yield — `log₂ λ = m·l − n − log₂ m!` — and rearrange, so
that the left-hand side is the first count read against the second:

```text
    (unknowns) − (equations)  =  m·l − n  =  log₂ λ  +  log₂ m!
```

Exact, up to the two sub-bit corrections folded into `≈`: `|F| ≈ 2^l` (the
background note measures subspace point densities `0.494` and `0.452` on the
challenge curve) and `#E = 2^n + O(2^{n/2})` by Hasse.

Read as a statement about surplus equations, `S := n − m·l`:

```text
    λ  =  2^{−S} / m!
```

> **One surplus equation is one halving of the expected number of
> decompositions.**  Overdetermination is not a property of the encoding
> that can be tuned against the yield; it *is* the yield, on a log scale,
> shifted by `log₂ m!`.

### 2.1 The surplus is invariant under chaining

The count above is the **direct** descent of one `S_{m+1}`.  The harness does
not use it: for `m ≥ 3` it chains `S₃` over `m − 2` intermediate points, which
is a different system — `(m−1)·n` equations in `m·l + (m−2)·n` unknowns
(`koblitz_bench::n_vars_for`).  Chaining adds `(m−2)·n` unknowns, and it adds
`(m−1)n − n = (m−2)n` equations.  **Equal counts, so the surplus is the same
number:**

```text
    [(m−1)n]  −  [m·l + (m−2)n]  =  n − m·l  =  S
```

Checked against every cell in the two notes that print both counts —
`n = 5, ℓ = 4` (10 eqs, 17 vars, `S = −7 = 5 − 12`), `n = 7, ℓ = 3`
(14, 16, `S = −2`), `n = 9, ℓ = 6` (18, 27, `S = −9`).

This has a consequence for the scaling target, whose primary metric is
`m·ℓ + (m−2)·n` and whose instruction is to drive it down
(`research/notes/ecc2k130/RESEARCH_KOBLITZ_SCALING_TARGET.md`).  That metric is
the **matrix size**, and lowering it is real and worth doing — it is what
makes an instance affordable.  It is not a move toward the attack's regime,
because it leaves `S` untouched.  Size and faithfulness are separate axes and
this repository has been optimising only the first.

**Cross-check against an independently derived number.**  §3.2 of the
background note computes the yield at the free-oracle optimum for `m = 4`
(`l* = 26.70`) as `λ = 2^{−28.78}`, by a route that never mentions the
descended system.  The identity gives
`S = 131 − 4(26.70) = 24.20`, hence `λ = 2^{−24.20}/24 = 2^{−28.79}`.  The
two agree to the last digit printed.  That is the only validation this note
has, and it is worth more than the derivation being short.

## 3. It holds for every arity — and the arity knob works against itself

Three dimensions per `m`, all at `n = 131`, all derived:

| `m` | saturating `l` (`λ = 1`) | square `l` (`m·l = n`) | gap `log₂(m!)/m` |
|---:|---:|---:|---:|
| 2 | 66.00 | 65.50 | 0.50 |
| 3 | 44.53 | 43.67 | 0.86 |
| 4 | 33.90 | 32.75 | 1.15 |
| 5 | 27.58 | 26.20 | 1.38 |
| 6 | 23.42 | 21.83 | 1.58 |
| 8 | 18.29 | 16.38 | 1.91 |

The saturating column is the background note's §2 table, reproduced here
from the identity rather than copied — which is the check that the identity
is the same statement.  The point of the third column: **the dimension at
which decompositions appear and the dimension at which the system stops
being overdetermined differ by one to two, for every arity worth running.**
They are not two regimes with room between them.

Now locate the attack's own operating points.  `l*` is the dimension that
minimises the free-oracle floor (`RESEARCH_ECC2K130_DECOMPOSITION.md` §5.3),
so these are the cells the method would actually run in:

| `m` | `l*` | `m·l*` | surplus `S` | `λ = 2^{−S}/m!` | floor | vs rho |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 43.33 | 86.66 | `+44.34` | `2^−45.34` | `2^89.25` | `2^+28.44` |
| 3 | 33.00 | 99.00 | `+32.00` | `2^−34.58` | `2^68.58` | `2^+7.77` |
| **4** | 26.83 | 107.32 | `+23.68` | `2^−28.26` | `2^56.40` | **`2^−4.41`** |
| 5 | 22.76 | 113.80 | `+17.20` | `2^−24.11` | `2^48.44` | `2^−12.37` |
| 6 | 19.89 | 119.34 | `+11.66` | `2^−21.15` | `2^42.85` | `2^−17.96` |
| 8 | 16.12 | 128.96 | `+2.04` | `2^−17.34` | `2^35.61` | `2^−25.20` |

Read the surplus column against the floor column.  They move in opposite
directions, and that is the note's one structural finding:

> **The arity that makes the floor beat rho is the arity that destroys the
> overdetermination the algebra was living on.**  Raising `m` from 2 to 8
> improves the floor by `2^53.6` and drives the descended system from 44
> surplus equations to 2 — from comfortably overdetermined to square, which
> is the regime in which a Boolean system is hardest for any generic
> algebraic method.

This is the precise sense in which "just reach `m = 4`" is not a way past
the barrier.  It is a real improvement to the floor — `m = 4` is the first
arity whose floor sits below rho — bought by moving the solver into a worse
regime by the same parameter.  Route 4 is worth running for its own
falsifier (does the first fall degree grow with `n` on `m = 4` systems?);
it is not a way around the identity, because the identity has no outside.

## 4. What this licenses: a scoping rule for extrapolation

Every cheap-refutation measurement on a descended Semaev system is a
measurement at a specific surplus, hence at a specific `λ`.  Quoting it as
evidence about the attack's regime is an extrapolation across the distance
between the two, and the identity prices that distance:

| measurement | `m` | `m·l` | surplus | `λ` there | distance to the `m = 4` operating point |
|---|---:|---:|---:|---:|---:|
| linear `NO`-certificate alive (`solver_15`, `d = 6`) | 3 | 18 | `+113` | `2^−115.6` | `2^87.3` in `λ` |
| certificate dead (`solver_15`, `d = 7`) | 3 | 21 | `+110` | `2^−112.6` | `2^84.3` |
| `D_refute = 6`, `FFD = 3` (`RESEARCH_DREG_MEASUREMENT.md`, `n = 5`) | 3 | 12 | `−7` | `2^4.4` | — (different `n`) |
| attack cost-optimum | 4 | 107.32 | `+23.68` | `2^−28.26` | — |
| saturating base | 4 | 135.58 | `−4.59` | `2^0` | `2^28.3` |

The rule that falls out, and the reason to write this down:

> **State the surplus with every refutation measurement.**  A refutation
> degree, a first fall degree, or a certificate observed at surplus `S` is
> an observation about a system in which a typical target has `2^{−S}/m!`
> decompositions.  It may be quoted about the attack only alongside the
> surplus the attack runs at, and the gap between them named.

Applied to what is already on the books: the `d ≤ 6` linear certificate is
an observation at `λ ≈ 2^{−116}`, about eighty-seven bits of yield away from
the cell the attack would use.  That is a sharper statement of the RR
panel's `2^38`-in-factor-base-size gap, in the unit that explains *why* the
certificate cannot be carried across.

The third row is the one that matters for what to run next, and it points
the other way: the `n = 5, m = 3` cell of the DREG note sits at surplus
`−7`, i.e. **already underdetermined**, which is why it refutes at degree 6
rather than cheaply, and why it is the right shape of instance to build a
ladder from.  Small `n` with `m·l ≳ n` is a faithful scale model of the
attack's regime in a way that large `n` with `m·l ≪ n` is not.

## 5. A correction to how this was first put

An earlier statement of this argument — in conversation, not in a committed
note — said the descended system "crosses from massively overdetermined to
underdetermined at precisely the dimension where targets start to
decompose," and offered `4 × 33.90 = 135.6 > 131` as the `m = 4` instance.

That is true at the **saturating** dimension and false at the dimension the
attack would actually choose.  The free-oracle optimum for `m = 4` is
`l* = 26.83`, where the system carries `+23.68` surplus equations and is
still firmly overdetermined; the attack does *not* run in the
underdetermined regime.  What the identity gives is a **gradient, not a sign
change**: the surplus falls from `+113` where the cheap algebra was measured
to `+24` where the attack runs, and only reaches zero at arities (`m = 8`)
whose oracle cost is worse still.

The conclusion is unchanged and now carries a number instead of an
adjective — the cheap regime and the attack's regime are separated by about
ninety bits of yield — but the mechanism was misstated, and per §3 of
`AGENTS.md` that is **accounting**: nothing measured changed, and no gain is
claimed from restating it.

## 6. What this does not settle

- **It is not a hardness result and bounds no algorithm.**  It is a relation
  between two counts.  A solver that refutes a square system cheaply is
  excluded by nothing here; §1's falsification target is exactly that
  solver, and the RR panel's §8 caveat still applies — `F` is algebraically
  structured, and genuinely sub-quadratic behaviour on structured instances
  is excluded by no measurement in this repository.
- **"Overdetermined" is not the same knob as "cheap".**  The two thresholds
  in play are *different mechanisms* and this note does not merge them.  The
  linear certificate dies at `d = 7` because the `S₄` value set's affine
  span saturates `F_2^131` (measured: span dimension `71, 97, 123, 131` at
  `d = 4…7`), which happens at `m·l = 21`, far above the square crossover at
  `m·l = 131`.  The identity explains why overdetermination must fall as the
  base grows; it does **not** explain the span law, and nothing here derives
  one from the other.
- **The `λ` column is derived at `n = 131`.**  The yield law behind it is
  measured at `n ≤ 19` (twelve cells, ratio in `[0.57, 1.09]`), and its
  conclusion is spot-checked at full size only by the certified `m = 2, 3, 4`
  witnesses of the background note's §4 — which confirm that decompositions
  exist, not the rate.
- **Only subspace factor bases.**  The equation count `n` and the unknown
  count `m·l` both assume `F = {P : x(P) ∈ V}` for an `F_2`-subspace `V`.
  An orbit-union base has no `V`-basis to descend in, so this identity says
  nothing about the `2^125.55` orbit-union rows; those are priced by
  materialisation, not by algebra.
- **`m!` is the ordered/unordered correction, and nothing more.**  Where the
  base or the arity makes repeated summands or inverse-cancelling tuples
  non-negligible, the `log₂ m!` shift is the wrong constant — that is the
  same over-counting the background note's §2 had to correct once already.

## What this owes the scoreboard

Nothing, and stating why is part of §7 of `AGENTS.md` rather than an
exemption from it.  This note introduces no variant, measures no cost, and
claims no ratio to either boundary; there is no row to add, and adding one
would put a derivation on a page whose standing constraint is that every
number on it comes from a frozen experiment file.  What it changes is which
existing rows may be quoted about which regime, which is prose, not an axis.

If §7 of this note is run, its results do carry a scoreboard obligation, and
it is recorded there.

## 7. What to run

Ordered by cost, all inside what this repository already has:

1. **Report the surplus wherever a descent system is tabulated.**  `n − m·ℓ`
   is three existing fields, so this is a column, not an experiment.  Two
   tables want it: `dreg_sweep`, beside `FFD` and `D_refute`; and
   `ic_descent_degrees`, which already prints `n`, `vars` and `eqs` next to
   each other and measures `D_av/D_sr` — the ratio of the degree reached to
   what a structureless system of the same shape would reach.  **That ratio
   against the surplus is the falsification target of §1, at the cost of
   plotting two columns this repository already computes.**  Semi-regularity
   is itself a statement about a system with no exploitable structure, and
   an overdetermined system is not generic; if `D_av/D_sr` is flat in the
   surplus, §4's scoping rule is wrong and should be withdrawn.

   It also lands on the comparison the DREG note asks for, and not gently.
   That note proposes `n = 9` against `n = 15` because they share an unknown
   count (27) and differ in field degree, which is meant to separate "the gap
   grows with `n`" from "the gap grows with the matrix".  By §2.1 their
   surpluses are `S = 9 − 18 = −9` and `S = 15 − 12 = +3` — **opposite
   signs**: one is underdetermined and has `2^6.4` expected decompositions per
   target, the other is overdetermined and has `2^{−5.6}`.  The pair is
   therefore confounded a third way, by twelve bits of yield, and reading a
   `D_refute` difference between them as a field-degree effect would be
   wrong.  Matching the surplus instead of the unknown count is the fix, and
   it costs nothing but a different choice of `ℓ`.
2. **Build the ladder at fixed surplus, not at fixed `ℓ`.**  Choose
   `(n, m, ℓ)` along `m·ℓ ≈ n + log₂ m!` so that every rung is a scale model
   of the attack's regime, and measure `D_refute` along it.  This is the
   direct test of §1's falsification target and the honest version of the
   extrapolation Route 4 needs.
3. **Then the `m = 4` cells**, which is Route 4's own falsifier, read with
   the surplus named.

## References

- **I. Semaev**, *Summation polynomials and the discrete logarithm problem on
  elliptic curves*, ePrint 2004/031.
- **C. Diem**, *On the discrete logarithm problem in elliptic curves*,
  Compositio Math. 147 (2011) — factor bases as subspaces, which is the
  shape the unknown count `m·l` assumes.
- **M. Kosters, S. L. Yeo**, *Notes on summation polynomials*,
  arXiv:1503.08001 — the first-fall-degree collapse this note generalises
  from a remark about `3d ≪ 131` to an identity in the surplus.
- **C. Petit, J.-J. Quisquater**, *On polynomial systems arising from a Weil
  descent*, ASIACRYPT 2012 — the complexity claim stated in the first fall
  degree, and therefore the claim the scoping rule in §4 applies to.
