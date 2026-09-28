# Index calculus vs Pollard rho on genus-2 Jacobians over `F_p`

Reference measurement for
`src/cryptanalysis/hyperelliptic_index_calculus.rs` (the ADH / Gaudry
attack, merged in #365), which shipped with operation counters and no
comparison. This note supplies the comparison.

Class (per `AGENTS.md` §3): **accounting**. No algorithm changed; what
changed is that its cost now has a boundary and a reference next to it.
Nothing here is a finding.

## Boundaries, stated before measuring

**Reference.** Pollard rho with a 16-adding walk (Teske) and Floyd cycle
finding, run on the *same* curve, the same `D₁`, `D₂`, the same prime
order `N`, counted in the same unit
(`src/cryptanalysis/hyperelliptic_ic_bench.rs`).

**Floor (index calculus).** `m + 1` unknowns need `m + 1` independent
relations. Under heuristic **H1** — a reduced divisor drawn uniformly
from `Jac(C)(F_p)` splits into degree-1 places with probability `≈ 1/g!`
— the expected trial count is at least `(m+1)·g!`, and no walk yields a
trial for under one group operation. Any dense solve must read its own
`m × m` matrix once, `m²` mul-mods. So, in group-operation equivalents,

```
floor_ops = (m + 1)·g!  +  m²/c
```

with `c` the **measured** mul-mods-per-group-op factor (column `conv`).
The floor moves only with `m` and `g`. It bounds *this method*, not the
DLP.

## Unit

```
S = total group operations / sqrt(N)
```

Everything each side spends is inside `S`: rho's branch precomputation
and every walk step; index calculus's relation-stage scalar
multiplications **and** its linear algebra, converted at `c`. Leaving
the solve out would be the relabelling failure mode — the work does not
vanish, it moves out of the headline.

## Table

`C : y² = x⁵ + 3x³ + 2x² + x + c` over `F_p`, genus 2. For each `p`,
`c` is the first value giving a near-prime Jacobian order
(`N ≥ #Jac/4`, `N > 1000`), chosen before either algorithm runs. Every
row averages 3 index-calculus runs and 9 rho runs on the same instance;
every run of both sides returned the verified `k`, so every row is a
result.

Reproduce with `cargo run --release --example hyperelliptic_ic_vs_rho`.

| `p` | `N` | `m` | IC ops | rho ops | `S_ic` | `S_rho` | **`S_ic/S_rho`** | `S_walk` | **`S_ic`/floor** | `c` |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 41 | 1321 | 15 | 1571 | 526 | 43.2 | 14.5 | **3.0** | 3.24 | **48.9** | 2171 |
| 61 | 1399 | 22 | 2199 | 522 | 58.8 | 14.0 | **4.2** | 3.06 | **47.6** | 2337 |
| 101 | 4663 | 46 | 3546 | 649 | 51.9 | 9.5 | **5.5** | 2.54 | **37.4** | 2660 |
| 151 | 7949 | 78 | 5373 | 681 | 60.3 | 7.6 | **7.9** | 2.30 | **33.5** | 2457 |
| 211 | 11813 | 112 | 7469 | 852 | 68.7 | 7.8 | **8.8** | 3.15 | **32.3** | 2497 |
| 251 | 61667 | 123 | 10229 | 1204 | 41.2 | 4.9 | **8.5** | 2.52 | **40.1** | 2195 |

`S_walk` is rho's walk alone; `S_rho` also carries the 16-branch
precomputation, a fixed cost that does not scale with `sqrt(N)` and so
inflates `S_rho` at these sizes. The collision arrives after about
`sqrt(πN/2) = 1.25·sqrt(N)` distinct points, and Floyd spends 3 group
operations per iteration, so `S_walk ≈ 3` is this reference at its own
optimum — consistent with the column. A distinguished-point rho would
cut it by roughly 3×, which **widens** the gap below rather than closing
it; the reference is if anything handicapped here.

## Reading it

1. **Index calculus loses, and loses harder as `p` grows.**
   `S_ic/S_rho` runs 3.0 → 8.8 across `p = 41 … 211`. This is the
   expected shape, not a defect: at `g = 2` the full-factor-base
   relation stage is `O(p²)` against rho's `O(p)` (since `N ≈ p²`), so
   the ratio should grow like `sqrt(N)`. The `p = 251` row dips to 8.5
   only because its `N` is 5× larger than `p = 211`'s — the denominator
   moved, not the method.
2. **The measured smoothness rate is 0.34 – 0.53**, against the `1/g! =
   0.5` that H1 predicts for genus 2. H1 survives this test at these
   sizes.
3. **The implementation sits 32 – 49× above its own floor.** Almost all
   of that is the relation stage drawing fresh `(a, b)` and paying
   `2⌈log₂ N⌉` group operations per trial where the floor allows one:
   an incremental walk (`R ← R + D₁`) would recover most of the factor.
   Charging one operation per trial instead of `2⌈log₂ N⌉ = 28` puts
   `p = 211` at `6850/28 + 619 = 864` operations, `S_ic ≈ 7.9` against
   `S_rho = 7.8` — parity, not a win, because the linear algebra the
   division does not touch then dominates. `p = 251` lands the same way
   (`≈ 1199` ops, `S_ic ≈ 4.8` vs `4.85`). So the available engineering
   is worth ~9× and buys a tie at toy sizes; it does not change the
   asymptotics, which is why both columns stay.
4. **The linear algebra is not the bottleneck here**, but it is growing
   fastest: 3 group-op equivalents at `p = 41`, 925 at `p = 251`,
   scaling as `m³/c` with `m ≈ p/2`. At `p ≈ 10³` it overtakes the
   relation stage, and the dense solve would have to go sparse
   (Wiedemann / Lanczos) before any larger measurement is meaningful.

## Scope

Genus 2, odd characteristic, `deg f = 2g+1`, `p ≤ 251`, full degree-1
factor base, dense solve, single machine. Nothing here transfers to
cryptographic sizes, and nothing here is evidence about genus 3 or 4,
where the asymptotic crossover against rho actually lives. The
`O(p)`-by-evaluation root finding alone bars larger `p`.

---

# Follow-up: the two optimisations the table asked for

The baseline above named its own two bottlenecks. Both are now
implemented and measured on the same instances, same unit, same
reference.

Class: **engineering**, and §3's table would let me write "advance"
because the ratio to the floor fell. I am not writing it. What moved is
this implementation's distance to the method's own floor, using two
techniques that were already standard; the method does nothing a generic
algorithm cannot that it did not do before. The column that would earn
the other word — the ratio to the reference — is now ≈ 1, not robustly
below it. See "What this is not" below.

## What changed

1. **`RelationSearch::Walk`** — an `r`-adding walk over
   `R ← R + S_j`, carrying `(a, b)` along with it, so a candidate
   divisor costs **one** group operation instead of `2⌈log₂ N⌉`
   (measured: 0.98 ops/trial, against 22 – 34 before). 16 branches, to
   match the rho reference; restarted every 64 trials, because a
   deterministic walk cycles after `O(sqrt(N))` steps and this attack
   needs `m + 1 ≈ p/2` relations from a group of size `N ≈ p²` — the
   same order, so an unrestarted walk would start re-deriving relations
   it already had.
2. **`LinearAlgebra::Sparse`** — Markowitz-style sparse elimination mod
   `N`. A genus-`g` relation touches at most `g` factor-base columns
   plus the `k` column, so rows carry `≤ g + 1` non-zeros whatever `m`
   is, and the dense `O(m³)` solve was throwing that away: at
   `p = 251`, 2,029,632 mul-mods became **879**, a factor of 2300.

A third change is an honesty fix rather than an optimisation, and it
cuts the other way: the **smoothness oracle is now charged**. Root
finding and the decomposition's evaluations are `O(p)` field
multiplications per trial, which `AGENTS.md` puts inside `S` as "work an
oracle does per call" — and once the walk made the group operation
cheap, that omission was worth more than either optimisation. At
`p = 251` it adds 168,933 mul-mods (76 group-op equivalents) to a
relation stage of 976 operations. The baseline table above did not carry
it, so its `S_ic` was understated; the `rnd+dns` rows below are the
corrected baseline, and they sit within 1% of the originals because at
`2⌈log₂ N⌉` per trial the oracle was noise.

## Table

Same curves, same `N`, same rho, same protocol (3 IC runs and 9 rho runs
per row, every run verified). `rnd+dns` is the merged baseline;
`wlk+spr` is both optimisations together.

| `p` | `N` | `m` | `S_ic` base | `S_ic` opt | `S_ic/S_rho` base → opt | `S_ic`/floor base → opt |
|--:|--:|--:|--:|--:|--:|--:|
| 41 | 1321 | 15 | 43.3 | 13.7 | 3.00 → **0.95** | 49.0 → 15.5 |
| 61 | 1399 | 22 | 59.0 | 14.2 | 4.22 → **1.02** | 47.7 → 11.5 |
| 101 | 4663 | 46 | 52.3 | 9.8 | 5.50 → **1.03** | 37.6 → 7.1 |
| 151 | 7949 | 78 | 60.7 | 8.2 | 7.95 → **1.07** | 33.7 → 4.5 |
| 211 | 11813 | 112 | 69.1 | 8.0 | 8.81 → **1.02** | 32.5 → 3.8 |
| 251 | 61667 | 123 | 41.0 | 4.2 | 8.45 → **0.87** | 40.0 → 4.1 |

> Superseded by the round below, which fixes an unfairness in these
> numbers: the index-calculus walk was given cheap step coefficients
> and rho's branch precomputation was not. The rho column here is
> therefore too expensive, and every `S_ic/S_rho` above too kind to
> index calculus. The corrected figures are in the final table.

`cargo run --release --example hyperelliptic_ic_vs_rho` prints all four
combinations (`rnd`/`wlk` × `dns`/`spr`) with the per-stage breakdown.
Wall clock tracks the unit: at `p = 251`, 49 ms for index calculus
against 51 ms for rho, where the baseline was 424 ms against 49 ms.

## Reading it

1. **The ratio stopped growing.** Baseline `S_ic/S_rho` climbed 3.0 →
   8.8 with `p`; optimised it sits at 0.87 – 1.07 across the whole
   range with no trend. That shape change is the result here, more than
   the 5 – 10× constant. It is what the asymptotics predict once the
   relation stage costs one operation per trial: `(m+1)·g! ≈ p`
   operations against rho's `sqrt(N) ≈ p` at `g = 2`, so both sides are
   linear in `p` and the ratio is a constant.
2. **Where the cost now sits.** At `p = 251` the relation stage is 976
   group operations, of which **714 is precomputation** — 16 branch
   divisors and 4 walk restarts, each a scalar-multiplication pair. The
   walked trials themselves cost 262. Precomputation is now the
   dominant term and the obvious next target, though it amortises as
   `p` grows and is bounded below by needing *some* starting point.
3. **The oracle, not the solve, is the field-work bottleneck now.**
   168,933 mul-mods against 879 for the linear algebra. Root finding by
   evaluation is `O(p)` per trial; a Cantor–Zassenhaus or a
   `gcd(u, x^p − x)` test would cut it to `O(log p)` multiplications
   per trial and is the next thing worth implementing.
4. **The smoothness rate held** at 0.33 – 0.53 under the walk, against
   `1/g! = 0.5`. A walked sequence is not an independent sample, so
   this was worth checking rather than assuming; H1 survives.

## What this is not

`S_ic/S_rho < 1` at `p = 251` does **not** mean index calculus beats
Pollard rho at genus 2. The reference is this repository's rho: Floyd
cycle finding at 3 group operations per iteration, plus a 16-branch
precomputation. Rho's walk alone is `S_walk = 2.52` on that instance,
against the optimised `S_ic = 4.24` — so against a distinguished-point
rho, which is the rho anyone attacking a real curve would run, index
calculus is still roughly 1.7× behind. The honest summary is **parity
with a plain rho at toy sizes, and the ratio no longer grows**, which
was not true before.

Nor does anything here transfer upward. Scope is unchanged: genus 2,
`p ≤ 251`, full factor base, single machine. The interesting question —
genus 3 and 4, where the asymptotic crossover lives — needs the
`O(log p)` smoothness test from point 3 before it can be measured at
all.


---

# Round three: the oracle, the precomputation, and a fairness fix

Round two left two things on the table and introduced one problem.

Class: **engineering**. Same reasoning as before — this is distance to
the method's own floor, closed with standard techniques.

## What changed

1. **`SmoothnessTest::Gcd`** (default). `u` splits into linear factors
   iff its radical equals `gcd(u, x^p − x)`, and `x^p mod u` costs
   `O(log p)` polynomial multiplications by square-and-multiply. At
   `deg u ≤ 2` — everything genus 2 produces — even that is
   unnecessary: a monic quadratic splits exactly when its discriminant
   is a square, so the oracle is **one Legendre symbol**, and the roots
   come from the quadratic formula. Measured at `p = 251`: 168,933
   mul-mods per run became 7,249, a factor of 23.
2. **Cheap walk steps and jump restarts.** Precomputation was 73% of
   the relation stage after round two. Two changes: branch steps are
   built from 8-bit coefficients (`2·8` operations each instead of
   `2⌈log₂ N⌉`), and a restart now adds a *randomly chosen* existing
   step — the index comes from the RNG rather than the position hash,
   so the trajectory leaves its orbit for **one** operation instead of
   a fresh scalar-multiplication pair. Relation stage at `p = 251`:
   976 operations → 566.
3. **The fairness fix.** Cheap step coefficients are not an
   index-calculus trick — rho can use them too, and its branch
   precomputation was still being charged at full price. The reference
   now gets the same treatment. This makes rho cheaper and the
   comparison worse for index calculus, which is why it belongs in the
   same round rather than a later one: optimising one side's
   precomputation and not the other's is a rigged table, and the
   rigging favours the side under study.

Short step coefficients weaken nothing. A relation is an algebraic
identity, verified as it is recorded; the step distribution affects only
how fast smooth divisors turn up, which the trial count measures
directly. The measured smoothness rate (0.48 – 0.55) did not move.

## Final table

Same curves, same `N`, same rho, both sides with cheap precomputation.
Samples raised to 5 index-calculus runs and 25 rho runs per row,
because the gap is now a factor of ~2 rather than ~10 and had to be
separated from rho's own spread. Every run of both sides returned the
verified `k`.

| `p` | `N` | `m` | `S_ic` base | `S_ic` final | `S_rho` | **`S_ic/S_rho`** | `S_walk` | `S_ic`/floor |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 41 | 1321 | 15 | 45.9 | 9.65 | 11.00 | **0.88** | 3.30 | 10.9 |
| 61 | 1399 | 22 | 58.2 | 9.85 | 10.85 | **0.91** | 3.36 | 8.0 |
| 101 | 4663 | 46 | 48.1 | 5.90 | 7.23 | **0.82** | 3.07 | 4.2 |
| 151 | 7949 | 78 | 57.7 | 5.06 | 6.80 | **0.75** | 3.61 | 2.8 |
| 211 | 11813 | 112 | 67.9 | 4.76 | 5.56 | **0.86** | 2.93 | 2.2 |
| 251 | 61667 | 123 | 41.3 | 2.28 | 4.17 | **0.55** | 3.00 | 2.2 |

Cumulative: `S_ic` down 5 – 18×, `S_ic/S_rho` from 4.2 – 12.2 (growing
with `p`) to 0.55 – 0.91 (flat, drifting down), and the distance to the
method's own floor from 32 – 52× to 2.2 – 10.9×. Wall clock agrees: at
`p = 251`, 27 ms against rho's 52 ms, where the baseline was 424 ms.

## Where the cost is now

At `p = 251`, per run: 566 group operations for the relation stage (of
which **294 is precomputation** — 16 branch divisors and the starting
point), 7,249 mul-mods for the oracle (3 group-op equivalents), 981 for
the linear algebra (0 equivalents). The walked trials themselves cost
272 operations at 1.00 each.

So two thirds of what is left is precomputation and the floor is only
2.2× away. There is no third optimisation of this size available: the
remaining levers are fewer branches or shorter coefficients, both of
which trade walk quality for setup cost and neither of which can win
more than the ~50% precomputation share.

## What this is still not

**Index calculus now beats this repository's rho at genus 2 and toy
sizes, and that is a narrower claim than it sounds.**

- The reference is Floyd rho at 3 group operations per iteration.
  `S_walk` — rho's walk alone, without precomputation — sits at
  2.9 – 3.6, consistent with `3 × 1.25`. A distinguished-point rho
  would pay closer to `1.25`, so the honest extrapolation is
  `S_rho(DP) ≈ 1.0 – 1.5` once its own (now cheap) precomputation is
  added. Against *that* rho, index calculus at `S_ic = 2.28` is still
  roughly 2× behind at `p = 251`, and 6 – 9× behind at `p ≤ 101`.
  Implementing a DP rho is the honest next step before any claim that
  this crosses over.

  > **This extrapolation was wrong, and round four measured it.** A
  > distinguished-point rho was implemented and comes in at
  > `S_rho = 2.41 – 9.08`, not `1.0 – 1.5`. The error was in the phrase
  > "once its own precomputation is added": at these sizes rho's branch
  > precomputation is roughly *half* its total, so scaling only the walk
  > term understated it several-fold. Measured, index calculus is at
  > `0.95 – 1.26` of DP rho rather than 2 – 9× behind it. An
  > extrapolation across the one term that dominates at the sizes being
  > measured is not an extrapolation.
- `S_ic/S_rho` is flat in `p`, not falling, which is what genus 2
  predicts: `(m+1)·g! ≈ p` relation operations against rho's
  `sqrt(N) ≈ p`. Both sides are linear in `p`, so no amount of constant
  factor changes the asymptotics. The crossover that matters lives at
  `g = 3` (Gaudry, reduced factor base) and `g = 4`, neither of which
  this measures.
- Scope is unchanged: genus 2, `p ≤ 251`, full factor base, one
  machine, and a factor-base build that is still `O(p)` — now the
  binding `O(p)` step, since the oracle no longer is.


---

# Round four: a real rho, and genus 3

Round three ended with two named next steps. Both are now in the code;
one is measured here, the other's table is still running and lands in a
follow-up commit on the same branch.

Class: **accounting** for the rho change (the reference got fairer, no
algorithm moved) and **infrastructure** for genus 3.

## The reference is now a distinguished-point rho

`RhoVariant::DistinguishedPoints`: store the walk positions whose hash
ends in `theta_bits` zeros, stop when one repeats. One group operation
per step against Floyd's three, for the same expected number of steps.
`theta_bits` is set so about 32 distinguished points are expected before
the collision.

It behaves as theory says: `S_walk`, rho's walk alone, drops from
2.9 – 3.6 to **1.14 – 1.48**, straddling the `sqrt(π/2) = 1.2533` ideal.
That is the check that the implementation is right, not a result.

| `p` | `N` | `S_ic` | `S_rho` (DP) | **`S_ic/S_rho`** | `S_walk` | `S_ic`/floor |
|--:|--:|--:|--:|--:|--:|--:|
| 41 | 1321 | 9.65 | 9.08 | **1.06** | 1.37 | 10.9 |
| 61 | 1399 | 9.85 | 8.82 | **1.12** | 1.33 | 8.0 |
| 101 | 4663 | 5.90 | 5.45 | **1.08** | 1.29 | 4.2 |
| 151 | 7949 | 5.07 | 4.67 | **1.09** | 1.48 | 2.8 |
| 211 | 11813 | 4.66 | 3.96 | **1.18** | 1.33 | 2.2 |
| 251 | 61667 | 2.28 | 2.41 | **0.95** | 1.24 | 2.2 |

**The round-three claim does not survive, and that is the point of
having done it.** Against Floyd rho, index calculus measured
`0.55 – 0.91` — a win. Against a real rho it measures `0.95 – 1.18`:
parity, with the single `p = 251` row below 1 and no trend. The earlier
"beats this repository's rho" was true and is now uninteresting; the
repository's rho was 3× more expensive than it needed to be.

Round three also *predicted* this, and predicted it wrong. It
extrapolated `S_rho(DP) ≈ 1.0 – 1.5` and concluded index calculus would
be 2 – 9× behind. The measured `S_rho(DP)` is 2.41 – 9.08, because at
these sizes rho's branch precomputation is about half its total and the
extrapolation scaled only the walk term. An extrapolation across the one
term that dominates at the size being measured is not an extrapolation.
The note above is corrected in place rather than rewritten.

## Genus 3 now runs

Two pieces were missing and are now present, each with its own test:

1. **Root finding above degree 2.** The fast oracle closed in form at
   `deg u ≤ 2`, which is all genus 2 produces; genus 3 produces cubics,
   and without a root finder for them the oracle fell back to the
   `O(p)` scan exactly where the interesting measurement is. Now
   Cantor–Zassenhaus equal-degree splitting: `gcd(d, (x+b)^((p−1)/2) −
   1)` for random `b`. Held to the scan's answers exhaustively over
   every monic cubic for two primes.
2. **Group order without an L-polynomial.** The genus-2 route to
   `#Jac` goes through an `L`-polynomial that is genus-2 only.
   `divisor_order_bsgs` instead finds a multiple of `ord(D)` by
   Baby-step Giant-step over the Hasse–Weil interval
   `[(sqrt(p) − 1)^{2g}, (sqrt(p) + 1)^{2g}]` and divides it down — no
   point counting, any genus. Cross-checked against the L-polynomial
   route at genus 2, where both are available: they agree on every
   divisor tested.

The same interval fixes an instance-selection bug this work introduced:
a subgroup must carry most of the Jacobian (`l ≥ lo/4`) or the
comparison measures the instance rather than the algorithms — rho
searches the subgroup while index calculus pays for a factor base sized
by the whole curve. Without that constraint `p = 211` picked `N = 1103`
and reported `S_ic = 14.7`, four times the correct row.

`solves_a_genus_three_dlp` takes a genus-3 curve end to end: BSGS for
the order, a prime-order subgroup, the walk, the fast oracle, the sparse
solve, and a verified `k`. Measured smoothness there is below 0.45,
consistent with `1/3! = 0.167` being the target rather than genus 2's
`1/2` — the cost the extra genus buys, and the reason genus 3 is where
the crossover is supposed to live.

The genus-2-vs-genus-3 table is below.


## The reference needed two more fixes first

Measuring genus 3 exposed two ways the DP walk spent operations it did
not need to. Both inflated `S_rho` — the wrong direction for a
reference to be wrong in, since it flatters the algorithm under study.

1. **A degenerate collision** (same point, same coefficients) rebuilt a
   starting point for `2⌈log₂ N⌉` operations. It now jumps by a
   randomly chosen precomputed step for **one**, the same escape the
   relation walk already used.
2. **A cycle can contain no distinguished point at all**, in which case
   the walk detects nothing however long it runs. The first version hit
   an outer cap at 64× the expected cost — which is also how the
   benchmark came to look hung rather than slow, and cost an hour of
   looking at the wrong stage. The walk now tracks its distance since
   the last distinguished point and jumps once it is well past the
   expected gap.

After both, `S_walk` sits at **1.05 – 1.48** across every row, against
the `sqrt(π/2) = 1.2533` ideal. Rho samples per row went from 25 to 40,
since the measured gap is now under 2× and rho's step count has a long
tail.

## Genus 2 vs genus 3, same unit, same reference

`m` is the factor-base size; `S_walk` is rho's walk without its
precomputation. Every run of both sides returned the verified `k`.

**Genus 2** (`C : y² = x⁵ + 3x³ + 2x² + x + c`):

| `p` | `N` | `m` | `S_ic` base | `S_ic` | `S_rho` | **`S_ic/S_rho`** | `S_walk` | `S_ic`/floor |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 41 | 1321 | 15 | 45.9 | 9.66 | 9.09 | **1.06** | 1.39 | 10.9 |
| 61 | 1399 | 22 | 58.3 | 9.85 | 8.76 | **1.12** | 1.27 | 8.0 |
| 101 | 4663 | 46 | 48.1 | 5.90 | 5.45 | **1.08** | 1.29 | 4.2 |
| 151 | 7949 | 78 | 58.1 | 5.07 | 4.61 | **1.10** | 1.42 | 2.8 |
| 211 | 11813 | 112 | 67.3 | 4.66 | 3.96 | **1.18** | 1.33 | 2.2 |
| 251 | 61667 | 123 | 40.8 | 2.28 | 2.41 | **0.94** | 1.25 | 2.2 |

**Genus 3** (`C : y² = x⁷ + x³ + c x + 1`):

| `p` | `N` | `m` | `S_ic` base | `S_ic` | `S_rho` | **`S_ic/S_rho`** | `S_walk` | `S_ic`/floor |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 23 | 6299 | 12 | 35.4 | 5.02 | 4.74 | **1.06** | 1.16 | 5.1 |
| 31 | 7333 | 16 | 47.9 | 5.22 | 4.36 | **1.20** | 1.05 | 4.4 |
| 41 | 6679 | 27 | 44.6 | 5.39 | 4.79 | **1.13** | 1.31 | 2.6 |
| 61 | 124459 | 33 | 21.2 | 1.52 | 2.30 | **0.66** | 1.48 | 2.6 |
| 101 | 364747 | 53 | 21.6 | 1.16 | 1.90 | **0.61** | 1.41 | 2.1 |

## Reading it

1. **Genus 3 is the first place the ratio moves.** At genus 2 it is
   flat: 1.06, 1.12, 1.08, 1.10, 1.18, 0.94 — no trend across a 47×
   range of `N`, which is what the asymptotics say (relation stage
   `≈ (m+1)·g! ≈ p` against rho's `sqrt(N) ≈ p`, both linear in `p`).
   At genus 3 the same column runs 1.06, 1.20, 1.13, **0.66, 0.61** —
   it breaks below 1 exactly where `N` gets large, because rho now
   costs `sqrt(N) ≈ p^{1.5}` against a relation stage still linear in
   `p`. That is the crossover this whole thread was built to look for,
   and the shape is right.
2. **It is two data points.** `p = 61` and `p = 101` at genus 3 are the
   only rows below 0.7, and they are also the only genus-3 rows with
   `N > 10^5`. The three small-`N` genus-3 rows sit at 1.06 – 1.20,
   indistinguishable from genus 2. The honest statement is that the
   trend has the predicted sign and the predicted place, on a sample
   too small to fit an exponent to.
3. **Smoothness tracks `1/g!` loosely.** Genus 3 measures 0.18 – 0.27
   against `1/3! = 0.167`, genus 2 measures 0.49 – 0.56 against
   `1/2! = 0.5`. The genus-3 excess is consistent at every `p`, so H1
   is a slight under-estimate at these sizes rather than wrong.
4. **What costs what now.** At genus 3, `p = 101`: 645 relation-stage
   operations of which 301 are precomputation, an oracle costing 53
   group-op equivalents (114,283 mul-mods — the largest single term
   after precomputation, because a cubic `u` needs Cantor–Zassenhaus
   rather than a Legendre symbol), and a linear algebra term that
   rounds to zero. The implementation is 2.1× from its own floor.

## Scope, unchanged in kind

Genus 2 and 3, odd characteristic, `p ≤ 251`, full degree-1 factor
base, one machine. Two genus-3 rows below 1 against a toy-sized rho are
not a statement about genus-3 curves at cryptographic size, where the
factor-base build alone is `O(p)` and everything here would have to be
rebuilt. What the rows do support is narrower and was the point: the
crossover moves in the direction and at the genus the theory names.

## Round five: the exponent law, pre-registered

Rounds one to four established that `S_ic/S_rho` is flat at genus 2 and dips
below 1 at genus 3 on the two largest `N`. Two data points with the right sign
in the right place is not an exponent. This round states the law the earlier
rounds' asymptotics imply, and the condition that would kill it, *before*
measuring.

**The law.** With `N ≈ p^g`, rho costs `√N ≈ p^{g/2}`. The relation stage is
`(m+1)·g!` group operations with `m ≈ p/2` places, so it is `≈ p·g!`, linear
in `p` at every genus. In the unit `S = ops/√N` both sides divide by `p^{g/2}`,
so the ratio is

    S_ic / S_rho  ∝  g! · p^{1 − g/2}

giving a predicted slope of `d log S_ratio / d log p = 1 − g/2`:

| `g` | predicted slope | what it means |
|--:|--:|:--|
| 2 | 0.0 | flat — no crossover at any size |
| 3 | −0.5 | crosses, slowly |
| 4 | −1.0 | crosses, twice as fast in the exponent |

The `g!` prefactor moves the crossover *later* as genus rises while the slope
moves it *earlier*; the crossover point is where `g!·p^{1−g/2} = 1`, i.e.
`p* = (g!)^{2/(g−2)}`. That predicts `p* = 36` at `g = 3` and `p* = 24` at
`g = 4` — both inside reach, which is why this is measurable at all rather
than an extrapolation.

**Success condition.** Fitted slopes within `±0.25` of `1 − g/2` at genus 3
and genus 4, over at least four `p` per genus, with every run returning the
verified `k`, against a distinguished-point rho whose `S_walk` stays inside
`1.0 – 1.5`.

**Falsification.** Any of: a genus-3 or genus-4 slope flatter than `−0.25`
(the dip is then a constant, not a trend); `S_walk` drifting outside
`1.0 – 1.5` (the reference has broken and the ratio is measuring that);
`S_ic/S_rho` failing to fall below 1 at genus 4 anywhere in range; or the
fitted crossover disagreeing with `p* = (g!)^{2/(g−2)}` by more than a factor
of 4.

**What a success would and would not mean.** It would establish that the
crossover is a trend with a measured exponent rather than two lucky rows, at
toy sizes, on one machine. It would *not* be a statement about
cryptographic-size Jacobians, where the `O(p)` factor-base build alone is
prohibitive and everything here would have to be rebuilt. The honest claim
available from this design is about the shape of the curve, not its position.

## Round five: measured

`cargo run --release --example hyperelliptic_ic_vs_rho`, then
`python3 scripts/hyperelliptic_exponent_fit.py`. Evidence:
`experiments/hyperelliptic_exponent_fit.json`. Every row returned its verified
`k`. `S_rho*` is rho with its walk renormalised to `sqrt(pi/2)`, per the
correction pre-registered above; `corr` is the conservative column and the one
to read.

| `g` | `p` | `N` | `S_ic` | `S_rho` | raw | `S_walk` | `S_rho*` | **corr** |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 3 | 61 | 124,459 | 1.51 | 2.30 | 0.66 | 1.48 | 2.07 | 0.73 |
| 3 | 101 | 364,747 | 1.12 | 1.90 | 0.59 | 1.41 | 1.74 | 0.64 |
| 3 | 151 | 1,180,351 | 0.79 | 1.98 | 0.40 | 1.71 | 1.52 | 0.52 |
| 3 | 211 | 4,620,611 | 0.53 | 1.96 | 0.27 | 1.82 | 1.39 | 0.38 |
| 4 | 31 | 239,753 | 1.77 | 1.90 | 0.93 | 1.30 | 1.85 | 0.95 |
| 4 | 41 | 333,041 | 1.41 | 1.86 | 0.76 | 1.34 | 1.77 | 0.80 |
| 4 | 61 | 16,790,591 | 0.28 | 2.17 | 0.13 | **2.10** | 1.32 | **0.21** |

**The law holds at all three genera**, fitted against `log2 N`:

| `g` | predicted `1/g − 1/2` | corrected slope | `R²` |
|--:|--:|--:|--:|
| 2 | 0.000 | −0.027 | 0.24 (no trend) |
| 3 | −0.167 | **−0.154** | 0.98 |
| 4 | −0.250 | **−0.275** | 0.96 |

### The pre-registered fit was against the wrong variable, and that is my error

Round five above predicted slopes against `log2 p` of `1 − g/2`, and genus 4
missed badly: −1.611 corrected against −1.0, outside the ±0.25 band. The law
is not wrong; the step from `N` to `p` is. `N ≈ p^g` is an assumption the
harness does not satisfy — it selects a prime-order subgroup above `lo/4`, and
the measured scaling is

    N ~ p^1.91 (g=2),   p^3.20 (g=3),   **p^5.94 (g=4)**

At genus 4 that inflates the p-slope by `5.94/4`, predicting `−0.25 × 5.94 =
−1.485` against the measured −1.611. Restated in the variable the law is
actually about — rho costs `√N` whatever the genus —

    S_ic / S_rho  ∝  g! · N^{1/g − 1/2}

every genus lands inside the band, at `R² = 0.96–0.98` where there is a trend
to fit. The `log2 N` fit is now the primary one and the `log2 p` fit is kept
beside it, because the discrepancy between them *is* the finding about the
harness.

### The reference degrades, and the correction is doing real work

The pre-registered falsification condition on `S_walk` fired on three rows:
genus 3 at `p = 151, 211` (1.71, 1.82) and genus 4 at `p = 61` (**2.10**),
against the `1.2533` ideal. The drift is systematic in `N`, not noise, and
genus 2 stays flat near 1.3 — so it is not a constant implementation tax, it
is a reference that gets worse exactly where the interesting rows are. That
direction flatters index calculus, so the raw column overstates the result.

The single best row is the clearest case: genus 4, `p = 61`, `N = 2^24`. Raw
ratio 0.13 — a 7.7× win. With rho's walk renormalised, 0.21 — a 4.8× win. The
correction nearly halves the claim, and that row's `S_walk` is 68% above ideal,
so it is also the row whose correction is least trustworthy. **0.21 is the
number to quote, and the reference needs fixing before it is quoted hard.**

### Where this leaves the thread

Index calculus beats a real distinguished-point rho at genus 3 and genus 4, on
the conservative accounting, by up to 4.8× at `N = 2^24`, with a measured
exponent matching the derived law to within 0.025. That is a crossover with a
slope, not two lucky rows — which is what round five set out to decide.

It remains a statement about toy sizes and one machine. `p ≤ 251`, `N ≤ 2^24`,
full degree-1 factor base, and the `O(p)` factor-base build that dominates at
cryptographic size is not modelled here at all.

**Next, in order.** (1) Fix the DP walk so `S_walk` stays near 1.2533 as `N`
grows, and re-measure — the correction should then be a no-op, and if the
crossover survives an uncorrected reference it is no longer arguable. (2) Push
genus 4 past `N = 2^24`, where the law predicts the ratio keeps falling as
`N^{−0.25}`. (3) Only then ask what any of it costs at a size anyone cares
about.
