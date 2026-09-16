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
