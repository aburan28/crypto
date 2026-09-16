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
