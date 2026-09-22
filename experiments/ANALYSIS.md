# Isogeny-Class ECDLP Experimental Analysis

**Suite:** `cargo run --release -- isogeny ...`
**Date:** 2026-05-16 — 2026-05-18 (two follow-up rounds appended)
**Module:** `src/isogeny/` (ported from `claude/focused-jennings-553ded` worktree)

## Update log

- **Round 1 (initial)**: rho hits cap ~60% of vertices; bit ceiling
  effectively 18-bit due to brute-force point counting and Vélu
  enumeration.
- **Round 2 (2026-05-17, Pollard-ρ sterile-collision fix in
  [`src/cryptanalysis/pollard_rho.rs`](src/cryptanalysis/pollard_rho.rs:189))**:
  recover partial answer via CRT-style brute force over the gcd
  coset when `gcd(rhs, n) > 1`.  Lifts success rate from 39 % to
  98–100 % across 14/16/18-bit.
- **Round 3 (2026-05-18, scale-up)**: replaced `O(p)` brute-force
  Frobenius trace with Shanks-Mestre `O(p^{1/4})` BSGS in
  [`src/isogeny/cm.rs`](src/isogeny/cm.rs); replaced the `O(p²)`
  brute-force point enumeration inside `velu_isogeny_odd` and
  `compose_isogenies` with `O(p · log² p)` Tonelli-Shanks point
  sampling.  Bit ceiling lifted from 18 to 22-bit at full 4-trial
  experiment cost, 26-bit (4 trials) feasible at ~1 hr.
- **Round 10 (2026-09-22, the 22–50-bit rows, and three tables round 9
  left behind)**: round 9 stopped at 14/16/18 bits because it recorded
  that "no frozen artifact exists" for the wider rows.  That was wrong —
  `13_22bit_dpoly`, `14b_30bit`, `15b_40bit` and `16c_50bit` were all
  committed, on the same `seed: 0xC0FFEE` config family, cited by nothing.
  All four are regenerated here from their own `config` blocks, starting
  curves verified unchanged, plus the two low-cap companions.  22-bit
  `30 → 9` vertices, 30-bit `5 → 5` (median `38 266 → 16 352`), 40-bit
  `15 → 8`, 50-bit `14 → 32`.
  Round 9 also updated §3.2, §3.4 and §3.4.1 and left §3.2.1, §3.3 and
  §1's parity table on pre-fix numbers — §3.2.1 claimed 24 successes at 14
  bits where §3.2 counts 16 vertices in total.  The cause was structural:
  `analyse_isogeny_vertices.py` computed exactly the tables round 9
  checked.  It computes all six now, and the pooled CSV covers all seven
  widths.
  **Three findings changed, not just their numbers.**  (i) The even-`n`
  catastrophic-failure story in §1 is gone: even-`n` curves succeed at
  88–100 %, because round 2 repaired the sterile-collision restart that
  caused it — and that table's headers had the parities backwards.
  (ii) §3.2.1's "within 15 % of the textbook prediction" does not survive;
  conditioning on odd `n` no longer does any work now that even-`n` curves
  succeed.  (iii) §3.4 gains its best evidence yet: a 16-sample 50-bit
  class at CoV `0.542` against the generic prediction of `0.52`.
  **Two limits recorded rather than papered over:** at 50 bits both classes
  saturate `max_graph_nodes = 16`, so class size there is a censored lower
  bound; and §1/§5 had described a `bits.clamp(8, 30)` that round 3 raised
  to `62`, which made the 40- and 50-bit rows read as secretly-30-bit runs.
- **Round 9 (2026-09-21, Vélu `v(Q)` correction — §3 data
  regenerated)**: [`src/isogeny/velu.rs`](src/isogeny/velu.rs) used
  `v(Q) = 4 x_Q g_y(Q)² + 2 g_x(Q)²` for non-2-torsion kernel points
  where Vélu's value is `2 g_x(Q)` (`g_x(Q)` on the T-row).  Codomains
  were wrong, and **wrong in a way that broke isogeny**: on `toy-a` the
  recorded 3-isogenous curve had `#E = 2059` against the start curve's
  `2019`, so by Tate it was not in the isogeny class at all.  Every §3
  figure below rests on graph vertices, so all of them were regenerated
  on the fixed code; superseded values are kept inline as "was".
  **What this does not license:** the regenerated runs are today's
  code, so the deltas conflate the Vélu fix with rounds 2–4 (the rho
  gcd-recovery in particular lifts success from 39 % to ~100 %).  They
  are a corrected baseline, not a measurement of the Vélu fix alone.
- **Round 8 (2026-05-20, Pohlig-Hellman + multi-target GLS)**: added
  [`pohlig_multi_target_rho`](src/isogeny/attack.rs) — wraps the
  round-7 `multi_target_dp_rho` in a Pohlig-Hellman outer loop.
  Factors `n` into prime powers, lifts the base & targets into each
  order-`q^e` subgroup via scalar-multiplying by `n/q^e`, runs GLS
  there, CRT-combines the per-prime residues to recover the secrets
  mod `n`.  When GLS stalls inside a subgroup (e.g. for the
  `q = 2^e` factor of a smooth `n`), the loop falls back to per-
  target round-6 dp-parallel rho in that subgroup.

  Tests:
  * `pohlig_multi_target_rho_solves_m4_robustly` — exercises the
    16-bit `n = 2 · prime` setup at **all** seeds that failed in
    round 7.  Round 8 solves them cleanly.

  Optional dispatch: setting `CRYPTO_ISOGENY_PH_GLS=1` routes
  `rho_on_curve` through `rho_via_pohlig_gls`, which picks `m=4`
  random secrets and reports the per-target iteration count.

  **Empirical result is mixed.**  At 22-bit (4 trials):

  | dispatch | mean iters / target | wall |
  |----------|---------------------:|------:|
  | round-6 single-target | 3 219 | 1.0 s |
  | round-8 PH + GLS m=4  | 14 818 | 0.5 s |

  Wall time halves, but the *per-target* iter count goes up 4.6×.
  Reason: each multi-target step updates m+1 coefficients
  `(a, b_0, …, b_{m-1})` instead of just `(a, b)`, so the per-step
  cost rises faster than the cycle-length saving at this scale.
  The √m speedup from GLS is real asymptotically but doesn't
  materialise yet at 22-bit.  Likely needs (a) a tighter inner
  loop with fewer BigUint allocations, or (b) larger `n` where the
  cycle-length saving dominates the coefficient overhead.

  Round 8 is therefore **opt-in only** — the default sweep still
  uses round-6 DP single-target rho, which delivered the 14 → 60
  bit progression already in the record.

- **Round 7 (2026-05-20, multi-target GLS amortisation, partial)**:
  added
  [`src/isogeny/attack.rs::multi_target_dp_rho`](src/isogeny/attack.rs) —
  a Galbraith-Lin-Scott multi-target rho.  K=4 parallel walkers
  track state `(a, b_0, …, b_{m-1})` such that the running point
  equals `a·g + Σ b_i·h_i`.  Collisions yield linear equations in
  `(d_0, …, d_{m-1})`; an incremental Gaussian elimination
  (`RrefState`) reduces each new equation against existing pivots
  in `O(m²)`.  Theoretical per-target speedup: **√m**.

  **The catch — rank stalling.**  RrefState requires every pivot
  to be invertible mod `n`.  When `n = 2 · p` (which is the case
  for *every* curve in our 16-bit experiment because
  `p ≡ 1 (mod 2) ⇒ n = p + 1 − t` and `p + 1` is even), some
  pivot column never gets a row with an odd entry, and rank stalls
  at `m − 1`.  This is *not* a coding bug — it's a Smith-normal-
  form issue: the random adders' span in `(Z/n)^{m+1}` happens to
  miss the right coset modulo 2 for some seeds.  Empirical sweep:

  | m | seeds that solve in < 10k iters | seeds that stall at rank m−1 |
  |---|--------------------------------:|-----------------------------:|
  | 1 | all                             | none                         |
  | 2 | most                            | rare                         |
  | 3 | most                            | rare                         |
  | 4 | ~30 %                           | ~70 %                        |

  Round-7 commits the infrastructure (multi-target rho + RrefState
  with full correctness tests) and a single 16-bit lucky-seed test
  showing the algorithm is *correct* when rank progression
  succeeds.  The rho dispatch in `rho_on_curve` is **unchanged** —
  the experiment still uses round-6 single-target DP rho.  Fixing
  the rank-stalling problem requires either (a) factoring `n` and
  solving mod each prime factor with CRT, or (b) generating
  adders constrained to span the full Smith-normal-form basis.
  Both are tractable but out of scope for this round.

- **Round 6 (2026-05-19, distinguished-points + K parallel walkers)**:
  added
  [`src/isogeny/attack.rs::dp_parallel_rho_on_curve`](src/isogeny/attack.rs)
  using K=4 OS threads sharing a `HashMap`-backed DP table.  Each
  walker runs the r-adding step from round 5; emits a DP whenever
  the low `dp_bits ≈ log₂(n)/4` bits of `x` are zero; checks the
  shared table for a cross-walker collision.  Closes the long-tail
  problem that bit round 5 at 50-bit.

  | bits | round 5 (r=20, single) | round 6 (DP + K=4)  | improvement |
  |------|-----------------------:|--------------------:|------------:|
  | 30   | 3.2 s                  | 3.1 s               | (noise)     |
  | 40   | 6 m 40 s               | **53 s**            | **7.5×**    |
  | 50   | killed at 1 h 30+      | **33 m 58 s**       | **\>2.6×**  |
  | 60   | infeasible             | **~14 hr**          | new scale (was unreachable) |

  Notes:
  * Median-ratio dropped a little (0.51–0.68 vs 0.66–0.93 in round 5).
    This is consistent with the K-walker cross-collision effect:
    walkers find each other's DPs faster than one walker finds its
    own.  Mean iterations is ~25 % below √(πn/2).
  * Wall-time at 40-bit is 53 s for 4 trials × 15 vertices = 60
    rho calls.  That's 0.88 s per rho on average, against
    `√(πn/2) ≈ 1.3 M` iters per rho — confirming the parallel
    walkers add up to roughly 4× CPU utilisation in addition to
    the algorithm-level √K saving.
  * **60-bit headline numbers** (2 starting curves, 11 vertices):
    * `n ≈ 1.15 × 10¹⁸`
    * `√(πn/2) ≈ 1.35 × 10⁹`
    * Observed median: `1.08 × 10⁹` (**ratio 0.81**)
    * Observed mean: `1.05 × 10⁹` (mean ratio 0.78)
    * CoV: 0.44 — *lower* than at any of the smaller scales,
      consistent with the K-walker variance reduction
    * Success rate: **11/11**
    * Wall: ~14 hr at 4 cores

- **Round 5 (2026-05-19, r-adding walk)**: replaced the 3-partition
  Pollard ρ (default in
  [`src/cryptanalysis/pollard_rho.rs`](src/cryptanalysis/pollard_rho.rs))
  with a 20-bucket r-adding walk written inline in
  [`src/isogeny/attack.rs::r_adding_rho_on_curve`](src/isogeny/attack.rs).
  Bucket choice uses a splitmix64-style mixer on `(x, y)` to avoid
  collisions with small-order torsion structure.  Round-2's gcd-
  recovery branch is preserved.

  **Result: net improvement is real but uneven.**

  | bits | r=3 wall (round 4) | r=20 wall (round 5) | improvement |
  |------|-------------------:|--------------------:|------------:|
  | 18   | 2 m 25 s           | 0.4 s               | 360×        |
  | 22   | 1 s                | 1.1 s               | (noise)     |
  | 30   | 2 s                | 3.2 s               | (noise; tail) |
  | 40   | 13 m 07 s          | **6 m 40 s**        | **2.0×**    |
  | 50   | ~1 hr              | killed at 1 h 30    | (regressed) |

  The 18-bit gain is mostly attributable to the dropped restart
  cost; 30/50-bit show a long-tail problem where a single starting
  curve's class has a bad walk and the harness has no parallel
  walker to share work.  The conclusion: **r-adding alone is not
  enough** to push past 50-bit cleanly — the next required upgrade
  is distinguished-points with K parallel walkers, which would
  amortise the long tail across cores.

- **Round 4 (2026-05-18, division-polynomial replacement)**: replaced
  the `O(p · log² p)` point enumeration in `velu_isogeny_odd` with
  random-point cofactor sampling (`m = #E / ℓ`, `Q = m·P`, repeat
  until all distinct cyclic subgroups found).  Same trick for
  `ℓ = 2` in `graph.rs::two_isogenies_via_cubic_roots`.  Per-isogeny
  cost is now `O(ℓ² · log² p)` instead of `O(p · log² p)`.  At
  `p = 2^22`, that's a `10⁶×` speed-up; 22-bit full sweep dropped
  from 5 min to **1 s**; 30-bit (was previously infeasible)
  completes in 2 s; 40-bit completes in 13 min with the rho cap
  scaled to `2³ · 2^{bits/2}` ≈ 8.4M iters.  This is the
  *asymptotically correct* implementation of `ℓ`-torsion kernel
  finding; the cofactor sampling approach is what Schoof's paper
  describes as the "natural" alternative to full division
  polynomials when `#E` is already known.

## 0. Module status & merge note

The isogeny module was uncommitted work in the `claude/focused-jennings-553ded`
worktree.  That branch was several hundred commits behind `main`, so a real
`git merge` was impractical — I ported the additions surgically:

* copied `src/isogeny/` (9 files: `mod`, `class_group`, `cm`, `velu`,
  `volcano`, `graph`, `attack`, `experiment`, `secp256k1_analysis`)
* added `pub mod isogeny;` to [src/lib.rs](src/lib.rs)
* added `rayon = "1.10"` to [Cargo.toml](Cargo.toml) (used by
  `experiment::run_experiment` for the per-curve parallel sweep)
* wired the `Isogeny { op }` subcommand and `cmd_isogeny` dispatcher
  into [src/main.rs](src/main.rs:96)

`cargo check` and `cargo build --release` both pass clean.

## 1. Experiment scale and a critical caveat

The CLI accepts `--bits 64`.  Round 3 lifted the ceiling: the clamp now
lives in `generate_random_curves` and reads `bits.clamp(8, 62)`, the bound
being `p < 2^62` so the `u128` products inside the `SmallCurve` helpers
cannot overflow, and `next_prime_near_2pow` no longer clamps at all.

*(This paragraph said `bits.clamp(8, 30)` and "every bit-width past 30 is
silently capped to 30" until round 10.  That was true of round 1's code and
was left standing through the round-3 scale-up.  It is worth naming because
it made the 40- and 50-bit rows below read as secretly-30-bit runs: they are
not, and the `actual bits` column is now measured from each artifact's own
prime rather than asserted.)*

With BSGS point counting the runs complete well past 2^18 — 50-bit takes
45 minutes — so the table covers every width that has a frozen artifact:

| bits requested | actual bits | wall time | n succ / n total | notes |
|----------------|-------------|-----------|------------------|-------|
| 14             | 14          | < 1 s     | 16 / 16          | cap=2^14 (was 35 / 89) |
| 14             | 14          | < 1 s     | 16 / 16          | cap=2^18 — still no extra successes (was 35 / 89) |
| 16             | 16          | 3 s       | 70 / 73          | cap=2^14 (was 30 / 80) |
| 18             | 18          | < 1 s     | 25 / 30          | cap=2^14 (was 41 / 106) |
| 22             | 22          | < 1 s     | 8 / 9            | cap=2^18 (was 30 / 30) |
| 30             | 30          | < 1 s     | 5 / 5            | cap=2^18 (was 5 / 5, median 38 266 → 16 352) |
| 40             | 40          | 12 s      | 8 / 8            | cap=2^23 (was 14 / 15, 13 m 07 s) |
| 50             | 50          | 45 m 02 s | 32 / 32          | cap=2^28 (was 14 / 14, ~1 hr) |
| 64             | clamped→62  | not run   | –                | clamp is 62 now, not 30; untested at this size |

**The even-order failure mode is gone, and this section used to say
otherwise.**  Round 1 found rho failing on even-`n` curves far more often
than on odd-`n` ones and traced it to
[src/cryptanalysis/pollard_rho.rs:189](src/cryptanalysis/pollard_rho.rs:189):
on a collision the code checks `rhs.gcd(n) != 1`, and where the gcd is
non-trivial it can only recover the log mod `n / gcd`, so it declared the
collision *sterile* and restarted.  With `n` even about half of all
differences `(h_b − t_b)` are even, so about half of collisions were
sterile, and after `max_restarts = 8` rho gave up well short of the cycle.

Round 2 fixed it — the gcd coset is brute-forced CRT-style instead of
discarded — and on the regenerated data the split has closed at every
width measured:

| bit-width | n odd (trace odd) | n even (trace even) |
|-----------|-------------------|---------------------|
| 14-bit    | 5 / 5 (100%)      | 11 / 11 (100%)      |
| 16-bit    | 64 / 67 (96%)     | 6 / 6 (100%)        |
| 18-bit    | 5 / 6 (83%)       | 20 / 24 (83%)       |
| 22-bit    | 1 / 1 (100%)      | 7 / 8 (88%)         |
| 30-bit    | 3 / 3 (100%)      | 2 / 2 (100%)        |
| 40-bit    | 4 / 4 (100%)      | 4 / 4 (100%)        |
| 50-bit    | 16 / 16 (100%)    | 16 / 16 (100%)      |

Even-`n` curves now succeed at 88–100 %, against the 7–18 % this section
reported.  Where the two parities differ at all (18-bit, 83 % both ways)
they differ together, which is the signature of the iteration cap rather
than of a parity effect.

*(was, on the buggy Vélu and pre-round-2 rho: 14-bit `11 / 60` and
`24 / 29`; 16-bit `5 / 54` and `25 / 26`; 18-bit `5 / 69` and `36 / 37`,
read as "EVEN-order curves fail catastrophically (7-18 % success)".  Two
things were wrong with it by round 10 and neither was caught in round 9:
the numbers predate both the rho fix and the Vélu fix, and the column
headers paired the parities backwards — for `p` odd, `n = p + 1 − t` with
`p + 1` even, so `n` parity **equals** `t` parity, which the paragraph
under the old table derived correctly while the header above it said the
opposite.)*

What survives from the original observation is the narrower claim it was
evidence for: the residual failures are not "rho cycle exceeded √(πn/2) by
orders of magnitude".  Raising the cap from `2^14` to `2^18` still converts
no extra 14-bit failure into a success.  The sterile-collision restart was
an implementation artifact of the generic `pollard_rho_dlp` module, never a
cryptographic property of the curves — and it is now repaired rather than
merely diagnosed.

## 2. Single-curve probes (positional experiments)

### `isogeny secp256k1` — full GLV / MOV / structural certificate

```
β:       7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
         β³ ≡ 1 (mod p): true
λ:       5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72
         λ² + λ + 1 ≡ 0 (mod n): true
embedding degree: ≥ 200  → MOV/Frey-Rück infeasible
quadratic twist:  not pairing-friendly for k ≤ 32

Small-ℓ structural survey over Z[ω], disc = -3:
  ℓ=2,5,11,17,23 inert
  ℓ=3            ramified
  ℓ=7,13,19      split
```

All four security indicators (β³≡1, λ²+λ+1≡0, MOV infeasibility, twist
safety) check out.  The split primes 7,13,19 in Z[ω] match the standard
GLV literature for secp256k1 — these are the exact ℓ at which the
endomorphism ring Z[(1+√-3)/2] has degree-1 prime ideals, giving rise to
the 3-isogeny endomorphism that powers the GLV decomposition.

### `isogeny volcano --curve toy-a --ell 2`

```
y² = x³ + 4x + 4  over F_2003
j(E) = 84,  trace = -15,  #E = 2019
End(E) disc = -7787  (fundamental, conductor 1)
Position: depth 0, crater_size 1, on_crater = true
```

The curve sits on a degenerate volcano: a single crater vertex, no descent
edges.  This is consistent with `Δ_E = -7787` being a *prime* fundamental
discriminant — there is no ℓ=2 prime above it (Legendre `(-7787/2)` shows
2 inert), so the 2-volcano is the trivial single-vertex graph at this
curve.

### `isogeny volcano --curve toy-j0 --ell 3`

```
y² = x³ + 1  over F_103   (j = 0)
trace = 20,  #E = 84,  End(E) disc = -12
fundamental disc = -3, conductor 2
Position: depth 0, crater_size 1, on_crater = true
```

j = 0 over F_103 has CM by Z[ζ₃] (disc -3 with conductor 2 → -12).
The 3-volcano collapses because 3 is *ramified* in Q(√-3): the unique
prime ideal above 3 acts trivially on isogenies of degree 3 from j=0.
This is the well-known j=0 obstruction (e.g. SIKE's avoidance of j=0
starting nodes).

### `isogeny graph --curve toy-a --ell-list 2,3`

```
Vertices: 2,  Edges: 1,  Components: 1,  Diameter: 1
[0] j=84   -> 1 (ℓ=3)
[1] j=676  -> 0 (ℓ=3)
```

The (2,3)-graph adds one 3-edge: at ℓ=3 there *is* a horizontal isogeny.
This matches the per-discriminant prediction — `class_above_prime(-7787, 3)`
returns a non-trivial form, so the 3-action is non-trivial and yields one
neighbour.

### `isogeny cm --discriminant -23`

```
h(-23) = 3                      (class number)
Reduced forms: [1,1,6], [2,1,3], [2,-1,3]
Cl(O) ≅ Z/3                     (group table is cyclic)
Above ℓ:  ℓ=2 → [2,1,3]         (split, non-principal)
          ℓ=3 → [2,-1,3]        (split, non-principal — the conjugate class)
          ℓ=13 → [2,-1,3]       (split, congruent to ℓ=3 in Cl(O))
          ℓ=5,7,11 inert/ramified
```

D = -23 is the smallest negative fundamental discriminant with h>1.  The
class group is Z/3 generated by either of the two non-principal forms
(which are inverse to each other).  Crucially the action on isogenies
of degree ℓ=2 and ℓ=3 are *different* generators of the same Z/3 — so the
2-volcano and 3-volcano share the same crater but differ in which edge is
travelled.  This is exactly the structure that powers CSIDH-style group
actions.

## 3. Aggregate rho-cycle statistics across the isogeny class

I ran the experimental harness at `bits ∈ {14, 16}` with 10 random starting
curves per run, primes={2,3}, max graph nodes per class = 16.  The harness
runs Pollard ρ at every vertex of the discovered class.

### 3.1 Per-curve table (14-bit, cap 2^18)

```
 #  j_E    |cls| trace  order   f_disc cond  succ/N  min  med   max  med/√(πn/2)
 0   7057    1    -59   16471    -6907   3  1/1      36    36    36  0.224
 1   1224    1   -185   16597    -3491   3  1/1      59    59    59  0.365
 2   8770    2    126   16286   -49768   1  2/2      53    90   128  0.566
 3  12331    2    107   16305   -54195   1  2/2      64    70    76  0.437
 4   4044    2    -42   16454   -63880   1  2/2      45    65    85  0.404
 5  10350    1    184   16228     -883   6  1/1     220   220   220  1.378
 6   2652    2     64   16348   -15387   2  2/2      48    66    84  0.412
 7  13652    2   -138   16550    -1864   5  2/2      29    34    38  0.208
 8  12300    2   -188   16600     -303  10  2/2      18    21    24  0.130
 9   7640    1     21   16391   -65203   1  1/1       7     7     7  0.044
```

### 3.2 Pooled statistics (successful rho runs only)

| bits | n_succ | n_fail | median iters | √(πn/2) | median ratio | mean ratio | CoV   |
|------|--------|--------|--------------|---------|--------------|------------|-------|
| 14   | 16     | 0      | 50           | 161     | **0.314**    | 0.40       | 0.81  |
| 16   | 70     | 3      | 196          | 320     | **0.613**    | 0.62       | 0.69  |
| 18   | 25     | 5      | 224          | 642     | **0.349**    | 0.51       | 0.98  |

*(was, on the buggy Vélu: 14 → `35 succ / 54 fail`, median ratio
`1.021`; 16 → `30 / 50`, `1.029`; 18 → `41 / 65`, `1.103`.)*

After the round-2/round-3/round-4 fixes (rho gcd-recovery + BSGS
point counting + cofactor-sampling ℓ-isogeny construction):

| bits | cap  | vertices | succ | median iters | √(πn/2)   | ratio | wall   |
|------|------|----------|------|--------------|-----------|-------|--------|
| 14   | 2^14 | 16       | 16   | 50           | 161       | 0.31  |  < 1 s |
| 16   | 2^14 | 73       | 70   | 196          | 320       | 0.61  |  3 s   |
| 18   | 2^14 | 30       | 25   | 224          | 642       | 0.35  |  < 1 s |
| 22   | 2^18 | 9        | 8    | 802          | 2 567     | 0.31  |  < 1 s |
| 30   | 2^18 | 5        | 5    | 16 352       | 41 068    | 0.40  |  < 1 s |
| 40   | 2^23 | 8        | 8    | 817 327      | 1 314 195 | 0.62  |  12 s  |
| 50   | 2^28 | 32       | 32   | 34 279 439   | 42 054 244 | 0.82 |  45 m 02 s |

*(was, on the buggy Vélu: 22 → `30 / 30`, median `1 384`, ratio `0.54`;
30 → `5 / 5`, median `38 266`, ratio `0.93`; 40 → `15 / 14`, median
`870 880`, `13 m 07 s`; 50 → `14 / 14`, median `30 200 000`, ratio `0.72`,
`~1 hr`.  And earlier still, at 14/16/18: `89 / 87`, `80 / 80`,
`106 / 105` vertices.)*

**All seven rows are now regenerated on the corrected Vélu.**  Round 9
recorded here that "no frozen artifact for them exists under
`experiments/`" and left the 22–50-bit rows stale on that basis.  That was
wrong: `13_22bit_dpoly.json`, `14b_30bit.json`, `15b_40bit.json` and
`16c_50bit.json` were all committed, carrying the same
`seed: 0xC0FFEE`, `primes: [2,3]`, `max_graph_nodes: 16` config family as
the 14/16/18 runs.  They were simply cited by nothing — §6 listed neither
them nor the rows they fed, which is how they went unnoticed.  They are
listed there now.

Each regeneration used its own artifact's `config` block verbatim, and the
`start` curve of every trial is unchanged, which is the check that only the
class walk moved.

**The cap column is new and the rows are not comparable without it.**  A
rho cap has to scale with √n to be reachable at all, so the 14/16/18 runs
cap at `2^14` where the 50-bit run caps at `2^28`.  A success rate only
means something against another row at the same cap; the old table hid
this by omitting the column.

**At 50 bits the class size is censored, not measured.**  Both starting
curves' classes hit `max_graph_nodes = 16` exactly (`3` and `11` before),
so "32 vertices" is two truncated walks, and the true classes are only
known to be `≥ 16`.  This is the same cap that round 9 found the *broken*
walk saturating at 14 and 18 bits — there it was the bug manufacturing
vertices, here it is a real class outgrowing the budget, which is what a
class number growing like √p should do.  Raising `max_graph_nodes` is the
obvious next round and is not done here.

The post-fix median ratio runs 0.31–0.82 of the textbook √(πn/2)
expectation, and **rises with `bits`**: 0.31, 0.61, 0.35, 0.31, 0.40,
0.62, 0.82 at 14 → 50.  The mechanism the earlier rounds proposed still
fits: at small `bits` the gcd-recovery branch reports tiny iteration
counts when the chosen generator sits in a 2- or 3-torsion subgroup,
pulling the median down, and that contamination thins as `n` grows.  The
50-bit row is the cleanest look at the asymptote and sits at 0.82.

*(was: "0.5–0.9 … the ratio approaches the theoretical 1.0 (0.93 at
30-bit)".  On the regenerated data 30-bit is 0.40, not 0.93 — the old
30-bit median of 38 266 came from a class walk that was partly outside the
isogeny class.  The trend with `bits` is the same; the width that
best shows it is now 50, not 30.)*

The 22→30 bit transition is no longer the interesting one: both complete
in under a second where the v1 brute-force implementation would have taken
months.  The cost wall is now at 50 bits, at 45 minutes.

**50-bit**: 32/32 vertices succeed.  The order is `n ≈ 1.13 · 10¹⁵`; the
median rho cycle of 34.3 M iters matches √(πn/2) ≈ 42 M to a ratio of
0.82, well within the geometric-distribution noise floor.  All 32
successes completed below the `8 · 2^25 ≈ 268 M` cap — which is exactly
the CLI's own `8·2^{bits/2}` default — confirming the cap-scaling formula
in `cmd_isogeny` is calibrated correctly for this regime.  What changed
from the previous run is the class structure, not the rho cost: the
`3`-vertex and `11`-vertex classes reported here were an artifact of the
broken walk, and both classes now saturate the 16-node budget.

### 3.2.1 Pooled stats CONDITIONAL ON n=odd (sterile-collision bug avoided)

After we exclude the even-n curves whose failure rate is a partition-bug
artifact:

| bits | n_succ | total odd-n | median iters | √(πn/2)    | median ratio | mean ratio | CoV   |
|------|--------|-------------|--------------|------------|--------------|------------|-------|
| 14   | 5      | 5 (100%)    | 59           | 160        | **0.368**    | 0.30       | 0.57  |
| 16   | 64     | 67 (96%)    | 196          | 320        | **0.613**    | 0.62       | 0.70  |
| 18   | 5      | 6 (83%)     | 273          | 642        | **0.425**    | 0.94       | 0.91  |
| 22   | 1      | 1 (100%)    | 794          | 2 566      | **0.309**    | 0.31       | n/a   |
| 30   | 3      | 3 (100%)    | 11 844       | 41 069     | **0.288**    | 0.35       | 0.66  |
| 40   | 4      | 4 (100%)    | 889 123      | 1 314 196  | **0.677**    | 0.66       | 0.67  |
| 50   | 16     | 16 (100%)   | 34 279 439   | 42 054 244 | **0.815**    | 0.89       | 0.54  |

*(was, on the buggy Vélu: 14 → `24` of `29` odd-n, median ratio `0.859`,
CoV `1.07`; 16 → `25 / 26`, `1.129`, `1.04`; 18 → `36 / 37`, `0.950`,
`0.75`.  **Those numbers were impossible against their own neighbour** —
§3.2 reports 16 vertices in total at 14 bits, so 24 odd-`n` successes
cannot exist.  Round 9 regenerated §3.2 and left this table, §3.3 and §1's
parity table untouched, because `analyse_isogeny_vertices.py` computed
only the tables round 9 checked.  It computes all six now.)*

**The headline this section used to carry does not survive.**  It read
"median rho-iteration count is within 15 % of the textbook prediction at
all three scales", on median ratios of 0.86 / 1.13 / 0.95.  On the
corrected data the median ratios are 0.29–0.82 and the *conditioning on
odd `n` no longer does any work*: §1 shows even-`n` curves succeeding at
88–100 % since round 2 fixed the sterile-collision restart, so there is no
longer a contaminated subpopulation to exclude.  This table is now close
to a sub-sample of §3.2 rather than a correction to it, and the honest
reading is §3.2's: the ratio is below 1 and rises with `bits`.

The CoV is 0.54–0.91 where sample size allows one at all, against the
theoretical 0.52 for an un-truncated Pollard ρ.  The 50-bit row — 16
samples, CoV 0.54 — is the closest to theory and the best-sampled, which
is the opposite of the old table, where the largest CoVs sat at the widths
with the most (spurious) vertices.

### 3.3 Conductor stratification (conditional on n=odd)

| bits | conductor | n  | median     | mean       | CoV  |
|------|-----------|----|------------|------------|------|
| 14   | 1 (max-order) | 3  | 64         | 49         | 0.75 |
| 14   | 3         | 2  | 48         | 48         | 0.34 |
| 16   | 1 (max-order) | 64 | 196        | 200        | 0.70 |
| 18   | 1 (max-order) | 1  | 1 500      | 1 500      | 0.00 |
| 18   | 3         | 4  | 262        | 382        | 0.70 |
| 22   | 3         | 1  | 794        | 794        | 0.00 |
| 30   | 1 (max-order) | 3  | 11 844     | 14 406     | 0.66 |
| 40   | 1 (max-order) | 2  | 1 361 463  | 1 361 463  | 0.07 |
| 40   | 3         | 1  | 260 628    | 260 628    | 0.00 |
| 40   | 5         | 1  | 488 481    | 488 481    | 0.00 |
| 50   | 1 (max-order) | 16 | 34 279 439 | 37 597 512 | 0.54 |

*(was, on the buggy Vélu: 14-bit cond 1/3/9 at n = 10/11/3; 16-bit cond
1/19 at n = 23/2; 18-bit cond 1/3/7 at n = 17/16/3.  Same staleness as
§3.2.1 — those totals are pre-Vélu-fix and exceed the corrected vertex
counts.)*

**The corrected data cannot stratify by conductor.**  Only two cells in
the whole table have more than four samples — 16-bit cond=1 (n = 64) and
50-bit cond=1 (n = 16) — and both are max-order, so there is no
cond=1-vs-cond>1 comparison at any width with the samples to support one.
Every suborder cell is n ≤ 4.

The original conclusion happens to survive, but it is now an absence of
data rather than a measured null: **no statistically defensible signal**
that suborder curves are easier or harder than max-order curves within an
isogeny class.  Testing it needs classes large enough to contain several
conductors, which at these widths means raising `max_graph_nodes` and
re-running — the same upgrade §3.2 flags for the censored 50-bit classes.

### 3.4 Within-class variation

For starting curves with at least two successful odd-`n` vertices.
**The framing changed with round 9:** this section used to read "curves
whose ℓ ∈ {2,3} graph reached the 16-vertex cap", but on the corrected
Vélu no 14- or 18-bit class reaches the cap at all — the saturated
classes were the broken walk manufacturing vertices.  Genuine
cap-reaching classes survive only at 16 bits, which is why that width
now carries the bulk of the evidence.

14-bit:

| curve | \|cls\| | conductor | n_succ | mean | stdev | CoV |
|-------|-------|-----------|--------|------|-------|-----|
| #3 | 2 | 1 | 2 | 70 | 8 | 0.121 |

16-bit:

| curve | \|cls\| | conductor | n_succ | mean | stdev | CoV |
|-------|-------|-----------|--------|------|-------|-----|
| #0 | 16 | 1 | 16 | 212 | 114 | 0.541 |
| #1 | 16 | 1 | 15 | 194 | 121 | 0.624 |
| #7 | 16 | 1 | 14 | 158 | 127 | 0.805 |
| #8 | 16 | 1 | 16 | 238 | 184 | 0.775 |

18-bit:

| curve | \|cls\| | conductor | n_succ | mean | stdev | CoV |
|-------|-------|-----------|--------|------|-------|-----|
| #1 | 2 | 3 | 2 | 237 | 18 | 0.078 |
| #2 | 2 | 3 | 2 | 526 | 359 | 0.681 |

40-bit:

| curve | \|cls\| | conductor | n_succ | mean | stdev | CoV |
|-------|-------|-----------|--------|------|-------|-----|
| #2 | 2 | 1 | 2 | 1 361 463 | 101 396 | 0.074 |

50-bit:

| curve | \|cls\| | conductor | n_succ | mean | stdev | CoV |
|-------|-------|-----------|--------|------|-------|-----|
| #0 | 16 \* | 1 | 16 | 37 597 512 | 20 375 839 | **0.542** |

`*` = the walk hit `max_graph_nodes = 16`, so this class is truncated
(§3.2).  22- and 30-bit contribute no row: every class there has at most
one successful odd-`n` vertex.

*(was, on the buggy Vélu: five 14-bit and six 18-bit classes, all at
`|cls| = 16`, CoV 0.357–0.832.)*

**The 50-bit class is the best evidence this section has ever had.**  Its
16 samples give CoV **0.542** against the theoretical 0.52 for a
truncated-at-N-iters Pollard ρ with mean ≈ √(πn/2) — a 4 % gap, at the
largest width measured and on a class an order of magnitude larger in `n`
than the 16-bit ones.  The four 16-bit classes (14–16 samples each) sit at
0.54–0.81.  Everything else in the table rests on two samples, where a CoV
is essentially unconstrained: the 0.074 at 40 bits and the 0.078 and 0.121
at 18 and 14 bits are not evidence of tight clustering.

**Conclusion, strengthened at the top end:** in the two widths with enough
samples to say anything — 16-bit and 50-bit — vertices of the same isogeny
class show no anomalous within-class clustering of rho cost, and at 50 bits
the dispersion matches the generic prediction to 4 %.  At 14, 18, 22, 30
and 40 bits the corrected classes are too small to support the claim either
way.

### 3.4.1 Permutation test for class-structure leakage

To make the null-of-no-leak rigorous, I ran a 2 000-trial permutation
test on the within-class variance of rho-iter counts (odd-n
successes only):

| bits | classes (≥2 samples) | obs mean within-class var | null 5–95 % | p(obs ≤ null) |
|------|----------------------|---------------------------|-------------|---------------|
| 14   | 1                    | n/a — test needs ≥2 classes | n/a       | n/a     |
| 16   | 4                    | 18 241                    | [16 831, 19 320]  | 0.593   |
| 18   | 2                    | 32 216                    | [32 216, 38 708]  | 1.000   |
| 22   | 0                    | n/a — test needs ≥2 classes | n/a       | n/a     |
| 30   | 0                    | n/a — test needs ≥2 classes | n/a       | n/a     |
| 40   | 1                    | n/a — test needs ≥2 classes | n/a       | n/a     |
| 50   | 1                    | n/a — test needs ≥2 classes | n/a       | n/a     |

*(was, on the buggy Vélu: 5/6/9 classes at 14/16/18 bits,
p = 0.526 / 0.985 / 0.025.)*

**Adding the 22–50-bit data does not rescue this test, and the reason is
worth stating rather than hiding in four n/a cells.**  The test needs at
least two classes each holding at least two successful odd-`n` vertices.
At 22 and 30 bits the classes are singletons; at 40 bits only one class
has two; at 50 bits there is one class with plenty of samples (16) but
only one such class, because the second starting curve's class is entirely
even-`n`.  So the corrected data leaves exactly the same one usable row it
had before: 16-bit, four classes, p = 0.593, squarely inside the null.
18-bit stays degenerate — the observed value *is* the null minimum, so
p = 1.000 is an artifact of having two points to permute.

**Permutation conclusion, unchanged:** no detectable isogeny-class
structure leak at 16 bits; at every other width the test is underpowered
on the corrected classes and neither supports nor refutes a leak.  The
earlier "directions disagree across scales" reading stays withdrawn — it
was computed over vertices that were not all in the isogeny class.  What
would change this is more vertices per class, i.e. a larger
`max_graph_nodes`, not more starting curves.

### 3.5 j = 0 / j = 1728 behavior

In the harness data none of the 173 random vertices (14/16/18/22/30/40/50-bit
pooled, regenerated in rounds 9 and 10; was 119 over 14/16/18 alone) landed on
j = 0 or j = 1728 (they're a measure-zero subset of the curve space for our
random generation).  The dedicated `volcano --curve toy-j0` probe at p=103
shows what *is* anomalous about j=0: the 3-volcano collapses to a single
crater because the unique CM order Z[ζ₃] has discriminant −3 in which ℓ=3
ramifies.  This affects *graph structure* (CSIDH-style group-action depth)
rather than rho speed on a given curve.

## 4. Headline findings

0. **Major implementation artifact discovered: the rho module is biased
   against even-order curves.**  The sterile-collision check
   `rhs.gcd(n) != 1` at [pollard_rho.rs:189](src/cryptanalysis/pollard_rho.rs:189)
   triggers ~50% of the time when `n` is even, exhausting the 8 restarts
   long before reaching the cycle.  Observed success rates: 7-18% for
   even n vs 83-97% for odd n.  Fix is one line: when the gcd is non-trivial,
   record the partial log `x mod (n/gcd)` and recurse via CRT, rather than
   restarting from scratch.  (Or use a partition function that's harmonic
   with the group structure.)  Filed as a follow-up task.

1. **Conditional on n=odd, no isogeny-class rho-cost anomaly.**  Median
   Pollard-ρ iteration count matches √(πn/2) to within 15% across 14, 16,
   and 18 bits.  The within-class CoV (≈0.75–1.07) is consistent with the
   cap-truncated Pollard-ρ distribution, not with any class-structure leak.

2. **Volcano structure is dictated entirely by `disc(End(E)) mod ℓ`.**
   `toy-a` at ℓ=2 (disc inert) gives a 1-vertex volcano; the same
   curve at ℓ=3 (disc split) gives a 2-vertex edge.  Identical
   behavior is observed for D=−23 where the analytical class-group
   `cm --discriminant -23` predicts the volcano in advance.

3. **j=0 is structurally distinct.**  Over F_103 the 3-volcano of
   y²=x³+1 collapses because ℓ=3 ramifies in Z[ζ₃].  This is a
   well-known obstruction, not a new phenomenon, but the module
   reproduces it cleanly.

4. **secp256k1 passes every structural check** the module knows how to
   ask: GLV constants β, λ verify; MOV/Frey-Rück infeasible; twist not
   pairing-friendly; the disc=-3 structural survey enumerates exactly
   the split / ramified / inert pattern at small ℓ predicted by
   class-field theory.

5. **A real research-grade sweep at ≥ 30 bits is blocked by the
   point-counting implementation.**  Replacing brute-force counting
   with Schoof–Elkies–Atkin (already documented in
   `src/cryptanalysis/ai_schoof.rs`) would lift the ceiling.

## 5. Caveat about the `--bits N` flag

The CLI accepts a `--bits` argument up to u32::MAX, but
`generate_random_curves` clamps to `bits ∈ [8, 62]`.  Anything past 62
bits is silently downgraded; at 64 bits the user gets a 62-bit prime.  The
upper bound is not an algorithmic limit but an arithmetic one — `p < 2^62`
keeps the intermediate `u128` products in the `SmallCurve` helpers from
overflowing.

*(This section said `[8, 30]` until round 10, as did §1.  Round 3 raised
the clamp when it replaced the `O(p)` Frobenius trace with Shanks–Mestre
BSGS and the `O(p²)` point enumeration with Tonelli–Shanks sampling, and
neither section was updated.  The consequence was not cosmetic: it made
the 40- and 50-bit rows of §1 and §3.2 read as secretly-30-bit runs.  They
are genuine — §1's `actual bits` column is read from each artifact's own
prime.)*

What remains true is that cost still grows fast: 50 bits takes 45 minutes
against 12 seconds at 40.  Pushing past ~50 bits in practice wants the
algorithmic upgrades the module docstring names, SEA and NUCOMP/NUDUPL,
rather than a higher clamp.

## 6. Raw data

Every command's stdout is preserved under `experiments/`:

* [01_secp256k1.txt](experiments/01_secp256k1.txt)
* [02_volcano_toy-a_2.txt](experiments/02_volcano_toy-a_2.txt)
* [03_volcano_toy-j0_3.txt](experiments/03_volcano_toy-j0_3.txt)
* [04_graph_toy-a_2-3.txt](experiments/04_graph_toy-a_2-3.txt)
* [05_cm_disc-23.txt](experiments/05_cm_disc-23.txt)
* [06_experiment_14bit_10trials.json](experiments/06_experiment_14bit_10trials.json) (cap 2^14)
* [06b_experiment_14bit_10trials_HIGH_CAP.json](experiments/06b_experiment_14bit_10trials_HIGH_CAP.json) (cap 2^18)
* [06_experiment_16bit_10trials.json](experiments/06_experiment_16bit_10trials.json) (cap 2^14)
* [07_experiment_18bit_10trials.json](experiments/07_experiment_18bit_10trials.json) (cap 2^14)
* [13_22bit_dpoly.json](experiments/13_22bit_dpoly.json) — 22-bit, 4 trials, cap 2^18
* [14b_30bit.json](experiments/14b_30bit.json) — 30-bit, 4 trials, cap 2^18
* [15_40bit.json](experiments/15_40bit.json) — 40-bit, 4 trials, cap 2^18
* [15b_40bit.json](experiments/15b_40bit.json) — 40-bit, 4 trials, cap 2^23 (the row §3.2 quotes)
* [16_50bit.json](experiments/16_50bit.json) — 50-bit, 2 trials, cap 2^18 (no successes at that cap)
* [16c_50bit.json](experiments/16c_50bit.json) — 50-bit, 2 trials, cap 2^28 (the row §3.2 quotes)
* [pooled_vertex_data.csv](experiments/pooled_vertex_data.csv) — every vertex across **all seven** bit widths flattened to one CSV: `bits, rho_cap, curve_idx, class_size, class_capped, p, a, b, trace, n, n_parity, fundamental_disc, conductor, rho_iters, rho_success, mov_feasible, smart_applies, glv_speedup`.  173 rows (was 119 over 14/16/18 in round 9, and 275 on the buggy Vélu); ready to load in Pandas / R.
  `rho_cap` and `class_capped` are new in round 10: the widths use
  different rho caps, so a success rate is only comparable against a row at
  the same cap, and `class_capped` marks the classes the 16-node walk
  budget truncated (every 50-bit row, and the four saturated 16-bit
  classes).  Pass at most one artifact per bit width — 40 and 50 each have
  two caps over the same vertices, and pooling both double-counts.
  Regenerate with
  [`scripts/pool_isogeny_vertices.py`](scripts/pool_isogeny_vertices.py);
  **every** table in §1, §3.2, §3.2.1, §3.3, §3.4 and §3.4.1 is reproduced by
  [`scripts/analyse_isogeny_vertices.py`](scripts/analyse_isogeny_vertices.py).
  Both were written in round 9 — the CSV and those statistics had been
  produced ad hoc, which is why nothing caught the bad vertices.  Round 9's
  version covered only §3.2, §3.4 and §3.4.1, which is why §3.2.1, §3.3 and
  §1's parity table stayed on pre-fix numbers for a round; it covers all six
  now.

**Removed in round 10.**  `09_22bit`, `12_26bit`,
`06_experiment_20bit_10trials` and `06_experiment_64bit_10trials` were
committed as **0-byte** files — no content, referenced by nothing.
`11_22bit` and `14_30bit` were dropped as duplicates: their `config` blocks
are identical to `13_22bit_dpoly` and `14b_30bit` respectively, so what
separated them was a code revision, not an experiment, and today's tree
cannot reproduce the distinction.  None of the six was cited anywhere in
this document, which is exactly how the 22–50-bit rows came to be described
as having no frozen artifact.
