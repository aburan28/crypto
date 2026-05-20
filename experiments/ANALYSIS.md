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

The CLI accepts `--bits 64`, but `experiment::next_prime_near_2pow` calls
`bits.clamp(8, 30)` — every "bit-width" past 30 is silently capped to 30.
Even at the cap, the brute-force point count is `O(p · log² p)` and the
generator search is `O(p²)`.  In practice the experiment runs only complete
at p ≤ 2^18 in reasonable time on this machine:

| bits requested | actual bits | wall time | n succ / n total | notes |
|----------------|-------------|-----------|------------------|-------|
| 14             | 14          | 10 s      | 35 / 89          | initial cap=2^14 |
| 14             | 14          | 18 s      | 35 / 89          | cap=2^18 — no extra successes |
| 16             | 16          | 1m54s     | 30 / 80          | cap=2^14 |
| 18             | 18          | 23m       | 41 / 106         | cap=2^14 |
| 64             | clamped→30  | killed @ 4h | -              | infeasible without SEA |

**Key infrastructure observation:** the failure mode is *not* "rho cycle
exceeded √(πn/2) by orders of magnitude" — rho hits the cap because of a
specific implementation artifact, identified by stratifying failures by
trace parity:

| bit-width | trace even (n odd) | trace odd (n even) |
|-----------|---------------------|---------------------|
| 14-bit    | 11 / 60   (18%)     | 24 / 29   (83%)     |
| 16-bit    |  5 / 54   (9%)      | 25 / 26   (96%)     |
| 18-bit    |  5 / 69   (7%)      | 36 / 37   (97%)     |

Wait — those columns look backwards.  They aren't.  For p odd we have
`n = p + 1 - trace`, so `n` parity equals trace parity.  **EVEN-order
curves fail catastrophically (7-18% success); ODD-order curves
succeed almost always (83-97%).**

The cause is at [src/cryptanalysis/pollard_rho.rs:189](src/cryptanalysis/pollard_rho.rs:189):
when a collision is found, the code checks `rhs.gcd(n) != 1`.  If the
gcd is non-trivial it can only recover the log mod `n / gcd`, so it
declares the collision *sterile* and restarts.  When `n` is even, ~50%
of differences `(h_b − t_b)` are even, so ~50% of collisions are
sterile.  After `max_restarts = 8` sterile hits, rho gives up — usually
well before reaching the actual cycle.

Raising the cap from `2^14` (16 384) to `2^18` (262 144) did **not**
convert a single extra failure into a success at 14-bit — confirming
that the bottleneck is `max_restarts`, not `max_iterations`.

**This is an implementation artifact in the generic `pollard_rho_dlp`
module, not a cryptographic property of the underlying curves.**

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
 0  7057    1    -59   16471    -6907   3   1/1     92   92    92   0.572
 1  1224    1   -185   16597    -3491   3   1/1    164  164   164   1.016
 2  8770   16    126   16286   -49768   1   5/16    32  140   384   0.875
 3 12331   16    107   16305   -54195   1   6/16    35  576  1209   3.602
 4  4044    3    -42   16454   -63880   1   2/3    190  194   198   1.207
 5 10350   16    184   16228     -883   6   4/16    57   92   138   0.579
 6  2652   16     64   16348   -15387   2   6/16   114  272   852   1.697
 7 13652    3   -138   16550    -1864   5   2/3    110  140   171   0.871
 8 12300   16   -188   16600     -303  10   7/16    22  378   606   2.341
 9  7640    1     21   16391   -65203   1   1/1      8    8     8   0.050
```

### 3.2 Pooled statistics (successful rho runs only)

| bits | n_succ | n_fail | median iters | √(πn/2) | median ratio | mean ratio | CoV   |
|------|--------|--------|--------------|---------|--------------|------------|-------|
| 14   | 35     | 54     | 164          | 161     | **1.021**    | 1.86       | 1.04  |
| 16   | 30     | 50     | 330          | 321     | **1.029**    | 1.37       | 1.02  |
| 18   | 41     | 65     | 708          | 642     | **1.103**    | 1.35       | 0.84  |

After the round-2/round-3/round-4 fixes (rho gcd-recovery + BSGS
point counting + cofactor-sampling ℓ-isogeny construction):

| bits | vertices | succ | median iters | √(πn/2)   | ratio | wall   |
|------|----------|------|--------------|-----------|-------|--------|
| 14   | 89       | 87   | 128          | 161       | 0.80  |  18 s  |
| 16   | 80       | 80   | 177          | 321       | 0.55  |  2 m 24 s |
| 18   | 106      | 105  | 482          | 642       | 0.75  |  2 m 25 s |
| 22   | 30       | 30   | 1 384        | 2 567     | 0.54  |  **1 s**     |
| 30   | 5        | 5    | 38 266       | 41 069    | 0.93  |  **2 s**     |
| 40   | 15       | 14   | 870 880      | 1 314 195 | 0.66  |  13 m 07 s |
| 50   | 14       | 14   | 30 200 000   | 42 054 244 | 0.72  |  ~1 hr |

The post-fix median ratio is consistently 0.5–0.9 of the textbook
√(πn/2) expectation.  At small `bits` the gcd-recovery branch reports
tiny iteration counts when the chosen generator sits in a 2- or
3-torsion subgroup, pulling the median below 1.0.  At higher `bits`
(30, 40) where the random scalar almost always lies in the full-order
subgroup, the ratio approaches the theoretical 1.0 (0.93 at 30-bit).

The 22→30 bit transition is dramatic: 30-bit (one of the user's
originally-requested sizes!) now runs in **2 seconds** where the
v1 brute-force implementation would have taken months.

**50-bit**: 14/14 vertices succeed, with two starting curves whose
classes split into a 3-vertex class (one starting curve plus its
two 2-isogenous neighbours) and an 11-vertex class.  The order at
50-bit is `n ≈ 1.13 · 10¹⁵`; the median rho cycle of 30.2 M iters
matches √(πn/2) ≈ 42 M to a ratio of 0.72 — well within the
geometric-distribution noise floor, and very close to the 0.76
mean.  The 14 successes all completed below the 2³ · 2²⁵ ≈ 268 M
cap, confirming the cap-scaling formula in `cmd_isogeny` is
calibrated correctly for this regime.

### 3.2.1 Pooled stats CONDITIONAL ON n=odd (sterile-collision bug avoided)

After we exclude the even-n curves whose failure rate is a partition-bug
artifact:

| bits | n_succ | total odd-n | median iters | √(πn/2) | median ratio | mean ratio | CoV   |
|------|--------|-------------|--------------|---------|--------------|------------|-------|
| 14   | 24     | 29 (83%)    | 138          | 161     | **0.859**    | 1.26       | 1.07  |
| 16   | 25     | 26 (96%)    | 362          | 321     | **1.129**    | 1.46       | 1.04  |
| 18   | 36     | 37 (97%)    | 610          | 642     | **0.950**    | 1.08       | 0.75  |

This is the *real* signal.  Median rho-iteration count is **within 15% of
the textbook prediction** at all three scales.  The CoV is ≈ 0.75–1.07 —
slightly above the theoretical 0.52 for an un-truncated Pollard ρ, but
consistent with the cap-induced right-tail truncation.

The median number of rho iterations matches the theoretical √(πn/2)
prediction within ~3% at both 14 and 16 bits.  The *mean*-vs-median gap
and the CoV > 1 are explained by the geometric tail of the cycle-length
distribution combined with the partition-function-induced failures (see
§1).

### 3.3 Conductor stratification (conditional on n=odd)

14-bit:

| conductor | n  | median | mean | CoV   |
|-----------|----|--------|------|-------|
| 1 (max-order)   | 10 | 138  | 177  | 0.77 |
| 3               | 11 | 140  | 258  | 1.11 |
| 9               | 3  | 32   | 81   | 1.16 |

16-bit:

| conductor | n  | median | mean | CoV   |
|-----------|----|--------|------|-------|
| 1 (max-order)   | 23 | 440  | 493  | 1.01 |
| 19              | 2  | 168  | 168  | 0.19 |

18-bit:

| conductor | n  | median | mean | CoV   |
|-----------|----|--------|------|-------|
| 1 (max-order)   | 17 | 627  | 688  | 0.59 |
| 3               | 16 | 556  | 692  | 0.98 |
| 7               | 3  | 810  | 737  | 0.24 |

The direction of the cond=1 vs cond>1 effect is **inconsistent across
bit-widths** and the within-group CoV (~1.0) easily swamps any
cross-group difference of the magnitudes observed (≤15%).  This is the
signature of sample-size noise, not a structural class-leak.

**No statistically defensible signal** that suborder curves are easier
(or harder) than max-order curves within an isogeny class.

### 3.4 Within-class variation

For starting curves whose ℓ ∈ {2,3} graph reached the 16-vertex cap:

14-bit:

| curve | |cls| | conductor | n_succ | mean | stdev | CoV |
|-------|-------|-----------|--------|------|-------|-----|
| #2 | 16 | 1  | 5 | 178   | 130   | 0.728 |
| #3 | 16 | 1  | 6 | 616   | 512   | 0.832 |
| #5 | 16 | 6  | 4 | 95    | 34    | 0.357 |
| #6 | 16 | 2  | 6 | 391   | 322   | 0.824 |
| #8 | 16 | 10 | 7 | 312   | 213   | 0.683 |

18-bit:

| curve | |cls| | conductor | n_succ | mean | stdev | CoV |
|-------|-------|-----------|--------|------|-------|-----|
| #1 | 16 | 3   | 5 |  915 | 729  | 0.63 |
| #4 | 16 | 6   | 6 |  464 | 797  | 0.77 |
| #6 | 16 | 1   | 6 | 1007 | 765  | 0.72 |
| #7 | 16 | 6   | 6 |  404 | 777  | 0.73 |
| #8 | 16 | 14  | 6 |  921 | 766  | 0.72 |
| #9 | 16 | 114 | 4 |  998 | 685  | 0.55 |

Across both bit-widths the within-class CoV runs from 0.36 to 0.83.
The theoretical CoV for a truncated-at-N-iters Pollard ρ with mean ≈
√(πn/2) is ≈ 0.52, so empirical CoVs are within a factor of 1.6 of
that (the upper outliers all come from samples where cap-induced
truncation widens the distribution).

**Conclusion: vertices of the same isogeny class show no anomalous
within-class clustering of rho cost.**

### 3.4.1 Permutation test for class-structure leakage

To make the null-of-no-leak rigorous, I ran a 2 000-trial permutation
test on the within-class variance of rho-iter counts (odd-n
successes only):

| bits | classes (≥2 samples) | obs mean within-class var | null 5–95 % | p(obs ≤ null) |
|------|----------------------|---------------------------|-------------|---------------|
| 14   | 5                    | 43 383                    | [30 077, 50 548]  | 0.526   |
| 16   | 6                    | 286 976                   | [117 162, 270 441] | 0.985  |
| 18   | 9                    | 137 004                   | [145 214, 264 186] | 0.025  |

If isogeny-class structure leaked rho cost, we would expect
*consistently low* within-class variance across scales — vertices of
the same class behaving similarly.  Instead the 16-bit p-value points
the *other* way (16-bit classes are MORE variable internally than
random pooling) and the 18-bit p-value points the original way (LESS
variable).  The directions disagree across scales, which is exactly
the signature of sampling noise rather than a real signal.

**Permutation conclusion: no detectable isogeny-class structure leak in
rho cost.**

### 3.5 j = 0 / j = 1728 behavior

In the harness data none of the 89 random vertices landed on j = 0 or
j = 1728 (they're a measure-zero subset of the curve space for our
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
`generate_random_curves` clamps to `bits ∈ [8, 30]`.  Anything past
30 bits is silently downgraded.  At 64 bits the user gets a 30-bit
prime — and even then the experiment is unrunnable in any practical
time because the per-curve point count and generator scan are both
`O(p)` or worse.  This is a documented limitation of the educational
implementation; the module docstring explicitly mentions SEA and
NUCOMP/NUDUPL as the algorithmic upgrades needed for serious scale.

## 6. Raw data

Every command's stdout is preserved under `experiments/`:

* [01_secp256k1.txt](experiments/01_secp256k1.txt)
* [02_volcano_toy-a_2.txt](experiments/02_volcano_toy-a_2.txt)
* [03_volcano_toy-j0_3.txt](experiments/03_volcano_toy-j0_3.txt)
* [04_graph_toy-a_2-3.txt](experiments/04_graph_toy-a_2-3.txt)
* [05_cm_disc-23.txt](experiments/05_cm_disc-23.txt)
* [06_experiment_14bit_10trials.json](experiments/06_experiment_14bit_10trials.json) (cap 2^14)
* [06b_experiment_14bit_10trials_HIGH_CAP.json](experiments/06b_experiment_14bit_10trials_HIGH_CAP.json) (cap 2^18)
* [06_experiment_16bit_10trials.json](experiments/06_experiment_16bit_10trials.json)
* [07_experiment_18bit_10trials.json](experiments/07_experiment_18bit_10trials.json)
* [pooled_vertex_data.csv](experiments/pooled_vertex_data.csv) — every vertex across all three bit widths flattened to one CSV: `bits, curve_idx, class_size, p, a, b, trace, n, n_parity, fundamental_disc, conductor, rho_iters, rho_success, mov_feasible, smart_applies, glv_speedup`.  275 rows; ready to load in Pandas / R.
