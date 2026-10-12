# Lead 16.5: class-wide multi-target rho, measured (round 4, 2026-10-11)

Source: `RESEARCH_ISOGENY_FIBERS_ECDLP.md` §14.5 (prediction), §16.5 (this
lead), §13.6 (one table serves the class). Code:
`examples/r4_16_5_class_multitarget_rho.rs` (native Rust, uses `crypto_lib`).
Toy instances only. No sub-rho claim, and nothing here is about a deployed
curve.

## Question, hypothesis, falsifier

- **Question.** Take L keys on L different curves of one `F_p`-isogeny class,
  where the curves are linked by explicit small-degree isogenies. Transport
  every key to one common curve. Does solving them cost what L-target rho
  costs on a single curve?
- **Hypothesis (§14.5).** After transport, total cost is about √(2Lr) group
  operations (Kuhn–Struik). Solving the L keys independently costs about
  L·√(πr/2). Transport adds a small, additive cost.
- **Falsifier (stated before measuring).** The transport argument is wrong
  somewhere if arm B (transport, then multi-target rho on E₀) is not within
  sampling error of arm C (L targets planted directly on E₀, same
  multi-target rho). Comparison: the ratio of mean totals B/C at each
  (r, L), with a 95% interval.
- **Boundary and unit.** The unit is S = total group operations / √r in the
  order-r subgroup. One group operation is one affine addition or doubling.
  Two references are derived, not measured:
  - independent rho: L·√(π/2);
  - Kuhn–Struik: Σ_{i<L} C(2i,i)·4^{−i}·√(π/2), which tends to √(2L).

## What ran

1. **Class construction (native).**
   - The base curve E₀: y² = x³ + ax + b over F_p with #E₀ = h·r, r prime,
     6 | h ≤ 48. The order is counted by a native Mestre-style BSGS and
     confirmed on extra points.
   - A breadth-first walk over rational 2-isogenies and 3-isogenies.
     - 2-isogeny kernels come from the rational roots of x³+ax+b.
     - 3-isogeny kernels come from the rational roots of ψ₃ = 3x⁴+6ax²+12bx−a².
     - Roots come from `crypto_lib`'s `isogeny_walk::poly::roots`.
     - Each isogeny is Vélu's, written in x only, with Y = y·dX/dx.
   - The walk stops once it has 65 distinct j-invariants.
   - Every edge is checked four ways:
     - images lie on the codomain;
     - φ(T+P) = φ(T)+φ(P) for a random point T;
     - every 3-isogeny codomain equals `crypto_lib`'s
       `isogeny_walk::kernel::velu_codomain`, and so does its j-invariant;
     - an explicit dual exists: a Vélu kernel on the codomain followed by an
       F_p-isomorphism (x,y) ↦ (s x, ±s√s y), and it is certified by
       dual∘φ = [ℓ] on the generator.
2. **Generators.** P_i = ψ_i(P₀) along the forward path. The dual path then
   sends P_i to [d_i]P₀, where d_i = deg ψ_i, so c_i = d_i mod r. For every
   node the run asserts transport(P_i) = [d_i]P₀.
3. **Rho.**
   - All arms use one distinguished-point rho:
     - a target-independent 32-adding walk X ← X + α_j P₀;
     - the partition and the DP test come from a salted hash of x;
     - DP probability 2^−dpbits;
     - walks start at a·P + b·Q with random b;
     - restarts cost one addition: start ← start + (σP + τQ);
     - a walk is abandoned after 20·2^dpbits steps.
   - Collisions are resolved for both y-equal and y-opposite points.
   - Targets are solved in sequence with the DP table kept across them
     (Kuhn–Struik).
   - Every recovered log is checked against its target.
4. **Arms.** Each trial plants fresh random keys k_i ∈ [1, r) on E₁..E_L.
   E₀ itself is excluded, so every arm-B target is really transported.
   - **A.** L independent single-target runs: target i on E_i with base P_i,
     a fresh DP table and a fresh step table for each run.
   - **B.** Each Q_i is mapped to E₀ by its dual path, then one multi-target
     run solves the transported targets. Each log becomes k_i = m_i·d_i^{−1}
     mod r, and the run checks [k_i]P_i = Q_i on E_i and k_i = planted key.
   - **C.** L fresh keys planted on E₀, then the same multi-target run.
5. **Accounting.**
   - The totals include every walk step plus setup:
     - 32 scalar multiples for the step table;
     - per target, two scalar multiples for the start and two for the
       restart shift;
     - one addition per restart.
   - Verification scalar multiples are recorded separately and are not
     included.
   - Transport is counted in field multiplications (M) and inversions (I).
     It is added to the B total at 1 group-op equivalent per inversion plus
     (M − 3I)/3. In this code one affine addition is 1I + 3M.

Commands (shared target dir, release, single-threaded):

```sh
export CARGO_TARGET_DIR=/tmp/claude-0/-home-user/717e4727-b5a7-5b25-924a-7b31edb43331/scratchpad/r4target
B=$CARGO_TARGET_DIR/release/examples/r4_16_5_class_multitarget_rho
$B --rbits 28 --trials 100 --dpbits 8  --seed 28 --out run_r28.json   # stderr -> run_r28.txt
$B --rbits 32 --trials 40  --dpbits 8  --seed 32 --out run_r32.json   # stderr -> run_r32.txt
$B --rbits 36 --trials 20  --dpbits 10 --seed 36 --out run_r36.json   # stderr -> run_r36.txt
$B --fit run_r28.json,run_r32.json,run_r36.json --out fit.json        # stderr -> fit.txt
venv/bin/python plot.py .                                              # SVGs from the JSON only
```

`--rbits` sets p ∈ [2^(rbits+3), 2^(rbits+4)). The r that came out is in the
table: 2^28.79, 2^31.40 and 2^37.18. The last is a little above the brief's
2^36 because the search took h = 6.

Replay check: `--rbits 28 --seed 28` was rerun with the final binary. Its
`summary`, `trials` and `class` came out identical to `run_r28.json`.

Base curves (toy). `curve_id::prime` computed each one's ICV1 slug, which is
recorded in the run JSON as `base_curve.icv1_slug`. The slugs are not yet in
`docs/curves/registry.json`, so this file leaves them out. The checker flags
unregistered slugs even inside code fences. Register them with
`scripts/build_curve_registry.py` before naming them in prose.

```
r=2^28.79: p=2790788831   a=666561205    b=2453577277   #E=6·465130349
r=2^31.40: p=68190058247  a=65137035636  b=33544543934  #E=24·2841254393
r=2^37.18: p=931154352071 a=194304939514 b=873553658046 #E=6·155192406833
```

Class shapes (measured):

| r | 3-isogeny edges | 2-isogeny edges | max depth |
|:--|--:|--:|--:|
| 2^28.79 | 63 | 1 | 17 |
| 2^31.40 | 49 | 15 | 5 |
| 2^37.18 | 63 | 1 | 17 |

At r = 2^28.79 and r = 2^37.18 the tree is mostly a 3-isogeny cycle walked
in two directions. All 64 curves in each class (65 with E₀) are listed in
the JSON under `class.curves`.

## Measured results

S = total group operations / √r, shown as mean ± standard error over trials.
The reference columns are derived. The transport column is the mean per
target. "Verified" means every recovered log equals the planted key and was
checked on its original curve.

| r | L | trials | A (indep.) | B (transport + multi) | C (control, E₀) | KS ref | L√(π/2) | B/C [95%] | A/B | transport / target | verified |
|:--|--:|--:|--:|--:|--:|--:|--:|:--|--:|:--|:--|
| 2^28.79 | 1 | 100 | 1.312 ± 0.068 | 1.414 ± 0.066 | 1.393 ± 0.070 | 1.253 | 1.253 | 1.015 [0.879, 1.152] | 0.93 | 9.0 M + 1.0 I | yes |
| 2^28.79 | 4 | 100 | 5.442 ± 0.119 | 2.880 ± 0.067 | 2.907 ± 0.074 | 2.742 | 5.013 | 0.991 [0.924, 1.057] | 1.89 | 11.2 M + 1.2 I | yes |
| 2^28.79 | 16 | 100 | 21.820 ± 0.224 | 6.056 ± 0.072 | 6.110 ± 0.073 | 5.613 | 20.053 | 0.991 [0.958, 1.024] | 3.60 | 24.8 M + 2.8 I | yes |
| 2^28.79 | 64 | 100 | 87.543 ± 0.488 | 12.909 ± 0.077 | 12.752 ± 0.065 | 11.292 | 80.212 | 1.012 [0.997, 1.028] | 6.78 | 78.8 M + 8.8 I | yes |
| 2^31.40 | 1 | 40 | 1.292 ± 0.106 | 1.347 ± 0.108 | 1.148 ± 0.105 | 1.253 | 1.253 | 1.173 [0.894, 1.452] | 0.96 | 9.0 M + 1.0 I | yes |
| 2^31.40 | 4 | 40 | 5.180 ± 0.176 | 2.897 ± 0.096 | 2.958 ± 0.112 | 2.742 | 5.013 | 0.979 [0.883, 1.076] | 1.79 | 11.2 M + 1.2 I | yes |
| 2^31.40 | 16 | 40 | 21.506 ± 0.474 | 5.670 ± 0.093 | 5.828 ± 0.100 | 5.613 | 20.053 | 0.973 [0.928, 1.018] | 3.79 | 20.2 M + 2.2 I | yes |
| 2^31.40 | 64 | 40 | 82.869 ± 0.777 | 12.062 ± 0.122 | 12.084 ± 0.108 | 11.292 | 80.212 | 0.998 [0.972, 1.025] | 6.87 | 34.6 M + 3.8 I | yes |
| 2^37.18 | 1 | 20 | 1.334 ± 0.123 | 1.170 ± 0.140 | 1.057 ± 0.117 | 1.253 | 1.253 | 1.107 [0.753, 1.461] | 1.14 | 9.0 M + 1.0 I | yes |
| 2^37.18 | 4 | 20 | 5.046 ± 0.360 | 2.859 ± 0.117 | 2.898 ± 0.150 | 2.742 | 5.013 | 0.986 [0.859, 1.114] | 1.77 | 11.2 M + 1.2 I | yes |
| 2^37.18 | 16 | 20 | 19.635 ± 0.518 | 5.811 ± 0.178 | 5.934 ± 0.177 | 5.613 | 20.053 | 0.979 [0.897, 1.061] | 3.38 | 24.8 M + 2.8 I | yes |
| 2^37.18 | 64 | 20 | 81.143 ± 1.109 | 11.648 ± 0.209 | 11.354 ± 0.168 | 11.292 | 80.212 | 1.026 [0.979, 1.073] | 6.97 | 78.8 M + 8.8 I | yes |

Counts across all recorded trials of B and C:

| r | abandoned walks | useless collisions |
|:--|--:|--:|
| 2^28.79 | 40 | 0 |
| 2^31.40 | 8 | 0 |
| 2^37.18 | 0 | 0 |

All 40,800 recovered logs were verified: (100 + 40 + 20) trials × (1 + 4 +
16 + 64) targets × 3 arms. The B logs were checked on their original curve
E_i.

**Fits (`fit.json`).** Least squares of ln S_mean on ln L over L ∈ {1, 4, 16, 64}.
The 95% intervals come from 2000 bootstrap resamples of the trials within
each L.

| r | β (A) | β (B) | β (C) | c_A (S ≈ c·L) | S_B/√L at L=64 | S_C/√L at L=64 |
|:--|:--|:--|:--|--:|--:|--:|
| 2^28.79 | 1.009 [0.987, 1.032] | 0.532 [0.513, 0.554] | 0.533 [0.512, 0.555] | 1.325 | 1.614 | 1.594 |
| 2^31.40 | 1.003 [0.970, 1.043] | 0.523 [0.489, 0.559] | 0.558 [0.522, 0.600] | 1.298 | 1.508 | 1.510 |
| 2^37.18 | 0.987 [0.949, 1.035] | 0.548 [0.504, 0.604] | 0.565 [0.521, 0.615] | 1.307 | 1.456 | 1.419 |

For comparison, the exponent derived from the Kuhn–Struik reference over the
same four L values is ln(11.292/1.253)/ln 64 = 0.53, an endpoint slope. Its
S/√L at L = 64 is 1.41.

## Verdict: confirmed, on toy classes

- **The falsifier did not trigger.** In all 12 (r, L) cells the 95%
  interval for B/C contains 1. At L = 64 the ratio is 1.012 [0.997, 1.028],
  0.998 [0.972, 1.025] and 1.026 [0.979, 1.073]. Arm B tracks the
  single-curve multi-target control target by target as well; see
  `fig_per_target.svg`.
- **Scaling.** Total cost grows as L^0.52–0.55 for B and C, against the
  derived 0.53. Arm A grows as L^0.99–1.01.
- **Saving.** At L = 64 the class-wide batch is 6.8–7.0× cheaper than
  solving each curve separately. The derived ratio is
  80.2/11.29 = 7.1 (measured A/B: 6.78, 6.87, 6.97).
- **Transport cost.** It is additive and negligible:
  - 9 M + 1 I per isogeny step; the depth is at most 17;
  - at most 78.8 M + 8.8 I per target, against about 0.2√r ≈ 4·10³ to
    8·10⁴ rho steps per target at L = 64.
- **Gap to the ideal.** B and C sit 3–14% above the Kuhn–Struik number, by
  the same amount in both arms. The cause is accounted setup and DP tails:
  the 32-entry step table, the four scalar multiples per target, and about
  2^dpbits steps to reach a DP after each collision. This overhead is
  largest at the smallest √r. It is engineering, and it does not separate
  B from C.

Read as protocol guidance: issuing per-user curves of one order inside one
isogeny class, connected by small-degree rational isogenies, buys no
multi-target protection over a single shared curve. This is the
multi-target face of §13.6.

## Limits and what was not done

- **Toy sizes only.**
  - r runs from 2^28.8 to 2^37.2, with L ≤ 64.
  - Each r uses one class, built from one seed.
  - Trial counts are 100, 40 and 20, so the L = 1 cells are noisy (±10–13%).
- **Rational ℓ-torsion only.** The classes come from rational 2- and
  3-torsion kernels, so every edge has degree 2 or 3. Elkies-type
  isogenies, larger ℓ and vertical steps with large conductor (Prop. 3.3)
  were not exercised. The cost of *finding* a path between curves was not
  measured: the paths are given by the construction, and §13.6 prices that
  phase at Õ(h(D)) once.
- **Not run through ecbench.** ecbench pools no multi-target runs (its spec
  says targets are "never pooled into a multi-target run"). The rho here is
  the example's own DP rho:
  - it is validated by checking every recovered log against the planted key;
  - its single-target S of 1.29–1.33 (arm A, per target) sits near the
    √(π/2) = 1.25 reference plus the accounted setup.
  The rho has no negation map, so a negation-map variant would shift A, B
  and C by the same factor of about √2. No ecbench session or replay
  certificate exists for these numbers.
- **Accounting.** The transport is converted with the conversion stated in
  §What ran, item 5. Because transport is under 0.1% of any B total, the
  conversion factor does not matter.
- **Wall time.** It was not measured as a metric. The contended wall times in
  the JSON (104 s, 111 s and 447 s) are practicality notes on a shared
  4-core host.

## Next step

- Repeat with an Elkies-isogeny class: walk Φ_ℓ through `isogeny_walk`, use
  ℓ up to about 20, and take the path cost from the walk. The aim is to price
  √élu-sized transports.
- Move the rho into ecbench as a pooled multi-target method so the numbers
  carry replay certificates.

## Files

- `run_r28.json`, `run_r32.json`, `run_r36.json`: the raw records, holding
  the class, every trial's per-arm counts, per-target steps and
  verification flags.
- `run_r*.txt`: the example's stderr.
- `fit.json` and `fit.txt`: the fits.
- `plot.py`: plots recorded values only.
- `fig_total_vs_L.svg`, `fig_per_target.svg` and `fig_class_tree.svg`: the
  graphs.
- `report.typ` and `report.pdf`: the report.

**What should move into `src/`.** Two pieces, if the class-wide question
comes up again:

- A u64 prime-field Vélu evaluator for ℓ ∈ {2, 3} with an explicit,
  certified dual. `isogeny_walk::kernel` currently gives codomains and odd
  kernels, but not 2-isogenies, point evaluation or duals.
- A shared-table multi-target DP rho in ecbench.
