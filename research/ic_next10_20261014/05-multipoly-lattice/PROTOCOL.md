# Item 05 protocol: multi-polynomial and shared-target Coppersmith reach (frozen before data)

Frozen 2026-10-14, before any optimiser or LLL output existed. Item "### 5." of
`reports/Index calculus leads next steps.md`.

## Prior-work check (done before freezing)

`grep -ril "multi-polynomial|shared-target|unravelled"` over `research/` and `benchmarks/` finds only
W10 (`research/lit2026_w10_s5_coppersmith_20261011/`), which lists these routes as *not* covered (§7).
No other run, example or journey entry analyses a multi-polynomial Coppersmith lattice for Semaev
polynomials. `conductor check` could not reach a control plane (connection refused); continued per RULES.md.

## Hypothesis

H1: some lattice family that W10 excluded (several polynomials over the same p; shared unknowns across
targets; sublattice / unravelled-linearisation strategies) has an asymptotic small-root reach δ* for the
`S₅` (m = 4) relation search that is at least 1.5x the single-polynomial δ*, and reaches the δ needed for
relation collection to beat rho.

H0: every such family stays below the W10 single-polynomial ceiling 1/32 for `S₅`.

## Arms

| arm | system (small unknowns |x| < p^δ; large unknowns are ~p) | role |
|:--|:--|:--|
| C1 | single `S₃`, `S₄`, `S₅` (W10 family) | **positive control**: must reproduce W10's exact rationals (S₃ ℓ=1..8: 9/100, 13/112, 29/225, 3/22, 90/637, 335/2312, 959/6498, 22/147; S₄ ℓ=1,2: 1331/70875, 79/3087) and the box limits 1/6, 1/24, 1/80 |
| C2 | Boneh–Durfee `1 + x(A+y)` | known-answer control for the ceiling argument (0.2929 must sit under the computed ceiling 1/2) |
| N1 | two `S₃` on disjoint unknowns `(x₁,x₂)`, `(x₃,x₄)` | null control: no shared unknown, so any "gain" over 1/6 asymptotically would be an artefact of the optimiser |
| V1 | `S_{m+1}` + curve equations `y_i² = x_i³ + a x_i + b` (y large) | multi-polynomial |
| V2 | Semaev chain: `S₃(x₁,x₂,t)`, `S₃(x₃,t,z)` (= S₄) and the S₅ chain with intermediate `t` large | multi-polynomial, shared intermediate unknowns |
| V3 | shared target: `S_{m+1}(x_S, x_A, z₁)`, `S_{m+1}(x_S, x_B, z₂)`, s = |S| shared summands, m = 2,3,4 | shared unknowns across targets |
| V4 | sublattices of any of the above (Boneh–Durfee style) and unravelled linearisation | strategy variants |

## Method

1. Generalise W10's family to `r` polynomials: shifts `p^{ℓ−|k|} x^γ Π f_j^{k_j}`, triangular under a
   common monomial order, full rank on a monomial set `M`. Success at `p → ∞`: `Σ_{β∈M} |k(β)| > Σ_{β∈M} λ·β`.
2. Exact finite-ℓ optimum over down-closed `M` (top-corner leading monomials) by min-cut + Dinkelbach,
   exact rationals, on symmetry orbits (Rust, `examples/ic_next10_05_multipoly_lattice.rs`).
3. Asymptotic δ* over per-group box extents (continuous LP limit, numeric quadrature, 1e-9).
4. A pointwise ceiling valid for every `M`, every ℓ and every sublattice, for top-corner leading monomials.
5. Exact MILP (HiGHS) over **all** monomial sets in a bounding box and **all** leading-vertex choices, at
   small ℓ, for the arms whose leading monomial is not the top corner (V1, V2).
6. Relation economics: a system with `n` small unknowns and `r` equations, solved by a poly-time oracle,
   costs `p^{max(r − (n−1)δ, 2δ)}` group operations up to polynomial factors; it beats rho iff
   `(r − 1/2)/(n − 1) < δ < 1/4`. For single `S₅` this is `1/6 < δ < 1/4` (best `δ = 1/5`).

## Sizes

Finite ℓ: S₃ systems ℓ ≤ 8 (≥ 6), S₄ systems ℓ ≤ 2, S₅ systems via asymptotic LP only (7 unknowns).
LLL validation (only for arms with any reach gain): S₃ shared-target vs single S₃, 40-bit p (56-bit if
the dimension allows), ℓ = 1 and 2 (3 if dimension ≤ ~250), 8 trials per (arm, δ), δ grid 0.08–0.24.
Success = at least n reduced rows vanish over Z at the planted root.

## Metric

δ* per arm (exact rational at finite ℓ, numeric asymptotic); ratio to the single-polynomial δ*; LLL
success fraction per δ.

## Kill criterion (from the report)

An LP δ* for `S₅` below 1/32 in every new family closes the lattice route to m = 4. Reaching ≥ 1/10
justifies a measured probe. LLL validation runs when a family beats its single-polynomial δ* (the item
asks for any gain; the report's trigger was ≥ 1.5x).

## Memory / time

700 MB cap via `prlimit --as`, single thread, `nice`, each job ≤ 45 min.
