# Index calculus on the GHS plane quartic (toy planted scalar)

2026-10-14. Conductor **T-2491** (successor to T-2280). Dependency-free solver:
`rustc --edition 2021 -O quartic_ic.rs`.

Frozen protocol: [PROTOCOL.md](PROTOCOL.md) (do not edit in this PR). This directory began as the plan-only PR on main (no solver, no run); this README now records the smoke runs that execute the toy-scalar part of that plan.

## What ran

Runs on the first split instance per `q`, planted scalar in `<G> ⊂ E0(K)[ℓ]`,
costs in unit `S = ops / √ℓ`. Recovery `ic` means the PROTOCOL-style path on
`F_q`-rational transfer divisors (Gal-orbit pairs), with no
`transfer_pushforward_dlog` fallback.

### Recorded run 20261014T0241Z (q ∈ {5, 7, 11}, 10 targets per variant)

| q | variant | trials | outcome |
|---|---------|--------|---------|
| 5 | `rho` / `full` / `2lp` | 10 each | **PASS** — all `correct:true`; IC `recovery:ic`; mean rho steps 9.0; S = 650 896.03 |
| 7 | `rho` / `full` / `2lp` | 10 each | **PASS** — all `correct:true`; IC `recovery:ic`; mean rho steps 11.5; S = 2 751 662.59 |
| 11 | `rho` / `full` / `2lp` | 10 each | **PASS** — all `correct:true`; IC `recovery:ic`; mean rho steps 30.7; S = 37 920 187.89 |
| 5, 7, 11 | `--check-lin-equiv` | 1 each | **PASS** — `pic_jacobian` (homomorphism, associativity, `#Jac`/`ℓ` annihilators, cubic reduction smoke) |

90/90 trials correct. Run report: [reports/run-20261014T0241Z/REPORT.md](reports/run-20261014T0241Z/REPORT.md).
The rho arm records step counts only (not priced operations), so this vehicle
supports no IC-vs-rho `S` ratio and no "below rho" claim. Before this run,
the prior smoke run `run-20261014T0200Z` was reproduced exactly.

### Prior smoke run 20261014T0200Z (superseded, kept)

| q | variant | trials | outcome |
|---|---------|--------|---------|
| 7 | `rho` | 3 | **PASS** — all `correct:true` |
| 7 | `full` | 3 | **PASS** — `recovery:ic`, all correct |
| 7 | `2lp` | 1 | **PASS** — `recovery:ic` |
| 7 | `--check-lin-equiv` | 1 | **PASS** — `pic_jacobian` (homomorphism, associativity, `#Jac`/`ℓ` annihilators, cubic reduction smoke) |
| 11 | `rho` / `full` | 1 each | **PASS** — IC `recovery:ic` |
| 11 | `--check-lin-equiv` | 1 | **PASS** |

Latest JSONL: [results/run-20261014T0241Z.jsonl](results/run-20261014T0241Z.jsonl).
Toolchain: [results/.env.txt](results/.env.txt).

## What works

- **E0** arithmetic on `y² = u·u''` over `K = F_{q³}`; ρ recovers planted scalars.
- **Transfer** `φ` to `F_q`-rational effective degree-6 divisors on the plane quartic
  (typically two Gal(`K/F_q`)-orbits of size 3).
- **Pic / lin-equiv checks** on those `F_q` divisors: φ homomorphism and
  associativity via divisor keys, `#Jac·D~0` / `ℓ·φ(P)~0`, plus cubic residual
  reduction to a degree-3 representative when the pencil admits one.
- **Pure IC target path**: class logs keyed by the deg-6 Gal-orbit transfer
  divisor; three-target recovery of `k` with `recovery:ic` (orbit residuals, no
  pushforward dlog fallback).

## What remains open

- Unique deg-3 Pic canonical forms for all classes (some distinct `φ(nG)` still
  collide after cubic reduction; deg-6 keys are used for dlog).
- Dense FB-only line IC: `C(F_q)` is too sparse in the `ℓ`-part on these toys
  (`λ` of rational points is often 4-primary); orbit-aware transfer classes carry
  the `ℓ`-dlog.
- PROTOCOL cost scaling vs ρ at larger `q` (no crossover claim).

No deployed-curve or sub-ρ security claims.

## Build and run

```bash
cd research/ghs_cubic_quartic_ic_20261012
rustc --edition 2021 -O quartic_ic.rs
./quartic_ic 7 10 full
./quartic_ic 7 1 --check-lin-equiv
python3 reports/plot.py results/run-20261014T0241Z.jsonl reports/run-20261014T0241Z
```

Links: PR #1761 (quartic build), #1763 (problem record), #1768 (this study).
