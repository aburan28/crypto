# No class arithmetic on the Jacobian of a plane quartic

- **Status:** `resolved` (toy smokes)
- **Kind:** implementation gap
- **Affected work:** `research/ghs_cubic_genus3_build_20261012/` (Conductor T-2267); follow-up IC study `research/ghs_cubic_quartic_ic_20261012/` (T-2280 / T-2491, PR #1768)
- **Next action owner:** unassigned
- **Last updated:** 2026-10-14

## Problem

We could not test numerically that the transfer `E0(F_{q³}) → Jac(C)(F_q)` is a homomorphism injective on the order-`ℓ` subgroup, nor run a discrete logarithm in `Jac(C)(F_q)`, because the repository had no arithmetic on divisor classes of a non-hyperelliptic genus-3 curve (a smooth plane quartic).

## What unblocked it (2026-10-14)

`research/ghs_cubic_quartic_ic_20261012/quartic_ic.rs` now provides:

1. **Pic / linear-equivalence checks** on `F_q`-rational transfer divisors (typically two Gal(`K/F_q`) orbits): `--check-lin-equiv` runs `pic_jacobian` (homomorphism and associativity via deg-6 divisor keys, `#Jac` / `ℓ` annihilators, cubic residual reduction smoke). Passes on `q ∈ {7,11}`.
2. **Pure IC target path** keyed by those Gal-orbit deg-6 divisors: planted-scalar recovery with `recovery:ic` and **no** `transfer_pushforward_dlog` fallback on `q=7` and `q=11` smokes.

Evidence: `research/ghs_cubic_quartic_ic_20261012/results/run-20261014T0200Z.jsonl`, README table.

## Remaining limitations

- Deg-3 Pic canonical forms are not yet unique for every class (some `φ(nG)` collide after reduction); dlog uses deg-6 orbit-pair keys.
- Dense FB-only line IC does not by itself reach the `ℓ`-part on these toys (`C(F_q)` is often 4-primary under `λ`).
- No deployed-curve or sub-ρ claims.

## Resolution and history

2026-10-12: opened with PR for T-2267.
2026-10-13: partial progress — pushforward identity / `transfer_pushforward_dlog` only (not full Pic; not pure IC).
2026-10-14: Pic lin-equiv checks + orbit-aware pure IC on toys (T-2491 / PR #1768). Status → `resolved` for the discriminating toy checks named in “What would unblock it.”
