# No class arithmetic on the Jacobian of a plane quartic

- **Status:** `open`
- **Kind:** implementation gap
- **Affected work:** `research/ghs_cubic_genus3_build_20261012/` (Conductor T-2267); the follow-up it names, a genus-3 index-calculus count against rho on a toy cubic-extension instance
- **Next action owner:** unassigned
- **Last updated:** 2026-10-13 (T-2463: hyperelliptic-route evidence added)

## Problem

We cannot yet test numerically that the transfer `E0(F_{q³}) → Jac(C)(F_q)` is a homomorphism injective on the order-`ℓ` subgroup, nor run a discrete logarithm in `Jac(C)(F_q)`, because the repository has no arithmetic on divisor classes of a non-hyperelliptic genus-3 curve (a smooth plane quartic).

The instances are the nine `F_q` plane quartics of run `20261013T0320Z` (`q ∈ {5, 7, 11}`, `#E0(K) = 4ℓ`). The study verifies the descent, smoothness, the L-polynomial `1 − t_E T³ + q³T⁶` and that every transfer divisor is an `F_q`-rational effective divisor of degree 6 on the quartic. The homomorphism and kernel statements rest on a derivation only. The existing Mumford/Cantor code (`src/prime_hyperelliptic/`, the study-local labs) is hyperelliptic.

## Evidence and uncertainty

- **Observed:** `research/ghs_cubic_genus3_build_20261012/results/run-20261013T0320Z.jsonl`: 9 of 9 cases descend, are smooth, and have `#Jac(C)(F_q) = #E0(K)` from point counts over `F_q`, `F_{q²}`, `F_{q³}`.
- **Derived or proved:** `π0_* ∘ φ = [2]`, so `ker φ ⊂ E0[2]` (README of the study; not formalized, not measured).
- **Suspected explanation:** not applicable; this is a missing implementation.
- **Missing capability:** equality test and addition for degree-0 divisor classes on a plane quartic over `F_q` (for example reduction to an effective divisor of degree 3 plus a fixed base, via interpolating cubics and residual intersection), and a canonical form.
- **What this evidence does not establish:** any cost. Matching group orders and a derivation do not replace an executed transfer followed by a solved logarithm.

## What would unblock it

- [ ] A native routine that decides linear equivalence of two `F_q`-rational degree-0 divisors on a smooth plane quartic, validated on the nine quartics by `#Jac · D ~ 0` for random `D` and by associativity.
- [ ] With it: `φ(P + Q) ~ φ(P) + φ(Q)` on random pairs and `ℓ · φ(P) ~ 0`, `φ(P) ≁ 0` for a generator `P` of `E0(K)[ℓ]`, on at least one case per field.

These are sufficient for an executed transfer. They are not sufficient for a cost statement, which needs relation collection (lines through two factor-base points) and linear algebra modulo `ℓ`, counted in the same unit as rho.

## Next discriminating step

- **Action or proof obligation:** implement the linear-equivalence test and run the two checks above on the `q = 7` split case `u_idx = [7, 110]`.
- **Inputs and controls:** the quartics as recorded; control: a divisor that is not in the image (a random `F_q`-rational class) must fail `ℓ · D ~ 0` when `ℓ ∤ ord(D)`.
- **Execution limits and stop condition:** toy fields only; stop when both checks pass or one fails on a recorded input.
- **Outcomes and next actions:** pass → run Diem-style relation collection on the same quartic and count against rho; fail → the transfer bookkeeping in the study is wrong and its README is corrected; infrastructure failure → record here.
- **Recheck trigger:** a plane-quartic or general-curve Jacobian routine landing in the tree.

## Partial progress (2026-10-13, T-2280)

`research/ghs_cubic_quartic_ic_20261012/quartic_ic.rs` executes the **pushforward identity** `π₀∘φ = [2]` on E0 for toy smokes (`--check-lin-equiv`) and recovers planted scalars via **`transfer_pushforward_dlog`** when the PROTOCOL cubic target path fails. That is not a full Pic linear-equivalence test on `F_q` divisors.

**Still open:** native `F_q` divisor-class equality/addition; PROTOCOL **target cubic** decomposition when the three residual points form a **Galois orbit not in the factor base** (`n_res3_all_fq = 0` → `no_decomposition` on the pure IC path).

## Different route: the hyperelliptic `m = 3, r = 5` type (2026-10-13, T-2463)

`research/ghs_cubic_hyperelliptic_genus3_20261013/` builds the other genus-3
type of the cover census, which is hyperelliptic, as an explicit
`w² = h(z)` over `F_q` (degree 7 when `{a, b} ⊂ F_q`), and executes the
transfer `φ = (1 + Fr + Fr²)∘π0^*` with Cantor arithmetic (run
`20261014T0029Z`, `q = 5, 7, 11, 13`, 12 instances): homomorphism on
400/400 random pairs per instance, `φ(G) ≠ 0` and `ℓ·φ(G) = 0` for a
generator of `E0(K)[ℓ]`, `π0_*∘φ = [4]` on every point of `E0(K)`, and
`ker φ` of order 2 or 4 (so `φ` is injective on the `ℓ`-part but not an
isomorphism). This gives an executed genus-3 transfer **for that type
only**. It does not supply class arithmetic on plane quartics, so this
record stays open for the `m = 2` covers.

## Resolution and history

2026-10-12: opened with PR for T-2267. The index row in `problems/README.md` is not added in that PR because the file is held by Conductor T-2240; add it when that reservation clears.

Unresolved (Pic arithmetic and `F_q` target cubics).
