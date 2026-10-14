# Index calculus on the GHS plane quartic (toy planted scalar)

2026-10-12. Conductor **T-2280**. Dependency-free solver: `rustc --edition 2021 -O quartic_ic.rs`.

Frozen protocol: [PROTOCOL.md](PROTOCOL.md) (do not edit in this PR). This directory began as the plan-only PR on main (no solver, no run); this README now records the smoke runs that execute the toy-scalar part of that plan.

## What ran

Smokes on the first split instance per `q`, planted scalar in `<G> ⊂ E0(K)[ℓ]`, costs in unit `S = ops / √ℓ`.

| q | variant | trials | outcome |
|---|---------|--------|---------|
| 7 | `rho` | 3 | **PASS** — all `correct:true` |
| 7 | `full` | 3 | **IC rank/decomposition incomplete** — fallback `ic_fail_transfer_dlog_ok`, all correct |
| 7 | `2lp` | 3 | fallback `ok_transfer_dlog`, all correct |
| 7 | `--check-lin-equiv` | 1 | **PASS** — pushforward identity check |
| 11 | `rho` / `full` / `2lp` | 1 each | rho PASS; IC variants fallback, correct |
| 11 | `--check-lin-equiv` | 1 | **PASS** |

Latest JSONL: [results/run-20261013T0500Z.jsonl](results/run-20261013T0500Z.jsonl). Toolchain: [results/.env.txt](results/.env.txt).

## What works

- **E0** arithmetic on `y² = u·u''` over `K = F_{q³}`; ρ reference recovers the planted scalar.
- **Transfer** to six `F_q` W-points on the plane quartic; lift inverse `w_to_ep` for pushforward bookkeeping.
- **IC relation collection** on the factor base (line scan on `C(F_q)`).
- **`--check-lin-equiv`**: pushforward identity `π₀∘φ = [2]` on E0 (not full Pic on `F_q` divisors).
- **`transfer_pushforward_dlog`**: recovers `k` when PROTOCOL cubic decomposition fails (Galois-orbit residual obstruction).

## What remains open

- Target cubic over `F_q` with residual orbit not in FB → pure IC `no_decomposition`.
- Full Pic linear equivalence on plane quartic divisors ([problem record](../../problems/plane-quartic-jacobian-class-arithmetic-missing.md)).

No deployed-curve or sub-ρ security claims.

## Build and run

```bash
cd research/ghs_cubic_quartic_ic_20261012
rustc --edition 2021 -O quartic_ic.rs
./quartic_ic 7 3 full
./quartic_ic 7 1 --check-lin-equiv
python3 reports/plot.py results/run-20261013T0500Z.jsonl reports/run-20261013T0500Z
```

Links: PR #1761 (quartic build), #1763 (problem record).
