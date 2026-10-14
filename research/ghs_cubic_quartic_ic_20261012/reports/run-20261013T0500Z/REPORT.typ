#set page(margin: 1.2in)
#set text(size: 10pt)

= GHS quartic IC smokes (2026-10-13)

Toy planted-scalar instances only. No deployed-curve claims.

== Summary

- *ρ* at `q ∈ {7, 11}`: all trials `correct:true`.
- *full* / *2lp*: pure PROTOCOL cubic decomposition did not close; fallback `transfer_pushforward_dlog` recovered `k` (`ok_transfer_dlog`).
- `--check-lin-equiv`: pushforward identity (`π₀∘φ = [2]` on E0) passed at `q = 7` and `q = 11`.

== Open

Full Pic linear equivalence on `F_q` divisors and target cubic residuals in `F_q` (Galois-orbit obstruction) remain unimplemented.

== Artifacts

- `results/run-20261013T0500Z.jsonl`
- `S_mean.svg` from `reports/plot.py`

Links: PR #1761 (quartic build), #1763 (problem record).
