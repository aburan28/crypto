# n=63 rung attempt: mathematically non-viable on this curve family (2026-10-02)

The ledger's next rung after the n=61 record was n=63, "the last u64-packing
rung" (`MAX_N=63`). The attempt fails closed for a mathematical reason, not an
implementation one: **neither Koblitz curve `K_a` over `F_{2^63}` has a usable
prime-order subgroup.**

- `KoblitzCurve::new(0, 63)` and `KoblitzCurve::new(1, 63)` both return `None`
  (fail-closed panics in the smoke receipts below, each under 0.1 s — the
  `r ≤ cofactor` guard fires long before any expensive work).
- `#E_0(F_{2^63}) = 9223372041104766164 = 2^2·29·43·127^2·421·757·359731` —
  largest prime 359,731 (~2^18.5), cofactor ~2.6·10^13.
- `#E_1(F_{2^63}) = 9223372032604785454 = 2·7^3·37·43·71·379·631·497701` —
  largest prime 497,701 (~2^18.9), cofactor ~1.9·10^13.

Composite `n = 3^2·7` splits both group orders into small cyclotomic-type
factors; there is no large prime-order subgroup to attack. The same holds at
`n = 62` (both `a`, largest prime 1,439,393 ~ 2^20.5). An irreducible field
polynomial does exist at degree 63 (the sparse search finds the trinomial
`x^63 + x + 1`; `is_irreducible_f2` and an independent sympy check agree) —
the obstruction is purely the subgroup structure.

**Conclusion:** `n = 61` is the u64-packing ceiling for this curve family, and
the landed n=61 compact-orbit record is the last u64 rung. The next rung is
beyond-63 multi-word field arithmetic; the natural target is **n=71, a=0**
(prime degree; `r = 5513228015079457` ~ 2^52.3 with cofactor 428,276),
which requires u128 packing.

## Files

- `curve_orders.json` — order/factorization/viability map for
  n ∈ {59, 61, 62, 63, 67, 71} × a ∈ {0, 1}; the trace recurrence is pinned
  against the landed n=61 fixture order (2305843009055095052).
- `smoke_time.txt`, `smoke_a1_time.txt` — fail-closed producer receipts
  (`construct:63:0:600` / `construct:63:1:600` panics).
- `smoke_out.jsonl`, `smoke_a1.jsonl` — empty (producers failed before
  emitting rows).
