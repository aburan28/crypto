# Crossbred kernel frontier, 2026-09-18

X2 of `RESEARCH_ECC2K130_ROUTE_TARGETS.md`: does a crossbred space exist
past the toy rungs? Frozen stdout of
`cargo run --release --example crossbred_bench`.

**Boundary, stated before the run.** `kernel_dim > 0` with `k < m·ℓ`.
The bench's own correctness gate is the `agree` column. A row that is
not `yes` is not a result.

**Unit.** Bit operations, one 64-bit word XOR = 64 bit ops, as the
example defines. This is **not** the field-product unit of the ECC2K-130
scoreboard. Do not put `xb/F4` into `S`.

**Host.** `ip-172-31-19-103`, rustc 1.98.1 (48a229cea 2026-09-01),
source `82c846b6`.

Files:

- `default_bench.txt` — documented default ladder `9:2 13:2 7:3 9:3`
  plus the `(D, k)` sweep at `K_0/F_2^9`, `m = 2`
- `n11_m3.txt` — `11:3` skipped (`KoblitzCurve::new` / factor base `None`)
- `n13_m3.txt` — `13:3`, space extracted, `agree = NO`
- `extra_rungs.txt` — `5:3` and `5:2` agree; `17:2` skipped
- `summary.json` — the tables cited, not recomputed
