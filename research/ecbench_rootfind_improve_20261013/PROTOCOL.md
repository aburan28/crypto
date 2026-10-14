# Protocol: ecbench rootfind improve loop (2026-10-13)

Written before the improve-loop sessions. Extends
`research/ecbench_rootfind_engine_20261012/` (PR #1765). Classification
target: **engineering**. No key-recovery, asymptotic sub-rho, or
deployed-curve-security claims. Rootfind/Gröbner `S` ratios are not
wall-time speedups.

## Diagnosis (from #1765 e2e)

- Against `mitm` / `subtract`, ops ratio is **1.000** (same query stream).
- Against inherited-F4, cold `S` is higher because rootfind **charges** the
  PairSumTable build (`oracle_setup` ≈ `|F|(|F|+1)/2` GAE) while F4's
  Boolean word XORs are unpriced (lower bound). Descriptive L0 wall still
  favours rootfind (especially `m = 3`).
- Highest-leverage levers that stay claim-hygienic: (1) smaller `|F|` via
  subspace `dimension`, (2) `m = 2` vs `m = 3`, (3) cyclic third-summand
  window (`engine=rootfind-window`), (4) secondary amortized reading of
  setup vs online phases / multi-target, (5) honest F4 wall baseline.

## Engine ids

| id | role |
|---|---|
| `pdp3-koblitz:m=2\|3,engine=rootfind` | committed full scan (#1765) |
| `pdp3-koblitz:m=3,engine=rootfind-window` | **new**: `⌈\|F\|/3⌉` cyclic window |
| `rootfind-subtract:m=2\|3` | unchanged prime arm |

## Cycles

1. **Select:** smoke `smoke-koblitz-rootfind-select.json` — dimensions
   3/4/5 at `m = 2`, full vs window at `m = 3` d=5, F4 m=3 baseline.
2. **Claim panel:** `claim-koblitz-rootfind-improve-n13.json` with the
   selected `|F|` and window arm; before/after vs #1765 tables; amortized
   secondary analysis from resources phases.

## Unit and reading rules

As #1765 PROTOCOL: floor `√(π/2A)` in `S`; reference
`rho.signed_frobenius_strong`; cold `S` with every phase charged.
Amortized / online figures that subtract `oracle_setup` are **secondary**
and labelled as such. Wall on macOS is L0 / descriptive only.
