# Protocol — rootfind faster e2e (2026-10-15)

**Classification.** Engineering until a ledger-complete `vs_rho` claim is rebuilt.
**Hygiene.** Public synthetic / known-answer only; fail-closed verify every solve;
no key-recovery, asymptotic sub-rho, or deployed-curve security claims. Multi-target
amortization is secondary to cold single-target.

## Goal

Improve **real CPU time** of Semaev root-finding decomposition and end-to-end IC,
not just accounted `S`. Chase the ~1.54× cold complete-IC vs rho-strong gap
honestly; document hard stops.

## Measurement

1. **Oracle / e2e:** `examples/semaev_rootfind_e2e_bench.rs` (CPU time via
   `getrusage`, not wall). Knobs: `E2E_SKIP_PRIME`, `E2E_SKIP_BINARY`,
   `E2E_SKIP_MULTI`, `E2E_MULTI_N`, `BENCH_OUT=…`.
2. **ecbench:** claim panels vs `rho.signed_frobenius_strong` and XOR-priced F4;
   complete costing; `--cpus none --allow-busy` → L0 on macOS (ops primary).
3. Evidence under `research/rootfind_faster_e2e_20261015/cycles/`.

## Non-goals / known dead ends (do not re-try)

`window_div`, `window@d4`, denser d=6/7, plain rootfind vs negfold, orbit subtract
(see `research/ecbench_rootfind_{improve,parity,faster,vs_rho}_*`).

## Cycles

Each cycle: change → measure → fail-closed verify → commit+push → update PR.
