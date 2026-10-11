# PDP follow-ups to the SAT/SMT/Nagao survey (2026-10-10/11)

The seven engineering items of
[`RESEARCH_PDP_SAT_SMT_NAGAO_SURVEY.md`](../notes/index-calculus/RESEARCH_PDP_SAT_SMT_NAGAO_SURVEY.md) §6,
attempted in one cloud session.  Conductor task T-2.  Every row is a
**stage diagnostic**: no `S`, no end-to-end cost, no `vs_rho` figure, and no
row moves a ratio to the floor.  Three items ran to a verdict; four are
blocked by the environment (no route to `github.com` or `crates.io` beyond
the attached repositories; no GPU) and carry frozen protocols instead.

| # | item | status | where | class |
|--:|:--|:--|:--|:--|
| 5 | complete the SMT ladder to `n = 23` | **done** — P1 pass, P2 fail, P3 fail (cvc5-bv1 5 vs native 4 found at `n = 19`), P4 unscorable; nothing decided at `n = 23` by any arm | [`smt_ladder/`](smt_ladder/), [`../pdp_smt_oracle_20261008/README.md`](../pdp_smt_oracle_20261008/README.md) | engineering candidate (constant), by the registered rule |
| 6 | yield-law census vs thesis Table 7.7 | **done** — mean yield law exact; satisfiable fraction `0.56–0.83×` the heuristic, `0.64×` Poisson at `λ ≈ 1.3`; thesis column is `#sat/#unsat` and its `2×` is not reproduced | [`yield_census/`](yield_census/) | accounting |
| 3 | Hamming-weight ONB base under CryptoMiniSat / CaDiCaL (native emitter) | **done** — 128 cells, 75 verified relations, zero disagreements; at `n = 23, m = 3` neither solver decides any of 8 targets in 180 s (P2 fails); every proved refutation costs `117 … 174×` the `m - 1` floor (P3 passes); PB-XOR arm blocked | [`hamming_nb_sat/RESULTS.md`](hamming_nb_sat/RESULTS.md), [`../../examples/hamming_nb_sat_panel.rs`](../../examples/hamming_nb_sat_panel.rs) | boundary |
| 7 | GPU pairs-and-solve leaf | **blocked** (no GPU); CPU reference arm measured to `ℓ = 11` (0.651 s per refuting target at `ℓ = 11`, `2^{21}` pairs); **scoped as PS-5 of [`../pdp_oracle_speedup_20261011/PLAN.md`](../pdp_oracle_speedup_20261011/PLAN.md)** | [`blocked/PROTOCOLS.md`](blocked/PROTOCOLS.md), [`blocked/pairs_and_solve_cpu_l11.txt`](blocked/pairs_and_solve_cpu_l11.txt) | pending (scoped) |
| 1 | XOR-learning CDCL (Xorricane) | **blocked** (source unreachable); protocol frozen | [`blocked/PROTOCOLS.md`](blocked/PROTOCOLS.md) | pending |
| 2 | WDSat `-x` audit, then SATIC | **blocked** (source unreachable; Python runner needs a native replacement first) | [`blocked/PROTOCOLS.md`](blocked/PROTOCOLS.md) | pending |
| 4 | Bosphorus preprocessing | **blocked** (source unreachable) | [`blocked/PROTOCOLS.md`](blocked/PROTOCOLS.md) | pending |

Tooling obtained in-session (Ubuntu 24.04 packages unless stated): cvc5
1.1.2, CryptoMiniSat 5.11.15, CaDiCaL 1.7.3, Z3 4.8.12; Z3 5.1.0 binary from
the `z3-solver` wheel on PyPI.  Kissat is not packaged.
