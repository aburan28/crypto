# Matrix F5 / GF(2) opt-in defaults (2026-10-17)

Three exact switches had passed their recorded incremental gates on x86 CI
and were left opt-in only because their protocols tied a default change to
a 2× gate. This study measured each one on a second host and on the other
users of `gf2_elim`, and decided each against a rule
(`scripts/DECISION_RULE.md`) that was committed before the deciding run.

Report: [`REPORT.pdf`](REPORT.pdf) (source `REPORT.typ`). Runs:
`reports/run-20261017T0331Z` (3 rounds, all arms, plus matrix and
inherited F4) and `reports/run-20261017T0341Z` (15 rounds on F5; matrix F4
10 rounds).

| switch | outcome |
|---|---|
| `KIC_F5_DIRECT_PACK` | **on by default** (`=0` is the control): 1.006–2.173× on every F5 cell in run 1 |
| `KIC_F5_UNPACK_DIRECT` | stays opt-in: 0.986× on `f5_n24_m24_d4`, SelectiveEchelon, 4 threads (rule needs ≥ 1.00) |
| `KIC_GF2_REUSE_TABLE` | stays opt-in: 0.955× on matrix-F4 `n16_m24` (rule needs ≥ 0.97); F5 cells pass |

Every arm returned identical rows, bases, decisions and word counts in both
runs, so none of this changes a recorded stage unit. With all three on, the
largest F5 cases ran 1.11–1.27× faster than the previous defaults (run 2).
Most of that comes from direct pack.

Host: one shared 4-vCPU Xeon @ 2.10 GHz (AVX-512F), L0, wall times. Round to
round the times vary by up to ±25 %. No IC-versus-rho statement.
