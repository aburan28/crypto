# The root-finding decomposition oracles inside `ic.pipeline`: an ecbench ladder against the inherited-F4 oracle and the strong rho

2026-10-13.  Study directory `research/pdp_root_finder_ecbench_20261013/`.
Protocol: [`PROTOCOL.md`](PROTOCOL.md) (registered before the runs),
[`AMENDMENT_1.md`](AMENDMENT_1.md) (the oracle's accounting corrected after
the first size, the ladder restarted).  Oracles: `ic_framework::rootfind`
(`rootfind-mitm`, `rootfind-pairs`), the plug-in form of the oracles of
[`../pdp_root_finder_20261012/`](../pdp_root_finder_20261012/README.md)
(PR #1752).  Run report: [`reports/`](reports/).  Digest: `REPORT.pdf`.
**Class: RESULT-CLASS.**

## 0. Summary

RESULT-SUMMARY

## 1. What was measured

Index calculus end to end — factor base, oracle set-up, relations, linear
algebra, verification — in ecbench's one unit, `S` = charged
group-addition equivalents / `√r`, on `E_0: y² + xy = x³ + 1` over
`F_{2^n}` with the standard subspace factor base of dimension `d`
(`koblitz-standard-subspace`, the base of the repository's F4/F5/F6
workers), eight public one-target workloads a size, arms interleaved:

| arm | role | method |
|:--|:--|:--|
| `rho-strong` | reference | `rho.signed_frobenius_strong` (lanes 32, `dp_bits` 4, `step_cap_factor` 2000) |
| `rho` | baseline | `rho.signed_frobenius` |
| `ic-f4` | baseline, `n ≤ 23` | `ic.pipeline`, `pdp3-koblitz:m=3,engine=inherited-f4,degree=3,node_budget=8192` |
| `ic-mitm` | baseline | `ic.pipeline`, `mitm:m=3,negation_folded=1` (the harness's table oracle) |
| `ic-rootfind` | candidate | `ic.pipeline`, `rootfind-mitm:m=3` |
| `ic-rootfind-pairs` | candidate, `n ≤ 31` | `ic.pipeline`, `rootfind-pairs:m=3` |
| `rho-strong-aa` | control | `rho.signed_frobenius_strong` (the A/A wall-noise floor) |

| `n` | `d` | `icv1` slug | `log₂ r` | rounds | runs an arm |
|--:|--:|:--|--:|--:|--:|
| 13 | 5 | `icv1-f2m13-t181-515ee569` | 10.97 | 5 | 40 |
| 19 | 7 | `icv1-f2m19-t797-b6cf2467` | 17.00 | 5 | 40 |
| 23 | 8 | `icv1-f2m23-t5197-69e76b73` | 21.00 | 3 | 24 |
| 31 | 10 | `icv1-f2m31-tm90707-c95f16f5` | 20.46 | 3 | 24 |
| 37 | 12 | `icv1-f2m37-tm534059-32aad96b` | 27.78 | 3 | 24 |

The F4 arm's `S` is a **lower bound** (its solver's word XORs are counted,
not charged; the F6-IC ladder read it the same way).  The table and
root-finding oracles charge every group operation, so their `S` is
complete.  `rootfind-mitm` charges the same `2A` additions a refuted
target as `mitm:m=3` (one block more at most on a hit), so the two table
oracles are read as equal in `S` by construction; what the new id changes
is the arithmetic behind those additions (one shared inversion a block of
32 abscissae), which wall time sees and `S` does not, and the memoryless
three-summand arm the harness had not had.

## 2. Results

RESULT-TABLES

## 3. Decision

RESULT-DECISION

## 4. Claims

RESULT-CLAIMS

## 5. Environment, isolation and audit

RESULT-ENV

## 6. Reproduction

```sh
cargo build --release --bin ecbench
D=research/pdp_root_finder_ecbench_20261013
$D/run_sessions.sh 13 19 23 31 37        # sessions, replay-all audits, comparisons, tables, resources
$D/after_sessions.sh 13 19 23 31 37      # vs_rho claim reports and the database SQL
/tmp/report-venv/bin/python $D/figures/plot.py
```

Sessions are claimed atomically under the isolation lock; a directory that
exists is refused, so a rerun writes a new one.

## 7. Files

- `PROTOCOL.md`, `AMENDMENT_1.md` — the pre-registration and its one amendment.
- `specs/SPEC-n{13,19,23,31,37}.json` — the frozen specs.
- `sessions/n*/` — the sealed sessions (`spec.json`, `plan.json`,
  `records.jsonl`, `session.json`, `host.json`, `comparisons/`).
- `audits/n*.audit.json` — the replay-all receipts; `audits/n*.verify.log`.
- `tables/` — `ecbench table`, the saved `ecbench compare` outputs, `ecbench resources`.
- `claims/` — the `vs_rho` claim reports and the checker's output per workload; `db/` — the database SQL per session.
- `superseded-binary-1/` — the first binary's `n = 13` session and interrupted `n = 19` session (Amendment 1); not in any table.
- `figures/` — `plot.py` and the SVG/PNG figures; `reports/run-*/` — the run report with its graph and PDF.
