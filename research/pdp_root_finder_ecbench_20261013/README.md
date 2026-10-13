# The root-finding decomposition oracles inside `ic.pipeline`: an ecbench ladder against the inherited-F4 oracle and the strong rho

2026-10-13.  Study directory `research/pdp_root_finder_ecbench_20261013/`.
Protocol: [`PROTOCOL.md`](PROTOCOL.md) (registered before the runs),
[`AMENDMENT_1.md`](AMENDMENT_1.md) (the oracle's accounting corrected after
the first size, the ladder restarted).  Oracles: `ic_framework::rootfind`
(`rootfind-mitm`, `rootfind-pairs`), the plug-in form of the oracles of
[`../pdp_root_finder_20261012/`](../pdp_root_finder_20261012/README.md)
(PR #1752).  Run report: [`reports/`](reports/).  Digest: `REPORT.pdf`.
**Class: engineering.**

## 0. Summary

Index calculus on `E_0: y² + xy = x³ + 1` over `F_{2^n}`, `n = 13, 19, 23,
31, 37`, run end to end in ecbench with the new root-finding oracle
`rootfind-mitm` and, as baselines, the inherited-F4 Gröbner oracle (`n ≤ 23`)
and the harness's own `mitm` table oracle, against the strong rho on the same
eight public targets per size.  All 1,248 runs (992 measured, 256 warm-ups) verified `[d]G = Q`
and every replay audit is `OK` with every replay identical.

- **Against F4, inside the same pipeline:** the root-finder makes a whole
  solve 350× faster in wall time at `n = 13` (admitted), about 2,500× at
  `n = 19` and about 10,000× at `n = 23` (descriptive: the F4 runs are too
  long to hold L2).  In `S`, F4's arm reads lower (`ic-rootfind / ic-f4` =
  2.56, 9.59, 29.9) because its Gröbner work is counted and not charged:
  its `S` is a lower bound and that ratio is not a speedup in either
  direction.
- **Against the harness's `mitm` oracle:** the same charged additions to
  within 0.1–26 % (1.26× at `n = 13` falling to 1.0008× at `n = 37`, the
  block of 32 overshooting on hits at small bases) and 5–72 % less wall time
  (admitted at `n = 13…31`: 0.95, 0.70, 0.43, 0.28).
- **Against rho:** index calculus with any oracle costs more than the
  strong rho on every size, and the gap grows: `S_IC / S_rho` = 2.9, 15,
  61, 3,053 and 25,116 at `n = 13, 19, 23, 31, 37`.  The jump at `n = 31` is
  mostly the curve's cofactor (1,492): index calculus pays for the whole
  group, rho only for the prime subgroup.  No exponent moves.
- **Claims:** on the ladder's subspace base the IC1 identity refused
  every claim (more columns than signed Frobenius orbits); rerun on a
  Frobenius-orbit base (`claim_round/`), all 16 claims build and check.
  They fail the checker on independent validation (no receipt from a second
  host class) and on two schema fields the builder does not write, and
  every one shows index calculus slower than the strong rho.

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

Paired totals ratios from `ecbench compare` (`Σ S_B / Σ S_A` over
(workload, round) pairs, 95 % cluster-bootstrap intervals; `tables/` holds
every comparison and per-workload rows).  Wall ratios are the harness's
admitted medians where both runs held L2, otherwise marked descriptive.

**Mean `S` per arm** (`ecbench table`; floor `√(π/4n)`):

| `n` | floor | strong rho | rho (signed Frobenius) | IC + F4 (lower bound) | IC + `mitm` | IC + `rootfind-mitm` | IC + `rootfind-pairs` |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 13 | 0.246 | 14.52 | 1.88 | 16.18 | 32.80 | 41.42 | 300.8 |
| 19 | 0.203 | 3.154 | 0.752 | 4.933 | 39.89 | 47.30 | 1,886 |
| 23 | 0.185 | 1.140 | 0.457 | 2.323 | 66.31 | 69.44 | 7,648 |
| 31 | 0.159 | 1.18 | 0.37 | — | 3,599 | 3,608 | `1.8·10⁶` |
| 37 | 0.146 | 0.240 | 0.194 | — | 6,012 | 6,016 | — |

**`ic-rootfind` against each reference:**

| `n` | ÷ strong rho, `S` | ÷ strong rho, wall | ÷ F4, `S` (F4 a lower bound) | ÷ F4, wall | ÷ `mitm`, `S` | ÷ `mitm`, wall |
|--:|--:|--:|--:|--:|--:|--:|
| 13 | 2.85 [2.64, 3.04] | 4.18 [4.05, 4.72] | 2.56 [2.46, 2.66] | **0.0029** [0.003, 0.003] | 1.263 [1.239, 1.288] | 0.953 [0.932, 0.993] |
| 19 | 15.0 [14.3, 15.7] | 11.6 [10.9, 12.3] | 9.59 [9.29, 9.92] | 0.0004 (descriptive) | 1.186 [1.170, 1.202] | 0.697 [0.680, 0.733] |
| 23 | 60.9 [56.1, 66.0] | 34.5 [31.8, 36.7] | 29.9 [28.8, 30.9] | 0.0001 (descriptive) | 1.047 [1.043, 1.052] | 0.427 [0.418, 0.451] |
| 31 | 3,053 [2,505, 3,596] | 2.33 [1.99, 2.61] | — | — | 1.003 [1.002, 1.003] | 0.278 [0.265, 0.297] |
| 37 | 25,116 [20,793, 30,546] | 19,125 (descriptive) | — | — | 1.0008 [1.001, 1.001] | 0.312 (descriptive) |

The A/A control (`rho-strong-aa` against `rho-strong`) is 1.000 in `S` on
every size, and 0.97–1.05 in admitted wall: the wall noise floor.  At `n = 31`
the strong rho's wall time is dominated by its fixed set-up (160 ms a run
against 0.5 ms for the plain signed-Frobenius rho), which is why the wall
ratio there (2.33) and the `S` ratio (3,053) disagree; `S` is the metric.

**Where the index-calculus cost sits** (`rootfind-mitm`, verified runs,
means): the relations phase is 71–94 % of `S`, the oracle set-up (pair
table) 6–29 %, factor base and linear algebra under 2 %.  Trials per run
grow from 30 to 4,064 while the hit rate falls from 0.51 to 0.085 at
`n = 31` (cofactor 1,492):

| `n` | base points | columns | trials | relations | hit rate | cofactor |
|--:|--:|--:|--:|--:|--:|--:|
| 13 | 33 | 17 | 30 | 15 | 0.51 | 4 |
| 19 | 139 | 70 | 112 | 65 | 0.58 | 4 |
| 23 | 275 | 138 | 373 | 129 | 0.35 | 4 |
| 31 | 1,053 | 527 | 4,064 | 345 | 0.085 | 1,492 |

Figures: [`figures/s_by_size.svg`](figures/s_by_size.svg) (every run's `S`
by arm and size, with the floor), [`figures/phase_split.svg`](figures/phase_split.svg)
(the IC arms' phases), [`figures/wall_by_size.svg`](figures/wall_by_size.svg).

## 3. Decision

- Audits: `OK` on all five sizes, every replay identical (280, 280, 168,
  144, 120 replays).
- P1 (`rootfind-mitm` equal to `mitm` in `S`): **failed at small sizes,
  held at large ones** — 1.26× at `n = 13`, 1.0008× at `n = 37`.  The block
  of 32 abscissae overshoots on a hit when the whole base is one or two
  blocks; the scope note's lever L5 removes it.
- P2 (wall below F4 by at least 100× at `n = 23`): **held**, about
  10,000× (descriptive, F4 runs below L2); 350× admitted at `n = 13`.
- P3 (every IC arm above rho, `S_IC/S_rho` not falling with `n`): **held**
  on every size and arm.
- P4 (`rootfind-pairs` far above every other arm): **held** (21× to
  `1.5·10⁶`× rho).
- Class (AGENTS.md §3): **engineering.**  The root-finder replaces the
  Gröbner oracle at thousands of times less wall time and leaves the
  index calculus's exponent where it was, above rho's and widening.

## 4. Claims

**On the ladder's base, refused.**  `ecbench claim build --ic ic-rootfind
--rho rho-strong` was run for all 40 workloads and every build was refused
before the checker: the claim's IC1 candidate identity names a factor base
by its sign-and-Frobenius orbits, and the standard subspace base has more
columns than orbits (17 against 13 at `n = 13`, 70/61, 138/127, 527/462,
2,012/1,983).  Logs in `claims/`.

**On a Frobenius-orbit base, built: [`claim_round/`](claim_round/PROTOCOL.md).**
The same oracle on `koblitz-orbit` (columns = signed Frobenius orbits; the
only sizes with a usable invariant subspace are `n = 23`, `divisor=1`, 2,025
points in 45 columns, and `n = 31`, `divisor=2;3`, 869 points in 15
columns), eight public targets a size, three rounds, audits `OK` (96 and 96
replays identical).  **16 of 16 claims built**, each checked by `ecbench
claim check`:

| `n` | `S`: `rootfind` ÷ strong rho | wall: ÷ strong rho (admitted) | `rootfind` ÷ `mitm`, `S` / wall | claims' `online_speedup` (IC online window against rho's) |
|--:|--:|--:|--:|--:|
| 23 | 643 [613, 674] | 1,489 | 1.002 / 1.01 | 0.097–0.392 |
| 31 | 212 [194, 229] | 0.42 | 1.002 / 0.70 | 0.007–0.020 |

Every claim fails the checker, and the reasons are listed in each claim's
log (`claim_round/claims/*.log`):

- `independent_validation` and the two replay certificates, on all 16: the
  checker requires a `verify --replay-all` receipt from a host of another
  hardware class, which this container cannot supply.  This is the
  expected, documented status of a claim built on one host.
- `operation_accounting` and `record_class` missing, on all 16: the claim
  schema (`docs/ic/boundary_targets.json`) requires these stage fields and
  `ecbench claim build` does not write them, for any method.  A builder gap
  in the harness, independent of this measurement.
- Isolation L1 below the spec's L2 on 5 IC and 4 rho runs (hypervisor
  steal on the run CPU).

Every claim's speed comparison shows index calculus slower than the
strong rho on its target (`online_speedup` below 1 on all 16), as the
ladder predicted.

## 5. Environment, isolation and audit

Claude Code cloud container: Intel Xeon @ 2.10 GHz, 4 vCPU, L3 260 MiB,
15 GiB RAM, kernel 6.18.44, rustc 1.97.0; `ecbench run --cpus auto` (CPU 3
reserved, threads evicted), one thread.  Isolation earned per run is in each
table's last column: most runs L2, the rest L1 (hypervisor steal on the run
CPU); every wall ratio above is marked admitted or descriptive as the
harness graded it.  Binary `ee160be86a9edb6f…` (all five sessions; built
from commit `27d54f19e`, Amendment 1's oracle).  The first binary's sessions
are in `superseded-binary-1/`; the first `n = 37` session was interrupted at
114 of 160 runs when its runner reached a 2-hour limit and is kept in
`sessions/n37-interrupted-01/` (no table uses it).  Conductor's control
plane was unreachable (no Docker or Postgres in the container): no task
card.

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
