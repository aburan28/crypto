# Yield sweep: results

**Decision: the ledger is populated and the predictions are scored. Four
of five hold; the one that fails is a finding.** On the same subspace
base and the same targets, every exact oracle (subtract, meet-in-the-middle,
and the algebraic descent solved by Buchberger, F4 or crossbred) yields
the same relations per trial, so yield is a property of the base and the
field, not of the oracle. The SAT solver does not: on 4 of its 32
verified cells it reports a different relation count from the exact
oracles on identical trials, and it times out on 11 of 32 runs at the
largest sizes. Solver cost per decomposition *falls* with the field
degree at fixed base dimension, the opposite of what was registered,
because the system gains an equation per bit of `n` while its variable
count stays `2d`.

This is a stage diagnostic (AGENTS.md §2, §6): it ranks bases, oracles
and solvers on identical targets and prices nothing new into `S`. Every
descent arm's `S` is a lower bound because its solver work is counted in
the solver's own unit and left unpriced; the tables say so.

| | |
|---|---|
| protocol | [`PROTOCOL.md`](PROTOCOL.md), committed before the sessions ran |
| specs | [`spec-koblitz.json`](spec-koblitz.json), [`spec-prime.json`](spec-prime.json) |
| binary | `ecbench`, release, built from `cd0228d4c` (the commit that adds the record's `solver` block); every record carries its SHA-256 |
| host | Intel Xeon @ 2.80 GHz, 4 vCPU virtual machine, Linux x86-64; `--cpus auto` as root; levels L0–L2 by run, wall time descriptive |
| sessions | [`sessions/koblitz`](sessions/koblitz) `ECBS1h0d1ce1d04e54` (208 records: 197 verified, 11 timeout), [`sessions/prime`](sessions/prime) `ECBS1h548d49a71881` (84 verified of 84) |
| correctness | every verified execution recovered the planted scalar; the 11 timeouts are `sat-cdcl` descent runs that exceeded 240 s, kept as their own status |
| replay certificates | `audit-koblitz.json` `7ee590cf4dc52ca8bc8ef76793f9981e3f84df60b436e3d653b50e9abefb695d`, `audit-prime.json` `9ec5c00e0eebe0d0f45bfd8deacea65806453a3787644c701e49b19b5cdb2d07`; each `ok`, 24 of 24 replays identical, algebraic arms included |
| tables | [`table-koblitz.md`](table-koblitz.md), [`table-prime.md`](table-prime.md), by `ecbench table` |
| ledger | the Yield view of the lab browser (`docs/browser/#yield`), and `ic_yield` in the `ecbench` database; `run_solver` holds the solver blocks |

## Koblitz: yield by base, field and oracle

Mean yield (relations per trial) over 4 targets, with the coefficient of
variation across targets; one figure per base and curve because every
exact oracle on that base reads the same. `d` is the subspace dimension;
the system an algebraic arm solves has `2d` variables and `n` equations
of degree 2.

| curve | log₂ r | d = 6 yield (CV) | d = 8 yield (CV) |
|---|---:|---:|---:|
| `icv1-f2m13-t181-515ee569` | 11.0 | 30.73 % (0.15) | 97.89 % (0.02) |
| `icv1-f2m17-tm101-00378d4e` | 16.0 | 1.549 % (0.27) | 19.57 % (0.09) |
| `icv1-f2m19-t797-b6cf2467` | 17.0 | 0.389 % (0.20) | 7.473 % (0.15) |
| `icv1-f2m23-t5197-69e76b73` | 21.0 | 0.040 % (0.15) | 0.443 % (0.07) |

The SAT arm differs where it verified: `d = 8` reads 96.35 % at n = 13
and 19.578 % at n = 17 (exact oracles 97.89 % and 19.57 %); see P1.

## Koblitz: what each oracle pays per decomposition, same base and targets

Per target, means over the four targets. Lookups are the table oracles'
unit; the solver arms report their own unit per call and the whole
pipeline's `S` as a lower bound (solver work unpriced). The reference is
the strong signed-Frobenius rho on the same point.

| curve | base | subtract `S` | mitm `S` (lookups / relation) | descent `S`, lower bound | sat-cdcl conflicts / call | buchberger-f2 monomial ops / call | f4-f2 word XORs / call (max Macaulay rows) | crossbred-f2 ops / call (max Macaulay rows) | rho `S` |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| f2m13 | d6 | 62.98 | 59.15 (3.3) | 13.9 | 899 | 137,759 | 10,564 (697) | 23,696 (8,788) | 14.16 |
| f2m13 | d8 | 44.58 | 737.2 (1.0) | 17.3 | 2,668 | 183,605,767 | 18,166,674 (5,137) | 39,852 (12,155) | 14.16 |
| f2m17 | d6 | 227.9 | 15.30 (69.1) | 7.84 | 1,111 | 5,566 | 297 (319) | 33,235 (403,734) | 3.85 |
| f2m17 | d8 | 294.2 | 117.8 (5.1) | 5.97 | 14,615 | 3,021,777 | 362,021 (2,441) | 63,956 (117,912) | 3.85 |
| f2m19 | d6 | 873.2 | 23.65 (261) | 18.1 | 1,020 | 4,718 | 156 (307) | 37,803 (1,444,397) | 3.02 |
| f2m19 | d8 | 717.4 | 112.7 (13.6) | 6.21 | 14,928 (3 of 4 timed out) | 925,957 | 48,906 (2,205) | 79,904 (421,311) | 3.02 |
| f2m23 | d6 | 1,894.5 | 28.74 (2,554) | 26.9 | timeout (4 of 4) | 4,525 | 118 (289) | 47,715 (15,160,652) | 1.17 |
| f2m23 | d8 | 2,422.2 | 36.26 (226) | 10.4 | timeout (4 of 4) | 18,341 | 1,008 (1,149) | 119,736 (5,486,624) | 1.17 |

Read the descent column with its caveat: `S` there counts the group
operations of the walk, the base, the matrix and the verification, and
not the solver's work, which has no pinned ratio to the unit. The solver
columns are that work, in each solver's unit; converting them is the
pricing question §12 of the boundary ledger leaves open. The f4 and
crossbred Macaulay figures are the largest matrix any call built (the
record's `macaulay_rows`); crossbred's "ops" are its partial count.

Solving degrees: F4 and Buchberger reached degree 3 on every base except
`d = 8, n = 13` (degree 4), against semi-regular degrees of 3 to 5 for
the shapes recorded; the solving degree never exceeded the semi-regular
bound, and at `d = 6, n = 19` it equalled it.

## Prime field: yield and table cost against base size

| curve | log₂ r | base | yield (CV) | mitm lookups / relation | subtract lookups / relation | mitm `S` | subtract `S` | rho `S` |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| `icv1-fp16-t295-8d3c3165` | 15.7 | 16 | 1.084 % (0.18) | 94 | 2,998 | 11.08 | 160.9 | 2.39 |
| | | 32 | 4.289 % (0.22) | 24 | 1,504 | 10.15 | 119.9 | |
| | | 64 | 15.58 % (0.12) | 6.5 | 743 | 22.63 | 118.9 | |
| `icv1-fp18-tm337-d28d5e09` | 17.7 | 16 | 0.207 % (0.37) | 546 | 17,443 | 15.48 | 376.9 | 1.39 |
| | | 32 | 0.939 % (0.12) | 108 | 6,845 | 10.33 | 330.5 | |
| | | 64 | 3.501 % (0.18) | 30 | 3,706 | 13.46 | 282.4 | |
| `icv1-fp20-t727-cd198a38` | 19.1 | 16 | 0.102 % (0.19) | 996 | 31,860 | 19.59 | 552.8 | 1.13 |
| | | 32 | 0.373 % (0.22) | 278 | 17,735 | 12.94 | 592.2 | |
| | | 64 | 1.527 % (0.13) | 67 | 8,437 | 10.80 | 453.7 | |

## Predictions, scored

| | statement | verdict |
|---|---|---|
| P1 | exact oracles agree on yield cell by cell; a budget-exceeded SAT verdict breaks it | **pass for subtract, mitm, Buchberger, F4 and crossbred**: identical trial and relation counts on every verified cell. **SAT disagrees on 4 of its 32 verified cells** (n = 13 targets 1–3, n = 17 target 1): on 53 identical trials it reports 49 relations against 50, and on the other three the walk diverged after the first disagreement (70 vs 55, 33 vs 29, 323 vs 410 trials). No call reported a budget-exceeded verdict, so this is a verdict the SAT oracle returned that the exact oracles contradict. Not explained here; a decomposition the SAT arm accepts or rejects differently from an exact oracle on the same point is what AGENTS.md §6 requires be cross-checked before SAT can stand in for an exact oracle. |
| P2 | yield falls with n at fixed d; CV across targets below 0.5 | **pass**: d = 6 falls 30.7 → 1.55 → 0.39 → 0.040 %, d = 8 falls 97.9 → 19.6 → 7.5 → 0.44 %; every CV is 0.02–0.27. |
| P3 | solver cost per call rises with n at fixed d; F4's solving degree stays within the semi-regular bound | **first clause fails**: Buchberger's operations per call fall with n at both dimensions (137,759 → 5,566 → 4,718 → 4,525 at d = 6; 184 M → 3.0 M → 0.93 M → 18 k at d = 8), F4's word XORs likewise, and SAT's conflicts per call are flat at d = 6 and rise then time out at d = 8. The system has `2d` variables and `n` equations, so a larger field makes it more overdetermined and cheaper per call; the n = 13 systems, nearly square, are the expensive ones. **Second clause passes**: the solving degree never exceeded the recorded semi-regular degree. |
| P4 | prime yield rises roughly in proportion to base size; folded-table lookups per relation fall | **direction passes, the rate does not**: yield rises about fourfold per doubling of the base (1.08 → 4.29 → 15.6 %, and likewise on the other curves), as two-summand pairs over a base of `B` points predict, not in proportion. Lookups per relation fall with base size on every curve, for both oracles. |
| P5 | every verified execution correct; every replay identical | **pass**: 281 of 281 verified executions recovered the planted scalar; 48 of 48 replays identical, algebraic arms among them. |

## What this establishes, and what it does not

- **Established.** Yield on a subspace base is set by the base and the
  field, not by the oracle: every exact oracle returns the same relations
  from the same trials. The cheapest way to decide decomposability on
  these bases, counted in the unit, is a table (meet-in-the-middle) at
  small sizes and the algebraic descent as the base grows, with the
  caveat that the descent's solver work is in its own unit. The
  algebraic systems get cheaper per call as `n` grows at fixed `d`.
  SAT's verdicts are not interchangeable with the exact oracles' on
  these systems, and SAT is the only solver that timed out.
- **Not established.** Any speed: the descent arms' `S` excludes the
  solver and the table arms' `S` is above rho on every curve but n = 13
  (table-koblitz.md). Anything about orbit bases (not in this session)
  or above `r = 2^21`. Why the SAT verdicts differ.
- **For a candidate combination.** The ledger now answers, per target,
  "which base decomposes this point, and what did each oracle pay to
  decide it". The next question it cannot answer alone is the pricing
  of solver work in group-addition equivalents; until a ratio is pinned,
  a descent arm's `S` stays a lower bound and no combination built on
  it can claim a speed.

## Reproduce

```bash
cargo build --release --bin ecbench
./target/release/ecbench verify --dir research/ecbench_yield_sweep_20261004/sessions/koblitz --replay 24 --exit-code
./target/release/ecbench db sql research/ecbench_yield_sweep_20261004/sessions/* | sqlite3 -bail yield.db
sqlite3 yield.db "SELECT curve_slug, fb_family, oracle, round(avg(yield),4) FROM ic_yield WHERE status='verified' GROUP BY 1,2,3"
```
