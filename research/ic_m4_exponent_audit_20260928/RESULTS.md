# m = 4 exponent audit — results

This file is written after phase 2 finished. [PREREGISTRATION.md](PREREGISTRATION.md) is
unchanged, and so are its Amendments 1 and 2. The raw readout is
[runs/readout.txt](runs/readout.txt), printed by `analyze.py` and byte-identical to the run
log.

## Run

| | |
|:--|:--|
| engine | `m4_exponent_audit-2809b498`, sha256 `79ccfdf3…acb55e` (Amendment 1), from `runs/binary.sha256` |
| branch head at launch | `b1133883` (`runs/branch_head.txt`) |
| wall time | 2026-09-29 07:17:29Z to 07:38:57Z |
| commands | `run_audit.sh`, the §8 list in `phase2-commands.txt` |
| cells | the nine registered cells; three are structurally excluded as registered (`K_0/2^11`, `K_0/2^17`, `K_1/2^13`: no such curve in the tooling) |

## Registered verdict

**Closed for this engine at these sizes.**

The primary fit, over Semaev `m = 4` word XORs on all nine cells, gives
`ĉ = 0.985` bits per unit `n`. The bootstrap 95% band is [0.983, 1.018], with
`B = 10,000` and no invalid replicates. The band's lower end is above `c* = 0.25`, so §6
reads *closed*. §5.2 lets a closed verdict stand whatever the null gives. The null passed
anyway.

| cell | targets | censored | lower median (word XORs) | log₂ |
|:--|--:|--:|--:|--:|
| `K_0/2^9`, ℓ = 2 | 16 | 0 | 68,235 | 16.06 |
| `K_0/2^13`, ℓ = 3 | 16 | 0 | 1,493,937 | 20.51 |
| `K_0/2^15`, ℓ = 4 | 16 | 0 | 5,259,157 | 22.33 |
| `K_0/2^19`, ℓ = 5 | 16 | 0 | 50,076,470 | 25.58 |
| `K_1/2^9`, ℓ = 2 | 16 | 0 | 39,941 | 15.29 |
| `K_1/2^11`, ℓ = 3 | 16 | 0 | 719,041 | 19.46 |
| `K_1/2^15`, ℓ = 4 | 16 | 0 | 5,079,215 | 22.28 |
| `K_1/2^17`, ℓ = 4 | 16 | 0 | 17,475,739 | 24.06 |
| `K_1/2^19`, ℓ = 5 | 16 | 0 | 80,892,472 | 26.27 |

No Semaev target censored. Every cell is in the fit.

## Controls and instrument checks

- **Baseline reproduction (§5.1).** Passed exactly in phase 1, before this run.
- **Enumeration null (§5.3).** `ĉ_enum = 0.828`, band [0.828, 0.845], matching the model's
  0.828 and inside the registered [0.60, 1.00]. **Pass.**
- **Random-system null (§5.2).** All eight targets censored at every cell. At each retained
  Semaev cell the null is therefore dearer than its budget, so the control **passes**. The
  null's growth was not measurable within its budget, which §5.2 says in advance is not a
  failure.
- **Oracle against enumeration truth (§5.4).** 0 disagreements.
- **Group re-addition (§5.4).** 0 returned decompositions failed it.

## Reported, not decisive

- **Refuted-only fit.** `ĉ = 0.984`, band [0.983, 0.986].
- **Prime `n` only** (`K_0/2^{13,19}`, `K_1/2^{11,17,19}`). `ĉ = 0.833`, band
  [0.830, 0.953]. It gives the same reading, *closed*, as the primary fit. The point estimates
  differ. On nine points that span two curves and both composite and prime `n`, this difference
  is one of §10's listed confounds. It is not explained further here.
- **Satisfiable-only.** Not fitted, as §4 expected, because some cells have fewer than 3
  satisfiable targets. There were 0–6 per cell.
- **Refutation degree (secondary arm, §6).**
  - The exact Boolean root count matched the engine's full-tree count on every target that
    ran: 0 faults.
  - At `ℓ = 2`, both cells, `K_0/2^9` and `K_1/2^9`, resolved their first unsatisfiable
    target at `D = 6`, the registered `d_max`, in 200 and 219 s.
  - Both `ℓ = 3` cells, `K_0/2^13` and `K_1/2^11`, hit the 300 s per-cell CPU limit
    before finishing a target. All four `degree` processes were killed by that limit, which is
    machine protection (§7), so these targets are **censored, not negative**.
  - With resolved degrees at only one `ℓ`, the slope `s` is not fitted.

## What this says, and what it does not

Per §10, **this closes `m = 4` for the frozen `2809b498` engine at `n ≤ 19`.**

- It says nothing asymptotic.
- It says nothing about another solver family, or about the branch-head engine.
- It does not close the algebraic-decomposition route.

One observation lies outside the decision rule and is recorded only as an observation. At
these sizes the Semaev cost grows about as fast as `2^n` (`ĉ ≈ 0.98`).

- That is faster than the enumeration null's 0.83, which counts in point additions rather
  than word XORs.
- It is also faster than the `n/2` exponent of generic rho.
- Within this engine and range, then, the Gröbner route is not only above the `c* = 0.25` bar
  but also grows faster than both simple alternatives.
- The units differ, and the fit is nine points on two curves.

## Successors (§11, unchanged)

- **`m = 5`.** The next rung, at `c* = 0.30`.
- **Branch-head engine.** A rerun of these nine cells on it would test whether the post-freeze
  engine changes move the slope. It needs its own reproduction baseline
  (`baseline/head-1722bad1`). The two-cell drift in §5.1 does not suggest a lower slope.
- **Reopening conditions.** Any of these reopens `m = 4`:
  - an engine that measures below 0.25 at `m = 4`;
  - a published degree bound for Weil-descended chains;
  - §3.2's symmetry lever, if it lowers the slope `s` rather than cutting `D` by a constant.
