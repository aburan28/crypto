# Disclosed-point F4/F5 standard-subspace recovery pilot

Status: **registered plan only; not dispatched**. This protocol freezes a
bounded diagnostic before any fresh competitive panel. It does not consume the
one-dispatch slot of
[generic-backend-qualification-v2](../generic-backend-qualification-v2/PROTOCOL.md)
(seed `2026092902`, Actions run
[36580669479](https://github.com/aburan28/crypto/actions/runs/36580669479)).
Never redispatch that seed.

## Dependencies

1. Merge and retain the source-bound encoder gate and dimension-6 inventory
   control from the feasibility PR
   ([#952](https://github.com/aburan28/crypto/pull/952) or its successor):
   `generic_solver_feasibility.py`, `STATIC-FEASIBILITY.md` / `.json`, and
   `standard-subspace-d6-inventory-control.json`.
2. Reuse the reviewed generic worker commit
   `765c3c5f19032bd852163805f257c56babef2040` with the same pinned source
   objects as v2. Build with `generic_build.py` and keep the receipt.
3. Run only on **previously disclosed** public points. Do not generate a new
   tournament seed or call `tournament.py prepare` for fresh fixtures.

## Hypothesis and boundaries

**Hypothesis.** On the five toy cells, replacing the ambient-length
`subgroup_orbits` base with `standard_subspace` dimension 6 makes the pinned
`m=3` Semaev template fit under `MAX_VARS=64`, so at least one of `f4` / `f5`
/ `inherited_f4` actually enters its algebraic engine (observed PDP dispatch
≠ `unsupported` from the encoder cap) on every disclosed pilot point. Separately,
record whether ordinary-query attempts produce accepted relations, novel rank,
and any complete verified one-target recovery under a declared trial budget.

**Floor / reference (diagnostic only).**

- Floor: encoder variable count `m·ℓ+(m−2)·n` with `m=3`, `ℓ=6` must be `≤ 64`
  on every cell (35–49 for the five degrees). This is the static gate already
  audited; the pilot does not move it.
- Reference: the admitted optimized `pair_table` / `tiny_gauss` /
  `subgroup_orbits` path on the **same disclosed points**, run only as a
  wiring control. Paired online ratios on these reused points are **not**
  promotion evidence and must not update the scoreboard's competitive rows.

**Success.** For each algebraic arm under test, every pilot cell shows
observed engine dispatch consistent with the named solver (not encoder
`unsupported`), retains the full ordinary-query outcome mix including failures,
and records usable-base census, folded columns, matrix/rank evidence when a
solve is attempted, and null-vs-complete recovery honestly. Zero-yield cells
stay in the table.

**Stop / negative.** Encoder rejection after the static gate claimed a pass;
source/build mismatch; changing dimension, summand count, or `MAX_VARS` after
seeing outcomes; treating a single recovered scalar without checked IC descent
as success; or promoting this pilot into a family qualification.

**Inadmissible.** Fresh public targets; redispatch of seeds `2026092901` or
`2026092902`; claiming an end-to-end speedup or ECC2K-130 transfer; isolating
"F4 vs pair_table" while the factor-base policy also changed without labelling
the comparison `factor-base-policy`; filling static `unsupported` as measured
zero yield for the ambient-orbit v2 arms.

## Frozen inputs

| Item | Binding |
| --- | --- |
| Worker commit | `765c3c5f19032bd852163805f257c56babef2040` |
| Factor base | `standard_subspace`, dimension **6** only |
| Solvers | `f4`, `f5` (required); optional `inherited_f4` |
| Relation LA | `dense` (match inventory control); optional sparse follow-up only after dense rows exist |
| Summands | 3 |
| Cells | n17a1, n19a0, n23a0, n23a1, n31a0 |
| Points | Exactly the five public targets already used in the dimension-6 inventory control (one per cell). Digests and coordinates are sealed by that JSON once #952 merges; do not substitute readiness or improvement-round points without amending this protocol. |
| Budget | Frozen for the measured series: `max_trials=64`, `batch_trials=1`, child wall 300s, memory 8 GiB soft, `RAYON_NUM_THREADS=1`. Exhaustion is incomplete, not UNSAT. An earlier env-error attempt with `KIC_F5_AVX512_UNPACK` set is retained under `runs-env-error/` (exclusive-phase workers reject any `KIC_*`). A 1024-trial attempt on n17a1/f4 hit the 300s cap with no receipt; that row is retained as operational incomplete and does not count as a solver verdict. |
| Host label | Record `rustc --version`, CPU, threads. Valgrind instruction counts stay out of scope unless 3.22.0 is present. |
| Inventory points SHA-256 | `be3053b4b637311e8255f294807fde412510fbda2464e678a6a16fdea397d403` (`standard-subspace-d6-inventory-control.json`) |
| Runtime env | Tournament-style child env only: `PATH`/`HOME`/`LANG`/`LC_ALL`/`TZ`, `RAYON_NUM_THREADS=1`, `IC_ARTIFACT_CACHE=off`, `IC_F2_BACKEND=cpu`. Do **not** set `KIC_F5_AVX512_UNPACK` for this exclusive-phase pin. |

Classification of any later cost movement against the incumbent on these points
is at most **engineering** or **factor-base-policy**; it cannot be an **advance**
against the encoder floor, which is already satisfied by construction at `ℓ=6`.

## Cost accounting

Primary diagnostic outputs: observed PDP dispatch, query terminal classes,
accepted rows, novel rank, descent presence, and complete/incomplete status.
Wall time may be logged as practicality only. Instruction-count `S` and
rho-ratio speedups stay **unset** on this pilot unless a separately frozen
calibrated profile path is added in an amendment. Missing phases stay null.

## After the pilot

If at least one F4/F5 arm shows real dispatch and non-degenerate recovery on
the disclosed set, register a **new** competitive protocol (new seed) that:

1. Runs `generic_solver_feasibility.py --require-pass` on the panel before
   fixture generation.
2. Excludes every point seed `2026092902` could have generated (full v2
   exposure corpus once run 36580669479 finishes or is censored), plus all
   prior sealed/supplemental exposures.
3. Labels the comparison `factor-base-policy` when `standard_subspace` replaces
   orbit sampling.
4. Keeps the three sealed improvement confirmation/replay sets closed.

If every algebraic arm fails after clearing the encoder cap, retain the
negative pilot as evidence and do not register a fresh F4/F5 complete-solve
campaign until the mechanism changes.
