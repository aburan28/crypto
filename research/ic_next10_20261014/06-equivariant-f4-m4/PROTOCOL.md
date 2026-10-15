# Protocol: item 6, the m = 4 Weil-descent degree cells and equivariant F4

Frozen 2026-10-14, before any main-arm data. Item "### 6." of
`reports/Index calculus leads next steps.md`.

## 0. What was already run (check-first, 5 minutes)

- `research/dreg_m4_chain_20261006/`: the `m = 4` chained system's **Macaulay
  refutation degree** (bounded Macaulay/XL linear algebra, `dreg_ladder`):
  `ℓ = 1` → 5, `ℓ = 2` → 6, `ℓ = 3` → ≥7 (degree 6 built in full at `(9,3)` and
  `(10,3)`). No `ℓ ≥ 4` cell was ever run.
- `research/dreg_degree7_ell5_20260930/`: `m = 3`, `(10,5)` refutes at exactly
  7; `(15,5)` is `≥7` (degree 6 built, `research/dreg_fixed_surplus_20260923`).
- `research/ic_m4_exponent_audit_20260928/`: `m = 4` `ℓ = 3` censored at 300 s.
- `RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md`: G2 measured at **rank level only**
  (n ≤ 23). No `ict-g2-equivariant*` branch or worktree exists
  (`git branch -a`, `git worktree list` show only `feat/brave-bardeen-9jxoaj`).
- No file measures an **F4** (as opposed to Macaulay/XL) solving degree on the
  chained system. Every committed degree is the XL refutation degree.

### Pilot (disclosed, engineering only)

Before this file was written, three msolve runs were made to check the tool
works on GF(2): the committed `m = 4` draw `(6,2)` u0 (XL degree 6) refuted with
F4 max step degree 4 in 0.06 s; the committed `(9,3)` u0 (XL `≥7`) refuted at F4
degree 4 in 2.7 s; one `(n, ℓ) = (11, 4)` draw at seed 20261014 refuted at F4
degree 4 in 133 s, peak under the 3 GB cap. These three pilots are re-run
inside the arms below and are **not** counted twice. They motivated adding the
sensitivity control D (an instrument that always reads 4 must be shown able
to read more).

## 1. Hypotheses

- **H1 (the item's question).** At `m = 4` the solving degree keeps growing
  with `ℓ` (Kill direction: ≥8 at `ℓ = 4`) or levels off.
- **H2 (instrument).** The F4 solving degree `D_F4` (msolve) is at most the XL
  refutation degree `D_XL` on every draw where both are measured (theory: any F4
  computation of max degree `D` lies in the degree-`D` ideal truncation, which
  contains the XL row space but can be larger through degree falls).
  Prediction: `D_F4 < D_XL` on most draws.

Predictions (written before data): `D_F4 = 4` at `ℓ = 4`, `n = 11` (pilot);
`D_F4 ∈ {4, 5}` at `ℓ = 5`; `D_F4 ≤ 5` at m = 3 `(15,5)`. I.e. levels off.

## 2. Instruments

- **A. F4 (primary).** `msolve` (the PATH binary), `-t 1 -g 2 -v 2`, DRL,
  characteristic 2, input = the chained system exported by
  `examples/ic_next10_06_chain_export.rs` plus the `N` field equations
  `x_i^2 + x_i`. Metric **`D_F4`** = maximum `deg` over all F4 rounds in the
  verbose log (the largest total degree of any matrix row), on a draw whose
  exact solution count is 0 and whose reduced GB is `[1]`. Recorded too: per
  round degrees, largest matrix, wall time, peak RSS.
- **B. XL (repo's frozen instrument).** `dreg_ladder` (Macaulay refutation
  degree), used only as the positive control C1.
- Draw identity: the exporter replays `dreg_ladder`'s cell seed and sampling
  and its exact count, so `--unsat-index K` is the committed draw `K`. Identity
  is checked by matching `v_basis` and `x_r` against the committed JSONL.

## 3. Arms and sizes (four unsatisfiable draws per cell unless stated)

- **C1 positive control, XL (reproduce known cells).** `dreg_ladder --m 4`,
  `(6,2)` u0–u3, expect `6 6 6 6` (committed). And `(9,3)` u0 to `d_max = 6`,
  expect `≥7` (committed; ~386 s, ~2.85 GB) if it fits the 3 GB cap.
- **C2 calibration, F4 on committed draws.** `D_F4` on the committed `m = 4`
  draws `(4,1) (5,1) (6,1) (5,2) (6,2) (7,2) (8,2) (9,3) (10,3)` (seed 20261006)
  and `m = 3` draws `(10,5)` (seed 20260925, XL 7), `(15,5)` (seed 20260928
  or the one that matches the committed basis, XL `≥7`). These give the paired
  `(D_XL, D_F4)` table.
- **Main M.** `m = 4`, seed **20261014**: `(11,4)`, `(13,4)`, `(13,5)`, and
  `ℓ = 6` at the smallest `n` where unsatisfiable draws occur within 400 draws
  and the run fits 3 GB / 45 min (try `n = 15, 16, 17`). The `m = 3` `(15,5)`
  cell is C2's `(15,5)`.
- **D sensitivity control (can the instrument read > 4?).** Random cubic
  Boolean systems with the same `N`, equation count and terms per equation as
  the `m = 4` `(6,2)` and `(9,3)` cells (4 systems each, Python generator,
  seeded). Generic systems of this shape have semi-regular `D_reg` 7 and 8;
  prediction `D_F4 ≥ 5`. If D reads 4 too, the instrument is declared blind
  and H1 is INCONCLUSIVE.
- **E instrument correctness.** On 2 satisfiable draws per `m = 4` main cell
  the reduced GB must not be `[1]`; on every unsatisfiable draw it must be
  `[1]`. Any mismatch voids the arm.
- **FGHR arm: not run.** FGHR 2-torsion symmetrisation acts on the symmetric
  `S_{m+1}` formulation (and in this repo, `fghr_full.rs`, on `F_{p^3}`
  curves); the cells are defined on the chained `S₃` system, which has no such
  symmetry to quotient. Recorded as a deviation from the design, not a result.

## 4. Metric and decision rules

- Cell reading = median of the four `D_F4` values.
- **Levels off** if the `m = 4` readings at `ℓ = 4` and `ℓ = 5` (and `ℓ = 6`
  if measured) equal the `ℓ = 3` reading. **Grows** if any `ℓ ≥ 4` reading is
  above the `ℓ = 3` reading. A cell that cannot be run is BLOCKED with its wall.
- **Kill criterion (report, applied to `D_F4`).** `≥8` at `(4,4)` together
  with `≥8` at `(15,5)` → growth confirmed, ECC2K-130 IC closed. `≤7` at both
  `(4,4)` and `(4,5)` → the `D = 6` row reopens, after which `K` and G2's
  constant decide. Because the report wrote its criterion in XL degrees, both
  instruments' readings are reported and the verdict says which it rests on.

## 5. Equivariant F4 (priority 2)

Extend `g2_wide_probe`'s block decomposition to a full linear-algebra pass:
RREF each character block over `F_{2^d}`, inverse-transform the reduced rows
to `F₂` rows, and check that their span equals the full `F₂` Macaulay row
space (rank and row-space equality), at `l = 3`, degree 3, `n = 11 … 23`, and
`n = 31` if it fits. Metric: agreement (pass/fail) and wall time of
blocks + inverse transform vs. full `F₂` elimination, same process.

## 6. Budget

Single thread, `nice`, `prlimit --as=3000000000`, 45 min per job.
