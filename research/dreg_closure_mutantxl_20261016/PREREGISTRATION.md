# Pre-registration: does MutantXL closure refute the descent systems below XL?

Registered 2026-10-16, before the closure was measured on any draw below.

## Why

Every refutation degree this repository records for the chained `S₃`
decomposition systems is **XL's**. Studies 2026-09-23 to 2026-10-06 define it
as the first degree `D` whose plain Macaulay row space contains `1`. The
ECC2K-130 cost model uses those degrees
(`research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md` §2 and §5), and
its notes say the degree is "neither an F4/F5 solving degree nor a
certified first fall degree".

An F4- or MutantXL-type solver does more at a fixed degree. It multiplies the
row space's **degree falls** (elements of degree `< D`) by monomials that keep
the degree `≤ D`, and repeats until saturated. Every product lies in the
ideal, so this **closure** refutes at a degree never above XL's, and possibly
below it. The first fall degree of these systems is 3, so falls are present
from degree 3. The coverage survey of 2026-10-16 found no XL-vs-MutantXL
experiment in the repository.

If the closure refutes lower, every `n = 131` cost row built on the XL degree
overstates what an F4-type oracle needs. Cost grows roughly as `C(N, ≤D)` in
the unknowns `N`, so one degree is worth many bits.

## Instrument

`closure_refutation_scan` (`src/cryptanalysis/koblitz_bench.rs`) builds each
degree `D` twice with `ReducedBasis` (`src/cryptanalysis/inherited_f4.rs`):

- plain, which is XL;
- `from_system_closed_full`: every fall row of the reduced row echelon basis,
  linear falls included, times every monomial `t` with `deg t ≤ D − deg(fall)`
  over the occurring variables, absorbed, repeated to saturation.

`examples/dreg_closure` replays the committed draws (same subspace basis, same
`x(R)`, same field polynomial), so XL and closure see identical systems.

Validation, before any registered draw:

- **Positive control.** f₁ = x₁x₂x₃ + x₄x₅x₆, f₂ = f₁ + x₇ + 1,
  f₃ = x₇x₈x₉ + x₈x₉ + 1. The closure refutes at 3, XL at 4, and the
  repository's XL engine `solving_degree` agrees on 4.
  (`closure_scan_refutes_below_xl_on_the_constructed_control`.)
  - My hand derivation said XL = 5. The test showed 4, because `x₇·f₃ = x₇`
    is already a degree-4 row. I corrected the expectation, not the code.
- **Never above XL; saturation; XL identity.** On six unsatisfiable chained
  draws, the closure is never above XL, every degree saturates, and the plain
  basis refutes where `solving_degree` does
  (`closure_scan_is_never_above_xl_on_chained_draws`).
- **XL identity on every registered draw.** At every degree it builds, the
  driver asserts that the plain basis refutes exactly at the committed XL
  degree and never below it. A violation aborts the run, and that run is
  reported as such.
- **Existing behaviour.** The splitting solver's closure is unchanged:
  falls of degree ≥ 2, seeded as before. All `inherited_f4` and
  `koblitz_bench` tests pass (44).

## Cells and draws (frozen)

All unsatisfiable draws in these committed files, nothing else. The
"committed XL" column comes from `dreg_score score` on the same files.

| m | cell | N | committed XL | files |
|--:|---|--:|---|---|
| 3 | (5,2) (7,2) (8,2) (10,2) (12,2) | 11–18 | 5 | `dreg_fixed_surplus_20260923/runs/cell-{5-2,7-2}-7.jsonl`, `dreg_ell_grid_20260925/runs/cell-{8-2,10-2,12-2}-7.jsonl` |
| 3 | (4,3) | 13 | 5 | `dreg_ell_grid_20260925/runs/cell-4-3-7.jsonl` |
| 3 | (5,3) (7,3) (9,3) | 14–18 | 6 | `dreg_ell_grid_20260925/runs/cell-5-3-7.jsonl`, `dreg_fixed_surplus_20260923/runs/cell-{7-3,9-3}-7.jsonl` |
| 3 | (7,4) (8,4) | 19, 20 | 7 7 6 7; 6 7 6 7 | `dreg_surplus_control_20260925/runs/cell-{7-4,8-4}-7.jsonl` |
| 3 | (11,4) (13,4) | 23, 25 | 6 | `dreg_fixed_surplus_20260923/runs/cell-{11-4,13-4}-6.jsonl` |
| 3 | (10,5) | 25 | 7 | `dreg_degree7_ell5_20260930/runs/cell-10-5-7.u{0..3}.jsonl` |
| 4 | (4,1) (5,1) (6,1) | 12–16 | 5 | `dreg_m4_chain_20261006/runs/cell-{4,5,6}-1-7.jsonl` |
| 4 | (5,2) (6,2) (7,2) (8,2) | 18–24 | 6 | `dreg_m4_chain_20261006/runs/cell-5-2-8.jsonl`, `cell-{6-2-8,7-2-7,8-2-7}.u{0..3}.jsonl` |
| 4 | (9,3) (10,3) | 30, 32 | ≥7 | `dreg_m4_chain_20261006/runs/cell-{9,10}-3-6.u{0..3}.jsonl` |

**Scan rule.**

- Degrees run from 3 up to the committed XL degree (for a `≥D` bound, up to
  `D − 1`). The scan stops once the closure refutes.
- A degree is built only if `C(N, ≤D) ≤ 70,000` columns. The dense basis costs
  about the cube of that.
- So these cells are measured only up to the stated degree:
  - m = 3 (10,5): up to 5, two below XL;
  - m = 4 (9,3), (10,3): up to 4.
- These are **exploratory**: they can show the closure refuting ≥ 2 below XL,
  but they cannot test "one below".
- Excluded as untestable under the cap: m = 3 (11,5), (13,5), (15,5), where
  `C(N, ≤5)` > 83k.

## Predictions

- **Q1 — closure lowers the refutation degree.** In each testable cell every
  draw refutes under the closure at least one degree below its committed XL
  degree.
  - A cell is testable when degree XL − 1 fits the cap.
  - Confidence: moderate, about 55%.
  - Reason: falls exist from degree 3 (FFD 3), and the inherited-F4 notes
    report that closure at degree 3 refutes or pins branches the plain
    degree-3 matrix must split.
- **Q2 — by one degree, not more**, wherever Q1 holds. Confidence: low.
- **Q3 — the growth with ℓ survives.** At m = 3, the closure degree is
  non-decreasing in ℓ (medians over the ℓ = 2, 3, 4 cells), and some ℓ step
  rises. Confidence: moderate to high.

## Verdict rules

Per draw:

- **closure degree** = the first scanned degree whose closure refutes;
- **"not below XL"** if the closure has not refuted by degree XL − 1;
- **unmeasured** if XL − 1 exceeds the cap or the run was killed.

**Q1**, per testable cell:

- **below** if every draw refutes below XL;
- **not below** if no draw does;
- **mixed** otherwise.

The study verdict is **"closure lowers"** if at least 2/3 of the testable cells
are below. It is **"no effect"** if at least 2/3 are not below. Otherwise it is
**mixed**.

**Q2:** among cells that are below, the per-draw gap XL − closure, reported
as counts of gaps of 1 and of 2 or more.

**Q3:** per-ℓ medians of the m = 3 closure degree. Q3 holds if they are
non-decreasing and at least one step strictly increases.

## Cost accounting, descriptive

For each refuting degree the record keeps:

- columns;
- XL rank and closure rank;
- closure products;
- the word operations `InheritCost` charged.

The report compares `C(N, ≤D)` at the closure degree with `C(N, ≤D)` at the
XL degree. It makes no wall-clock claim. Wall times are practicality notes for
the four-core, 16 GB container, which is not isolated (`AGENTS.md` §10).

## Stop conditions

- A per-file wall limit of 3 hours (`run.sh`). A killed file is recorded as
  killed and its draws as unmeasured. A kill is a resource limit, not
  evidence.
- An XL-identity assertion failure aborts that file. It is reported as an
  instrument failure, and the study's affected cells are void.

## Scope and class

- **Systems.** The chained `S₃` decomposition systems (`b = 1`) of the
  `E_0` family's descent studies, random subspaces, `n ≤ 13`. Field subfields
  are as disclosed in the source studies (`AGENTS.md` §8b). None is chosen.
- **Class: stage diagnostic.** No `S`, no rho ratio, no speed claim.
  - A lower closure degree changes the input of the `n = 131` cost model. It
    is not a measured ECC2K-130 cost.
  - The verdict against rho in the descent note is degree-independent
    (`RESEARCH_ECC2K130_DESCENT_DEGREE.md` §2): the m = 4 free-oracle floor
    leaves an oracle budget of `2^5.71` per target, below any Macaulay matrix.
    A closure result cannot reverse it, and the report will say so.
- **Conductor.** `CLAUDE.md` asks for a Conductor task before code changes.
  The `conductor` CLI is not installed in this container and its hooks fail
  open (`7e7464881`), so no reservation was taken. This is disclosed, not
  hidden.
