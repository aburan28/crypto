# Split and hybrid point decomposition against matched rho (stage diagnostic)

**Question.** Does cutting one point decomposition into smaller
summation-polynomial pieces ("fanning out" `R = P₁ + … + P_m` through
auxiliary abscissae) lower the solving degree, and does it make the
decomposition cheaper? Does guessing `g` Boolean unknowns? What does each
variant cost against matched rho on the same curve?

**Classification.** This is a stage diagnostic: one PDP solve, plus a
relation-phase projection priced at a pinned ratio. No end-to-end IC run,
speedup or `vs_rho` claim follows from it.

- **Inputs:** [`PROTOCOL.md`](PROTOCOL.md) (preregistered before the grid)
  and [`AMENDMENT_1.md`](AMENDMENT_1.md) (solution read-out defect, rerun).
- **Analysis:** [`ANALYSIS.json`](ANALYSIS.json), with its table in
  [`ANALYSIS.md`](ANALYSIS.md), written by the native `pdp_split_analyze`.
- **Report PDF:** [`REPORT.pdf`](REPORT.pdf).
- **Run reports:** under [`reports/`](reports/).
- **Registration:** [`cairn/`](cairn/README.md).

## What was measured

**Systems.** Weil-descended summation-polynomial systems on Koblitz curves
(`b = 1`), with factor base `V = span(1, z, …, z^{ℓ−1})`. The native
implementation is [`src/cryptanalysis/pdp_split.rs`](../../src/cryptanalysis/pdp_split.rs).

| plan | system | unknowns |
|---|---|---|
| `unsplit` (m = 3) | `S₄(x₁,x₂,x₃,x_R)` | `3ℓ` |
| `split-2+1` (m = 3, the S₃ chain) | `S₃(x₁,x₂,u)`, `S₃(u,x₃,x_R)` | `3ℓ + n` |
| `split-2+2` (m = 4, two-way) | `S₃(x₁,x₂,u)`, `S₄(u,x₃,x₄,x_R)` | `4ℓ + n` |
| `split-2+1+1` (m = 4, three-way) | S₃ chain with two cuts | `4ℓ + 2n` |

**Hybrid and solver.**

- Hybrid `g ∈ {0, 2, 4}`: every one of the `2^g` branches is solved and
  charged.
- Solver: the repository's F4 over `F₂[v]/(v² − v)`.
- Solving degree means F4's (the highest step degree at which the basis
  grew). Cost means F4's deterministic 64-bit word XORs.

**Grid.** 8 cells × 4 seeded targets × the plans × 3 guess counts gives
192 runs. Every run's verified decomposition set was checked against an
exhaustive enumeration.

**Matched rho.** ecbench session `ECBS1he23809badf5b` covers the same four
curves: strong and plain signed-Frobenius rho, 8 public targets × 5 rounds,
L2 on most runs. The audit passed with `--replay-all` (receipt SHA-256
`0b56a9b7…cc0c`).

## Correctness and replication

- **First grid (`runs/`): one cell failed gate 2.** In `3:0:13:5`, the S₃
  chain at `g = 0` found none of the decompositions on 3 of 4 targets. The
  cause was the solution read-out, not F4: the library's
  `solutions_from_reduced_basis` silently returned no points on a basis
  whose linear elements are not inter-reduced. See
  [`AMENDMENT_1.md`](AMENDMENT_1.md).
- **Fix.** It was fixed in `pdp_split::solutions_from_basis`, with a
  regression test. The grid was rerun into `runs-v2/`.
- **Amended grid (`runs-v2/`):** all 192 runs are complete and exact. No
  timeouts, no false solutions, and no decomposition missed.
- **Replication.** All 192 runs of the amended grid were compared with the
  first grid. Solving degree and word XORs are identical on every one. The
  only 3 differences are the verified counts of the three amended runs
  (0 → 1, 3, 1).
- **cairn.** The `3:0:13:4` cell and the rho audit were re-executed through
  cairn replay objectives and accepted. This is the same host replaying
  itself, not an independent replication.

## Results (amended grid; mean over 4 targets)

F4 word XORs per target, all branches summed. Solving degrees are in
`ANALYSIS.md` and Figure 1.

| cell `m:a:n:ℓ` | plan | g = 0 | g = 2 | g = 4 | solving degree g = 0 → 4 |
|---|---|--:|--:|--:|---|
| 3:1:11:4 | unsplit | 3.08e7 | 6.74e6 | **4.05e5** | 7 → 5 |
| | S₃ chain | 1.19e8 | 2.23e8 | 1.34e7 | 4 → 3 |
| 3:0:13:4 | unsplit | 4.94e7 | 5.28e6 | **6.31e5** | 7 → 5 |
| | S₃ chain | 3.95e8 | 9.45e8 | 2.49e7 | 4 → 3 |
| 3:1:17:4 | unsplit | 3.95e7 | 4.66e6 | **1.64e5** | 7 → 4 |
| | S₃ chain | 1.12e9 | 2.90e9 | 7.11e7 | 4 → 3 |
| 3:1:19:4 | unsplit | 4.99e7 | 2.54e6 | **1.78e5** | 7 → 4 |
| | S₃ chain | 3.44e9 | 9.19e9 | 1.62e8 | 4 → 3 |
| 3:0:13:5 | unsplit | 1.93e10 | 1.45e9 | **6.48e7** | 7 → 6 |
| | S₃ chain | 9.28e8 | 3.66e9 | 1.66e9 | 4 → 4 |
| 3:1:17:5 | unsplit | 2.43e10 | 3.51e8 | **1.03e8** | 7 → 6 |
| | S₃ chain | 3.68e9 | 1.25e10 | 7.65e9 | 4 → 4 |
| 4:1:11:3 | two-way 2+2 | 2.80e10 | 9.65e9 | 1.65e8 | 7 → 6 |
| | three-way 2+1+1 | 2.49e8 | 2.18e9 | **3.32e7** | 4 → 3–4 |
| 4:0:13:3 | two-way 2+2 | 2.95e10 | 1.93e10 | **1.13e8** | 7 → 6 |
| | three-way 2+1+1 | 3.82e10 | 1.00e10 | 4.79e8 | 4 → 4 |

Bold marks the cheapest variant measured in each cell. The unsplit `m = 4`
system needs a symbolic binary `S₅`, which the repository lacks, so it was
not measured.

### Against matched rho (projection, relation phase only)

`S_relations = (F + 1)/λ · X̄ · ratio / √r`. Here `λ` is the
decompositions per target from 256 exhaustively counted targets, and
`ratio` is `docs/ic/calibration.json`'s pinned `ns_per_word_xor` for the
curve, calibrated on another host. Linear algebra, descent and
factor-base costs are omitted, so this is a lower bound on full IC.

| cell | cheapest variant | projected `S_relations` | strong rho `S` (measured) | ÷ strong | signed-Frobenius rho `S` (measured) | ÷ signed-Frob. |
|---|---|--:|--:|--:|--:|--:|
| 3:1:11:4 | unsplit, g = 4 | 956 | 42.0 | 22.8 | 2.45 | 391 |
| 3:0:13:4 | unsplit, g = 4 | 2,286 | 16.6 | 137 | 2.09 | 1,091 |
| 3:1:19:4 | unsplit, g = 4 | 1,266 | 2.38 | 532 | 0.577 | 2,195 |
| 3:0:13:5 | unsplit, g = 4 | 78,065 | 16.6 | 4,689 | 2.09 | 37,267 |
| 3:1:17:5 | unsplit, g = 4 | 357,347 | 4.05 | 88,144 | 1.09 | 326,874 |
| 4:0:13:3 | two-way, g = 4 | 494,714 | 16.6 | 29,713 | 2.09 | 236,169 |
| 3:1:17:4, 4:1:11:3 | — | no projection: `λ = 0` in 256 targets | | | | |

Strong rho's `S` at these sizes is dominated by its 32-lane setup. That is
why the cheaper `rho.signed_frobenius` is shown beside it.

## Predictions, scored

| | prediction | outcome |
|---|---|---|
| P1 | the chain has a lower solving degree than unsplit at every `m = 3` cell | **held**: 4 vs 7 at `g = 0` in all six cells |
| P2 | the chain costs more XORs than unsplit at `g = 0`, with the ratio growing in `n` at `ℓ = 4` | **held at ℓ = 4** (3.9×, 8.0×, 28×, 69× at n = 11, 13, 17, 19); **falsified at ℓ = 5**, where the chain is 20.8× (n = 13) and 6.6× (n = 17) *cheaper* at `g = 0` |
| P3 | three-way has a lower solving degree than two-way at `m = 4` | **held**: 4 vs 7 in both cells. Cost goes both ways: three-way is 113× cheaper at n = 11 and 1.3× dearer at n = 13 (`g = 0`) |
| P4 | guessing lowers the per-branch degree; unsplit is cheaper at `g ≥ 2`; the chain is dearer at `g = 4` than `g = 0` | **unsplit held** everywhere. The **chain part held only at ℓ = 5**: at ℓ = 4, `g = 4` is cheaper than `g = 0` for the chain. The chain's degree drops only at ℓ = 4, `g = 4` |
| P5 | every complete run passes the exhaustive check | **failed in the first grid** (read-out defect, Amendment 1); **held in the rerun** |
| P6 | no variant's projected `S` is below strong rho's on the same curve | **held**: the closest is 22.8× strong rho, or 391× plain signed-Frobenius rho |

## Reading (interpretation, labelled as such)

- **Splitting traded degree for unknowns.** The S₃ chain always reached a
  lower solving degree than the unsplit `S₄` (4 against 7). It paid with `n`
  extra unknowns per cut. At `ℓ = 4` that cost more, and more so as `n` grew.
- **At ℓ = 5 the unsplit system's degree-7 matrices became the larger cost**,
  and the split won at `g = 0`. That crossover is the measurable version of
  "splitting can help". It is a constant-factor effect on these sizes, not an
  exponent.
- **Guessing beat splitting in every cell.** In every `m = 3` cell the
  cheapest variant was the unsplit system with `g = 4` guessed unknowns.
  Guessing removes unknowns, while splitting adds them. That matches the
  original argument: the auxiliary abscissa `u` carries all `n` bits of
  coupling.
- **No variant came within an order of magnitude of rho on these curves.**
  The projection omits linear algebra and descent, so full IC would be
  dearer still. Nothing here bears on deployed curve sizes or ECC2K-130.

## Limits

- **Small curves.** `r ≤ 2^18`, four targets per cell. Four of the eight
  cells have a decomposition yield too low to project from 256 targets, or
  a yield of zero.
- **Pinned pricing.** The XOR→addition ratio was calibrated on another
  host, so the projection is labelled a projection throughout.
- **Wall time is descriptive.** The PDP runs had no isolation (cloud
  container, 4 vCPU). Only counts are compared.
- **Not measured.** The unsplit `m = 4` system (no binary `S₅`), and
  Frobenius-invariant factor bases. The standard subspace was used, as in
  ecbench's `koblitz-standard-subspace`.
- **Library defect.** The defect in
  `pq_f4_f2::solutions_from_reduced_basis` is worked around here, not fixed
  (see the PR).

## Reproduce

```sh
cargo build --release --bin pdp_split_dreg --bin pdp_split_analyze --bin ecbench
./research/pdp_split_dreg_20261011/run_grid.sh          # writes runs-v2/ (must not exist)
target/release/ecbench run --wait --spec research/pdp_split_dreg_20261011/SPEC-rho.json --out NEW_DIR
target/release/pdp_split_analyze --runs research/pdp_split_dreg_20261011/runs-v2 \
  --rho-session research/pdp_split_dreg_20261011/sessions/rho-01 \
  --calibration docs/ic/calibration.json --replica research/pdp_split_dreg_20261011/runs \
  --out ANALYSIS.json
```

## Dashboard panels checked

Every panel listed below was checked and left unchanged:

- `docs/index-calculus-scoreboard.html` (`#lab-progress`, `#lab-best`,
  `#lab-ecbench-all`), `docs/ic/progress-timeline.json` and
  `docs/ic-leaderboard.html`. This study adds no measured IC/rho ratio and
  no cross-method ecbench session (its session is rho against rho). Its IC
  figures are relation-phase projections at a pinned ratio, which the
  dashboard rules keep off the measured series.
- `docs/curves/` and `figures/`: no affected graph.
