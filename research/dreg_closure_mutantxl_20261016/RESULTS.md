# Results: MutantXL closure refutes the descent systems at degree 3–4, where XL needs 5–7

Scored natively by `dreg_score closure runs/*.jsonl` (`score-output.txt`)
against `PREREGISTRATION.md`. Every value comes from a committed record in
`runs/`.

- Run: 2026-10-16 02:35–02:39 UTC. All 41 files exit 0, and no file hit
  the 3-hour limit.
- Binary: `runs/binary.sha256`.
- Run report with graph and PDF: `reports/run-20261016T0235Z/`.

## The measurement

There are 92 unsatisfiable draws, four per cell. "Committed XL" is the
refutation degree the source studies recorded. The driver re-established it
on every degree it built: the plain basis refuted exactly there and never
below. This XL-identity assertion held on all 92 draws, and every closure
saturated.

| m | cell (n, ℓ) | N | committed XL | closure degree | gap | Q1 |
|--:|---|--:|---|---|---|---|
| 3 | (5, 2) | 11 | 5 5 5 5 | 4 3 4 4 | 1 2 1 1 | below |
| 3 | (7, 2) | 13 | 5 5 5 5 | 4 4 4 3 | 1 1 1 2 | below |
| 3 | (8, 2) | 14 | 5 5 5 5 | 4 4 4 4 | 1 1 1 1 | below |
| 3 | (10, 2) | 16 | 5 5 5 5 | 4 4 4 4 | 1 1 1 1 | below |
| 3 | (12, 2) | 18 | 5 5 5 5 | 3 3 3 3 | 2 2 2 2 | below |
| 3 | (4, 3) | 13 | 5 5 5 5 | 4 4 4 4 | 1 1 1 1 | below |
| 3 | (5, 3) | 14 | 6 6 6 6 | 4 4 4 4 | 2 2 2 2 | below |
| 3 | (7, 3) | 16 | 6 6 6 6 | 4 4 4 4 | 2 2 2 2 | below |
| 3 | (9, 3) | 18 | 6 6 6 6 | 4 4 4 4 | 2 2 2 2 | below |
| 3 | (7, 4) | 19 | 7 7 6 7 | 4 4 4 4 | 3 3 2 3 | below |
| 3 | (8, 4) | 20 | 6 7 6 7 | 4 4 4 4 | 2 3 2 3 | below |
| 3 | (11, 4) | 23 | 6 6 6 6 | 4 4 4 4 | 2 2 2 2 | below |
| 3 | (13, 4) | 25 | 6 6 6 6 | 4 4 4 4 | 2 2 2 2 | below |
| 3 | (10, 5) | 25 | 7 7 7 7 | 4 4 4 4 | 3 3 3 3 | exploratory: below |
| 4 | (4, 1) | 12 | 5 5 5 5 | 4 3 4 4 | 1 2 1 1 | below |
| 4 | (5, 1) | 14 | 5 5 5 5 | 4 4 4 3 | 1 1 1 2 | below |
| 4 | (6, 1) | 16 | 5 5 5 5 | 3 4 4 3 | 2 1 1 2 | below |
| 4 | (5, 2) | 18 | 6 6 6 6 | 4 4 3 3 | 2 2 3 3 | below |
| 4 | (6, 2) | 20 | 6 6 6 6 | 4 4 3 3 | 2 2 3 3 | below |
| 4 | (7, 2) | 22 | 6 6 6 6 | 4 3 4 4 | 2 3 2 2 | below |
| 4 | (8, 2) | 24 | 6 6 6 6 | 4 3 4 4 | 2 3 2 2 | below |
| 4 | (9, 3) | 30 | ≥7 ≥7 ≥7 ≥7 | 4 4 4 4 | ≥3 | exploratory: below |
| 4 | (10, 3) | 32 | ≥7 ≥7 ≥7 ≥7 | 4 4 4 4 | ≥3 | exploratory: below |

**Q1: closure lowers the refutation degree. Holds, more strongly than
predicted.**

- All 20 testable cells are below: every draw refutes under the closure
  below its XL degree, and no cell is mixed.
- The three exploratory cells, where `C(N, ≤ XL−1)` exceeds the 70k cap,
  are below too, at 4.
- Registered confidence was about 55%.

**Q2: by one degree, not more. Refuted.**

- Of the 84 resolved-XL draws, 26 refute one degree below XL and 58 refute
  two or more below.
- The gap grows with the XL degree:
  - one degree (sometimes two) where XL is 5;
  - two (sometimes three) where XL is 6;
  - three where XL is 7.

**Q3: the growth with ℓ survives. Refuted.**

- The m = 3 closure-degree medians are 4 at every ℓ:
  - ℓ = 2: 20 draws;
  - ℓ = 3: 16 draws;
  - ℓ = 4: 16 draws;
  - ℓ = 5: 4 draws, exploratory.
- None rises. Over the same cells XL rises from 5 to 7.
- At m = 4 the closure also reads 3–4 at ℓ = 1, 2 and 3, against XL's
  5, 6, ≥7.

**What the data say.** On every measured cell, the closure refutes at a degree
that does not grow with ℓ, with m (3 or 4), with n, or with the unknowns
(11 to 32). The degree-ℓ growth that every earlier study in this line
measured is a property of plain XL, not of the ideal.

- The degree-3 refutations are concentrated in some cells. One example is
  m = 3 (12, 2), where all four draws refute at 3.
- No predictor for them was registered, and none is claimed.

## Soundness control (added after measurement, disclosed)

The pre-registration validated "closure ≤ XL" only on unsatisfiable draws.
A closure that wrongly refuted would also have produced these numbers, so a
negative control was added after the scored run, before this write-up.

- **The test.** Every committed *satisfiable* draw in the same 41 input
  files is replayed at degrees 3 to 5 (as far as the cap allows), with the
  same driver and a `--satisfiable` flag.
- **The pass condition.** Neither XL nor the closure may refute a system
  that has solutions.

**Result:** **0 of 345** satisfiable draws refuted by XL or the closure, every one built to degree 5 (one above every closure refutation in the study). Split across two invocations: the first (`runs/satisfiable-control.{jsonl,log}`) was killed by the session tool's 10-minute background limit after 315 complete draws (it ran file by file, so `(8, 4)` and everything before it were complete); the remaining 30 draws, `(7, 4)`, `(11, 4)` and `(13, 4)`, ran in `runs/satisfiable-control-2.{jsonl,log}` (exit 0, peak 12.3 GB at `(13, 4)`). The kill is a tool limit, not evidence.

- Records: `runs/satisfiable-control.{jsonl,log}` and `runs/satisfiable-control-2.{jsonl,log}`; tallied in `score-output.txt`.
- Binary SHA-256: `dcffc6b685b54162b01e517c1174924ac916adaf56f2878c714bd3f0037cc094`.
  It is the frozen binary's source plus the control flag. The closure code
  is identical.

## Cost, descriptive

Cell medians, log₂. The columns are `C(N, ≤D)`, the width of the
degree-`D` matrix. The word ops are the counted 64-bit word operations
(`InheritCost`) of the whole closure build at the refuting degree: base
elimination plus every closure round. The XL matrix at its own degree was
not built at the larger cells, so it has no word-op count. Its width is
exact.

| m | cell | N | XL degree | cols at XL | cols at closure | closure products | closure word ops |
|--:|---|--:|--:|--:|--:|--:|--:|
| 3 | (5, 2) | 11 | 5 | 10.0 | 9.1 | 12.0 | 19.0 |
| 3 | (12, 2) | 18 | 5 | 13.6 | 9.9 | 12.3 | 20.5 |
| 3 | (9, 3) | 18 | 6 | 14.9 | 12.0 | 14.6 | 26.3 |
| 3 | (7, 4) | 19 | 7 | 16.5 | 12.3 | 14.9 | 26.4 |
| 3 | (13, 4) | 25 | 6 | 17.9 | 13.9 | 16.6 | 31.0 |
| 3 | (10, 5) | 25 | 7 | 19.5 | 13.9 | 16.5 | 30.1 |
| 4 | (8, 2) | 24 | 6 | 17.5 | 13.7 | 16.5 | 28.8 |
| 4 | (10, 3) | 32 | ≥7 | ≥22.1 | 15.3 | 17.8 | 33.8 |

The full table is in `score-output.txt`.

- **Matrix width.** The closure's matrix is 2^1 to 2^7 narrower than XL's.
- **Closure work.**
  - Closure products and word ops grow roughly like
    columns^1.2 and columns^2.2 respectively. For example, (10, 3) at
    m = 4 has 2^15.3 columns, 2^17.8 products and 2^33.8 word ops.
  - So the closure's work is close to the `w = 2` pricing of its own
    matrix width.
  - For comparison, a dense elimination of XL's wider matrix would be about
    `cols³/64`. That is an estimate, not a measurement.
- **Wall time.** The slowest single degree took 39 s. Peak memory was
  1.0 GB on the 4-core, 16 GB, non-isolated container. These are
  practicality notes, not claims (`AGENTS.md` §10).

## What this means for ECC2K-130 (derived, extrapolated)

`examples/descent_degree_cost` now prices a "4, constant: the MutantXL
closure degree" row. Its output is
`research/notes/ecc2k130/descent_degree_cost-output-20261016.txt`, and every
earlier row reproduces unchanged.

| m | degree | w | ℓ | N | total | vs rho |
|--:|---|--:|--:|--:|--:|--:|
| 3 | XL fit `⌈ℓ/2⌉+4` | 1 | 4 | 143 | 2^158.96 | +98.15 |
| 3 | **4, closure** | 1 | 40 | 251 | **2^82.27** | **+21.46** |
| 3 | 4, closure | 2 | 45 | 266 | 2^100.25 | +39.44 |
| 4 | XL fit `ℓ+4` | 1 | 4 | 278 | 2^173.13 | +112.32 |
| 4 | **4, closure** | 1 | 33 | 394 | **2^68.43** | **+7.62** |
| 4 | 4, closure | 2 | 34 | 398 | 2^93.91 | +33.10 |

- **The model input.** If the closure degree stays at 4 to `n = 131`, the
  best descent row moves from 2^173.13 to 2^68.43 (`m = 4`, `w = 1`).
  Against rho's 2^60.81, the remaining distance falls from 112 bits to
  7.6 bits.
  - At `w = 2`, the pricing the measured closure work tracks, it falls to
    33.1 bits.
  - Units are as in the descent note: column touches against group
    operations, unconverted.
- **This is an extrapolation, not a measurement.** The closure degree is
  measured only to `N = 32` and `ℓ = 5`. The priced rows sit at
  `N ≈ 251` to `398` and `ℓ ≈ 33` to `45`.
  - A flat degree over 11 to 32 unknowns is evidence. It is not a proof
    that the degree stays flat over 400.
  - The m = 4 free-oracle floor is still 2^56.46 per relation set. It
    leaves an oracle budget of 2^5.71 per target, below any degree-4
    matrix at `N ≈ 370`.
  - So **no row beats rho**, and the degree-independent verdict of the
    descent note §2 stands.
- **Against the success criterion.** This is a measured improvement of the
  instrument the IC cost model uses for the descent route: a lower
  refutation degree on all 92 committed draws, with a narrower matrix.
  - It is not a measured ECC2K-130 cost.
  - It is not a speedup claim against rho.
  - The `n = 131` rows are derived.
- **§5a pricing.** This study ranks no two arms.
  - The descent rows' units are column touches, and rho's are group
    operations.
  - The touches-per-group-operation factor `K` is unmeasured
    (descent note §2).
  - So a 0 / 0.1 / 1 GAE band cannot order them until `K` is measured. The
    `+7.62` row is a distance in unconverted units, not a separation in
    either direction.

## Successor lead

The open question with the largest effect on the derived rows is whether
the closure degree stays at 4 as `N` grows.

- **Test.** Measure the closure at larger unknown counts with a sparse
  closure (the dense basis caps `N` at about 32 for degree 4).
- **Cells first.** m = 3 (13, 5), (15, 5) and m = 4 (12, 3), at degree 4
  only.
- **What it would show.** A single degree-4 refutation at `N ≈ 40` to 50
  would extend the flat range. A draw that needs degree 5 would bound it.
