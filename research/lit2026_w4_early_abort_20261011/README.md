# lit2026 W4: early abort at the refutation degree (2026-10-11)

Results report: [`REPORT.pdf`](REPORT.pdf) (source `REPORT.typ`).
Run report: [`reports/run-20261011T2025Z/report.pdf`](reports/run-20261011T2025Z/report.pdf).
Raw probe output is in `reports/run-20261011T2025Z/raw/`, the Rust-derived summary in
`summary.json`, and the graphs `degrees.svg` and `arms.svg`, drawn from it by `plot.py`.

**Outcome.**

- Mode 1 is `no-change`. Mode 2 is `negative`.
- Both switches stay off by default. The trace stays as a diagnostic.

## Question

Horak–Semaev, arXiv [2603.24890](https://arxiv.org/abs/2603.24890), study the
inconsistency probability of sparse equations over F2. The program read this as
follows: most decomposition attempts have no solution, so proving 1 ∈ I dominates,
and F4/F5 should stop at the degree where 1 appears.

W4 asks two questions:

- At which degree does the production decomposition oracle refute, and at which does
  it decide?
- Does an exact early abort save work?

## Production path

`SolverEngine::InheritedF4 { max_degree: 3 }` is the default decomposition oracle,
used by `groebner_decompose` and, with the geometry gate, by F6-IC
(`docs/ic/GROEBNER_F4F5F6_AUDIT_20261010.md`).

- It is a splitting tree. Each node does linear elimination, then reads a
  support-local degree-3 Macaulay basis, either inherited from its parent or built
  fresh.
- On these descended systems (24 or more variables) the ladder is degree 3 alone.
- No F5 engine is on this path.
- A node whose basis holds 1 already stops at once. An inherited child stops
  specialising at the refuting row.

## Change (behind switches, default off)

All code is in `src/cryptanalysis/koblitz_groebner.rs`, in `reduce_inherited` and
`solve_rec`. `inherited_f4.rs` is untouched, so W5 is not affected.

- **`RefutationDegreeTrace`**: thread-local, enabled with
  `set_refutation_degree_trace`. It records:
  - for every node, the first from-scratch degree, 2 or 3, holding 1 or a forced
    variable;
  - refutations before any matrix;
  - stage-unit word operations spent by nodes that end refuted against nodes that
    end in a split;
  - the root system's first refutation degree, from 2 to 4.
- **`KIC_F4_EARLY_ABORT`**:
  - `1`: a degree-2 probe before a fresh degree-3 build. It preserves the tree.
  - `2`: a degree-2 probe before a child specialises its parent's bases. It is
    exact.
  - `3`: both.
  - `set_early_abort_mode` overrides the variable per thread.
- **Tests**: `early_abort_tests` checks that roots, planted roots and tree counters
  are identical across modes 0–3 on 12 wide sparse systems, plus the trace
  accounting.
- **Probe**: `examples/f4_refutation_degree_probe.rs` writes per-attempt traces, and
  its `--summarize` mode tabulates a run directory.

## Commands

The `CARGO_TARGET_DIR` values were:

- candidate: `/Volumes/SSD990-2/crypto/targets/lit2026-w4`
- original (worktree `/Volumes/SSD990-2/crypto/wt/lit2026-baseline` at `4b6e09b78`):
  `/Volumes/SSD990-2/crypto/targets/lit2026-baseline`

```bash
cargo build --release --example f4_refutation_degree_probe --example f6_ic_probe
# traces (k = 0, 1, 2)
KIC_F4_EARLY_ABORT=$k f4_refutation_degree_probe --a 0 --n 23 --dim 8 --targets 4 --engines f4,f6p
KIC_F4_EARLY_ABORT=$k f4_refutation_degree_probe --a 0 --n 31 --dim 9 --targets 1 --budget 2000 --engines f4,f6p
KIC_F4_EARLY_ABORT=$k f4_refutation_degree_probe --a 1 --n 17 --dim 6 --targets 8 --engines f4,f6p
KIC_F4_EARLY_ABORT=$k f4_refutation_degree_probe --a 1 --n 23 --dim 7 --targets 4 --engines f4,f6p
# timing: arms orig (4b6e09b78 binary), prev (=0), cand1 (=1), cand2 (=2); 5 rotated rounds
KIC_F4_EARLY_ABORT=$k f6_ic_probe <same cell flags> --engines f4,f6p --reps 3
f4_refutation_degree_probe --summarize reports/run-20261011T2025Z/raw > summary.json
```

## Run conditions

- **Host**: a shared Apple M4 Pro with a 1-minute load average of 63.6–110.9.
  Isolation was L0, unpinned.
- **Commit**: the run used the working tree of this PR's commit. The code is the
  same except `--summarize`, which was added after the runs.

## Results

Curves are named by ICV1 slug. The subgroup is the prime-order subgroup that
`KoblitzCurve::new` selects. The base is a standard subspace, m = 3.

In the first table, "F4 refutations" counts nodes refuted at degree 2 against those
refuted only at degree 3. "Refuted-node share" is the share of the stage unit spent
in nodes that end refuted.

| cell | attempts (hits) | F4 refutations, deg 2 / deg 3 only | forced (degree) | refuted-node share |
|:--|--:|--:|--:|--:|
| K1/2^17 dim 6, `icv1-f2m17-tm101-00378d4e` | 8 (1) | 203 / 305 | 189 (3) | 0.173 |
| K1/2^23 dim 7, `icv1-f2m23-tm5197-1f85e9e1` | 4 (0) | 387 / 115 | 178 (3) | 0.038 |
| K0/2^23 dim 8, `icv1-f2m23-t5197-69e76b73` | 4 (0) | 425 / 770 | 409 (3) | 0.203 |
| K0/2^31 dim 9, `icv1-f2m31-tm90707-c95f16f5`, budget 2000 | 1 (0) | 536 / 67 | 205 (3) | 0.0076 |

The second table gives the stage unit, `word_xors`. It is deterministic and was
identical in all 5 rounds of every arm. The original and previous arms are equal.

| cell, engine | original = previous | cand1 (=1) | cand2 (=2) | original/cand1 | original/cand2 |
|:--|--:|--:|--:|--:|--:|
| K1/2^17, f4 | 15 015 618 | 14 959 088 | 15 058 684 | 1.00378 | 0.99714 |
| K0/2^23, f4 | 158 512 775 | 158 321 408 | 158 846 162 | 1.00121 | 0.99790 |
| K0/2^31, f4 | 281 583 026 | 281 350 163 | 281 787 203 | 1.00083 | 0.99928 |
| K1/2^17, f6p | 12 635 914 | 12 636 832 | 12 665 650 | 0.99993 | 0.99765 |
| K0/2^23, f6p | 124 949 801 | 124 951 833 | 125 135 729 | 0.99998 | 0.99851 |
| K0/2^31, f6p | 278 237 634 | 278 238 729 | 278 409 768 | 1.00000 | 0.99938 |

- **Identity**: every per-target verdict and counter is identical across the 4 arms
  and 5 rounds: hit, exhausted, reductions, splits and infeasible. The one relation
  found, K1/2^17 target 3 with base indices [21, 48, 57], is identical across modes 0,
  1 and 2.
- **Wall time**: the L0 medians are in the run report. The original-against-previous
  null control ranges 0.14–3.15 per round, so no wall claim is made.

### Reading

- Refuted nodes make up about 42–43 % of the reductions, but they are cheap. The
  inherited basis already stops at the refuting row.
- The cost is in the nodes that propagate and split. Those are decided only at
  degree 3:
  - no node was forced at degree 2;
  - no root system was refuted at degree 2 or 3.
- The attempt-level "refutation degree" is the depth of the tree, not a Macaulay
  degree.

### Lead (outside W4)

In F6-IC, 67–74 % of node word operations are spent in nodes that the geometry gate
eventually refutes. See REPORT.pdf.

## Canonical graphs checked and left unchanged

These were left unchanged because no speedup was recorded:

- `docs/ic/progress-timeline.json`
- `docs/index-calculus-scoreboard.html`
- `docs/performance-gains/`

## Limits

- Isolation was L0 only.
- `word_xors` omits Macaulay row generation. Mode 1 also skips the 2–10 degree-3
  builds per cell that it replaces, and those are not counted.
- Only the four cells shown, with m = 3 and `max_degree` 3.
