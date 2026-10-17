# IC improvement round 2: shared rules (2026-10-15)

Three follow-ups to [`../ic_improve_20261014/README.md`](../ic_improve_20261014/README.md), run at the same time.
Each lives in its own directory, `research/ic_improve_20261015/<NN>-<slug>/`.

## What is measured

Every item measures an **improvement over the incumbent IC component it replaces**, as a paired comparison:

- the same curves and targets;
- the same seeds where they apply;
- the same instrument and counted unit on both arms.

The incumbent is the IC configuration the repository uses today. Your task names it.

Rho is a reference column only. **No verdict depends on beating rho.**

**Verdict vocabulary** (from ecbench-bounds):

| Verdict | Meaning |
|---|---|
| ADVANCES | better, CI excludes 1 or 0 |
| TRADES | better on one phase or level, worse on another |
| MATCHES | CI covers no change |
| REGRESSES | worse |

State the level of the change: **exponent** (fitted slope with CI over ≥4 sizes), **constant** (paired ratio with CI) or **primitive** (cost per operation of one stage).

**Correctness** gates every arm: each logarithm and each decomposition is verified, and any arm whose outputs disagree with the incumbent fails.

## Machine budget (shared by three items)

- 4 cores, 15 GB RAM and about 6 GB of free disk.
- Every compute job runs single-threaded under `nice`. Set `RAYON_NUM_THREADS=1` and `OMP_NUM_THREADS=1`, and pass `-t 1` to msolve.
- Cap memory with `prlimit --as=<bytes> -- <cmd>`, using your item's cap.
- No single job may run longer than 60 minutes of wall time. Aim to finish the whole item within about 4 hours.
- One shared `target/`. Never run `cargo clean`. Cargo's lock serialises builds, so wait for it rather than fighting it.
- **Rust for research results** (AGENTS.md "no Python"). Use Python only for figures, unless your task says otherwise.
- Tools available: `msolve`, `gp`, `pandoc`, and headless Chromium at `/opt/pw-browsers/chromium-1194/chrome-linux/chrome`.

## Files and the shared build

- Write inside your directory.
- New Rust code goes in new files: `examples/ic_improve2_<NN>_*.rs`, or new modules under `src/` if your task allows it.
- When you must touch an existing `src/` file:
  - Write the new module completely first.
  - Then register it with the smallest possible edit, made in one go.
  - Run `cargo build --release --lib` straight away.
  - Other items build from the same tree, so the tree must never be left broken for more than a few minutes.
  - Stay inside the source files your task lists.
- Do not edit `benchmarks/README.md`. Propose the row in your final message.
- No `git commit`, `git push` or `git checkout`. The coordinator commits.
- Run `conductor check --summary "<item>" --scope path:<dir>` first. If it fails, continue.

## Evidence

1. Write `PROTOCOL.md` **before** collecting data. It must state:
   - the incumbent and the candidate;
   - the paired design;
   - the sizes;
   - the metric;
   - the decision rule for each verdict;
   - a positive control (a known improvement the instrument must detect) or a null A/A arm.
2. Write `README.md` as the report. Include:
   - the verdict and level;
   - exact numbers with CIs;
   - deviations from the protocol;
   - reproduce commands.
3. Produce an SVG graph with a PNG render, and a PDF of the README.
4. For performance results, also copy the PDF and figures into `benchmarks/20261015-ic-improve2-<NN>-<slug>/`.
5. Keep raw data under 20 MB.
6. When the machine is the limit, say exactly where.
