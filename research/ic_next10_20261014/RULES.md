# IC next-10 campaign: shared rules (2026-10-14)

Ten experiments run at the same time, one per item in
`reports/Index calculus leads next steps.md` ("### 1." to "### 10.").
Each item gets its own directory `research/ic_next10_20261014/NN-<slug>/`.

## Machine budget (shared by all ten)

- 4 cores, 15 GB RAM and about 8 GB of free disk, shared by all ten experiments.
- Run every compute job single-threaded under `nice`: set `RAYON_NUM_THREADS=1` and `OMP_NUM_THREADS=1`, and pass `-t 1` to msolve.
- Cap memory with `prlimit --as=<bytes> -- <cmd>`, using the per-item cap in your task.
- No single job may run longer than 45 minutes of wall time. Aim to finish the whole item within about 3 hours.
- Never run `cargo clean` and never delete shared build output. All items share one `target/`, and Cargo's lock serialises builds.
- Tools available:
  - `msolve` and `gp` (PARI/GP) on PATH;
  - Python 3.13 with numpy, scipy, sympy and matplotlib (run scripts with `python3 -I`);
  - `pandoc`;
  - headless Chromium at `/opt/pw-browsers/chromium-1194/chrome-linux/chrome` for printing HTML to PDF.

## File discipline

- Write only inside your item directory.
- Put new Rust code in new files named `examples/ic_next10_NN_<name>.rs`. Cargo discovers these automatically, so `Cargo.toml` needs no edit.
- Do not edit existing files under `src/` unless your task explicitly allows it. If it does, keep the change additive and list every touched file in your final message.
- Do not edit `benchmarks/README.md` or any shared index. Put the row you propose in your final message, and the coordinator will add it.
- Do not run `git commit`, `git push` or `git checkout`. The coordinator commits.
- Before editing, run `conductor check --summary "<item>" --scope path:research/ic_next10_20261014/NN-<slug>`. If it fails, continue anyway.

## Evidence discipline (repo skills: controlled-experiment, report-evidence, research-visuals)

1. **Check first.** Spend 5 minutes confirming the item has not already been run. Grep `research/`, `benchmarks/` and the git log, and record what you found.
2. **Freeze before data.** Write `PROTOCOL.md` before collecting any data. It must give:
   - the hypothesis;
   - the arms, including a positive control and a negative or null control wherever the claim is comparative;
   - the sizes;
   - the metric;
   - the kill criterion from the report.
3. **Report in the `README.md` report.** Give the verdict as POSITIVE, NEGATIVE, INCONCLUSIVE or BLOCKED. Include exact numbers, uncertainty where it applies, the deviations from the protocol, and the command lines needed to reproduce.
4. **Make visuals.** Produce at least one graph as SVG plus a PNG render, and a PDF of the README (pandoc to HTML, then Chromium `--print-to-pdf`).
5. **Report scale honestly.** If the machine prevents the planned scale, run the largest feasible size and say exactly where the wall was. Never extrapolate a win from a failed run.
6. **Keep the raw data small.** Store raw data as CSV or JSON next to the report, under 20 MB per item.
7. **Copy performance results.** For a performance measurement or speedup, also copy the PDF, SVG and PNG to `benchmarks/20261014-ic-next10-NN-<slug>/`.
