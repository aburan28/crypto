# lit2026 workstream brief (shared by all workstream agents)

You carry out one workstream (Wn) of the lit2026 speedup program in `aburan28/crypto`.
Read these files first:

- `research/lit2026_speedup_program/README.md`
- `research/lit2026_speedup_program/workstreams.json`
- `research/lit2026_speedup_program/journey/SCHEMA.md`
- `research/literature_search_20261011/README.md` and `papers.csv`

## Setup

1. You run in your own git worktree. Put your branch on the program commit, which is
   stacked on PR aburan28/crypto#1681:

   ```
   git fetch origin feat/elliptic-curve-crypto-research-88c7cc
   git checkout -b lit2026/<wn>-<slug> origin/feat/elliptic-curve-crypto-research-88c7cc
   ```

   - Set `export CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/lit2026-<wn>` for
     every cargo command, because SSD990 is filling up.
   - Never use bare `git stash`.
2. Read `AGENTS.md`, especially these sections:
   - "Preserve scope and report evidence"
   - "Implementation language: no Python"
   - "Research searches must leave visual reports"
   - the claim-hygiene rules
3. Read `.claude/skills/research-visuals/SKILL.md`. Also read
   `.claude/skills/ecbench-measure/SKILL.md` if you compare ECDLP methods end to end.
4. Conductor:
   - Run `conductor check --summary "<wn>: …" --scope dir:<paths>` before editing.
   - Then `conductor task create --title … --scope …` and `conductor task claim <ref>`.
   - If a check names a holder of your files, coordinate rather than overwrite.
   - Report scope expansion before editing outside your reserved paths.

## Implementation and measurement

- All algorithm, harness, diagnostic and measurement code is Rust, integrated with the
  existing library, examples or bins and their tests. Python is allowed only to plot
  recorded values and render PDFs.
- **Build contention.** The host is a shared M4 Pro carrying load 60–290 from other
  sessions. A release build of `crypto_lib` takes about 20+ minutes.
  - Issue one cargo command at a time, run it in the background, and wait for it.
  - Prefer `--example`/`--bin` targets while iterating.
  - Before pushing, run the touched examples' tests and the relevant
    `cargo test --release --lib <filter>`.
- **Original (program-baseline) arm.**
  - A shared detached worktree at commit `4b6e09b78` (main when the program started)
    is at `/Volumes/SSD990-2/crypto/wt/lit2026-baseline`. Build baseline binaries from
    it with
    `(cd /Volumes/SSD990-2/crypto/wt/lit2026-baseline && CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/lit2026-baseline cargo build --release --example X)`
    (or `--bin X`).
  - Other agents share that target dir, so cargo may block on its lock; that is
    expected, so wait. Never edit files there.
  - Copy built binaries to your scratch directory before timing.
- **Arms.** Make your change switchable in one binary (an env var or CLI flag, like
  `KIC_F4_EAGER` in `research/inherited_f4_eager_rewrite_20261011/`). That gives the
  `previous` and `candidate` arms. The baseline binary gives `original`.
- **Interleaving.**
  - Run the arms interleaved (ABAB…) for at least 5 rounds and record `uptime` with
    each sample.
  - Prefer deterministic counters (rows, reductions, pairs, Ir, counted ops) as the
    gate. Report wall time only as L0 interleaved medians with min–max spread.
  - Check that outputs (solutions, relations, verdicts) are identical across arms.
- **Workloads.** Use the exact workloads the repo already measures where they exist:
  Weil-descended S3/S4 over GF(2^n) for Koblitz and binary curves, the inherited-F4
  probe cells, ecbench specs. Name curves by registered ids and state the subgroup.
- **Results that do not pan out.** If the idea gives no gain, record it as
  `no-change` or `negative` with numbers. That counts as a completed workstream.
  Never invent, extrapolate or round numbers in your favour.

## Deliverables (one PR per workstream)

1. **Code.**
   - The change, behind its switch.
   - Tests that pin output identity, or the diagnostic's correctness on a small
     known case.
2. **Study directory.** Put it at `research/lit2026_<wn>_<slug>_<YYYYMMDD>/`, containing:
   - `README.md`: question, papers, method, exact commands, commits, host, arms,
     results table, limits;
   - `reports/run-<UTC>/`: run report, raw outputs, graph SVG and its source, `report.pdf`;
   - a results PDF built with Typst. There is a venv with typst and matplotlib at
     `/private/tmp/claude-501/-Volumes-SSD990-crypto-elliptic-curve-crypto-research-88c7cc/ab9482d3-a14b-4bea-8cfa-130f6184453a/scratchpad/venv/bin/python`.

   Inspect every PDF page, for example by rendering PNGs with typst and reading them.
3. **Journey entry.**
   - Add `research/lit2026_speedup_program/journey/entries/<YYYYMMDD>-<wn>-<slug>.json`
     per `journey/SCHEMA.md`. Add one entry per workload measured.
   - Set your workstream's `status` in `workstreams.json` to "PR open" and add the PR
     reference to `prs`. Edit only your own workstream's object.
   - Re-render with `PY=<venv python> sh research/lit2026_speedup_program/render.sh`.
4. **Canonical graphs.**
   - If the result changes a canonical graph (`docs/ic/progress-timeline.json`,
     `docs/index-calculus-scoreboard.html`, `docs/performance-gains/`), update its
     source data and regenerate per AGENTS.md section 7.
   - Otherwise, list the graphs you checked and left unchanged in the README.
5. **Commit and PR.**
   - End each commit message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
   - Push, then run `gh pr create --base main` to open the PR **ready for review,
     never a draft**.
   - The body says it is stacked on aburan28/crypto#1681. It gives measured ratios
     with intervals and units, and identity checks. It ends with:
     `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.
6. **Stop.**
   - Opening the PR ends your work. Do not watch CI, poll checks or schedule
     check-ins.
   - Your final message reports:
     - the PR URL;
     - each measured result (arms, values, units, ratio, interval, samples,
       isolation);
     - identity checks;
     - each requirement of your workstream as verified complete / implemented but
       unverified / partial / blocked / not attempted, with the exact gap;
     - anything surprising or wrong you found in the repo.

If a step is truly blocked (a build failure you cannot fix, missing data), land what
you have with the blocker stated in the README and PR. Do not silently shrink the
workstream.
