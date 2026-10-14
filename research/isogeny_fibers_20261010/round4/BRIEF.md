# Round 4 brief: pursuing the §16 leads (shared by every lead)

Context: `RESEARCH_ISOGENY_FIBERS_ECDLP.md` in this worktree. §13 (fiber
invisibility), §14 (image-vs-fiber law, linearity barrier, Prop. 14.2 4E
test), §15 (Thm 15.1 isogeny-augmented GGM, Prop. 15.2 affine saturation
law, §15.4 table) and §16 (the leads) are the background. Read them first.

## Hard rules

- **Native code only.** AGENTS.md "Implementation language: no Python":
  every computation that produces a result is Rust, as a new file
  `examples/r4_<lead>_<name>.rs` that uses the `crypto_lib` crate
  (see existing examples for the import style, e.g. `use crypto_lib::...`).
  Python is allowed only to plot already-recorded numbers and render the
  PDF, with the shared venv
  `/tmp/claude-0/-home-user/717e4727-b5a7-5b25-924a-7b31edb43331/scratchpad/venv/bin/python`
  (has `matplotlib` and `typst`).
- **Do not edit anything under `src/`, `Cargo.toml`, `Cargo.lock`, or any
  existing file.** Other agents build in this same worktree at the same
  time; a `src/` edit rebuilds the library for everyone and can break their
  builds. If a lead genuinely needs a library change, write it as code in
  your example and say in your report what should move into `src/`.
- **Do not run git commands that write** (no add/commit/stash/checkout).
  The coordinator commits.
- Build and run with the shared target dir and release mode:
  `export CARGO_TARGET_DIR=/tmp/claude-0/-home-user/717e4727-b5a7-5b25-924a-7b31edb43331/scratchpad/r4target`
  then `cargo run --release --example r4_<...> -- <args>`. Cargo locks the
  target dir, so a build may wait for another agent's; that is expected.
- The machine has 4 cores shared by ~6 agents. Run single-threaded, keep
  each run under ~10 minutes, and **prefer counted operations** (group
  operations, field multiplications, relations, solver calls) over wall
  time. If you report wall time, label it as contended wall time.
- Before finishing, your example must pass
  `cargo clippy --release --example r4_<...> -- -D warnings` and
  `rustfmt --edition 2021 --check examples/r4_<...>.rs`.
- Claim hygiene: toy and synthetic instances only. No sub-rho claim, no
  claim about a deployed curve's security. Name curves with ledger slugs
  (`python3 scripts/check_curve_names.py` must stay at 0 problems for your
  files; for Koblitz curves write e.g. the `icv1-...` slug or describe
  them as "the a=0 Koblitz curve over F_2^n" in prose, never retired K-subscript
  style names). ECC2K-130 and sect233k1 may be named as standards.
- Report evidence faithfully (`.claude/skills/report-evidence/SKILL.md`):
  label each number measured or derived, state negative and null results
  plainly, and say what was not done.

## Deliverables (in `research/isogeny_fibers_20261010/round4/<lead>/`)

1. `RESULTS.md`: the question, hypothesis and falsifier from §16, what ran
   (exact commands), measured tables, the verdict (confirmed, refuted,
   inconclusive, with the reason), limits, and the next step.
2. The raw outputs your example printed (`*.json` or `*.txt`).
3. At least one graph (`*.svg` plus the small plotting script and its data)
   of recorded values only, following `.claude/skills/research-visuals/SKILL.md`.
4. `report.typ` and the compiled `report.pdf` with the graph embedded.
   Inspect the PDF (page count, that the graph is present).
5. A short final message to the coordinator: verdict in 3-5 lines, the key
   numbers, files written, anything that should move into `src/`.
