# lit2026 speedup program: from the 2025–26 literature to measured speedups

**Started.** 2026-10-11. **Program baseline commit.** `4b6e09b78` (main at start).
**Source.** [`../literature_search_20261011/README.md`](../literature_search_20261011/README.md).

The program turns every finding of the 2025–26 literature search into one of three
things:

- a measured change, against an interleaved baseline;
- a diagnostic with numbers;
- a recorded negative.

Workstreams W1–W12 are listed in [`workstreams.json`](workstreams.json). Each
workstream lands its own PR and study directory under `research/`, with a run report,
a graph and a PDF (research-visuals skill). Each PR also adds one or more entries here:
[`journey/entries/`](journey/entries/), schema in [`journey/SCHEMA.md`](journey/SCHEMA.md).

## The journey graph

- [`journey/journey.svg`](journey/journey.svg) shows, for each workload, the measured
  ratio of program baseline to candidate.
- Every point comes from one interleaved session in which the baseline binary at
  `4b6e09b78` ran beside the candidate. Ratios are never chained across sessions or
  hosts, so a rising line is a measured fact about that workload, not a product of
  step ratios.
- [`journey/board.svg`](journey/board.svg) shows each workstream's status and its
  recorded steps by kind: speedup, no-change, regression, diagnostic, negative.

Regenerate after adding an entry or changing a status:

```bash
sh research/lit2026_speedup_program/render.sh
```

When two PRs both regenerate, the rendered SVG and PDF files can conflict. Resolve
the conflict by re-running the script; never by hand-editing them.

## Rules every workstream follows

1. **Order.** Conductor claim first, then a baseline, then the change, then an
   interleaved measurement, then a report, then a PR. Opening the PR ends the
   workstream's turn.
2. **Measurement.**
   - Deterministic counters (rows, reductions, Ir, counted ops) are the correctness
     and speed gate.
   - Wall time is reported only as the harness admits it. On this shared Mac that
     means L0 with interleaved rounds and `uptime` recorded per sample.
   - Outputs (solutions, relations, verdicts) must be identical between arms.
3. **Arms.** Each speed claim compares at least the `original` arm (`4b6e09b78`) and
   the candidate in the same session. Add `previous` once an earlier step has landed
   on the same workload. A one-binary env switch, like `KIC_F4_EAGER` in
   `research/inherited_f4_eager_rewrite_20261011/`, gives the `previous`/`candidate`
   pair. A separately built `4b6e09b78` binary gives `original`.
4. **Status words.**
   - A step that does not speed anything up is recorded as `no-change` or `negative`,
     with its numbers.
   - No exponent or deployed-curve claims; AGENTS.md claim hygiene applies.
   - A constant-factor change to an IC stage is not an IC-versus-rho result unless an
     ecbench comparison against the reference arm says so.
5. **Host.**
   - One cargo command at a time per worktree.
   - `CARGO_TARGET_DIR` and worktrees go on `/Volumes/SSD990-2`.

## Status

See [`journey/board.svg`](journey/board.svg); the data is `workstreams.json`. The PDF
is [`PROGRAM.pdf`](PROGRAM.pdf).
