# End-to-end: root-finding vs F4 decomposition inside complete DLP solves

Produced by `examples/semaev_rootfind_e2e_bench.rs`. All times are **CPU
time** (process + reaped children via `getrusage`), not wall-clock. Every
solve recovered and verified the planted logarithm (`[d]G = Q`); the bench
aborts otherwise.

Two runs, because the msolve accounting method changed between them:

| file | what | status |
|---|---|---|
| `run1_console.txt`, `run1_results.json` | complete run, all sections (`EXIT=0`) | **binary rows are the evidence of record.** Its msolve cells used an up-front spawn-floor subtraction whose ~1 ms drift is visible as a flat ~1 ms at `\|F\| = 4–8`; superseded for every msolve cell by run 2. Root-finding and Buchberger cells agree with run 2 within noise. |
| `run1b_killed_console.txt` | rerun with paired-control msolve accounting, killed at the 20-bit row (`EXIT=143`) | kept for provenance only; no row is taken from it |
| `run2/console.txt`, `run2/results.json` | prime section, msolve only (`E2E_SKIP_BINARY=1 E2E_SKIP_BUCHBERGER=1 E2E_BEST_REPS_DIV=4 E2E_FULL_BUDGET_S=300`), paired-control msolve accounting | **root-finding and msolve prime cells, and the best-`\|F\|` block, are the evidence of record.** Buchberger cells are taken from run 1 (same seed, same target stream — root target counts match exactly — same in-process CPU clock); it was skipped here only because the machine was saturated (load > 100) and the 155 CPU-s solve is unaffected by the msolve accounting change. |

Paired-control accounting: each msolve decide is immediately followed (or
preceded, alternating) by an msolve run on the trivial ideal `{X, Y}`; the
mean paired difference in child CPU is what is charged, so process start-up
and file I/O are excluded. A per-decide value is reported only if it clears
two standard errors of the paired mean; otherwise the 2·SE bound is charged
(in the baseline's favour).

Budget rule (`E2E_FULL_BUDGET_S`): a baseline is measured in full when its
projected time is under the budget; otherwise it is COMPOSED (sampled
per-target cost — 100 targets for Buchberger, 400 for msolve — × the root
finder's measured target count on the same stream) and labelled COMPOSED in
the console and `*_measured=false` in the JSON. Run 2 used 300 s so that the
16-bit msolve row is fully measured and the 20-bit one composed, and divided
the best-`|F|` sample sizes by 4 (500/125/25/10 reps), under a machine load
average above 100. The 2·SE resolvability test still applies per row.

Not a security claim. Public synthetic known-answer instances; the speed-up
is a constant on one stage of an algorithm whose asymptotics are unchanged.
