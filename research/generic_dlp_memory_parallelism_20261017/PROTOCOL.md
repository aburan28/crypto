# Protocol: memory-level parallelism in BSGS, inversion-free rho start-up

Frozen 2026-10-17, before any timed run. Study:
`research/generic_dlp_memory_parallelism_20261017/`. Conductor task T-1.

## Question

Two engineering changes to the single-word generic-DLP engines. Neither
changes what is computed:

- **H1 (BSGS).** `bsgs_fast` probes and fills an open-addressing table that
  outgrows the cache once the interval width passes about `2^36`.
  - The kernels now tell the table each visited point's `x` before the
    probe: a round ahead for `chain`, and one pass ahead for `paired` and
    `coz`. The table then issues a software prefetch of the slot.
  - The hypothesis: one-target, cold, full-group solves get faster in wall
    time on tables larger than L2.
  - The baby steps, giant steps, candidates, false candidates, table
    entries and verified answers stay identical.
- **H2 (rho).** `prime_fast::FastCurve::{mul, combine}` used affine
  double-and-add, paying one binary-GCD inversion per group operation. They
  now run in Jacobian coordinates, and `combine` uses a single Straus–Shamir
  chain.
  - `rho_parallel` calls `combine` for its 1024-entry `r`-adding table and
    for every walk start. None of that work is in `total_steps`.
  - The hypothesis: `rho_parallel` gets faster in wall time, most of all
    where this start-up dominates the walk.
  - `total_steps`, `distinguished_points` and the recovered log stay
    identical at one thread.

Both are **engineering** in the sense of AGENTS.md §3. The counted operation
units do not move by construction. The claim, if any, is wall time on the
host class measured.

## Code under test

- **Baseline A:** `main` at `b1bf5b8935f68257dfdd3fc542cb2fef41db741c`.
  - The binaries were built before any edit and kept unmodified:
    `lit2026_w11_bsgs_time` and `lit2026_w11_bsgs_opcount`.
  - `generic_dlp_rho_time` is built from the same baseline tree with only
    that example file added. It calls nothing that main lacks.
- **Baseline A′:** the same binary as A, run as a separate arm. This is the
  A/A control.
- **Candidate B:** the branch `claude/generic-dlp-memory-parallelism`. The
  commit is recorded in `manifest.json` at run time.
- **Build:** `cargo build --release`, the default profile, rustc 1.97.0. No
  `target-cpu` flags; the prefetch is the baseline-x86-64 `PREFETCHT0`.

## Frozen inputs

`inputs/curves.json` holds every curve with `p`, `a`, `b`, the group order and
the generator.

- **BSGS.** These are the W11 study's curves
  (`research/lit2026_w11_generic_engineering_20261011/reports/run-20261011T2049Z/raw/bsgs_curves.jsonl`):
  `y² = x³ + 7` at 40, 44, 48 and 52 bits, with generator `(Qx, Qy)` and
  order `n`.
  - The driver is `lit2026_w11_bsgs_time`: `targets = 1`, target seed
    `bits`, full group, negation map.
  - Kernels: `BSGS_STEP=chain`, and `paired` at the module default `K = 16`.
  - Threads: `RAYON_NUM_THREADS=1` is the primary arm. A secondary arm uses
    4 threads at 48 bits.
- **Rho.** The curve is `find_a3_curve(bits, 3)` (`a = −3`, prime order) at
  32, 36, 40, 44, 48 and 52 bits.
  - The driver is `generic_dlp_rho_time <bits> <threads> <seed>`, seed 1.
  - Threads: 1 is the primary arm (deterministic). A secondary arm uses 4
    threads.
  - The configuration is `RhoConfig::for_bits(bits, true)` with `threads`
    overridden.
- **External reference (AGENTS.md §15).** PARI/GP 2.15.4 (Ubuntu noble
  `pari-gp`) runs `elllog` on the same curves and the same targets,
  single-threaded. It is timed with gp's own `getabstime()` around the
  `elllog` call alone.

## Arms, order, rounds, isolation

- **Arms.** A, A′ and B per configuration.
- **Order.** Rotated each round: `A B A′`, then `B A′ A`, then `A′ A B`, …
- **Rounds.** 7 at sizes up to 48 bits, 5 at 52 bits.
- **Isolation.** Every timed process runs under
  `tools/isolated_bench.py run --cpus 3` (exclusive lock, quiet-machine
  preflight, reserved CPU 3). The 4-thread arms use `--cpus 0-3`.
  - The isolation record of every run is kept in `raw/isolation.jsonl`.
  - The host is a 4-vCPU KVM guest. Expect ±5–10% residual noise (AGENTS.md
    §10); the A/A arm measures it.
- **Concurrent work.** No builds, tests or other agents run during a timed
  stage.

## Metrics

- **Primary:** each binary's internal timers.
  - BSGS: `build_s`, `search_s` and `total_s`.
  - Rho: `solve_s`, the whole `rho_parallel` call.
  - All of them exclude process start-up, curve and target construction,
    and the final verification print.
- **Secondary:** the isolation record's process `wall_seconds`, plus
  contention and switches.
- **Derived, not a phase:** rho's `setup_s`, the same number of `combine`
  calls timed alone.
- **Pinned outputs.** The following must be identical between A and B at
  one thread:
  - BSGS: `answers_digest`, `table_entries`, `baby_steps`, `giant_steps`,
    `candidates` and `false_candidates`.
  - Rho: `log`, `verified`, `total_steps` and `distinguished_points`.
  - At 4 threads the giant early exit and the DP arrival order depend on
    scheduling, so only the answer is pinned there.
  - `lit2026_w11_bsgs_opcount` must also report identical field-operation
    counts per visited point (`M`, `S`, `I`) for A and B. This is the
    deterministic measure beside the wall time.

## Analysis

For each configuration:

- the median and minimum of each arm;
- the ratio of medians, A/B and A/A′;
- the minimum and maximum of the per-round paired ratios A/B and A/A′.

A configuration **improves** only when every per-round A/B ratio is above
the largest per-round A/A′ ratio (or below the smallest, for a regression).
Otherwise it is reported as **no measurable change**. Contended runs are
reported and never pooled with uncontended ones.

## Success and stop conditions

- **Reject the candidate** if any pinned output differs, or if any answer
  fails verification. A faster run that decides anything differently is a
  different algorithm.
- **Success (engineering)** is H1 or H2 met at one or more measured
  configurations, under the rule above.
- **A regression** at any configuration is reported as one. It is not
  averaged away.

## What this cannot establish

- No change to `S` or to any counted operation, so no IC/rho ratio moves,
  and nothing here touches the index-calculus scoreboard, leaderboard or
  progress chart (AGENTS.md §7a).
- No ranking of ECDLP methods. This is a same-method before/after study.
  The PARI/GP numbers are an external wall-time reference (§15), not a
  counted-unit comparison, so §5a's pricing band does not apply.
- Results hold for this host class (x86-64, a Sapphire/Emerald Rapids–class
  Xeon under KVM) only. Arm64, GPUs and other x86 parts were not measured.
