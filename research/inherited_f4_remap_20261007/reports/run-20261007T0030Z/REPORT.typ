#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set heading(numbering: none)

= Inherited-F4 remap round: run report `run-1`

*Status: measured, passed its preregistered rule against its own baseline, superseded by `main` before merge, not merged.*
Started 2026-10-07T00:30Z (UTC). Study directory: `research/inherited_f4_remap_20261007/`.

== What ran

- *Protocol:* `PROTOCOL.md`, written before any confirmatory run. *Runner:* `run.sh` (thin shell around `ecbench exec`, Valgrind Callgrind, `taskset`, `jq`).
- *Baseline:* `ecbench` from `main` at `a3652eddafb949073bf0f96b210a75b671cfb346`, SHA-256 `7e6ad0d2…`. A byte-identical copy served as the A/A arm.
- *Candidate:* the baseline plus `results/run-1/candidate.diff` (SHA-256 `3779a2b0…`), SHA-256 `7bfdd8c9…`. It was built from a working tree and was never committed, because a Conductor reservation required by `CLAUDE.md` could not be taken. A lint-only revision is kept as `results/run-1/final.diff`.
- *Host:* Intel Xeon \@ 2.80 GHz cloud VM, 4 logical CPUs, 16 GiB, Linux 6.18 (kernel build `fc-v77`), BMI2 present; `rustc` 1.97.0; Valgrind 3.22.0 (`results/run-1/host.txt`).
- *Inputs:* 60 exact `ecbench` child inputs exported from frozen specs (`results/run-1/inputs.sha256`):
  - ladder: m = 13 (16), m = 19 (16), m = 23 (4);
  - round-1 holdouts: m = 13 (8), m = 19 (4);
  - fresh holdouts: m = 13 (8, target seed 202610071013), m = 19 (4, target seed 202610071019).
- *Arms:* the inherited-F4 and F6-IC oracles of the F6-IC ladder, unchanged.

== Runs requested and completed

#table(columns: 4, align: (left, right, right, left),
  [*measurement*], [*requested*], [*completed*], [*note*],
  [identical-output gate], [60 inputs × 2 binaries], [60 × 2], [all 60 `IDENTICAL` (`gate.tsv`)],
  [Callgrind total Ir], [44 inputs × 2], [44 × 2], [`ir.tsv`, `ir.raw.tsv`],
  [pinned CPU time, ABAB with A/A], [3 inputs × 5 rounds × 4], [55 of 60 rows], [worker restart killed the run in L23 seq 25 round 4 (`wall.tsv`)],
  [L23 seq 25 re-run, new kernel build `fc-v80`], [5 rounds × 4], [20 of 20], [`wall-L23-rerun.tsv`, same binaries],
)

== Measured: instructions (Callgrind Ir, whole process)

One run per binary and input; Callgrind counts are deterministic. The ratio is candidate / baseline.

#figure(image("ir_ratio_per_input.svg", width: 100%),
  caption: [Per-input ratio from `results/run-1/ir.tsv`, run-1, 44 inputs. Blue: inherited F4; orange: F6-IC. Dashed: the preregistered threshold on the pooled confirmatory m = 19 ratio.])

#table(columns: 4, align: (left, right, right, right),
  [*m, arm*], [*inputs*], [*pooled ratio (derived: sum cand / sum base)*], [*per-input range*],
  [13, inherited F4], [16], [0.984], [0.983–0.986],
  [13, F6-IC], [16], [0.978], [0.976–0.981],
  [19, inherited F4], [6], [0.936], [0.935–0.938],
  [19, F6-IC], [6], [0.907], [0.905–0.909],
)

Decision statistic (derived): pooled over the confirmatory m = 19 inputs (L19 seq 24, 27, 28, H19, H19b), 423,612,106,635 / 458,788,722,509 = 0.9233. The preregistered threshold was ≤ 0.97. The largest per-input ratio is 0.986 (H13-4).

== Measured: pinned process CPU time (practicality note)

`taskset -c 3`, run alone on the host. Each round runs baseline, candidate, the A/A copy, then the candidate again. Scheduling delay ≤ 7.3 ms per run.

#table(columns: 4, align: (left, right, right, right),
  [*input*], [*baseline median (s)*], [*candidate median (s)*], [*A/B median (range); A/A median (range)*],
  [L19 seq 25, inherited F4], [7.13], [6.44], [0.908 (0.889–0.919); 0.994 (0.982–1.008)],
  [L19 seq 24, F6-IC], [4.68], [4.05], [0.866 (0.838–0.877); 1.006 (0.976–1.030)],
  [L23 seq 25, rounds 1–4, kernel `fc-v77`], [90.0], [78.1], [0.866 (0.860–0.870); 1.000 (0.995–1.008)],
  [L23 seq 25, re-run, kernel `fc-v80`], [91.7], [78.8], [0.860 (0.841–0.873); 0.993 (0.979–1.016)],
)

== Superseding check (2026-10-11)

`ecbench` built from `main` at `7666b708e` (SHA-256 `6ae1a2a1…`), which contains #1646 (a word-wise `pext`/`pdep` rewrite for every map) and #1671 (eager row rewrites). It was measured on three of the same inputs (`results/main-7666b708-check.tsv`). Its counted outputs differ from the round-2 baseline (`solver_ops` higher on all three; F6-IC gate statistics differ on L19 seq 24), so this is not a gated A/B.

#figure(image("main_check_ir.svg", width: 90%),
  caption: [Raw Callgrind Ir of the round-2 baseline, the round-2 candidate and `main` `7666b708e` on L13 seq 25, L19 seq 24 and L19 seq 25.])

#table(columns: 4, align: (left, right, right, right),
  [*input*], [*round-2 candidate Ir*], [*main Ir*], [*main / candidate (derived)*],
  [L13 seq 25], [2,185,880,355], [2,286,390,878], [1.046],
  [L19 seq 24], [30,866,373,759], [27,799,877,737], [0.901],
  [L19 seq 25], [48,778,989,164], [46,233,505,868], [0.948],
)

== Raw records

`results/run-1/`: `gate.tsv`, `ir.tsv`, `ir.raw.tsv`, `ir/*.err`, `wall.tsv`, `wall-L23-rerun.tsv`, `inputs/`, `host.txt`, `binaries.sha256`, `candidate.diff`, `final.diff`, `lint-followup.tsv`; `results/main-7666b708-check.tsv`. The graph source is `plot.py` in this directory, which reads only those files.
