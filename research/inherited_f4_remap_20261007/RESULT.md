# Inherited-F4 remap round: result

**Status: measured, passed its preregistered rule, superseded, not merged.**
The candidate code is not on `main`. While it was held, `main` gained a
more general word-wise rewrite (#1646: `pext`/`pdep` segment plans for every
map, folds and composed maps included) and eager row rewrites (#1671). On the
three frozen inputs re-measured against `main` at `7666b708e`, `main` runs
fewer instructions at m = 19 than this candidate. See *Superseded* below and
`results/main-7666b708-check.tsv`.

**Outcome under the rule, against its own baseline (`a3652edd`): success.
Class: engineering.**

- **Gate.** On 60 frozen inputs every output is identical.
- **Instructions.** On the confirmatory m = 19 inputs the candidate executes
  **0.923×** the baseline's instructions, against a threshold of 0.97, and no
  input regresses.
- **What did not move.** The ecbench `S` is unchanged by construction: no
  row, order, layout or charge moved, so no counted unit did either.
- **Per counted word XOR, at m = 19.** Inherited F4 goes from about 95 to
  89 instructions; F6-IC from about 84 to 76. These figures are on top of
  round 1, whose baseline was about 120 and 92.

Protocol: [`PROTOCOL.md`](PROTOCOL.md), written before any confirmatory run.
Raw results: `results/run-1/`.

**Provenance.** This round was measured on a working tree, not a commit.
`main`'s `CLAUDE.md` now requires a Conductor reservation before code
changes, and the Conductor tools were missing from the measuring container.
By the user's decision the candidate was measured locally and held back from
commit. Its exact source is `results/run-1/candidate.diff`, a diff against
the baseline revision, with SHA-256 `3779a2b0…` (`candidate.diff.sha256`).
Applying it to `a3652edd` rebuilds the candidate binary.

## Host and binaries

From `results/run-1/host.txt` and `binaries.sha256`:

- **Host.** Intel Xeon @ 2.80 GHz (cloud VM), 4 logical CPUs, 16 GiB, Linux
  6.18 x86-64, with `popcnt avx2 avx512f bmi2 pclmulqdq`. `rustc` 1.97.0,
  Valgrind 3.22.0.
- **Baseline.** `main` at `a3652edd`, which includes round 1. `ecbench`
  SHA-256 `7e6ad0d2…`; the A/A copy is byte-identical.
- **Candidate.** The baseline plus `candidate.diff`. `ecbench` SHA-256
  `7bfdd8c9…`.

## The gate: identical outputs

`results/run-1/gate.tsv` has 60 inputs and 60 `IDENTICAL`:

| set | inputs | what they are |
| --- | --: | --- |
| L13 | 16 | ladder, m = 13 |
| L19 | 16 | ladder, m = 19 |
| L23 | 4 | ladder, m = 23 |
| H13 | 8 | round-1 holdout targets, m = 13 |
| H19 | 4 | round-1 holdout targets, m = 19 |
| H13b | 8 | fresh holdout targets, m = 13 |
| H19b | 4 | fresh holdout targets, m = 19 |

## The table: instructions (Callgrind `Ir`, whole process)

The ratio is candidate / baseline. `S`, `S / floor` and `S / rho` are
unchanged and are those of the
[F6-IC ladder](../f6_ic_ecbench_ladder_20261005/RESULT.md).

| m | arm | set | inputs | baseline `Ir` | candidate `Ir` | ratio (pooled) | per-input range | `S` | class | correct |
| --: | --- | --- | --: | --: | --: | --: | --: | --- | --- | --- |
| 13 | inherited F4 | ladder | 8 | 8.933 × 10⁹ | 8.790 × 10⁹ | **0.984** | 0.983–0.985 | unchanged | engineering | 8/8 identical |
| 13 | inherited F4 | round-1 holdout | 4 | 4.529 × 10⁹ | 4.459 × 10⁹ | **0.985** | 0.983–0.986 | unchanged | engineering | 4/4 identical |
| 13 | inherited F4 | fresh holdout | 4 | 6.069 × 10⁹ | 5.973 × 10⁹ | **0.984** | 0.983–0.986 | unchanged | engineering | 4/4 identical |
| 13 | F6-IC | ladder | 8 | 5.523 × 10⁹ | 5.400 × 10⁹ | **0.978** | 0.977–0.979 | unchanged | engineering | 8/8 identical |
| 13 | F6-IC | round-1 holdout | 4 | 2.814 × 10⁹ | 2.749 × 10⁹ | **0.977** | 0.976–0.979 | unchanged | engineering | 4/4 identical |
| 13 | F6-IC | fresh holdout | 4 | 3.736 × 10⁹ | 3.652 × 10⁹ | **0.978** | 0.976–0.981 | unchanged | engineering | 4/4 identical |
| 19 | inherited F4 | ladder | 2 | 9.704 × 10¹⁰ | 9.093 × 10¹⁰ | **0.937** | 0.936–0.938 | unchanged | engineering | 2/2 identical |
| 19 | inherited F4 | round-1 holdout | 2 | 1.100 × 10¹¹ | 1.029 × 10¹¹ | **0.936** | 0.935–0.936 | unchanged | engineering | 2/2 identical |
| 19 | inherited F4 | fresh holdout | 2 | 1.043 × 10¹¹ | 9.769 × 10¹⁰ | **0.936** | 0.936–0.936 | unchanged | engineering | 2/2 identical |
| 19 | F6-IC | ladder | 2 | 6.354 × 10¹⁰ | 5.767 × 10¹⁰ | **0.908** | 0.906–0.909 | unchanged | engineering | 2/2 identical |
| 19 | F6-IC | round-1 holdout | 2 | 7.213 × 10¹⁰ | 6.535 × 10¹⁰ | **0.906** | 0.905–0.907 | unchanged | engineering | 2/2 identical |
| 19 | F6-IC | fresh holdout | 2 | 6.379 × 10¹⁰ | 5.784 × 10¹⁰ | **0.907** | 0.906–0.907 | unchanged | engineering | 2/2 identical |

**Decision statistic.** The confirmatory m = 19 inputs are L19 seq 24, 27
and 28, H19 and H19b; L19 seq 25 was used in exploration. Their pooled
ratio is 423,612,106,635 / 458,788,722,509 = **0.9233**, which is ≤ 0.97.
The worst single input is H13-4 (inherited F4) at 0.986, so no input
regresses. Pooled over all 44 `Ir` inputs the ratio is 0.928.

**Reading it.** The remap gains grow with m: the row layout widens, and a
deletion step moves whole words. F6-IC gains more than inherited F4 here,
the reverse of round 1, because more of its remaining instructions sit in
the remap. At m = 13 the gain is about 2%.

**Rounds 1 and 2 combined.** Against `fa80835a` on L19 seq 25 (inherited F4,
the exploration input of both rounds), the two rounds take
65,786,679,518 instructions to 52,095,119,421 to about 48.8 × 10⁹, about
**0.74×**. That is a chained figure from two baselines, not a single paired
run, and round 2 never shipped (see *Superseded*): `main` reached its own
figure for this input by another route.

## Time (practicality note)

`results/run-1/wall.tsv` has process CPU time from the child output, with
each run pinned by `taskset -c 3` and run alone on the host:

- **Design.** Five rounds. Each runs the baseline, the candidate, the
  byte-identical baseline copy (A/A), and the candidate again.
- **Ratio.** A/B is the round's mean candidate time over its baseline.
- **Scheduling delay.** At most 5.5 ms per run (7.3 ms in the re-run).

**Interruption.** A worker restart killed the run during L23 seq 25, after
four rounds; round 4 is missing its second candidate run. Those rows are
kept as they are. The restart also changed the VM's kernel build (same
CPU model), so L23 seq 25 was re-run in full on the new kernel into
`wall-L23-rerun.tsv`, with the same binaries, rather than appending to the
old rows.

| input | arm | baseline CPU, median (min) | candidate CPU, median (min) | A/B, median (range) | A/A, median (range) |
| --- | --- | --: | --: | --: | --: |
| L19 seq 25 | inherited F4 | 7.13 s (7.01) | 6.44 s (6.40) | **0.908** (0.889–0.919) | 0.994 (0.982–1.008) |
| L19 seq 24 | F6-IC | 4.68 s (4.64) | 4.05 s (3.97) | **0.866** (0.838–0.877) | 1.006 (0.976–1.030) |
| L23 seq 25, rounds 1–4 (interrupted) | inherited F4 | 90.0 s (89.6) | 78.1 s (77.3) | **0.866** (0.860–0.870) | 1.000 (0.995–1.008) |
| L23 seq 25, re-run (new kernel) | inherited F4 | 91.7 s (90.1) | 78.8 s (77.6) | **0.860** (0.841–0.873) | 0.993 (0.979–1.016) |

Every A/B round falls outside its A/A spread. Time falls further than
instructions: at m = 19 the instruction ratios are 0.937 for inherited F4
and 0.907 for F6-IC. One `pext` replaces a loop of dependent loads and
stores, so it saves more cycles than instructions. Wall time is a
practicality note here and decides nothing.

## What changed

- **`ColMap`.** Each layout map carries per-word `live` masks, plus the
  first target where a word's kept columns land on consecutive targets.
- **Deletion steps.** They build both in their own layout loop, and are
  consecutive on every word.
- **Folds and layout extensions.** They delete nothing, get all-ones masks,
  and are marked scattered.
- **Composed maps.** They take only the live masks, in the gather's pass.
- **`remap_row`.** It skips deleted bits with one AND per word. On a CPU
  with BMI2 it moves a consecutive word with `pext` and a shift; other words
  go bit by bit.
- **Tests.** Unit tests pin both paths against a plain bit-by-bit scatter
  on deletion-like, fold-like and shuffled maps, and pin the step-built
  masks against masks read off the map.

**Rejected, with their numbers, from L19 seq 25** (`PROTOCOL.md`, exploratory
record):

- a branch-free generic mask pass, 0.961 against 0.947;
- consecutive testing of composed words, which cost more than `pext` saved
  on them;
- drafts held as unique `Rc` slices, 0.944 against 0.936.

## What this does not establish

- **No `S`, ratio or crossover moves.** The F6-IC ladder's verdict
  stands: IC is slower than rho at m = 13, 19 and 23.
- **One host class.** The `pext` path was measured on an Intel x86-64 VM
  with BMI2. The masked bit-by-bit path, which is what x86-64 without BMI2
  and Arm64 run, is pinned by a test but not timed. AMD cores before Zen 3,
  which run `pext` in microcode, are not measured, and no gain is claimed
  on them.
- **No asymptotic claim.**

## Follow-up after measurement: lint only

After the suite ran, `cargo clippy -D warnings` flagged five manual
`% k == 0` tests in the candidate:

- one in `ColMap::without_deletions`;
- four in the new unit tests.

They were rewritten as `is_multiple_of`. The intended source for merging was therefore
`results/run-1/final.diff` (SHA-256 in `final.diff.sha256`), not the
measured `candidate.diff`. The only difference outside tests is that one
expression.

A build of that source gives byte-identical outputs to the measured
candidate on L13 seq 25 and L19 seq 24. Its instruction counts on L13 seq
25 and 24 are 0.9983× and 1.0016× the measured candidate's
(`lint-followup.tsv`), which is code-layout noise from the rebuild. The
suite was not re-run on it.

## Superseded (2026-10-11)

The candidate was held from commit because the Conductor reservation that
`CLAUDE.md` requires could not be taken in the measuring container. Before it
could land, `main` merged two changes to the same code:

- **#1646** ("word-wise inherited-F4 rewrite"). Every layout map, composed
  maps and folds included, gets `pext`/`pdep` segment plans. That is a
  superset of this round's change, which moved only consecutive deletion
  words with `pext` and masked the rest.
- **#1671** ("rewrite kept rows eagerly through the step plan").

On 2026-10-11 `ecbench` was built from `main` at `7666b708e` (SHA-256
`6ae1a2a1…`) and measured on three of this round's frozen inputs, with
Callgrind total `Ir`, whole process (`results/main-7666b708-check.tsv`):

| input | round-2 baseline (`a3652edd`) | round-2 candidate | `main` `7666b708e` | main / candidate |
| --- | --: | --: | --: | --: |
| L13 seq 25 (inherited F4) | 2,222,337,308 | 2,185,880,355 | 2,286,390,878 | 1.046 |
| L19 seq 24 (F6-IC) | 33,959,343,192 | 30,866,373,759 | 27,799,877,737 | 0.901 |
| L19 seq 25 (inherited F4) | 52,004,620,860 | 48,778,989,164 | 46,233,505,868 | 0.948 |

Observations:

- **At m = 19, `main` runs fewer instructions than this candidate** on both
  inputs. Merging the candidate's remap over `main`'s would replace a more
  general mechanism with a narrower one, so the candidate is not merged. Its
  source is kept as `results/run-1/final.diff`.
- **At m = 13, on L13 seq 25, `main` runs more instructions** than both this
  round's baseline and its candidate. That is one input only.
- **`main`'s outputs are not identical to this round's baseline.** The
  counted `solver_ops` rise on all three inputs (L13 seq 25: 11,575,434 to
  15,715,620; L19 seq 24: 405,966,170 to 520,513,284; L19 seq 25:
  549,394,694 to 668,516,039). On L19 seq 24 the F6-IC gate's statistics also
  change (`geometric_refutations` 4,219 to 4,795; `infeasible_branches` 576
  to 0). Recovered scalars agree. So the `main` column is a comparison of
  instructions between different counted algorithms, not a gated A/B.
  Whether those count changes are intended is for #1646 and #1671's
  evidence to say; it is not checked here.

What this round still establishes: on its own frozen baseline, moving
consecutive deletion words with `pext` and masking deleted bits cut
instructions by the amounts above, with every output identical. The masked
fallback is pinned by tests in the unmerged diff. The rejected variants and
their numbers stand.

**Dashboards.** No scoreboard, progress-chart, leaderboard or browser
figure changes, since nothing ships and no IC/rho ratio moves. Graphs
checked and left unchanged:

- `docs/index-calculus-scoreboard.html`, including the
  `inherited-f4-overhead-20261005` panel;
- `docs/ic/progress-timeline.json`;
- `docs/ic-leaderboard.html`.

The run report with its graphs is
[`reports/run-20261007T0030Z/`](reports/run-20261007T0030Z/REPORT.pdf).
