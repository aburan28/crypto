# Round 0025's wall-clock anomaly, re-timed in isolation

Registered before the measured run. §5 lists what was already seen.

## 1. Question

Round 0025 ([RESULTS](../RESULTS.md#round-0025-a-lean-rho-and-the-pairtriple-switch))
measured both IC arms **above** the lean rho in instructions at every cell.
It also measured them **below** it in native process time at the four
smallest cells, where IC's time was 0.87–0.98 of rho's.

Two readings compete:

- **Contention.** The round's confirmation stage ran pinned to CPU 3 from
  06:29 to 06:50. From 06:32 to 06:38 the m = 4 audit binary was rebuilt,
  linted and smoke-tested on the same 4-CPU container.
- **A real per-instruction difference.** Rho executes fewer instructions per
  unit time than the IC arm, so the instruction ratio overstates rho's
  advantage in time.

## 2. What was already known before this registration

Round 0025's own replay stage, which ran from 06:50 to 07:11 and overlapped
none of that build work, shows the same pattern. Incumbent/rho medians at
`n13a0, n17a1, n19a0, n19a1` are 0.91, 0.91, 0.95 and 0.95, and 1.04 at
`n23a0`.

On one `n13a0` fixture, the worker's own timer reads 0.23 ms for the IC arm
and 0.45 ms for rho. By callgrind count, IC executes 599k instructions and rho
418k.

## 3. Method

- **Runner.** [`run.py`](run.py), through
  [`tools/isolated_bench.py`](../../../../tools/isolated_bench.py) (AGENTS.md
  §10):
  - one lock and one preflight, which refuses to start while other processes
    use more than 0.1 CPU or PSI avg10 is above 5;
  - every movable thread moved off CPU 3 for the whole batch, and the driver
    itself kept off it;
  - each run pinned to CPU 3;
  - per run: wall time, user and system time, context switches, faults, PSI,
    and the CPU time other processes used. A run is `contended` if other
    processes used more than 0.1 CPU on average during it. Contended runs are
    excluded and counted.
- **Fixtures.** Round 0025's confirmation cases at `n13a0, n17a1, n19a0,
  n19a1` (the anomaly) and `n23a0` (the control, above one in round 0025).
  That is 12 cases per cell, using each arm's own recorded `job.json`.
- **Arms.**
  - `incumbent` and `rho` run the round's `worker` in its two modes;
  - `switch` runs `source_candidates/switch/worker`;
  - `rho_aa` is rho a second time, as the A/A noise floor.
- **Checks.** Every run's answer is checked against the round's recorded
  native answer, and a mismatch aborts the run.
- **Order.** Per round and case, the arms follow a Williams balanced Latin
  square, so over four rows every arm occupies every position once and
  follows every other arm exactly once. One untimed warm-up of every arm per
  cell comes first.
- **Size.** 15 rounds: 5 cells × 12 cases × 4 arms × 15 = 3,600 timed runs.
- **Timers.**
  - `wall`: the whole process, spawn to reap.
  - `cpu`: user plus system time.
  - `worker`: the worker's own `elapsed_seconds`. It runs from curve
    construction to report assembly, and excludes process start, job input
    and exit.

## 4. Metric and registered readout ([`analyze.py`](analyze.py))

- **Per case:** the median over rounds of arm X, divided by the median over
  rounds of `rho`.
- **Per cell:** the geometric mean over the 12 cases, with a 95% bootstrap
  over cases (10,000 draws, seed 20260929).

The predictions:

- **P1 (the anomaly survives isolation).** Incumbent/rho `wall` has its upper
  95% limit below 1 at all four small cells.
- **P2 (it is in the arms' own work, diluted by common startup).** At all
  four small cells, the incumbent/rho point estimate is lower on the `worker`
  timer than on `wall`.
- **P3 (control).** At `n23a0`, incumbent/rho `wall` is not below rho: its
  upper limit is at or above 1.
- **P4 (noise floor).** The `rho_aa`/`rho` `wall` interval lies inside
  [0.92, 1.09] at every cell. A difference inside that band is not a result.
- **P5.** At most 5% of timed runs are contended.

How each outcome is read:

- **P1 fails and P4 passes:** the round-0025 anomaly was a measurement
  artefact of the conditions, and the round's native column is annotated as
  such.
- **P1 and P2 hold:** the anomaly is real. The rho walk runs fewer
  instructions per unit time than the IC arm at these sizes, so on this host
  instruction counts are not proportional to time across the two arms.
  - This does not change round 0025's decision. Its `rho_gate` fails on
    instructions at every cell.
  - It also cannot reach the larger cells. There the instruction gap
    (1.4–13.7×) exceeds any such IPC factor, and round 0025's native times
    are above one as well.

## 5. Disclosed before registration

- **Smoke runs, not data.** Two were run into scratch and are not part of
  the record:
  - one round with a plain rotation order, where `rho_aa`/`rho` read
    0.84–0.90 because rho_aa usually ran right after an identical rho run.
    That is why the Williams order is used;
  - four rounds with the Williams order, reading incumbent/rho `wall` 0.87,
    0.87, 0.87 and 0.89 at the four small cells, 1.015 [0.971, 1.063] at
    `n23a0`, and A/A within 0.95–1.09. P3 is worded after this smoke.
- **The registered run is a fresh one.** It uses 15 rounds and a new output
  directory, and no smoke run is pooled with it.

## 6. Scope

- **Host.** One host class: a 4-vCPU Intel Xeon @ 2.80 GHz cloud container.
  Its frequency and host neighbours are not controlled, and the A/A band
  measures that residue.
- **Round.** One round's binaries and fixtures.
- **Not a new tournament result.** Round 0025's decision does not change
  either way.
