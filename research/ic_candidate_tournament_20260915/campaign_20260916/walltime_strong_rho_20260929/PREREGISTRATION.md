# Round 0025's small-cell time lead, against the strongest rho on main

Registered before the measured run. §3 discloses the smoke runs.

## 1. Question

The [isolated re-timing](../walltime_isolation_20260929/RESULTS.md) found the IC arm
0.87–0.94 of the lean rho's native time at `n13a0, n17a1, n19a0, n19a1`. The reason is that
the IC arm executes 1.65–2.03× more instructions per second.

One hypothesis says the rho walk is latency-bound. If so, a rho that advances many walks
in lockstep and shares one inversion across them would recover that throughput and erase
the lead. Main's strongest such rho is
[`examples/koblitz_rho_batch_ks_strong.rs`](../../../../examples/koblitz_rho_batch_ks_strong.rs),
rung 3. It overturned the compact-orbit batch lead in #943 and #966.

**Does IC's small-cell lead survive against the strongest verified rho?**

## 2. Arms, all on CPU 3 under `tools/isolated_bench.py`, one reservation

| arm | binary | targets |
|:--|:--|:--|
| `incumbent` | round 0025's `worker`, IC mode | the case's recorded `job.json` |
| `lean_rho` | round 0025's `worker`, rho mode | the same |
| `strong32` | strong rho, rung 3, `KIC_RHO_LANES=32`, `KIC_RHO_DP_BITS=2` | one derived target per process, batch seed = case index + 1 |
| `strong1` | the same with `KIC_RHO_LANES=1` | the same |
| `strong32_aa` | `strong32` again (A/A) | the same |

- **Strong-rho binary.** Built at `main` `69173367`, with the first audit's pinned lock,
  under the benchmark lock. Its sha256 is
  `46d5b93fcaae38008e324456823310f819ad2fb1bad0e747a3291044c907ee24`.
- **Cases and order.** 12 confirmation cases per cell at the four small cells, plus
  `n23a0` as a control. Ten rounds follow a Williams square for five arms: each arm follows
  every other exactly twice and holds every position twice.
- **Checks.** Answers are checked on every run:
  - the worker's against the recorded solutions;
  - the strong rho's by its own `all_verified`.
- **Target pairing.** The strong rho draws its own targets, so its cases are paired with
  the worker's by index, not by target. The geometric mean over 12 cases averages over
  target luck.

## 3. Disclosed before registration (smoke, not data)

- **Default settings abort.** At `n13a0` and `n17a1`, the strong rho aborted at its
  per-target step cap with its default distinguished-point bits (8 at `n = 19`), and also
  at 8 or 1 lanes and at rung 2. With `KIC_RHO_DP_BITS=2` it completes; that is why 2 is
  registered for every cell.
- **In-process medians of 9 runs, strong rho:**
  - lanes 1: 0.40, 0.60, 0.66, 0.63 and 0.90 ms;
  - lanes 32: 0.49, 0.45, 0.72, 0.75 and 0.79 ms;

  at `n13a0, n17a1, n19a0, n19a1, n23a0`. The lean rho's worker-timed medians from the
  registered isolation run are 0.35, 0.39, 0.42, 0.45 and 0.49 ms.
- **Whole process at `n19a0`.** It was about 4.5–4.9 ms for the strong rho and 1.7–2.0 ms
  for the lean rho. The strong rho is a standard `std` binary. The tournament worker was
  built for minimal process start in round 0010.
- **One-round plumbing run of `run.py`.** It went into scratch, is not data, and was not
  pooled.
  - 300 timed runs, 8 contended.
  - The strongest rho in-process was `strong1` at `n13a0` and `strong32` at `n17a1`, by
    about 7%, and the lean rho at the other three cells.
  - IC/strongest in-process read 0.80, 0.85, 0.78, 0.78 and 0.98.
  - The A/A intervals were 0.90–1.20 on one round.
  - That is why P1 below is reported, not predicted.

## 4. Metric and registered readout ([`analyze.py`](analyze.py))

- **Primary metric: in-process time.**
  - For the worker, `elapsed_seconds`, which runs from curve construction to report
    assembly.
  - For the strong rho, `in_process_ms`, which includes its setup, normal-basis
    construction, walks and verification.
  - Whole-process wall is reported beside it. It is not primary, because at these sizes it
    measures process start-up, which differs by binary packaging, not by algorithm.
- **Strongest rho.** Per cell, the rho arm with the lowest geometric-mean in-process time.
- **Ratio.** IC over the strongest rho: per case, a median over rounds; per cell, the
  geometric mean over cases with a 95% bootstrap (10,000 draws, seed 20260929).

The predictions:

- **P1 (reported, not predicted).** Which rho arm is strongest at each cell. The plumbing
  run in §3 was seen before registration.
- **P2.** IC/strongest rho, in-process, has its upper limit below 1 at all four small cells.
  That is, the lead survives.
- **P3.** The `strong32` A/A interval lies inside [0.9, 1.1] at every cell.
- **P4.** At most 5% of timed runs are contended.

**Reading.**
- **If P2 fails** at a cell, IC's lead there does not survive the strongest rho.
- **If a strong-rho arm is strongest** at a cell, the round-0025 instruction and native
  comparisons, both made against the lean rho, overstate IC there by that margin.

## 5. Scope

- **Where.** One host class, one round's IC binary, five cells at `n ≤ 23`.
- **Not a new rho.** The strong rho is main's, unchanged, used at non-default
  distinguished-point bits because its default fails at these sizes.
- **No decision change.** Round 0025's decision is unaffected either way.
