# Round 0025's wall-clock anomaly, re-timed in isolation: results

This was run after [PREREGISTRATION.md](PREREGISTRATION.md) was committed, and that
file is unchanged. Raw records are in [runs/registered/runs.jsonl](runs/registered/runs.jsonl),
with the readout in [readout.txt](runs/registered/readout.txt) and throughput in
[throughput.txt](runs/registered/throughput.txt).

## Answer

**The anomaly is real, not contention.** With CPU 3 reserved for the whole batch, the
isolated re-timing reproduces it, and all five predictions hold. The cause is
throughput: the IC arm executes **1.65–2.03× more instructions per second than rho**
on this host, so the instruction ratio overstates rho's lead in time.

## Conditions

- 3,600 timed runs; every answer matched round 0025's recorded answer.
- 89 runs (2.5%) were marked `contended` and excluded. In 70 of them, the CPU was
  used by the agent's own `claude` process; the rest were migration and RCU kernel
  threads.
- Every kept run saw 0.000 s of other-process CPU. The median run had one
  involuntary context switch.
- The A/A control (`rho_aa`/`rho`) came in at 0.979–1.025, with every interval inside
  [0.949, 1.051].

## Wall time, arm / rho

Geometric mean over 12 cases, 95% bootstrap over cases:

| cell | incumbent / rho | switch / rho | A/A | round 0025 replay, incumbent / rho |
|:--|--:|--:|--:|--:|
| `n13a0` | **0.872** [0.862, 0.883] | 0.890 | 0.985 | 0.907 |
| `n17a1` | **0.921** [0.901, 0.941] | 0.926 | 1.025 | 0.910 |
| `n19a0` | **0.912** [0.890, 0.936] | 0.939 | 1.016 | 0.947 |
| `n19a1` | **0.940** [0.896, 0.989] | 0.934 | 1.025 | 0.951 |
| `n23a0` (control) | 0.985 [0.954, 1.021] | 1.033 | 0.979 | 1.044 |

- **P1 holds.** At all four small cells, the upper 95% limit is below 1.
- **P3 holds.** At `n23a0` the IC arm is not below rho.
- **Agreement with round 0025.** The replay-stage medians, which ran clear of the
  build overlap, are within about 4% of these.

## Why: throughput, not contention

The worker's own timer runs from curve construction to report assembly. On it,
incumbent/rho is 0.706, 0.738, 0.797 and 0.840 at the four small cells, below the
whole-process ratios, so **P2 holds**. That is the arms' own work. Process start and
exit are common to both arms, and they dilute it.

Over the same span, round 0025's callgrind receipts count **more** instructions for
the IC arm:

| cell | instructions, inc / rho | worker time, inc / rho | instructions per second, inc / rho |
|:--|--:|--:|--:|
| `n13a0` | 1.437 | 0.706 | 2.03 |
| `n17a1` | 1.385 | 0.738 | 1.88 |
| `n19a0` | 1.361 | 0.797 | 1.71 |
| `n19a1` | 1.386 | 0.840 | 1.65 |
| `n23a0` | 1.867 | 1.015 | 1.84 |

Rho also takes about 30 more minor page faults per process (89–96 against 56–71).

The walk being latency-bound is **a hypothesis this run does not test.** The walk is
one dependent chain of field operations per step, while the collector's scans are
independent iterations. No cycle counters are available in this container, so
instructions per cycle were not measured directly.

## What changes

- **Round 0025's decision does not change.** `rho_gate` fails on instructions at
  every cell. At the large cells, the instruction gap (up to 13.7×) exceeds any
  factor of about 2, and those cells' native times are above one too. At `n23a0` the
  ratio is already at parity.
- **Instruction counts, the tournament's primary unit, favour rho** by roughly 1.7–2×
  against wall time, between these two arms on this host. They remain the right
  contention-proof unit for *within-arm* changes. Across arms with different
  instruction mixes they are not a time proxy, and a claim about time needs the
  isolated native measurement beside them.
- **Round 0025's native column** at the four small cells is confirmed, not an
  artefact. It was measured without the isolation now required (AGENTS.md §10). This
  run supersedes it as the time evidence for those cells.

## Scope

- One host class: a 4-vCPU Intel Xeon @ 2.80 GHz cloud container.
- One round's binaries and fixtures.
- Five cells, `n ≤ 23`.
