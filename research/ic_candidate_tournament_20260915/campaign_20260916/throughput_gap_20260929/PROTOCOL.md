# Why does the IC arm run more instructions per second than rho? A cache and branch simulation

Written before the grid below runs. The one smoke run is disclosed in §4.

## Question

Two isolated re-timings measured the IC arm executing 1.65–2.03× more instructions per
second than the lean rho, at `n13a0, n17a1, n19a0, n19a1, n23a0`
([walltime_isolation_20260929](../walltime_isolation_20260929/RESULTS.md),
[walltime_strong_rho_20260929](../walltime_strong_rho_20260929/RESULTS.md)). A batched rho
did not close the gap. This container exposes no hardware counters, so this diagnostic uses
valgrind's simulators, whose counts are deterministic and unaffected by contention.

Three explanations compete:

- **H-branch.** Rho mispredicts many more branches per instruction.
- **H-memory.** Rho misses in D1 or the last-level cache much more per instruction.
- **H-latency.** Neither. Rho's instructions sit on longer dependency chains, which the
  simulators cannot see; the gap is what H-branch and H-memory leave unexplained.

## Method

- **Tool.** `valgrind --tool=cachegrind --cache-sim=yes --branch-sim=yes`, valgrind 3.22.0.
- **Binary.** Round 0025's `worker`.
- **Jobs.** Round 0025's confirmation `job.json` for cases `000`–`003` at each of the five
  cells, in both the `incumbent` (IC) and `rho` modes. That is 40 runs.
- **Recorded per run.** Instructions, I1/LL instruction misses, D1 and LL data misses, and
  conditional and indirect branches with their mispredicts. The top functions by
  mispredicts are also recorded.
- **Analysis.** For each arm and cell, rates per 1,000 instructions, summed over the four
  cases. A penalty model is fixed now and applied without refitting:
  - 15 cycles per mispredicted branch;
  - 12 cycles per D1 miss that hits LL;
  - 150 cycles per LL miss.

  From it we get the arms' predicted **excess cycles per instruction**. We compare the
  resulting predicted time ratio with the measured in-process ratio, assuming both arms
  retire their other instructions at one common rate. The model is a stated
  approximation for an attribution, not a measurement of cycles.

## Reading, fixed now

- **H-branch or H-memory explains the gap** if, at every cell, the penalty model with a
  common base rate reproduces the measured in-process IC/rho time ratio to within 15%.
- **H-latency is left standing** if the model predicts no more than half of the measured
  throughput gap.
- **Anything between** is reported as partial.
- **The penalties are not tuned to fit.** This file fixes them. A sensitivity line at half
  and double the penalties is reported beside the result.

## Disclosed before the grid

One smoke run: `rho` mode on `n13a0-000`.

- 419,441 instructions;
- 5,728 branch mispredicts, 8.6% of 66,275 branches;
- 3,664 D1 misses and 3,453 LL data misses, most of them writes.

The IC arm had not yet been simulated.

## Scope

- **Simulated, not measured.** The cache and branch predictor are valgrind's model of this
  host's cache sizes, not the Xeon's real predictor.
- **What it can find.** It attributes. It cannot prove H-latency; that is what is left over.
- **One host class, one binary, five cells.**
