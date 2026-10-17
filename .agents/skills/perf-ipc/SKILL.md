---
name: perf-ipc
description: Measure cycles, instructions and instructions per cycle (IPC) of a benchmark with Linux kernel `perf stat` hardware counters, with the run's host, isolation and repeat spread recorded as JSON, and per-unit costs (cycles or instructions per walk step, relation, trial). Use whenever a speed result, a rho or index-calculus step cost, or a "where does the time go" question needs cycles or IPC, and before quoting any IPC figure. Detects hosts without a hardware PMU (most cloud VMs) and refuses to estimate there.
---

# Cycles, instructions and IPC with `perf stat`

`tools/perf_ipc.sh` runs one command N times under `perf stat`. It counts
user-space `cycles`, `instructions`, `branches` and `branch-misses`, plus
`task-clock`, and writes a `perf-ipc/v1` JSON record per run and a summary
(median, min, max). Read `AGENTS.md` §10 first: IPC is a property of a run on
one host, and it is evidence only under the same isolation, identity and A/A
rules as wall time.

## 1. Can this host count?

```sh
tools/perf_ipc.sh --repeat 1 -- true; echo $?
```

- **Exit 0.** It counted; continue.
- **Exit 3, `"status": "unsupported"`.** The host exposes no hardware PMU.
  The record says why, with the host's event sources. Typical case: a cloud VM
  whose `/sys/bus/event_source/devices/` lists no `cpu`, `cpu_core`,
  `cpu_atom` or `armv8_pmuv3*` source, so `perf stat -e cycles:u` prints
  `<not supported>`. Claude Code's cloud containers are such hosts. **Do not
  estimate IPC there**, and do not derive "cycles" from wall time × nominal
  GHz. Either:
  - run on a host with counters: bare metal, or the isolab runner
    (`ecbench-independent-runner` skill) where its PMU is exposed; or
  - report instruction counts from Callgrind (`valgrind --tool=callgrind`),
    which are deterministic but are not cycles, and state that IPC is
    unavailable on this host.
- **"No working perf binary."** On Ubuntu, install with
  `apt-get install linux-tools-common linux-tools-generic`. When no package
  matches the running kernel, `/usr/bin/perf` refuses to run; the script then
  uses `/usr/lib/linux-tools-*/perf` directly, which works across kernel
  versions. Set `PERF=` to choose a binary.
- **`perf_event_paranoid` errors.** Counting user-space events of your own
  process needs `/proc/sys/kernel/perf_event_paranoid` ≤ 2. The script counts
  `:u` events only, so it needs no root at 2. Levels above 2 forbid it:
  Debian's 3, and the 4 that GitHub's Ubuntu runners set, where even software
  events fail. On a machine you administer, `sudo sysctl -w
  kernel.perf_event_paranoid=2` (as `.github/workflows/rho-formal.yml` does);
  elsewhere ask the host owner. Do not run benchmarks as root to get around
  it.

`tools/perf_ipc.sh --selftest -- CMD` exercises the script's parsing and JSON
with software events on any host. It labels its records `selftest`, and its
IPC is always null. Use it to check the tool, never as a measurement.

## 2. Measure

Build the release binary first (never time a debug build), then:

```sh
tools/perf_ipc.sh --isolate --repeat 5 --label "strong rho n53 4 fixtures" \
  --steps 3034474 --out research/<study>/runs/perf/strong-rho-n53.jsonl \
  -- env RAYON_NUM_THREADS=1 target/release/examples/koblitz_rho_fixture 53 0 signed_frobenius 4 strong
```

- **`--isolate`.** Runs every repeat under `isolated_bench run --cpus 3`: the
  lock, the reserved core and the contention record (AGENTS.md §10; build
  `target/release/isolated_bench` or set `ISOLATED_BENCH`). The isolation
  records land beside `--out` as `….isolation.jsonl`. Use it for anything you
  report. Unisolated counts are exploratory.
- **`--steps K`.** The command's unit of work: walk steps, relations, trials,
  or field operations. The summary then gives cycles and instructions per
  unit. Take K from the command's own output (for the rho fixture, the sum of
  `walk_steps`), never from an estimate. If K varies by run, measure one run
  per invocation.
- **`--events`.** Adds events, e.g. `L1-dcache-load-misses,LLC-load-misses`
  or `cycle_activity.stalls_total`. More events than the PMU has counters
  forces multiplexing. Then perf scales the counts, the record lists them
  under `multiplexed`, and the script exits 5. Reported numbers must come
  from runs with no multiplexing: drop events, or measure them in a separate
  pass.
- **Fixed workload.** Pin the thread count (`RAYON_NUM_THREADS=1`) and the
  seeds. Check that the command's outputs (digests, recovered scalars,
  counted units) are identical across repeats and arms before reading any
  counter.

### Under Claude Code

With the host checked (step 1) and the command, repeats, `--steps` and
`--out` fixed, the main session may dispatch `experiment-director` to write a
command card for the `tools/perf_ipc.sh` invocations, relay it verbatim to
`executor-mechanical`, and dispatch `experiment-director` again to assess the
run against the JSON lines on disk (AGENTS.md "Model tiers"). Reading the
numbers (step 4) stays with the reasoning tier.

## 3. Compare two builds

Same as AGENTS.md §10, with counters as the measure:

1. A/A: the baseline binary twice, with `--repeat 5` each. The spread of
   `cycles` and `ipc` is the noise floor.
2. A/B interleaved: baseline, candidate, baseline, … for at least five rounds,
   one `--repeat 1` invocation per slot, each with its own `--label`.
3. Report medians and minima of cycles per unit, instructions per unit and IPC
   per arm, with the A/A spread. A difference inside the A/A spread is not a
   result.

## 4. Read the numbers honestly

- **Speed is cycles per unit, not IPC.** A change can raise IPC by adding cheap
  instructions while running slower. Report cycles per unit and instructions
  per unit; IPC explains them, it does not replace them.
- **`instructions:u` vs Callgrind `Ir`.** Both count user-space instructions.
  They differ slightly (perf counts what the core retires, including in
  libraries, with a small start-up share). Say which one a number is.
- **Frequency.** `cycles` counts core clock cycles, so cycles per unit does not
  depend on turbo or frequency scaling, but wall time does. On hybrid CPUs
  (`cpu_core` and `cpu_atom` sources), pin to one core type and say which.
- **Virtual machines.** A guest with a virtual PMU counts its own vCPU's
  cycles, but neighbours share caches and memory bandwidth. Record the host
  class, and do not compare IPC across hosts or core types.
- **What a value means.**
  - Below about 1: usually memory-, latency- or branch-bound.
  - Above about 3: rare outside tight vector or ALU loops.
  - The branch-miss rate and the extra cache events show which; profile with
    `perf record` (same host requirements) before guessing.

## 5. Record it

- Keep the JSON lines with the study's runs: `perf-ipc/v1` run and summary
  records, plus the isolation records.
- The report states:
  - the host block from the summary: CPU model, kernel, event sources,
    `perf_event_paranoid` and perf version;
  - the binary hash and commit;
  - the exact command, repeats and isolation;
  - the per-unit values with spread;
  - any `unsupported`, `multiplexed` or failed runs.
- Performance results also go to `benchmarks/` (AGENTS.md "Performance results
  also land in `benchmarks/`", `research-visuals` skill).
- An `unsupported` record is a result too: commit it, so the report shows
  which host could not count and why.
