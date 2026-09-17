# Round 0010 pre-registration: a lean job and a faithful clock

Written before any measured stage of this round was run.

## Why a further round

Round 0009 left the single-target panel with the winner at 0.55 of rho's
instructions and 0.93–0.94 of its native process time, but one cell, `n13a0`,
measured 1.017 on native time on confirmation (0.960 on replay), so the strict
per-cell native gate was not met. The archived receipts show why. On that cell
the two arms differ by about 0.2 ms inside the process (the incumbent's solve
0.47 ms against rho's 0.67 ms, `elapsed_seconds` in the worker reports), while
the measured cold process took 5.8 ms with a spread from 4.9 to 9.2 ms. The
per-case ratio was above one in 7 of 12 confirmation cases and in 3 of 12 on
replay: the cell's verdict was decided by noise that neither arm produces.

Where that time goes, measured on this host before this round: the evaluator
created every child with a `preexec_fn`, which makes CPython `fork()` the
evaluator process instead of using `vfork`/`posix_spawn`; the wall then charged
a copy of the evaluator's page tables to every job. Re-creating the evaluator's
`execute` with a 400 MiB parent gave a median of 23.7 ms per job (spread
6 ms) against 2.9 ms for the same child spawned without the fork. The worker
itself then spent about 0.4 ms in the dynamic loader (three shared libraries,
197 relocations, 18 `mmap`s) and the standard runtime's start-up (signal
handlers, an alternate signal stack, a read of `/proc/self/maps`), and the
IC path alone paid 73–97 µs for the standard library's CPU-feature cache on
its first use, against 10 µs for one CPUID instruction on this virtual
machine.

## What changes, and for whom

**Evaluator** (`tournament.py`, frozen into this round): `execute` spawns the
child without a `preexec_fn`. The calling thread is pinned to the measured CPU
before the spawn, so the child inherits the affinity; the memory and core caps
are applied with `prlimit` while the child is still blocked on `stdin`, that
is, before the job is delivered and before the child can allocate anything for
it; the watchdog thread is started before the timing window opens. The
protocol remains a blocking reap of the complete cold process wall. The
contract records the spawn method under `native_timing_protocol` and every
receipt under `native_process.spawn`. `snapshot_build` now also snapshots and
seals a `.cargo/config.toml` found in the source root, and locates the built
worker under either cargo layout. Tests cover the caps, the pin, the session,
the timeout and the artifact lookup.

**Baseline source** (`--source-root`, [patch](round10-lean.patch) on the
round-0009 incumbent source), shared by every arm including rho:

- `.cargo/config.toml`: `target-feature=+crt-static` with an explicit build
  target, so the worker is a static-pie executable with no dynamic loader work.
- The worker uses the C entry point (`#![no_main]`) instead of the standard
  runtime's start-up; it reads one JSON document and writes one, and exits with
  the same codes as before.
- `koblitz_tiny_ic::has_pclmul` issues one `CPUID` (leaf 1, ECX bit 1) instead
  of the standard library's feature cache. This is the only change that is
  paid by the IC arm alone.

Nothing about the algorithms changes: the same factor base, table rule,
arithmetic, linear algebra, descent and certificate; the same rho; the same
general-arithmetic final check.

Development evidence before freezing (not a claim): the new executable
reproduces, byte for byte apart from `elapsed_seconds`, all 304 frozen
fixture outputs of round 0009 (152 IC jobs and 152 rho jobs over every stage)
and every output passes the independent checker; under Callgrind on
`n13a0-000` the shared `startup_and_input` phase fell from 385K to 84K
instructions for both arms (IC total 1.87M against 2.18M, rho 3.09M against
3.40M), with every other phase within 1%. Under the new spawn, on two
confirmation fixtures per cell with ten repetitions each, the round-0009
worker's cold process took 2.3–3.2 ms (IC) and 2.5–4.0 ms (rho) with spreads
of 0.1–0.2 ms, and the new executable 1.9–2.8 ms and 2.1–3.5 ms; the IC/rho
ratio was 0.81–0.93 on every cell for both executables.

## Objective, arms and gates

`--objective rho`, as rounds 0007–0009. Incumbent: the round-0009 winner
configuration (`batch_trials: 1`) on the baseline above. Challengers:
`lean_batch4`, the configuration control (`batch_trials: 4`); and
`lean_stdprobe`, an ablation control built from the baseline with the standard
library's feature cache restored ([patch](round10-lean_stdprobe.patch)),
identical in everything else. Both are expected to be *retained*: the
configuration control as in every round since 0006, and the ablation because
its instruction count can differ from the incumbent's only by the probe, far
inside the no-regression gate; its native time against the incumbent
attributes the probe's cost. The decision then records the incumbent's own
`winner_over_rho`, `beats_rho_strict` and `rho_parity`, which is the
measurement this round exists for: every upper paired 95% limit and every
cell below one in both instructions and native process wall, on confirmation
and replay. A promotion of either control would be reported as such.

## Parent, seed, budget

Fresh seed 2026091610; target count 1; pilot profile; 1,800 paired-job
budget; one pinned CPU; 8 GiB cap; 60-second watchdog; Valgrind 3.22.0 `Ir`;
native progress recorded. Rho is the shipped per-target signed-Frobenius
solver on the incumbent's executable.

## Boundary, floor, class, honesty

Unit and boundary unchanged: the whole process, from `execve` to exit, in
guest instructions and in cold wall. Base support, `m = 3`, no direct
relations, every verification obligation, the worker's phase dumps and
report fields are unchanged. What moved is shared process cost and the
instrument that clocks it; the rho solve itself is the shipped implementation,
untouched. Class: engineering. Absolute times and instruction counts are not
comparable with earlier rounds, whose native walls included the evaluator's
fork and whose instruction counts included the dynamic loader; only
within-round ratios are claimed, as before. The native ratio is still bounded
below by the process creation that both arms pay; no arithmetic-complexity,
family-wide or cryptographic-size claim follows. Fresh fixtures from the new
seed; every failure retained; panels stay separate.
