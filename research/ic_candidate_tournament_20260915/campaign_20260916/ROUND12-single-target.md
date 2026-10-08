# Round 0012 pre-registration: the C runtime and the heap out of the way

Written before any measured stage of this round was run.

## Where the job's cost sat after round 0011

With the general arithmetic in words, the archived round-0011 receipts put
the IC arm's median complete job at 2.14M instructions and 1.98 ms against
rho's 4.71M and 2.30 ms. Two costs that are not the job's own were now a
large share of both arms:

- **The C runtime's start-up.** The executable is a static-pie glibc image,
  and on this virtual machine a Rust program that reads stdin and exits
  takes 1.17 ms from `posix_spawn` to reap, against 0.26 ms for a static
  binary whose whole body is the `exit` system call: about 0.9 ms of every
  job, for both arms, is glibc relocating itself
  (`_dl_relocate_static_pie`, 59% of the `startup_and_input` phase),
  reading its tunables and probing the cache hierarchy with `CPUID`, which
  traps to the hypervisor. The same program linked statically against musl
  takes 0.40 ms as a static-pie and 0.36 ms without relocations.
- **The heap.** By function, the `final_verification` phase spent about 45%
  of its instructions in `malloc`, `free` and `calloc` (the general field
  element allocates on every operation), and the construction another
  share. The whole process allocates very little: on every round-0011
  fixture the worker asks for between 115 KB in 1,353 calls (IC, `n13a0`)
  and 321 KB in 5,700 calls (rho, `n23a0`), with at most 208 KB live at
  once, and the process lives a few milliseconds. A general-purpose
  allocator's bookkeeping is pure overhead for that profile; musl's is
  dearer still per call than glibc's, so moving to musl alone would
  *raise* the instruction count.

## What changes, and for whom

**Baseline source** `arena` (`--source-root`; [round12-musl.patch](round12-musl.patch)
and [round12-arena.patch](round12-arena.patch) on the round-0011 source),
shared by every arm including rho:

- `.cargo/config.toml` selects the `x86_64-unknown-linux-musl` target with
  `crt-static` and `relocation-model=static`: a statically linked executable
  with no relocations to apply at start-up.
- The worker installs a global allocator: a bump pointer over a 64 MiB
  zero region in the executable's `.bss` (nothing is touched until it is
  used), which frees only the most recent block, grows the most recent
  block in place and falls back to the system allocator if the region is
  ever exhausted. Three unit tests cover in-place growth and give-back,
  alignment and moving a block that cannot grow, and the fallback.

The library source (`src/`, `Cargo.toml`, `Cargo.lock`, `build.rs`) is byte
for byte the round-0011 library; only the worker and the cargo config
differ. Rho is the shipped solver, untouched; it runs in the same process
and gets the same start-up and the same allocator.

**Challenger** `arena_fastio` ([round12-fastio.patch](round12-fastio.patch)):
the IC-only serialisation change of round 0011, rebased on the arena
worker. Round 0011 measured it at 0.9651 [0.9614, 0.9692] of the incumbent's
instructions on both final stages but could not resolve it on the native
clock (upper limit 1.0086 on confirmation) under the 0.9 ms start-up floor
that the baseline has now removed.

**Ablation controls**, expected eliminated at selection, each carrying one
of the two baseline changes without the other: `arena_glibc` (the arena
worker with the round-0011 static-pie glibc link; isolates the link mode)
and `musl_sysalloc` (the musl executable on the C library's own allocator;
isolates the arena).

Development evidence before freezing (not a claim). Every one of the four
executables reproduces, byte for byte apart from `elapsed_seconds`, all 304
frozen fixture outputs of round 0011 (152 IC and 152 rho jobs), and every
output passes the independent checker. Under Callgrind on one confirmation
case per cell, the baseline's IC job fell to 0.72–0.83 of its round-0011
instructions and rho's to 0.77–0.88: `startup_and_input` 0.085M → 0.012M
on every cell, `final_verification` 0.427M → 0.248M on `n13a0` and
0.870M → 0.480M on `n23a0`, `curve_and_targets` 0.190M → 0.142M and
0.540M → 0.476M, `reporting_and_cleanup` 0.021M → 0.008M. Separately, the
arena on the glibc link gave 0.78–0.89 and musl on its own allocator
1.06–1.11 (1.51 on the `n19a0` IC case), which is why both are in the round
as controls. Native, on two confirmation fixtures per cell with twelve
repetitions under the evaluator's spawn, the baseline's IC job took
0.86–1.20 ms (from 1.68–2.08) and rho 1.03–1.73 ms (from 1.90–2.58), an
IC/rho ratio of 0.70–0.83 on every cell (from 0.81–0.88); the arena on the
glibc link measured the same wall as round 0011 (1.68–2.04 ms), and the musl
static-pie without the arena 1.10–1.56 ms. `arena_fastio` measured 0.689–0.797 of the
round-0011 IC instructions (0.95–0.96 of the baseline's) and, natively,
0.84–1.14 ms against the baseline's 0.87–1.17 ms on the same fixtures, a
difference the tournament's paired clock may or may not resolve. The
library's own test suite is unchanged by construction (identical `src/`), so the round-0011 record
stands: 2,425 pass and the same five pre-existing pair-table failures. The
worker's own tests pass on every tree (three arena tests; four with the
serialisation test on `arena_fastio`).

## Objective and gates

`--objective rho`, as rounds 0007–0011. Incumbent: the round-0011 winner
configuration (`batch_trials: 1`) on the baseline above. `arena_fastio` is
promoted only if it passes the no-regression gate against the incumbent
(instruction ratio at most 0.98 with the paired upper limit below one,
native upper limit below one, no cell more than 10% worse) and the strict
rho gate on confirmation and replay; a *retained* verdict is a fully
reported outcome. Both controls are expected to fail the no-regression gate
at selection; a control that survives it is reported as such. Either way the
decision records the winner's `winner_over_rho`, `beats_rho_strict` and
`rho_parity`.

## Parent, seed, budget

Fresh seed 2026091612; target count 1; pilot profile; 2,400 paired-job
budget (four arms where round 0011 had three); one pinned CPU; 8 GiB cap;
60-second watchdog; Valgrind 3.22.0 `Ir`; native progress recorded; the
round-0010 evaluator (spawn without a fork, caps before job delivery).

## Boundary, floor, class, honesty

Unit and boundary unchanged. Base support, `m = 3`, no direct relations,
every verification obligation, the worker's phase dumps and report fields
are unchanged; the final check is still `[d]G = Q` in the general
arithmetic on every recovered scalar. Class: engineering, and of the
process rather than the algorithm: the round moves the cost of starting a
process and of serving its heap, which both arms pay equally, so the native
ratio moves more than either arm's own solve. Absolute times and
instruction counts are comparable with rounds 0010 and 0011 (same
evaluator, same host, same compiler) and not with earlier rounds. No
arithmetic-complexity, family-wide or cryptographic-size claim follows.
Fresh fixtures from the new seed; every failure retained; panels stay
separate.
