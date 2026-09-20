# gpu/ecc — GPU kernels for 256-bit elliptic-curve arithmetic

CUDA kernels for prime-field ECC, aimed at the three things a GPU is good
for here: batch scalar multiplication, massively parallel Pollard rho for
the ECDLP, and a parallel baby-step giant-step engine for full-group and
interval logs. Written for secp256k1 with a fast path for its special prime,
but the field layer is generic over any odd modulus below 2^256.

Companion documents:

- [`OPTIMIZATION_BLACKWELL.md`](./OPTIMIZATION_BLACKWELL.md) — the cost
  model, static measurements on `sm_90`/`sm_100`/`sm_120`, and what to tune.
- [`../../docs/ecc_fpga_cost_model.md`](../../docs/ecc_fpga_cost_model.md) —
  how this compares with the VHDL engine in `hdl/ecc/`.

**No kernel here has been run on a GPU.** The environment this was written
in had no NVIDIA device. What *is* verified is stronger than it sounds
though — see "How this is tested" below.

## Layout

| File | What it is |
|---|---|
| `ecref.py` | Pure-Python oracle: curve arithmetic, the walk definition, and the generator for every header and vector file |
| `fp256.cuh` | 256-bit field: Montgomery CIOS, secp256k1 special reduction, inversion, batch inversion |
| `point.cuh` | Short-Weierstrass points: Jacobian add/double/mixed-add, scalar multiplication, batch normalisation |
| `rho.cuh` | The r-adding walk: partitioning, negation map, fruitless-cycle escape, batched stepping |
| `rho_host.hpp` | Host side: jump-table construction, walk replay, collision to discrete log |
| `kernels.cuh` | The kernels and their launch structure |
| `bsgs.cuh` | Baby-step giant-step: parallel chains, batched stepping, the x-keyed lock-free table, candidate emission |
| `bsgs_host.hpp` | Host side of BSGS: the plan (m, stride, table, chain lengths), seed constants, candidate verification, CPU driver |
| `kernels_bsgs.cuh` | The BSGS kernels |
| `bench.cu` | Device driver: self-test against the host, microbenchmarks, rho and BSGS runners |
| `test_cpu.cpp` | Verification harness — compiles the `.cuh` headers with g++ |
| `test_bsgs.cpp` | BSGS verification harness — same idea, for `bsgs.cuh` |
| `ptx_stats.sh` | Static instruction/occupancy analysis with no GPU present |
| `ptx_asm_check.py` | Interprets the `FP_PTX` inline assembly and checks it against the portable path |
| `modal_app.py` | Runs `bench` on a rented GPU via Modal — the selftest, throughput, the launch-bounds sweep and the BSGS solve |

## Build and test

```bash
make test          # six CPU test suites, plus the inline-asm check
make ptxcheck      # just the inline-asm check
make bench         # CUDA benchmark binary (needs nvcc)
make bench ARCH=sm_100    # datacenter Blackwell; sm_120 for RTX 50-series
```

`make test` needs only Python 3 and a C++17 compiler. It runs six
configurations:

| Suite | Curve | Reduction | What it covers |
|---|---|---|---|
| `test_secp_fast` | secp256k1 | special | the production path |
| `test_secp_mont` | secp256k1 | Montgomery | the generic path on the same curve |
| `test_toy_mont` | 40-bit toy, a ≠ 0 | Montgomery | generic doubling, and an end-to-end DLP solve |
| `test_bsgs_secp_fast` | secp256k1 | special | BSGS table, steppers, a 2^22-wide interval log |
| `test_bsgs_secp_mont` | secp256k1 | Montgomery | the same on the generic path |
| `test_bsgs_toy` | 40-bit toy | Montgomery | the same, plus whole-group logs with the cost table below |

`make test` also runs `ptx_asm_check.py`, which covers the one thing the
C++ suites structurally cannot: the `FP_PTX` inline assembly is guarded on
`__CUDA_ARCH__`, so the host never executes it. The script parses those
`asm(...)` blocks out of `fp256.cuh`, interprets them, and checks all seven
against the portable branch of the same function. It establishes that the
carry chains and operand numbering are right; it says nothing about
register allocation or real device behaviour.

On a GPU, start with `./bench selftest`: it runs every kernel and compares
the results against the same host code the CPU suites verify, including the
full rho walk state after 64 batched iterations. Built with `-DFP_PTX=1`
that is also the differential test the assembly ultimately needs, and the
only one that closes the gap the script leaves open.

On a GPU, `./bench bsgs --wbits 44` builds the table for a 2^44-wide
interval on the device and solves a planted log in it, and `--targets 16`
amortises that table over sixteen of them; `selftest` also
covers the BSGS kernels against the host driver.

### Without a GPU of your own

`modal_app.py` rents one:

```bash
pip install modal
modal setup                       # interactive; opens a browser
# or, headless (CI, a container, an agent session):
export MODAL_TOKEN_ID=...  MODAL_TOKEN_SECRET=...   # modal.com/settings/tokens

ECC_GPU=H100 modal run modal_app.py::selftest   # both configurations vs the host
ECC_GPU=H100 modal run modal_app.py::bench      # does the 55% become throughput?
ECC_GPU=H100 modal run modal_app.py::tune       # sweep RHO_MIN_BLOCKS
ECC_GPU=H100 modal run modal_app.py::bsgs --wbits 48   # BSGS: is the giant phase memory-bound?
```

`selftest` is the one that matters: it builds `FP_PTX=0` and `FP_PTX=1`,
runs each against the host reference, and exits non-zero if either
disagrees. A green run is what licenses turning `FP_PTX` on by default for
that architecture. `ECC_GPU` accepts any Modal type — `T4`, `L4`, `L40S`,
`A100`, `H100`, `H200`, `B200`, `RTX-PRO-6000` — and the build targets the
matching `sm_`, so Hopper and both Blackwell variants can each be checked.

**This has not been run against a GPU** — no Modal credentials were
available where it was written. The app definition itself has been checked:
it imports cleanly under `modal` 1.5.5, the image chain and all three
functions construct, `selftest`/`bench`/`tune` register as entry points,
and the `make` lines it issues were verified with `make -n`. So the first
authenticated run is a shakedown of the kernels more than of the harness —
but it is still a first run, and the image build (which compiles the CPU
suites and `ptxcheck` before any GPU is touched) is where a surprise would
surface.

## How this is tested

The kernels are header-only and every function is `__host__ __device__`, so
the *same source* compiles for the CPU. The test harness therefore checks
the real arithmetic, not a model of it:

1. **Against Python.** `ecref.py` generates field vectors (including 0,
   p-1, and 1·(p-1) edge cases), point vectors (P + (−P), P + O, n·P, and
   scalar 2^256−1), and 200-step walk traces with a per-step record of which
   table entry was used and whether the negation map fired. The C++ must
   match all of it exactly.
2. **Against itself, three ways.** `rho_step_batch<W>` (batched inversion),
   `rho_step_batch_lowmem<W>` (batched, minimal per-thread scratch) and
   `rho_step_thread_ref` (one inversion per walk) must produce bit-identical
   state for 400 iterations, and the same multiset of distinguished points.
3. **End to end.** On the toy curve the harness runs the complete pipeline —
   walks, distinguished points, host replay, collision, linear solve — and
   recovers a known discrete logarithm, verifying `k·P == Q` before
   declaring success. It reports the step count against the theoretical
   sqrt(pi·n/2) or sqrt(pi·n/4), which is how the negation map's benefit was
   confirmed rather than assumed:

```
[solve] toy DLP, neg_map=0, 128 walks
  solved in 1052288 steps (1.13x sqrt(pi n / 2)), 1042 DPs
[solve] toy DLP, neg_map=1, 128 walks
  solved in 630784 steps (0.96x sqrt(pi n / 4)), 588 DPs
```

## The walk

The rho engine is a distinguished-point r-adding walk. Its definition lives
in one place conceptually and three places in code (`ecref.py`, the batched
kernel, the host replay), which is exactly the sort of thing that drifts, so
the three are cross-checked against each other by the tests above.

```
partition(P)   = limb0(x) & (R - 1)
step           P <- P + M[partition(P)]        M[j] = c_j P + d_j Q
negation map   if y > (p-1)/2:  y <- -y
distinguished  ((limb0(x) >> 8) & dp_mask) == 0
```

Two details worth knowing:

**The hash is over the internal representation.** `x` above is whatever the
field layer stores — Montgomery form for a generic curve, a canonical
integer for secp256k1. Hashing that directly saves a conversion per step and
is still a deterministic function of the point, which is all a random walk
needs.

**Fruitless cycles are handled explicitly.** The negation map buys a
sqrt(2) reduction in steps, but it also lets the walk fall into a 2-cycle
with probability about 1/(2R) per step — with a small table, essentially
every walk is trapped within a few hundred steps and the search stalls
completely. The kernel detects the cycle (the new point's x matches the
point from two steps back) and escapes by *doubling the cycle's canonical
element*, the member with the lexicographically smaller x. Escaping from a
canonical element rather than from wherever the walk noticed is what keeps
the walk a deterministic function of the point, so two walks that enter the
same cycle at different members still leave it identically and their
collision survives. Longer cycles are left to the `max_steps` abort, so keep
`r_bits >= 8` with the negation map on.

Walks are identified by (index, restart counter) and their start point is
derived from that pair by a fixed PRNG, so a distinguished point is reported
as just `{x, walk, restart, steps}` — 48 bytes — and the host recovers the
coefficients by replaying only the two walks that actually collided.

## Kernels

| Kernel | Shape |
|---|---|
| `k_scalar_mul` | one thread per scalar, 4-bit window, optional branch-free table select |
| `k_scalar_mul_base` | fixed base point, for key generation and table building |
| `k_point_add` | one thread per pair, mixed Jacobian + affine |
| `k_to_affine<PER>` | `PER` points normalised per thread with one inversion |
| `k_rho_walk<W>` | `W` walks per thread, all state live across the step |
| `k_rho_walk_lowmem<W>` | same walk, keeps only the `W` prefix products |
| `k_rho_walk_ref` | one inversion per walk — the baseline to beat |
| `k_bench_*` | isolated field operations, dependent and independent chains |

The two batched variants differ only in memory strategy and are required by
the tests to produce identical state. `k_rho_walk_lowmem` re-reads each
walk's point in the backward pass of Montgomery's trick instead of holding
it, trading one coalesced load for a 5x smaller stack frame — measured 792
bytes against 4008 before the frame optimisations described in the Blackwell
document, and 792 against 1616 after.

## Baby-step giant-step

`bsgs.cuh` solves `Q = xG` for `x ∈ [x0, x0 + width)`: the whole group when
`x0 = 0, width = n`, an interval otherwise (a known-range key, or one
Pohlig–Hellman sub-problem). It is the deterministic `√width` method, and
what it buys over rho is paid in memory: the baby table costs 16 bytes per
entry, so a 2^44-wide interval needs 2^21 entries (32 MB) and a 2^66-wide
one needs 2^32 (64 GB, the top of one 80 GB device). Past that the method
is out of memory, not out of time, and rho or kangaroo take over. This is
the boundary the table below is measured against.

### The layout

```
baby table    { hash(x(jG)) -> j : 1 <= j < m }
giant walk    P_i = Q' - i*S,  Q' = Q - x0*G,  S = M*G,  i = 0, 1, ...
hit           x(P_i) == x(jG)   =>   x = x0 + i*M +- j   (host verifies both signs)
```

The table is keyed by the x-coordinate alone, so `jG` and `−jG` share one
entry: with `neg_map` the giant stride is `M = 2m − 1` and `m` entries cover
`2m − 1` residues. Because the giant phase stops at the first verified hit,
its expected cost is half its stride count, which sets the balance:

| layout | baby `m` | stride `M` | expected cost (random target, cold table) |
|---|---|---|---|
| `neg_map=1` | `√width / 2` | `2m − 1` | `m + width/(2M) ≈ 1.00 √width` |
| `neg_map=0` (textbook, the baseline) | `√(width/2)` | `m` | `m + width/(2m) ≈ 1.41 √width` |

Both are rows of the Galbraith–Wang–Zhang table that the Rust
`cryptanalysis::ecdlp_variants::bsgs` module implements sequentially (#7 and
#2 there); this is the same arithmetic laid out for a machine with a hundred
thousand threads and one memory system.

The same structure is also ported to the CPU, in Rust, as
[`src/cryptanalysis/bsgs_fast.rs`](../../src/cryptanalysis/bsgs_fast.rs):
single-word Montgomery arithmetic, the same x-keyed flat table, the same
batched inversion, chains handed to `rayon` instead of to warps.  It runs
on the same toy40 curve as the tests here, so its numbers and these sit in
one table.  Measured there: 37-49 Msteps/s on four threads, `S = 1.03`
cold and `0.55` amortised, and **6.3x to 10.9x** the general `BigUint`
implementation on the same instances (`cargo run --release --example
bsgs_fast_bench`).

### What maps onto the GPU

- **Independent chains.** Both phases are sets of chains that each add one
  constant point per step, `G` for the baby chains and `−S` for the giant
  chains. Chain `c` owns the index range `[cL, (c+1)L)`, is seeded with one
  short scalar multiplication (`c · LG`, a double-and-add over the bit
  length of `c`, not 256 doublings), and from then on every step is one
  affine addition. A thread runs `W` chains and shares one field inversion
  across them with Montgomery's trick, so a step costs `~6 + 270/W`
  multiplications — the same trade the rho kernel makes.
- **Register-resident state.** Nothing but its own thread touches a chain,
  so `k_bsgs_run<W>` loads the `W` points once, runs the whole launch's
  `iters` steps, and stores once; the rho walk, whose state is shared with
  the host's replay, must write every step. `ptxas` reports 94 registers
  and no spills at `W = 8` on `sm_90` and `sm_100` (168 for the rho walk).
- **A lock-free x-keyed table.** Open addressing with linear probing over
  8-byte slots: a 32-bit tag (the high half of a 64-bit hash of `x` in its
  internal representation, so no Montgomery conversion per step) and the
  32-bit index `j`. Insertion is one 64-bit `atomicCAS` per probe; lookups
  probe until an empty slot, so nothing is missed, and at load ≤ 1/2 a
  probe touches 2.5 slots on average — one 32-byte sector. A tag match is a
  *candidate*, and a false one occurs with probability `~load / 2^32` per
  probe; the host verifies every candidate with a scalar multiplication,
  so a false positive costs time and never correctness. (Zero false
  candidates in every run below.)
- **Grouped probes.** A table read is a dependent random access to global
  memory, hundreds of cycles against a few dozen multiplications of
  arithmetic between reads. The giant stepper issues all `W` first-slot
  loads of a thread before examining any of them so the memory system
  overlaps them; this is what the kernel's throughput will turn on, and it
  is the one thing the CPU suites cannot measure.
- **Exceptional cases handled, not assumed away.** A chain that starts at
  `O` (baby chain 0), one whose point equals the addend (`1G + G`, or a giant
  chain that lands on `S`), or its negative (the chain passes through `O`,
  which *is* the `j = 0` hit and is reported as one): the phase-A/phase-B
  split substitutes a unit denominator for those steps so the shared
  inversion stays valid, and the test suite drives chains through every
  one of them.
- **Early exit.** One flag in device memory, raised by the chain that
  emits a candidate and polled by every chain once per iteration (an
  L2-resident load).  Without it a launch always runs all its iterations,
  so a hit in the first one still costs `iters x chains` steps.
- **Table reuse, as a mode.** The table depends only on `G`, `x0`'s width
  and `m`, so many targets in the same interval share one build:
  `./bench bsgs --targets N` pays for it once and every further target
  costs only its giant phase.  Per-target cost falls from `~1.0` to the
  giant phase's own `~0.5` -- measured at **0.56 over 12 targets**
  against **2.13** for a single cold solve on the toy curve.
- **`W` is a free parameter.** `ptxas` gives every `W` from 4 to 32 the
  same 94 registers and the same 640 threads/SM, so the trade is purely
  arithmetic against per-thread stack: 39.8 multiplies per step at
  `W = 8` with a 1320 B frame, 22.9 at `W = 16` with 2640 B, 14.4 at
  `W = 32` with 5280 B.  The default is 16; which point wins depends on
  whether local-memory traffic or arithmetic binds, and only a device
  settles that (`--w`).

### How it is tested

`test_bsgs.cpp` compiles the same headers with g++ and checks, in order:
the hash table (insert, look-up, false positives on absent keys, refusal
when full); `bsgs_run_batch<W>` against the one-inversion-per-chain
reference for bit-identical chain state, table contents and candidate
lists, through the exceptional starts above; that every `j ∈ [1, m)` is
found from `x(jG)` computed independently by `scalar_mul`; that *every*
`x` in intervals of width 1, 2, 3, 7, 16, 61 and 200 is recovered under
both layouts (the index arithmetic at the seams); whole-group logs on the
toy curve at the seams of the layout (`x = 0, 1, m−1, m, M−1, M, M+1, n−1,
n−m`) and at random; and a 2^22-wide interval at a random 256-bit offset
on the compiled curve, so the interval path runs on secp256k1 itself. Every
recovered `k` is checked as `kG == Q`.

`./bench selftest` then compares the device against the same host driver:
the baby table as a set (insertion order across threads is not
deterministic under `atomicCAS`, so slot positions may differ; the entries
may not), the giant chain state after a launch bit for bit, the candidate
list, the reference kernel against the batched one, and finally a planted
interval log recovered end to end. **This has not been run on a GPU.** The
device code has been compiled to PTX with clang and assembled by `ptxas`
for `sm_90` and `sm_100`, in both the portable and `FP_PTX=1` builds, which
rules out build errors and nothing more; `modal_app.py::selftest` and
`::bsgs` are the commands that close the gap.

### Cost, measured

Per the repository rule (`AGENTS.md`): boundary first, one table, one
unit. The unit is `S = group additions / √n` with every phase charged —
the baby steps, the giant steps until the verified hit, both seeds, and one
addition per candidate verification. The **reference** is Pollard rho as
measured by `test_toy_mont` on the same curve with the same accounting:
`1.13 × √(πn/2) = 1.42` without the negation map and `0.96 × √(πn/4) =
0.85` with it. The **floor** for a deterministic table method is the
`√width` of its own balance, which is not a bound rho respects: rho is
cheaper in operations *with* the negation map and pays nothing in memory,
so what BSGS offers is determinism, table reuse across targets, and a
`1.0` that is a mean, not a tail. All rows are whole-group logs on the
40-bit toy curve (`n = 649 523 094 257`) from `test_bsgs_toy`, CPU
emulation of the kernels, 128 chains, every answer verified as `kG == Q`.

| row | class | table (entries / √n) | inversions | mean `S`, 8 random targets | ratio to rho ref | ratio to own floor | correct |
|---|---|---|---|---|---|---|---|
| rho, negation map (reference, `test_toy_mont`) | — | 0 | ops / W | **0.85** | 1.00 | — | 1/1 |
| rho, no negation map (`test_toy_mont`) | — | 0 | ops / W | 1.42 | 1.67 | — | 1/1 |
| BSGS textbook layout, `neg_map=0`, reference stepper | baseline | 0.707 | = ops | 1.44 (predicted 1.41) | 1.69 | 1.02 vs `√2` | 21/21 |
| BSGS textbook layout, `neg_map=0`, batched `W = 16` | engineering | 0.707 | ops / 16 | 1.44 (bit-identical) | 1.69 | 1.02 vs `√2` | 21/21 |
| BSGS negation layout, `neg_map=1`, batched `W = 16` | `√2` of the layout | 0.500 | ops / 16 | **1.07** (predicted 1.00) | 1.26 | 1.07 vs `1.0` | 21/21 |
| … same, giant phase only — the limit a batch tends to | multi-target | 0.500 | ops / 16 | **0.56** (predicted 0.50) | 0.66 | — | 8/8 |
| … textbook layout, giant phase only | multi-target | 0.707 | ops / 16 | 0.72 (predicted 0.71) | 0.85 | — | 8/8 |

The last two rows are the giant phase alone, averaged over the same 8
random targets: that is what a further target costs once a table exists,
and so the limit per-target cost falls to as a batch grows.  `test_multi_target`
runs that as an actual batch — one table, 12 targets on a `2^26` interval —
and measures **0.560 per target against 2.131 for a single cold solve**,
which is the same limit reached from the other side.

"Correct" counts the 13 seam targets plus the 8 random ones per layout, each
verified as `kG == Q`; the 8 random rows are the only ones in the mean,
since a target at `x = 0` or one that wraps mod `n` is found in the first
launch and says nothing about cost. The spread on the random rows is the
uniform giant index: `S` ran from 0.55 to 1.46 under the negation layout
(0.83 to 2.06 textbook), which is the tail a deterministic method has and
rho's distribution does not. The giant phase used to overshoot by up to a full
launch (`iters × chains`) because the host only verifies candidates
between launches.  It no longer does: a device-side flag, raised the
moment a candidate is emitted and polled by every chain once per
iteration, bounds the work past the hit by one iteration across the grid.
A candidate is only *probably* the answer, so the host clears the flag and
relaunches when verification rejects one, and each chain resumes from its
own stored index -- nothing skipped, nothing repeated.  The test suite
pins this down by solving one instance at `iters` of 4, 64 and 512 and
requiring the same step count: **2944 at all three**.  The launch size is
therefore a free knob, set for launch overhead and the host round-trip
alone, which is why it now defaults to 256.

The chips follow §3 of the rule: the batched stepper is *engineering*
against the reference stepper (`S` unchanged, inversions cut by `W`), the
negation layout against the textbook one is the `√2` the table predicts,
and nothing here is an advance on rho's floor — the point is the memory
boundary, drawn above. Wall-clock is a practicality note only:
about 380 k steps/s per core on the CPU emulation, which says nothing
about the device.

## Curve support

`ecref.py params --curve <name>` emits the parameter header. `secp256k1` and
`toyNN` (a deterministically generated NN-bit curve with prime order and
a ≠ 0) are built in. Any other short-Weierstrass curve over a prime below
2^256 works by adding its parameters to `ecref.py`; only the `FP_FAST` path
is secp256k1-specific, and it falls back to Montgomery automatically.
