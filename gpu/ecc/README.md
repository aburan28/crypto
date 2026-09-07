# gpu/ecc — GPU kernels for 256-bit elliptic-curve arithmetic

CUDA kernels for prime-field ECC, aimed at the two things a GPU is good for
here: batch scalar multiplication, and massively parallel Pollard rho for
the ECDLP. Written for secp256k1 with a fast path for its special prime,
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
| `bench.cu` | Device driver: self-test against the host, microbenchmarks, rho runner |
| `test_cpu.cpp` | Verification harness — compiles the `.cuh` headers with g++ |
| `ptx_stats.sh` | Static instruction/occupancy analysis with no GPU present |
| `ptx_asm_check.py` | Interprets the `FP_PTX` inline assembly and checks it against the portable path |
| `modal_app.py` | Runs `bench` on a rented GPU via Modal — the selftest, throughput, and the launch-bounds sweep |

## Build and test

```bash
make test          # three CPU test suites, plus the inline-asm check
make ptxcheck      # just the inline-asm check
make bench         # CUDA benchmark binary (needs nvcc)
make bench ARCH=sm_100    # datacenter Blackwell; sm_120 for RTX 50-series
```

`make test` needs only Python 3 and a C++17 compiler. It runs three
configurations:

| Suite | Curve | Reduction | What it covers |
|---|---|---|---|
| `test_secp_fast` | secp256k1 | special | the production path |
| `test_secp_mont` | secp256k1 | Montgomery | the generic path on the same curve |
| `test_toy_mont` | 40-bit toy, a ≠ 0 | Montgomery | generic doubling, and an end-to-end DLP solve |

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

### Without a GPU of your own

`modal_app.py` rents one:

```bash
pip install modal && modal setup
ECC_GPU=H100 modal run modal_app.py::selftest   # both configurations vs the host
ECC_GPU=H100 modal run modal_app.py::bench      # does the 55% become throughput?
ECC_GPU=H100 modal run modal_app.py::tune       # sweep RHO_MIN_BLOCKS
```

`selftest` is the one that matters: it builds `FP_PTX=0` and `FP_PTX=1`,
runs each against the host reference, and exits non-zero if either
disagrees. A green run is what licenses turning `FP_PTX` on by default for
that architecture. `ECC_GPU` accepts any Modal type — `T4`, `L4`, `L40S`,
`A100`, `H100`, `H200`, `B200`, `RTX-PRO-6000` — and the build targets the
matching `sm_`, so Hopper and both Blackwell variants can each be checked.

**This has not been run.** It is written against the same Modal conventions
as `ecc2k130/modal_app.py`, which has, but no Modal credentials existed in
the environment where it was written, so treat the first run as a shakedown
of the harness as much as of the kernels.

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

## Curve support

`ecref.py params --curve <name>` emits the parameter header. `secp256k1` and
`toyNN` (a deterministically generated NN-bit curve with prime order and
a ≠ 0) are built in. Any other short-Weierstrass curve over a prime below
2^256 works by adding its parameters to `ecref.py`; only the `FP_FAST` path
is secp256k1-specific, and it falls back to Montgomery automatically.
