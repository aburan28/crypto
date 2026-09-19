# gpu/ecc — GPU kernels for 256-bit elliptic-curve arithmetic

CUDA kernels for prime-field ECC, aimed at the two things a GPU is good for
here: batch scalar multiplication, and massively parallel Pollard rho for
the ECDLP. Written for secp256k1 — a fast path for its special prime, and
the walk folded by its full automorphism group Aut(E) = Z/6 — but the field
layer is generic over any odd modulus below 2^256 and the walk runs on any
short-Weierstrass curve.

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
| `ecref.py` | Pure-Python oracle: curve arithmetic, the walk and its folding, the Z/6 constants, and the generator for every header and vector file |
| `fp256.cuh` | 256-bit field: Montgomery CIOS, secp256k1 special reduction, inversion, batch inversion |
| `point.cuh` | Short-Weierstrass points: Jacobian add/double/mixed-add, scalar multiplication, batch normalisation, β |
| `rho.cuh` | The r-adding walk: partitioning, automorphism folding (negation map, Z/6 orbit), fruitless-cycle detection and escape, batched stepping |
| `rho_host.hpp` | Host side: jump-table construction, walk replay with λ bookkeeping, collision to discrete log |
| `kernels.cuh` | The kernels and their launch structure |
| `bench.cu` | Device driver: self-test against the host, microbenchmarks, rho runner |
| `test_cpu.cpp` | Verification harness — compiles the `.cuh` headers with g++ |
| `ptx_stats.sh` | Static instruction/occupancy analysis with no GPU present |
| `ptx_asm_check.py` | Interprets the `FP_PTX` inline assembly and checks it against the portable path |
| `modal_app.py` | Runs `bench` on a rented GPU via Modal — the selftest, throughput, and the launch-bounds sweep |

## Build and test

```bash
make test          # four CPU test suites, plus the inline-asm check
make ptxcheck      # just the inline-asm check
make bench         # CUDA benchmark binary (needs nvcc)
make bench ARCH=sm_100    # datacenter Blackwell; sm_120 for RTX 50-series
```

`make test` needs only Python 3 and a C++17 compiler. It runs four
configurations:

| Suite | Curve | Reduction | What it covers |
|---|---|---|---|
| `test_secp_fast` | secp256k1 | special | the production path, fold 6 included |
| `test_secp_mont` | secp256k1 | Montgomery | the generic path on the same curve |
| `test_toy_mont` | 40-bit toy, a ≠ 0 | Montgomery | generic doubling, and an end-to-end DLP solve |
| `test_toyj0_mont` | 36-bit toy, j = 0 | Montgomery | the Z/6-folded walk end to end, step counts averaged over 8 solves per fold |

`make test` also runs `ptx_asm_check.py`, which covers the one thing the
C++ suites structurally cannot: the `FP_PTX` inline assembly is guarded on
`__CUDA_ARCH__`, so the host never executes it. The script parses those
`asm(...)` blocks out of `fp256.cuh`, interprets them, and checks all seven
against the portable branch of the same function. It establishes that the
carry chains and operand numbering are right; it says nothing about
register allocation or real device behaviour.

On a GPU, start with `./bench selftest`: it runs every kernel and compares
the results against the same host code the CPU suites verify, including the
full rho walk state after 64 batched iterations of the fold-6 walk. Built
with `-DFP_PTX=1` that is also the differential test the assembly
ultimately needs, and the only one that closes the gap the script leaves
open.

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
```

`selftest` is the one that matters: it builds `FP_PTX=0` and `FP_PTX=1`,
runs each against the host reference, and exits non-zero if either
disagrees. A green run is what licenses turning `FP_PTX` on by default for
that architecture. `ECC_GPU` accepts any Modal type — `T4`, `L4`, `L40S`,
`A100`, `H100`, `H200`, `B200`, `RTX-PRO-6000` — and the build targets the
matching `sm_`, so Hopper and both Blackwell variants can each be checked.
`bench` takes `--fold 1|2|6` to compare the foldings on the device (the
binary defaults to 6 on secp256k1).

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
   scalar 2^256−1), and 200-step walk traces — one per folding level — with
   a per-step record of which table entry was used and which automorphism
   the canonicalisation applied. The C++ must match all of it exactly.
2. **The automorphism group.** On a j = 0 curve the harness checks β³ = 1,
   1 + β + β² = 0, λ³ = 1 mod n, and that `(βx, y)` really is `λ·P` for
   every point vector, by scalar multiplication. It then applies all six
   automorphisms to each point and requires the same representative back
   from every one, with an aut code that maps the point onto it.
3. **Against itself, three ways.** `rho_step_batch<W>` (batched inversion),
   `rho_step_batch_lowmem<W>` (batched, minimal per-thread scratch) and
   `rho_step_thread_ref` (one inversion per walk) must produce bit-identical
   state for hundreds of iterations, and the same multiset of distinguished
   points; every distinguished point must replay on the host to the reported
   x, with coefficients (a, b) that reproduce it as a·P + b·Q — which is what
   proves the λ bookkeeping on secp256k1 itself, where no solve can run.
4. **Fruitless cycles.** With a deliberately tiny table (R = 8), every
   cycle class is frequent: the detector's per-length counts are checked
   against the rates derived below (they agree to within about 10%, the
   accounting and counting noise at those sizes, on both folds),
   and a brute-force search finds cycles of every length up to 6 and starts
   a walk at *each* of their members, requiring every entry point to leave
   the cycle to the same exit point.
5. **End to end.** On the toy curves the harness runs the complete pipeline —
   walks, distinguished points, host replay, collision, linear solve — and
   recovers planted discrete logarithms, verifying `k·P == Q` before
   declaring success. It reports the step count against the theoretical
   sqrt(πn / (2·fold)), which is how the folding's benefit was confirmed
   rather than assumed. On the 36-bit j = 0 curve, 8 solves per fold with 16
   walks, R = 256 and one distinguished point per 2^7 steps (same code as
   the kernel, run on the CPU; steps, not time, are the unit):

| fold | folding | expected sqrt(πn / 2·fold) | mean steps, 8 solves | ratio | correct |
|---|---|---|---|---|---|
| 1 | none | 249,454 | 305,740 | 1.23× | 8/8 |
| 2 | negation map | 176,391 | 199,042 | 1.13× | 8/8 |
| 6 | full Aut(E) = Z/6 | 101,839 | 113,522 | 1.11× | 8/8 |

   Fold 1 : fold 6 is 2.69× against a predicted sqrt(6) = 2.45×, and
   fold 2 : fold 6 is 1.75× against sqrt(3) = 1.73×. A single rho run has a
   standard deviation of about half its mean, so the mean of 8 carries
   ±0.18× of noise; the ~10–20% excess over theory on all three rows is the
   usual r-adding-walk inefficiency plus the distinguished-point tail, and
   cancels in the ratios.

## The walk

The rho engine is a distinguished-point r-adding walk, folded by an
automorphism subgroup of the curve. Its definition lives in one place
conceptually and three places in code (`ecref.py`, the batched kernel, the
host replay), which is exactly the sort of thing that drifts, so the three
are cross-checked against each other by the tests above.

```
partition(P)   = limb0(x) & (R - 1)
step           P <- canonical(P + M[partition(P)])      M[j] = c_j P + d_j Q
distinguished  ((limb0(x) >> 8) & dp_mask) == 0

fold 1         canonical(P) = P
fold 2         negation map:  if y > (p-1)/2:  y <- -y
fold 6         beta orbit, then the negation map:
               x <- beta^k x for the k in {0, 1, 2} that minimises it
```

`fold` is the order of the automorphism subgroup the walk quotients out.
The walk becomes a function on E/⟨aut⟩ instead of E, so the expected number
of steps to a collision drops from sqrt(πn/2) to sqrt(πn / (2·fold)):
sqrt(2) for the negation map on any curve, sqrt(6) for the full Aut(E) of a
j = 0 curve such as secp256k1, whose endomorphism `(x, y) -> (βx, y)` is
multiplication by a cube root of unity λ mod n (Wiener–Zuccherato 1998;
Duursma–Gaudry–Morain 1999). Six is all the automorphisms a j = 0 curve
has, so sqrt(6) is the ceiling for this kind of folding. `ecref.py` derives
β and λ for any j = 0 curve of prime order (for secp256k1 they are the
constants libsecp256k1 ships, and the generator asserts that) and emits
them into the curve header; `CURVE_HAS_AUT6` is 0 on other curves and
`rho_fold_supported()` refuses fold 6 there.

Three details worth knowing:

**The hash is over the internal representation.** `x` above is whatever the
field layer stores — Montgomery form for a generic curve, a canonical
integer for secp256k1. Hashing and comparing that directly saves a
conversion per step and is still a deterministic function of the point,
which is all a random walk needs. Fold 6 costs one field multiplication per
step (βx; then β²x = −(x + βx), since 1 + β + β² = 0) plus two 256-bit
comparisons — about 10% on top of a batched step — for 1.73× fewer steps
than the negation map alone.

**The kernel forgets which automorphism it applied; the host recovers it.**
`rho_canonical` returns an *aut code* `(k << 1) | negated` naming the
automorphism, and the walk state carries nothing else. The host replay,
which re-runs the two colliding walks from their seeds, multiplies the
tracked coefficients (a, b) by `(−1)^negated · λ^k` at every step, so a
distinguished point is still reported as just `{x, walk, restart, steps}` —
48 bytes.

**Fruitless cycles are detected up to length 6 and escaped
deterministically.** Folding lets the walk close short cycles that carry
no information. Per step, with R table entries:

| length | mechanism | probability per step | fold 2, R = 256 | fold 6, R = 256 |
|---|---|---|---|---|
| 2 | same index twice, canonicalisation applied [−1] in between | 1 / (fold·R) | 2.0e-3 | 6.5e-4 |
| 3 | fold 6 only: same index three times, [ω] applied twice (1 + ω + ω² = 0) | 1 / (18 R²) | — | 8.5e-7 |
| 4 | indices j j' j j' with automorphisms (a, −1/a, a, −1/a) | (R−1) / (fold² R³) | 3.8e-6 | 4.2e-7 |
| 5, 6 | fold 6: an index pair under [−1] plus a triple under [ω]; fold 2: three pairs (even lengths only) | O(1/R³) | — | — |

A trapped walk emits no distinguished point until `max_steps` aborts it,
so with R = 256 and one DP per 2^20 steps, 4-cycles alone would trap most
negation-map walks — and 3-cycles most fold-6 walks — before their first
DP. Detection therefore reaches past length 2: each walk keeps the hashes
(low limb of x) of its last `RHO_CYCLE_DEPTH` = 5 points, and when the new
point matches one of them it has closed a cycle of length 2..6 whose
members it has all just seen. It escapes by *doubling the member with the
smallest hash*. That member is a function of the cycle alone, not of where
the walk entered it, so two walks trapped in the same cycle leave it
identically and their collision survives; reaching it costs at most
`length − 1` more steps around the cycle, which are neither counted nor
tested for distinguished points. A 32-bit hash coincidence with a
non-member simply lapses after `RHO_CYCLE_DEPTH + 1` checks. Longer cycles
are still left to `max_steps` — the bench defaults it to eight DP periods —
and `bench rho` reports the count of detected cycles by length and the
number of walks that hit `max_steps`, which is the number to watch: if it
is not near zero, raise `r_bits` or `RHO_CYCLE_DEPTH`.

Walks are identified by (index, restart counter) and their start point is
derived from that pair by a fixed PRNG, so a distinguished point is reported
as just `{x, walk, restart, steps}` and the host recovers the coefficients
by replaying only the two walks that actually collided.

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
document, and 792 against 1616 after. (Those figures predate the fold-6
multiplication and the five-word cycle ring in the walk state; re-measure
with `ptx_stats.sh`.)

## Curve support

`ecref.py params --curve <name>` emits the parameter header. `secp256k1`,
`toyNN` (a deterministically generated NN-bit curve with prime order and
a ≠ 0) and `toyNNj0` (an NN-bit prime-order curve with a = 0 and
p ≡ 1 mod 3, so that Aut(E) = Z/6 is available) are built in. Any other
short-Weierstrass curve over a prime below 2^256 works by adding its
parameters to `ecref.py`; only the `FP_FAST` path is secp256k1-specific,
and it falls back to Montgomery automatically. Fold 6 needs a j = 0 curve;
the generator derives β and λ itself whenever the curve has them.
