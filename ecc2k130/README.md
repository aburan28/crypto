# ECC2K-130 and ECC2K-95

A GPU-oriented client for the Certicom binary-curve challenges: Pollard rho with
a Frobenius-based iteration function, bitsliced, with the field arithmetic
emitted by a code generator.

Two targets are configured:

* **ECC2K-130**, still open, `GF(2^131)`, about `2^60.9` iterations. Uses the
  permuted type-II optimal normal basis of Bailey et al.
* **ECC2K-95**, solved by Harley's group in 1998, `GF(2^97)`, about `2^44`
  iterations. `2*97+1 = 195` is composite so that field has no type-II optimal
  normal basis; this one runs on a polynomial-basis backend with the
  orbit-invariant weight taken through a single linear map into a normal basis,
  which is the arrangement Harley's own client used. Its answer is public, so a
  recovered logarithm can be checked rather than merely believed.

The whole arithmetic layer is produced by a code generator that verifies every
routine it emits against an independent model of the field. The same source
compiles for CUDA (32-bit lanes) and for the CPU (64-bit lanes), so the exact
code that would run on a GPU is what the test suite exercises.

## Status

| | |
|---|---|
| Field arithmetic, iteration function, solver | implemented and tested |
| End-to-end discrete logarithms | recovered on `GF(2^23)` and `GF(2^41)` |
| CPU client | measured, 12.6 M iterations/s per core |
| CUDA client | host and device both compile; **never run on a GPU here** |
| Modal integration | validate, benchmark, autotune, search, fan out |
| ECC2K-95 instance | parameters recovered and independently verified |

No GPU was available while this was written, so the CUDA path is verified by
compiling it — host and device halves, plus `ptxas` register allocation — and
by running the identical arithmetic on the CPU. `modal_app.py` exists to close
that gap: it builds the client on a Modal GPU, runs the same validation there,
recovers discrete logarithms on the device, and autotunes the build knobs
against real hardware.

## The problems

Solve `Q = [k]P` on the Koblitz curve `y^2 + xy = x^3 + 1`.

**ECC2K-130**, over `F_{2^131}`, `|E| = 4l`:

```
l = 680564733841876926932320129493409985129        (129-bit prime)
```

Expected work is `2^60.9` iterations, producing about `2^35.6` distinguished
points. Open since 1997.

**ECC2K-95**, over `F_{2^97} = F_2[t]/(t^97 + t^6 + 1)`, `|E| = 4l`:

```
l = 39614081257132074233778707191                  (95-bit prime)
XP = 08A84FB02034F7771DC940097   YP = 1D2F10A471D48A720F18F6339
XQ = 0E0BC08AC5818F303E2B05E90   YQ = 134C028FC3393124D673E6F8E
```

Expected work is `2^44` iterations; Harley's group measured 2.16e13 in 1998,
against the `1.79e13` this repository predicts from `sqrt(pi*l/(4m))`. The
parameters were recovered from Harley's surviving client source and are checked
by the generator every time it runs: both points lie on the curve, both have
order `l`, and the published answer

```
k = 37837308472231540269443981458
```

satisfies `[k]P = Q`. That makes it a real end-to-end test of the whole pipeline
at a scale of a few GPU-hours rather than a few GPU-centuries.

## Design

### Word width

One word operation advances `W` walks at once. The device uses 32-bit words,
and the host picks the widest it has: AVX-512 gives 512 lanes and, through
`vpternlogd`, the same single-instruction three-input logic that `LOP3.LUT`
gives on the GPU. The arithmetic source is identical for all three widths;
only `bitslice.h` changes.

### Two field backends

`fieldbs.h` implements the permuted type-II optimal normal basis, where squaring
is an index permutation and the orbit weight is immediate. It needs `2m+1`
prime, which holds for `m = 131` but not for `m = 97` or `m = 109`.

`fieldpb.h` implements a polynomial basis `F_2[z]/(F)` for the fields that have
no such basis. A multiplication is the same Karatsuba tree followed by reduction
modulo a trinomial, and squaring costs 51 instructions rather than nothing. The
weight, which has to stay invariant under Frobenius for the walk to run on
orbits, comes from one generated linear map into a normal basis; the generator
picks the sparsest normal element it finds and checks that squaring really does
rotate the resulting coordinates.

Which backend a curve uses is a property of its config, so the walk, the kernel,
the solver and the reference are all written once:

| | multiply | square | weight |
|---|---|---|---|
| `GF(2^131)` normal basis | 8859 | free | 268 |
| `GF(2^97)` polynomial basis | 5627 | 51 | 1716 |

### Field representation

An element is `m` words in the permuted type-II optimal normal basis
`gamma_i = zeta^i + zeta^-i`, `i = 1..m`, where `zeta` is a primitive 263rd root
of unity. Two properties drive the whole design:

* squaring is the index permutation `i -> fold(2i)`, so it costs nothing when
  loops are unrolled, and
* the Hamming weight is invariant under squaring and under negation, so the
  iteration function is well defined on orbits of size `2m` without ever
  computing a canonical representative inside the loop.

The generator models this field as the symmetric vectors over `Z/263` modulo the
all-ones vector, with multiplication as cyclic convolution. That is a complete
and obviously correct model of `GF(2^131)`, and it is what every generated
routine is checked against.

Multiplication follows Bernstein-Lange: convert both operands to the optimal
polynomial basis `{1, c, c^2, ...}` with `c = gamma_1`, multiply as polynomials,
convert the 261-coefficient product back. Both conversions are generated from
the substitution `c = z + 1/z`, which splits by parity as
`H(c) = Heven(c^2) + c * Hodd(c^2)`; since squaring is an index permutation this
gives an `O(K log K)` conversion instead of the `O(K*m)` of a dense matrix.

### Iteration function

```
R_{i+1} = sigma^j(R_i) + R_i,    j = 3 + ((HW(x_{R_i}) / 2) mod 8)
```

with a point distinguished when `HW(x) <= 34`. The `sigma^j` selection is
branch-free: three conditional squarings, each one `LOP3` per word, and `j` is
never materialised.

### Batching and memory traffic

Each thread owns `ECC_BATCH` slots; a slot is 32 walks in the bit lanes of one
word. An iteration is two passes: the first computes weights, the
distinguished-point test and the running product of the addition denominators;
the second does one field inversion for the whole batch and then the affine
additions.

Pass 2 recomputes `sigma^j(x)` and `sigma^j(y)` from three stored weight bits
instead of reading back stored denominators. That trades 786 instructions for
1048 bytes of traffic per slot, which is worth it here — see the memory
analysis below.

### Distinguished points and restarts

A walk that reports is restarted in place from `R = Q + sum c_i sigma^i(P)`,
with `c` a 128-bit string from a PRF of the walk seed, computed with the same
bitsliced arithmetic and merged into the finished lanes under a mask. No linear
combination of `P` and `Q` is tracked in the loop; the server recomputes both
walks from their seeds when two of them collide, which is what keeps the inner
loop free of conditional counter updates.

Reports are `(seed, endpoint)`. Collision resolution recomputes each walk
counting how often each `sigma^j + 1` was applied, giving
`endpoint = [mu](alpha_0 P + Q)` with `mu = prod_j (1 + s^j)^{n_j}`, matches the
two endpoints up to Frobenius and negation, and solves for `k`.

## Results

### Validation

`make test` runs 52 checks across `GF(2^23)`, `GF(2^41)`, `GF(2^83)` and
`GF(2^131)`; all pass.

* every generated routine is checked against the independent field model inside
  the generator, on 64 random inputs at a time, before it is written out;
* the bitsliced multiply, square, inverse, `sigma^7` and Hamming weight are
  differential-tested against the reference implementation, which shares no code
  with them;
* the iteration function commutes with Frobenius and with negation, weight is
  constant on an orbit, and subgroup x-coordinates have even weight;
* the bitsliced start point matches the reference, and the tracked scalar
  reproduces it;
* the collision solver recovers a planted discrete logarithm;
* the distinguished-point comparison is checked exhaustively over all weights
  and cutoffs;
* the challenge point converts out of the normal basis back to Certicom's
  published polynomial-basis hex, `05 1C99BFA6 F18DE467 C80C23B9 8C7994AA`.

Longer runs check the reporting path rather than the algebra. On `GF(2^83)`,
12.6 M iterations produced 172 distinguished points at a rate within 30% of the
predicted one, and 60 of them were recomputed from their seeds by the reference
implementation with no mismatch, including walks that had been restarted. On
the challenge curve itself, run with a deliberately loose cutoff of 34 -> 50 so
that points actually appear, 3.1 M iterations produced 12368 reports, 25 of
which were recomputed and matched exactly.

### ECC2K-95

The generator refuses to emit the curve unless the recovered parameters check
out: both points on the curve, both of order `l`, and the published answer
satisfying `[k]P = Q`. The validation suite repeats the last of those against
the reference implementation.

Cost, measured per iteration rather than guessed:

| | GF(2^131), normal basis | GF(2^97), polynomial basis |
|---|---|---|
| multiply | 8859 | 5627 |
| square | free | 51 |
| weight | 268 | 1716 |
| whole iteration | about 60000 | about 32000 |

So ECC2K-95 runs about 1.9x faster per iteration than ECC2K-130 and needs
`2^44` of them instead of `2^60.9`, which is the difference between GPU-hours
and GPU-centuries.

One operational point matters more than the arithmetic. The number of reported
points comes out at roughly four times the number of parallel walks, because
each walk only gets `total/walks` steps and the cutoff has to be loose enough
that it reports within them. Filling a large GPU with 50M walks would therefore
demand 200M reports; the Modal search sizes the walk count first and picks the
weight cutoff from it, defaulting to 4M walks and about 16M reports.

### End-to-end

`make break-small` recovers planted discrete logarithms through the complete
pipeline — bitsliced walk, distinguished points, store, collision, recomputation
with exponent tracking, linear algebra mod `l`, and verification that
`[k]P = Q`. On `GF(2^41)` (`l` is 39 bits, about `2^16.6` expected iterations)
each instance takes well under a second on four cores.

### Throughput

Measured on one core of a 2.1 GHz Xeon, `GF(2^131)`, median of five runs:

| word | lanes | batch | iterations/s | cycles/iteration |
|---|---|---|---|---|
| `uint64_t` | 64 | 16 | 4.42 M | 475 |
| AVX2 | 256 | 32 | 10.7 M | 196 |
| AVX-512 | 512 | 32 | 12.6 M | 167 |

Four threads reach 40 M iterations/s. For comparison, the 2009 hand-written
qhasm implementation reached 533 cycles/iteration on a Core 2 with 128-bit
vectors.

On Apple Silicon (M4 Pro, 10P+4E cores, `GF(2^131)`, 64-bit lanes, sustained,
`--steps 8192 --launches 8`):

| threads | BATCH=32 | BATCH=64 | BATCH=96 |
|---|---|---|---|
| 4  |  ~10 M  |  ~10 M  |  ~10 M  |
| 8  |  ~19 M  |  ~19 M  |  ~19 M  |
| 14 |  ~33 M  |  ~34 M  |  ~34 M  |

The 64-bit lane path autovectorizes into NEON 64-bit pairs under `-O3`, which
is worth about 25% on this hardware. Throughput peaks near the DRAM ceiling
(roughly 270 GB/s on M4 Pro, 8 MB/iteration at 14 threads x BATCH=64). The
default `BATCH=32` is near-optimal; `OMP_PROC_BIND=close OMP_PLACES=cores`
recovers another ~5% on the host OpenMP runtime.

The batch size matters more than cache pressure would suggest: going from 4 to
32 walks per inversion is worth 60% because it drives the amortised inversion
from 25% of the multiplication budget down to 3%.

| Karatsuba leaf | 9 | 17 (default) | 33 |
|---|---|---|---|
| instructions per 131-bit multiply | 8859 | 8859 | 10701 |

### Device resources (offline, `ptxas` 12.9, no GPU)

For `sm_120` at `ECC_BATCH=16`, leaf 17:

| | walk kernel | reseed kernel |
|---|---|---|
| registers | 255 | 255 |
| stack frame | 16440 B | 15400 B |
| spill stores / loads | 5160 B / 7024 B | 1556 B / 1564 B |

`mulLeaf`, `multPrep`, `toOnb`, `hamming` and `sigmaJ` each compile with zero
spills. The remaining spills are in the Karatsuba glue and the kernel body,
which hold several 131-word intermediates at once. The per-thread stack frame
is the number to watch on real hardware: it is local memory, so occupancy
depends on how much of it stays in L1.

### Occupancy versus spilling

A block-size bound alone lets ptxas give every thread all 255 registers, which
caps an SM at about 256 threads. Asking for more resident blocks
(`MINBLOCKS`) buys warps by taking registers away:

| blocks/SM asked | registers | spill stores | spill loads | warps/SM |
|---|---|---|---|---|
| 1 | 255 | 5436 | 7256 | 4 |
| 2 (default) | 255 | 5436 | 7256 | 8 |
| 4 | 128 | 7520 | 9648 | 16 |
| 6 | 80 | 9540 | 12096 | 24 |
| 8 | 64 | 10356 | 12960 | 32 |

Two is free: 128 threads at 255 registers each is 65280 of the 65536 an SM has,
so occupancy doubles without a single extra spill, and that is now the default.
Past that the trade is real and only a device can settle it, which is what
`autotune` sweeps.

## Findings that update the implementation guide

**LOP3 fusion is worth 1.65x, and it changes the Karatsuba tradeoff.** A
131-bit multiplication costs 14554 two-input bit operations but only 8859
instructions once `xor(xor(a,b),c)` and `xor(a, and(b,c))` are fused into single
three-input operations. Because the fused form charges nothing for the XOR that
accompanies an AND, *schoolbook leaves beat Karatsuba leaves*: recursing all the
way down costs 18819 bit operations and 14770 instructions, while stopping at a
schoolbook leaf costs more bit operations and fewer instructions. The guide's
§5.1 question — whether to re-evaluate Karatsuba depth — resolves toward
shallower recursion than the bit-operation count suggests.

**The conversions are not a problem.** Generated from the parity recursion and
fused, `multPrep` is 331 instructions and the double-size inverse conversion is
903, against 325 and 909 for the hand-designed chains in the papers. The
Hamming weight tree is 268 instructions against 654 bit operations, because a
full adder is two `LOP3`s.

**Code size, not arithmetic, is the wall — and it is avoidable.** Compiling the
batch with everything inlined produced 285,096 PTX instructions for one kernel;
`ptxas` had not finished register allocation after five minutes and 5 GB of
resident memory. Marking the multiplication, inversion, conversions, leaf,
weight tree and `sigma^j` as real subroutines and refusing to unroll the batch
loop brings this to 36,553 PTX lines, a 3.5 s front-end compile and a 2.4 s
`ptxas` run. This is the §9 toolchain risk, and it is entirely a consequence of
inlining policy: no reverse-engineered assembler is needed.

**Keep walk restarts out of the hot kernel.** Building a fresh start point
costs 128 point additions, and inlining that call chain into the walk kernel
cost 5.8 KB of per-thread stack frame (22264 down to 16440 bytes) for a path
taken roughly once in 2^25 iterations per lane. Reporting a point now only
marks its lane, and a separate kernel — launched between walk kernels, and
skipped entirely when nothing reported — revives marked lanes. On the CPU the
same change was worth about 15% of throughput.

**Every generated routine register-allocates cleanly.** For `sm_120`, `mulLeaf`,
`multPrep`, `toOnb`, `hamming` and `sigmaJ` all compile to zero spill stores and
zero spill loads. What spills is the Karatsuba glue inside `mul` and the kernel
body, which hold several 131-word intermediates at once.

**On current hardware this loop is close to memory bound, which it was not in
2009.** With denominators recomputed rather than stored, one walk-step moves
about 98 bytes (state in and out, plus the Montgomery product chain) and costs
about 1700 word instructions. A 5090-class part is balanced at roughly 29
instructions per byte and this design sits at about 18, so memory is the tighter
of the two constraints by about 1.6x — where the 2009 GTX 295 implementation
spent only 11.9% of its cycles on DRAM. Compute grew about 100x since then and
bandwidth about 8x. Anything that reduces bytes per walk-step is now worth more
than anything that reduces bit operations, which is the opposite of the tradeoff
the original design faced.

## Running on Modal

`modal_app.py` builds the client into a CUDA image (a fat binary covering
Ampere through Blackwell) and exposes five entry points:

```
export ECC_GPU=RTX-PRO-6000
modal run modal_app.py::validate      # correctness, on the device
modal run modal_app.py::bench         # throughput
modal run modal_app.py::autotune      # sweep the build knobs on real hardware
modal run modal_app.py::search --hours 4
modal run modal_app.py::fanout --count 8 --hours 4
modal run modal_app.py::merge         # collisions across every run
```

The GPU comes from `ECC_GPU`, read when the file is imported and baked into the
function definitions, which works on every Modal version. Clients from 0.72
also accept `--gpu` on the entry points to override it per call.

`validate` runs the whole test suite on the GPU machine and then recovers
planted discrete logarithms with the CUDA engine itself, so a GPU run proves the
same thing the CPU run does. `autotune` rebuilds for the local compute
capability only and sweeps batch size, block size and Karatsuba leaf, writing
the ranking to a Modal Volume. `search` collects distinguished points into that
same volume, checkpointing its live walks alongside them, so a stopped run
resumes where it left off and several containers contribute to one corpus; each
gets its own run id, which keeps their seed spaces disjoint, and each loads its
siblings' corpora so a cross-container collision is caught as it happens.
`merge` scans the corpus for repeated hashes, which are the candidate
collisions.

`RTX-PRO-6000` is the interesting target: it is the GB202 part the guide
identifies as the best value for this workload, since none of a datacenter
GPU's tensor silicon is reachable from binary-field arithmetic.

The image builds a fat binary for sm_80 through sm_120 with nvcc, which takes
about three and a half minutes and is cached thereafter. `autotune` rebuilds
for the local compute capability alone, which takes seconds.

## Running it

`run.sh` wraps the Modal entry points. A container has a finite life, so a real
search is a loop: each pass runs for `HOURS`, is stopped with SIGTERM so the
client checkpoints, and the next pass resumes from that checkpoint.

```
./run.sh validate                      # prove the GPU engine before spending anything
./run.sh bench                         # throughput on the challenge curve
CURVE=97 HOURS=4 PASSES=6 ./run.sh search
COUNT=8 ./run.sh fanout                # the same across eight GPUs
./run.sh merge                         # scan the corpus for collisions
```

`CURVE=131` collects on ECC2K-130 itself. That run does not finish: 2^60.9
iterations is decades of GPU time, so treat it as collection, not as a solve.
Two things follow from a run that never ends, and both are handled rather than
left to bite:

* The distinguished-point cutoff cannot be sized against the whole expected run,
  or no walk ever reaches it. At the full-run choice for m=131 a four-hour pass
  on a fast GPU reports about **two hundred points in total**. The cutoff is
  therefore sized from the iterations the pass will actually do, measured by the
  bench it already runs, which puts it back at a few reports per walk -- roughly
  30M points and a gigabyte of corpus per pass. `DPW` overrides it.
* The corpus outgrows memory. Every pass reloads it to catch collisions in
  process, and a store entry costs far more than the 32 bytes it occupies on
  disk, so an uncapped reload is what eventually ends a long collection run.
  `LOADMAX` (default 50M points) caps it, newest file first. What that gives up
  is finding a collision in process; the corpus is still complete on disk and
  `./run.sh merge` still finds it there.

Every pass is launched from the same variables, because resume only works when
the shape matches: the checkpoint header records curve, run id, worker count and
batch, and changing one between passes makes the client refuse the checkpoint
and start over. Run `validate` first -- it recovers planted discrete logarithms
on the GPU itself, so a broken kernel fails in a minute rather than quietly
burning a day of credits.

The searcher reports every 60 seconds:

```
[1h04m] 842.1 M it/s  3.24T iters (18.412% of 2^44.0)  14.80M dp  14.79M distinct  corpus 473.6 MB  2h56m left
```

Rate, iterations against the expected total for the curve, points reported
against distinct orbits stored (the gap is walks that re-reported), corpus size
on disk and time left in the pass. Anything that is not a progress line -- a
collision, a verification failure, a checkpoint that could not be written -- is
printed as it happens.

## Persistence

A rho search is long enough that a run will be interrupted -- a spot instance
goes away, a time budget expires, a container is recycled -- so nothing the
client computes is allowed to exist only in memory.

At any instant the walks in flight hold about a quarter of everything the run
has computed. A walk reports after roughly 1/theta steps and the cutoff is set
so a walk reports a few times within the steps it gets, which puts the mean
unreported prefix at a quarter of the run. Killing a client without warning
therefore throws away 25% of the work, and that fraction does not shrink as the
run gets longer.

Three things are saved.

**The corpus.** `--dp-file F` appends 32-byte records of (seed, canonical orbit
hash): fixed width, so a file can be counted with a stat, appended to by several
writers and truncated by a dying container without becoming unparseable. A
short trailing record is ignored rather than misread.

**Old points as live state.** `--load F` reads a corpus back into the store at
startup, so a collision between today's walk and one from last week is found the
moment the new point arrives, not by an offline merge. Reload collisions are
solved immediately. Paths are canonicalised and deduplicated, so naming the same
corpus twice costs nothing.

**The walks themselves.** `--checkpoint F` writes every live walk -- point,
seed and step counter -- every `--checkpoint-every` seconds (default 300) and
again on exit, to a temporary file that is then renamed, so an interrupted write
leaves the previous checkpoint intact rather than a half-written one. The header
records m, thread count, batch size, lane width and run id; a client that does
not match refuses the file and starts fresh instead of misreading it. Cost is
about 49 bytes per walk, so 195 MB for the four million walks the Modal search
uses by default.

`SIGINT` and `SIGTERM` stop the client between launches rather than killing it:
it finishes the launch in flight, drains the reports, checkpoints and exits 0.
The Modal driver signals rather than terminates for exactly this reason, waits
for the client to finish writing, and only then commits the volume, so the
snapshot contains the checkpoint just written rather than the one before it.

Resuming needs the same shape it saved, so pass the same curve, run id, worker
count and build-time batch size; the header check turns a mismatch into a fresh
start rather than corruption.

## Build

```
make generate     # regenerate the arithmetic headers (needs python3 only)
make cpu          # CPU client
make test         # validation suite
make break-small  # end-to-end discrete logarithms on the small curves
make bench        # throughput on the challenge curve
make gpu          # CUDA client (needs nvcc)
make check-cuda   # type-check host and device with clang, no GPU or nvcc needed
make ptx          # device compile + ptxas report, no GPU needed
```

Build-time knobs: `BATCH` (walks batched per inversion, default 32), `THREADS`
(CUDA block size), `MARCH` (host architecture; `native` picks the widest word),
and `--leaf` to the generator (Karatsuba leaf size).

```
./ecc2k130-cpu --curve 131 --bench
./ecc2k130-cpu --curve 41 --instance 3      # recover a planted discrete log
./ecc2k130 --curve 131 --dp-file dps.bin --checkpoint state.ck   # collect
```

### Building on macOS

`make cpu` defaults to `CXX=g++`, but on macOS that resolves to Apple's
clang, which has no `-fopenmp`. Use the Homebrew toolchain:

```
make cpu CXX="/opt/homebrew/bin/g++-16 -isysroot $(xcrun --show-sdk-path)"
```

The CUDA targets (`make gpu`, `make ptx`, `make check-cuda`) need a CUDA
toolchain, which is no longer published for macOS — see the Modal section
above for the supported path.

## Index calculus: SAT point decomposition

An experiment, not a second attack. Index calculus does not beat rho on a
curve anyone cares about, and nothing here changes that -- the growth measured
below is the reason. What was missing from the public record is a measured
curve for a Frobenius-reduced factor base on a Koblitz curve, and this produces
one. Needs `pip install python-sat pycryptosat`; `make test-decomp` checks it.

```
cd codegen && python3 indexcalc.py --m 11 --points 3 --weight 3 --trials 5
```

**Factor base.** The Koblitz saving is that factor base columns are Frobenius
orbits, so a run needs |F|/m relations rather than |F|. The usual construction
takes the points whose x lies in a Frobenius-stable subspace -- but in a normal
basis the Frobenius is a coordinate permutation, so a stable subspace is an
F_2[t]/(t^m-1) submodule, which is to say a binary cyclic code of length m. Its
available dimensions are the degrees of the divisors of t^m-1, and at m=131 and
m=163 two is primitive mod m, t^m-1 = (t-1)(irreducible), and **the only stable
subspaces are the trivial ones**. The construction is unavailable at exactly
the interesting sizes; m=127 (ord 7) and m=233 (ord 29) are the lucky cases.

So use `{P : HW(x(P)) <= w}` in the normal basis instead. Weight is invariant
under a coordinate permutation, so this is Frobenius stable for every m. It is
not a low-degree algebraic condition, so a Groebner descent cannot use it -- but
a solver that does not need linearity can, and a cardinality constraint costs
O(mw) clauses. It is also the predicate the walk kernel already computes.

**Decomposition.** Chain Semaev's third summation polynomial rather than
resolving the (k+1)-th, whose degree is 2^(k-1) per variable: introduce the
intermediate sums and constrain each link with `S_3(u,v,w) = (uv+uw+vw)^2 +
uvw + b`. Squaring is sigma^1, a permutation of coordinates, so the square
costs nothing -- in a polynomial basis it would be m^2 XORs per link.

**Encoding.** `Prog.emitCnf` is `Prog.emit` with clauses instead of C, so a
multiplication costs the generated Karatsuba circuit's bit operations rather
than the m^3 AND terms its definition would need: 15588 gates at m=131 against
2.2 million. The circuit is XOR-dominated, so XOR clauses go to the solver
natively. Applying sigma to a decomposition of R gives a decomposition of
sigma(R), and the relation it yields is lambda times the original, so there is
no Frobenius symmetry to break inside one instance -- what is breakable is the
k! orderings of the same multiset of points.

**Measured.** Median over solved instances, k=3 points, w=3, CryptoMiniSat:

| m | gates | orbits | median solve |
|---|---|---|---|
| 9  | 1587 | 7  | 0.11 s |
| 11 | 2311 | 9  | 2.64 s |
| 23 | 8029 | 52 | none finished; one instance ran 25 min without returning |

Roughly 2^2.3 per bit of field over that range. Three points is not a fit, and
the encoding has not been tuned, but the shape is already clear enough to say
what the instrument is for: finding where the curve bends, if it does.

## What is not done

* Throughput on real hardware is unmeasured. The client has now run on an
  RTX PRO 6000 Blackwell (sm_120): the whole validation suite passes there and
  every planted discrete logarithm is recovered on the device itself, through
  both backends. What has not been measured is how fast it walks. `bench` is
  the number that settles that, and `autotune` decides the build knobs —
  including whether the register-budget leaf, which is chosen from static ptxas
  analysis and measures *worse* on the host, is right on a GPU. The §6 layout
  question — one thread per bitsliced multiply versus 32 threads cooperating —
  is likewise only partly answered: the register data says the leaf fits, but
  the 16 KB per-thread stack frame means occupancy needs measurement.
* The multiplier is optimal only within the Karatsuba family. Toom-3 over
  GF(2), which is where Bernstein's 11961-bit-operation chain comes from, is not
  implemented; a search over balanced and unbalanced Karatsuba splits and the
  guide's 128+3 decomposition found nothing better than 8859 instructions.
* No server. Corpora are files that can be reloaded and merged, but there is no
  UDP protocol, no hash-routed sharding, and no live multi-machine merging.
* The start-point PRF is a mixing function rather than AES. Any client that
  wants to interoperate with a different implementation must agree on it.
* Multi-GPU is not implemented: one process drives one device.
