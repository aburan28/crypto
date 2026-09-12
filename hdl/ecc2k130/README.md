# hdl/ecc2k130 — FPGA engine for the ECC2K-130 walk

The ECC2K-130 iteration function in portable VHDL-2008: a pipelined
Karatsuba GF(2^131) multiplier in the permuted type-II optimal normal
basis, a step unit that computes `R + sigma^j(R)` for a batch of walks
around it with one shared inversion, and a sequencer that runs many walks
through the step unit and reports distinguished points. It is the hardware
counterpart of the GPU client in [`ecc2k130/`](../../ecc2k130/): same
field, same basis, same iteration function, same distinguished-point rule,
and the same oracle — every vector the testbenches check against comes from
the client's own field model in `ecc2k130/codegen/`.

Reference part is the AMD Virtex UltraScale+ VU47P (AWS `f2.6xlarge`), the
same one `hdl/ecc/` and `hdl/sha1/` target, but there are no vendor
primitives: plain `ieee.std_logic_1164` and `numeric_std`, simulates under
GHDL, reads into any synthesis flow. Nothing here has been synthesised or
run on an FPGA; every number below is labelled measured (in simulation),
derived, or estimated.

## Files

| File | What it is |
|---|---|
| `gf131_pkg.vhd` | The field: coordinates, `sigma^k`, `sigma^j` selection, two-clock weight, the two constant basis-change maps built at elaboration, and an independent direct-product oracle |
| `gf2_kmul.vhd` | Recursive Karatsuba polynomial multiplier over GF(2), `LEVELS` deep, schoolbook leaf, latency `2 LEVELS + 1` |
| `gf131_mul.vhd` | GF(2^131) multiplier around it, II = 1, latency 10, no DSPs |
| `ec2k_batch_pipe.vhd` | **The step unit.** W walks per batch, their inversions shared through a product tree: `5 + 5/W` multiplies per step, bursts of independent multiplies streamed from a ready queue |
| `ec2k_step_pipe.vhd` | The simple step unit, one inversion per walk, 10 multiplies per step; kept as the readable reference |
| `ec2k_walker.vhd` | The rho sequencer: N walks, host load port, distinguished-point output, around one `ec2k_batch_pipe` |
| `ec2k_axil.vhd` | **The host interface.** AXI4-Lite register block around `NENG` walkers: load registers, a distinguished-point queue, step/point/drop counters, geometry readback |
| `gf131_mul_tb.vhd`, `ec2k_batch_tb.vhd`, `ec2k_step_tb.vhd`, `ec2k_walker_tb.vhd`, `ec2k_axil_tb.vhd` | Self-checking testbenches |
| `gf131_tb_pkg.vhd` | Hex helpers shared by the testbenches |
| `ecc2k_ref.py` | Reference model of every hardware algorithm, proved against `ecc2k130/codegen/field.py`; generates the vectors; gate counts |
| `vectors_ecc2k130.txt` | Generated: 256 products, 200 steps, 64 walks to a distinguished point |
| `host/` | `ecc2k130-fpga`, the host program: drives the register block over PCIe (or a behavioural model with `--sim`), speaks `worker.py`'s client contract, writes the campaign's dp records and checkpoints |
| `aws/` | Out-of-context synthesis script; the F2 custom logic `cl_ecc2k130` around `ec2k_axil`; the AWS CLI that builds the image and runs an F2 fleet — see [`aws/README.md`](aws/README.md) |

## Run

```bash
make            # analyse, elaborate and run the five testbenches
make rate       # steady-state clocks per step of the batched unit, 6144 steps
make check      # prove the hardware algorithms against the client's field model
make vectors    # regenerate the vector file
make cost       # gate counts of the multiplier from its constant matrices
make -C host test   # build the host program and run it against its own model
```

Needs GHDL and Python 3 (the generator imports `ecc2k130/codegen`). Current
output:

```
=== gf131_mul_tb ===
gf131_mul_tb: PASS -- 256 multiplications, II=1, latency 10 clk
=== ec2k_batch_tb ===
ec2k_batch_tb: 200 steps in 1292 clk = 646/100 clk per step, W = 16, 8 batches in flight (85 multiplies per batch)
ec2k_batch_tb: PASS
=== ec2k_step_tb ===
ec2k_step_tb: 200 steps in 2051 clk with 16 slots (10 multiplies per step)
ec2k_step_tb: PASS
=== ec2k_walker_tb ===
ec2k_walker_tb: 64 distinguished points from 797 steps in 17113 clk (32 walks, batches of 8, 4 in flight)
ec2k_walker_tb: PASS
=== ec2k_axil_tb ===
ec2k_axil_tb: 64 distinguished points, 797 steps through 2 engine(s) of 16 walks, 7011 status polls
ec2k_axil_tb: PASS

$ make rate
ec2k_batch_tb: 6144 steps in 33511 clk = 545/100 clk per step, W = 16, 8 batches in flight (85 multiplies per batch)
```

The multiplier vectors include zero, one (which is the all-ones vector in
this basis), single basis elements and their squares, and the top and
bottom coordinates alone; each product is checked against the file and
against `gf_mul_ref`, the direct normal-basis convolution computed in
simulation, which shares nothing with the multiplier's route. The step and
walk vectors are random points of order `l` on the challenge curve itself.
The 200 step vectors are twelve full batches of sixteen and a partial one,
so the flush path that pads a short batch is exercised on every run, and
`ec2k_batch_tb` passes for every `W` from 2 to 32 with 2 to 16 batches in
flight.

## The field

An element is 131 bits, bit `i-1` the coefficient of
`gamma_i = zeta^i + zeta^-i`, `zeta` a primitive 263rd root of unity. This
is the basis the client uses (`ecc2k130/README.md`, "Field representation"),
and it has the two properties the whole design is built on:

- **Squaring is wiring.** `a -> a^2` sends coordinate `i` to `fold(2i)`,
  `fold(k) = min(k mod 263, 263 - k mod 263)`. `gf_frob(a, k)` is a
  permutation for every `k`, so `sigma^k` costs no logic and no clock.
- **The weight is squaring-invariant**, so the iteration function is
  well defined on Frobenius orbits and `j` can be read off the weight of
  `x` without canonicalising anything.

Multiplication does not happen in the normal basis. Following
Bernstein–Lange, both operands are converted to the optimal polynomial
basis `{c, c^2, ..., c^131}`, `c = gamma_1`, multiplied as ordinary GF(2)
polynomials, and the 261-coefficient product is converted back. Both
conversions are fixed GF(2)-linear maps and the package computes them at
elaboration:

```
gamma_i = T_i(c)        T_0 = 0, T_1 = c, T_i = c T_{i-1} + T_{i-2}     (Dickson)
c^k     = sum_j C(k,j) zeta^(k-2j)      C(k,j) odd  iff  j & (k-j) = 0    (Lucas)
zeta^e + zeta^-e = gamma_fold(e)
```

There is no reduction step: the back-conversion takes every `c^k`,
`k = 2..262`, straight to coordinates. `ecc2k_ref.py check` proves both
maps, the multiplier route, `sigma^j`, the inversion chain and a full step
against the client's field model before any vector is written.

## The multiplier

```
stage 0        latch a, b
stage 1        a' = prep(a), b' = prep(b)          gamma -> c-powers
stage 2        Karatsuba level 1 pre-add           a0, a1, a0+a1 (66 bits)
stage 3        Karatsuba level 2 pre-add           (33 bits)
stage 4        Karatsuba level 3 pre-add           (17 bits)
stage 5        27 products of 17 x 17, schoolbook
stage 6        level 3 post-combine
stage 7        level 2 post-combine
stage 8        level 1 post-combine                261-bit product
stage 9        r = to_onb(h)                       c-powers -> gamma
```

`gf2_kmul` is a recursive entity: each level is a registered pre-add, three
instances of itself on the halves, and a registered post-combine. Over
GF(2) Karatsuba is exact with no subtractions, so the unequal halves
(66 + 65) cost nothing — the high half is zero-extended and the top product
bits are provably zero and dropped. Every stage is one or a few LUT levels:
pre-add and post-combine are one XOR each per bit, the leaf is 17 AND terms
per output bit (three per LUT6 with their XORs, two to three levels), and
the widest thing anywhere is the 66-input column of `prep`.

Gate count per multiplier, read off the constant matrices and the recurrence
(`make cost`), **derived**:

| Levels | Leaf | AND | XOR2 (product) | XOR2 (per multiply, with prep and to_onb) | Latency |
|---|---|---|---|---|---|
| 0 (schoolbook) | 131 | 17161 | 16900 | 22049 | 4 |
| 1 | 66 | 13068 | 13199 | 18348 | 6 |
| 2 | 33 | 9801 | 10520 | 15669 | 8 |
| **3 (default)** | **17** | **7803** | **9404** | **14553** | **10** |
| 4 | 9 | 6561 | 9512 | 14661 | 12 |

`prep` is 991 XOR2 per operand, `to_onb` 3167. No DSP48 is involved
anywhere — binary-field arithmetic has no carries, so the DSP's adder is
useless to it and the multiplier is pure LUT fabric.

Where the LUT optimum sits was a synthesis question; Vivado 2025.2 on the
VU47P answers it (`gf131_mul` out of context, 4.0 ns clock), **measured**:

| Levels | LUTs | FFs | worst slack |
|---|---|---|---|
| 2 | 5 547 | 2 892 | +1.79 ns |
| **3** | **4 855** | 4 647 | **+2.68 ns** |
| 4 | 5 019 | 7 262 | +2.79 ns |

Three levels are 12% fewer LUTs than two and the leaf stage is a LUT level
shallower; four adds back more combine XORs than it removes ANDs. The two
extra clocks of latency cost 1% of the step unit's throughput (5.45
against 5.39 clocks per step, `make rate`), and flip-flops are not what
this device runs out of. Compare `hdl/ecc/`: the secp256k1 multiplier is
256 DSPs plus ~10k LUTs for the same II = 1.

The `tag` rides alongside and the datapath never looks at it, so one
multiplier serves any number of independent contexts back-to-back.

## The step

```
j   = 3 + ((HW(x) >> 1) & 7)
x2  = sigma^j(x),  y2 = sigma^j(y)          permutations, free
d   = x + x2
lam = (y + y2) / d
x3  = lam^2 + lam + d                       lam^2 = sigma(lam), free
y3  = lam (x + x3) + x3 + y
```

The division is Itoh–Tsujii, and here the free squaring pays for the second
time. `1/d = (d^(2^130 - 1))^2`, and the addition chain
`1, 2, 4, 8, 16, 32, 64, 128, 130` reaches `2^130 - 1` in **eight
multiplies**, each `beta_{2k} = beta_k * sigma^k(beta_k)`. On a prime field
the same inversion is about 270 multiplies (`docs/ecc_fpga_cost_model.md`),
which is why `hdl/ecc/` has to batch inversions across walks and leaves the
inverter outside. Here a whole step — inversion included — is **ten
multiplies and nothing else**, and that is what `ec2k_step_pipe` does: up to
`2^SLOT_W` independent steps in flight around one multiplier, dataflow
scheduled through a ready queue, 10.4 clocks per step with 16 slots
(measured, saturated).

## The batched step unit

Eight of those ten multiplies are the inversion, and inversions batch:
with `W` walks, `1/d_1 .. 1/d_W` cost one inversion plus `3(W-1)`
multiplies (Montgomery's trick). `ec2k_batch_pipe` takes `W` walks as the
leaves of a binary product tree:

```
forward    t(n) = t(2n) t(2n+1)                      W-1 multiplies, log W levels
invert     t(1) = 1/t(1)                              8 multiplies
backward   t(2n) = t(n) t(2n+1),  t(2n+1) = t(n) t(2n)  2(W-1) multiplies, in place
lam        lam_i = (y_i + y2_i) t(W+i)                 W
y3         y3_i  = lam_i (x_i + x3_i) + x3_i + y_i     W
```

**`5W + 5` multiplies per batch, `5 + 5/W` per step: 5.31 at `W = 16`**
against 10 — the same trick the GPU client plays with 32 walks per word,
and the single largest lever on throughput.

Every tree level, every inversion step, and each of the `lam` and `y3`
passes is a *burst* of independent multiplies that issue on consecutive
clocks. A batch is a chain of `2 log W + 10` bursts. When the last product
of a burst retires — products retire in issue order, the multiplier being a
fixed-latency pipe — the batch is pushed on a ready queue with its next
burst; the issue engine streams bursts from the head of the queue with the
following batch prefetched, so no clock is lost between bursts. Bursts from
different batches interleave freely and a handful of batches in flight
(`2^LOG_NB`) keeps the multiplier saturated. The tree is updated in place:
every read of a value a burst overwrites happens at issue, the overwrite at
retire eight clocks later.

The unit is built to make the clock the multiplier's, not the
scheduler's:

- every memory has exactly one writer, so it maps to distributed RAM, not
  flip-flops: `x, y, j, tag, valid` are written only by the fill, the tree
  only by retire, `z` (holding `d`, then `x3`) only by the issue side;
- operands are read one clock before they are formed, so no path runs
  RAM → `sigma^j` → XOR in one clock, and the retire path is a memory
  write plus one XOR — there is no `sigma` and no popcount in it;
- the Hamming weight, 131 bits wide, is taken over two clocks: 22 groups
  of six bits (one LUT6 per output bit) and then the sum, both on input
  and on output;
- batch and burst state ride in the multiplier tag, so results route
  back with no matching logic.

A batch shorter than `W` — the tail of a run, a slow host, a testbench —
would wait forever for leaves that never come, so a partly filled batch
with nothing arriving for `FLUSH_CLK` clocks is padded with dummy leaves
(`x = gamma_1`, whose `d = gamma_1 + gamma_8` is nonzero and keeps the
product tree invertible) that produce no output.

Measured in simulation (`ec2k_batch_tb`, 200 steps including a partial
batch and its flush wait, then `make rate` with 6144 steps in full batches):

| W | Batches in flight | 200 steps | 6144 steps | Bound `5 + 5/W` |
|---|---|---|---|---|
| 8 | 8 | 7.44 | 5.95 | 5.63 |
| 16 | 4 | 8.61 | 7.16 | 5.31 |
| 16 | 8 | 6.46 | **5.45** | 5.31 |
| 16 | 16 | 5.91 | 5.32 | 5.31 |
| 32 | 8 | 6.51 | 5.21 | 5.16 |

Clocks per step, multiplier latency 10. At the default `W = 16` with 8
batches in flight the long run sits 0.14 above the bound; about 500 clocks
of the 33511 are the fill of the first batch and the drain of the last,
the rest is the inversion's eight dependent single-multiply bursts, whose
latency eight batches cover almost but not quite (sixteen do). Four
batches are not enough. Memory per batch of `W` is `6W` field elements
(four leaf arrays and a tree of `2W`), 768 × 131 bits at the default.

Degenerate inputs (`d = 0`, i.e. `sigma^j(x) = x`) are not special-cased,
matching the client: the chain returns `1/0 = 0`, the product tree zeroes
and every step in the batch produces garbage, as the client's batched
inversion does. On the challenge curve the case has probability about
`2^-131` per step.

## The walker

`ec2k_walker` is the sequencer `hdl/ecc/` leaves as "not here". Walks live
in a ready FIFO of `(id, x, y)`: a walk is either inside the step unit or
waiting in the FIFO, and a completed step goes straight back into the FIFO
unless its `x` has weight at most `DP_WEIGHT`, in which case it is reported
on the `dp_*` port with its id and step count and the walk goes idle until
the host loads a new start point into that id. That is the same division
of labour as the GPU client, whose kernel reports `(seed, endpoint)` and
whose host owns restarts, the corpus and collision resolution
(`ecc2k130/README.md`, "Distinguished points and restarts"). `NWALK` should
comfortably exceed the `W · 2^LOG_NB` walks the step unit holds, so a full
batch is always forming; when it does not (start-up, a host slow to reload)
the flush keeps things moving at reduced efficiency, never a stall.

`ec2k_walker_tb` plays host: it loads 64 start points 32 at a time, checks
every report's id, step count and point against the trajectory the client's
reference model walked from the same start (under a loose cutoff of 56 so
walks report within a few dozen steps), and reloads. Walks interleave
arbitrarily inside the engine; the per-id bookkeeping makes the check
independent of that. The clocks-per-step figure the walker testbench prints
is not a throughput measurement: the walk population dwindles as records
run out and the tail is spent in partial batches waiting for the flush;
`make rate` is the one that measures the saturated rate.

## The host interface

`ec2k_axil` is what a host sees: an AXI4-Lite slave in front of `NENG`
walkers, each with `2^ID_W` walks, addressed as one flat space of
`gid = engine · 2^ID_W + id`. It does three things and nothing else.

- **Load.** The host writes `LD_ID`, the five 32-bit words of `x` and `y`
  (word `k` is bits `32k+31..32k`, word 4 the top three bits — the same
  limb order as the client's `Elem`), then `LD_GO`; `STATUS.LD_BUSY` clears
  when the engine has taken it. Loading an id that is still walking
  duplicates the walk, since the walker has no kill; the host only reloads
  an id after that id has reported.
- **Report.** A distinguished step goes back into the walker's ready FIFO
  flagged `dp`; when it reaches the head it is moved into the engine's
  report register instead of the step unit and held there, valid, until
  the register block acknowledges it (a one-clock `dp_ack`; both ends are
  registers, so engine-to-block wires have a full clock). The block takes
  each engine's report into a holding register and a rotating pointer
  drains those into a `2^DP_FIFO_W`-deep queue the host reads through
  `DP_ID`, `DP_STEPS`, `DP_X*`, `DP_Y*` and pops with `DP_POP`. Nothing in
  this chain can lose a report — a full queue or a slow host only delays
  them, with the walker's FIFO as the backlog — so `DROPPED` and
  `STATUS.DP_OVERFLOW` read 0 and exist for register-map compatibility.
  (The first version had a one-deep holding register that a second report
  on the next clock overwrote, which two distinguished leaves in one batch
  do; a lost report idles that walk until the host restarts.)
- **Count.** `STEPS` (64-bit, high word latched on the low read) is the
  total steps of every engine; `DPS` the reports queued. `GEOM` reads back
  `ID_W`, `LOG_W`, `LOG_NB`, `DP_WEIGHT` and `NENG`, so the host can refuse
  an image built for a different campaign; `MAGIC` is `0x2C130001`.

`CTRL.RUN` holds the engines in reset while clear, `CTRL.CLEAR` zeroes the
counters and the queue. The register map with byte offsets is at the top of
`ec2k_axil.vhd`. `ec2k_axil_tb` replays the walker testbench's 64 walks
through the registers alone, spread over two engines, and checks the same
per-id bookkeeping.

The step tests weight *after* each step, never the point the host loaded,
so a distinguished start point would walk on forever; the host tests start
points itself and reports them with zero iterations, exactly as the GPU
client does.

`host/ecc2k130_fpga.cpp` is that host. It is a drop-in for the CUDA client
under `ecc2k130/aws/worker.py`: same flags, same progress line, same
32-byte `(seed, canonical x)` dp records, a checkpoint the supervisor
resumes from the same way, verification of every report by rewalking it
from its seed with the packed CPU arithmetic, and collision resolution on
`--load` through the same `Solver`. Its seeds are `(runId, epoch, gid,
counter)`, so a resumed worker never repeats a seed and every seed
rewalks to its point without the FPGA. `--sim` swaps the PCIe bus for a
behavioural model of the register block, which is what `host/test.sh`
runs the whole `worker.py` contract against, down to feeding the resulting
`dp.bin` to `merge.py`.

## Capacity — synthesised

Vivado 2025.2, `xcvu47p-fsvh2892-2-e`, out of context, one `ec2k_walker`
(256 walks, W = 16, 8 batches in flight):

| | LUTs | of which LUTRAM | FFs |
|---|---|---|---|
| `ec2k_walker` (whole engine) | **13 146** | 4 812 | 7 388 |
| ├ walker body (FIFO with held reports, counters) | 2 311 | 1 736 | 384 |
| └ `ec2k_batch_pipe` | 10 835 | 3 076 | 7 004 |
| &nbsp;&nbsp; ├ step unit body (memories, scheduler, operand stage) | 6 012 | 3 076 | 2 361 |
| &nbsp;&nbsp; └ `gf131_mul` (three-level Karatsuba) | 4 823 | 0 | 4 643 |

The register block around them (`ec2k_axil`, 4 engines) adds 1 449 LUTs
and 1 773 FFs in total, not per engine.

Worst slack at a 4.0 ns clock is +1.98 ns (engine alone and with the
register block around four of them), so the fabric is not what limits the
clock; the shell's `clk_main_a0` recipe is. The multiplier is 37% of an
engine and its leaf products (`gf2_kmul`'s `leaf.r`, 3.9k LUTs) are the
single largest item; the memories are the next (3.1k LUTRAM in the step
unit, 1.7k in the walker). No DSPs or block RAM are used.

The VU47P has **1 303 680 LUTs** (the "2.85M" in the marketing sheet is
logic cells), so an engine is 1.01% of the device: **48 engines is half
the device**, 64 is 65%, and ~70% is where routing of a design this
regular typically starts to fail timing. At 5.3 clocks per step and the
250 MHz shell clock that is 47 M steps/s per engine and **2.3 G steps/s
at 48 engines, 3.0 G at 64**; a 400 MHz clock recipe, which the slack
above supports if implementation agrees, would give 1.6x that.

For scale, the measured client rate on an RTX PRO 6000 Blackwell is
6.9 G iterations/s (`ecc2k130/RTX-PRO6000.md`), reached by bitslicing 32
walks per word and spending 8859 instructions per multiply. The GPU has no
carryless multiply; the FPGA is nothing but carryless multiplies. That is
the structural reason a binary Koblitz curve suits an FPGA where the prime
field secp256k1 (`docs/ecc_fpga_cost_model.md`) does not: there the FPGA
advantage was "roughly 2x, fragile"; here the multiplier costs no DSPs, the
inversion costs eight multiplies, and a step costs 5.3.

Two things the first synthesis taught, both fixed:

- **One writer per array, really.** The walker's step counters were
  written from the load path and from the retire path; Vivado built them
  from 8 192 flip-flops behind a 256:1 read mux, 25k LUTs per engine, and
  the hierarchical report blamed the step unit. Loads now wait for a clock
  with no retiring step so one `if/elsif` chain writes every array.
- **One address per read port.** Each distinct `array(expr)` in the
  operand stage became its own LUTRAM copy: the tree had eight, `x` five.
  Stage A now forms one address per port and the phase selects which
  port feeds an operand; `d = x + sigma^j(x)` is formed once at fill
  instead of at each of its three uses. Step unit LUTRAM went from 6 748
  to 3 076 and three of five `sigma^j` instances disappeared.

## What is not here

- **A wider inversion share.** `W = 32` takes the bound from 5.31 to 5.16
  multiplies per step for twice the memory per batch; the generic is
  there and passes, the default stays at 16 because the last 3% costs
  as much RAM as the first 47% did.
- **A second multiplier per step unit.** The scheduler issues one product
  per clock; doubling that means two multipliers behind one ready queue
  and a two-port tree. The same throughput comes for free from
  instantiating two step units, which is the plan.
- **A routed image.** The custom logic has been synthesised in the F2
  CL flow (the numbers above) but no build has yet completed place and
  route or run on an F2 instance; post-route Fmax and the `NENG` that
  closes timing are what `aws/build_afi.sh` is for.
- **Reading a walk back.** The engine's walk state is write-only from the
  host, so walks in flight are lost when a worker restarts; each of them
  is at most `2^DP_WEIGHT` steps of work and the fresh seeds make up for
  it. It also means no per-walk `--max-iters`: a walk that never reaches a
  distinguished point (a cycle before the first one — at a cutoff of
  `2^34` on a walk whose expected cycle length is near `2^65`, that is
  astronomically rare) is only cleared by dropping `CTRL.RUN`.
- **Restart on the FPGA.** Building a fresh start point costs 128 point
  additions with `sigma^i(P)`, once per report — roughly one in `2^25`
  steps at the real cutoff. The client keeps that off the hot path too.
- **The trace-zero check on load.** The walker trusts the host to load
  points of order `l`.
