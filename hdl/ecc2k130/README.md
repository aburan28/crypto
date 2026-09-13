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
GHDL, reads into any synthesis flow (the DSP leaves of the multiplier are
unsigned integer products the synthesiser puts in DSP48E2 blocks). Every
number below is labelled measured (in simulation or on the device),
derived, or estimated. **The 96-engine image runs on an `f2.6xlarge` at
6.02 G steps/s**, the 80-engine one at 5.02 G, the 64-engine one at 4.01
G and the 48-engine one at 3.01 G (`aws/README.md`, "What came back"):
333 MHz / 5.31 clocks per step per engine, with the distinguished points
sampled from each checked against the client's reference walk. The
current revision keeps the product tree in UltraRAM, counts eleven of the
multiplier's 27 leaves in DSPs (the engine is 7.3k LUTs from 8.1k) and
batches 32 walks 16 deep (5.22 clocks per step); 112-, 120-, 128- and
136-engine builds were in flight when this was written.

## Files

| File | What it is |
|---|---|
| `gf131_pkg.vhd` | The field: coordinates, `sigma^k`, `sigma^j` selection, two-clock weight, the two constant basis-change maps built at elaboration, and an independent direct-product oracle |
| `gf2_kmul.vhd` | Recursive Karatsuba polynomial multiplier over GF(2), `LEVELS` deep, schoolbook leaf, the first `DSP_LEAVES` leaves in DSPs, latency `2 LEVELS + LEAF_LAT` |
| `gf2_dsp_leaf.vhd` | The DSP leaf: a 17-bit GF(2) product as six 24 × 15 integer products with the coefficients three bits apart, plus the corner in LUTs; four clocks |
| `gf131_mul.vhd` | GF(2^131) multiplier around it, II = 1, latency 13 (10 with LUT leaves only) |
| `ec2k_batch_pipe.vhd` | **The step unit.** W walks per batch, their inversions shared through a product tree: `5 + 5/W` multiplies per step, bursts of independent multiplies streamed from a ready queue |
| `ec2k_step_pipe.vhd` | The simple step unit, one inversion per walk, 10 multiplies per step; kept as the readable reference |
| `ec2k_walker.vhd` | The rho sequencer: N walks, host load port, distinguished-point output, around one `ec2k_batch_pipe` |
| `ec2k_axil.vhd` | **The host interface.** AXI4-Lite register block around `NENG` walkers: load registers, a distinguished-point queue, step/point counters, geometry and clock readback |
| `ec2k_axil_cdc.vhd` | AXI4-Lite clock-domain bridge, so the block and the engines run on a faster clock than the host port's |
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
ec2k_batch_tb: 200 steps in 1330 clk = 665/100 clk per step, W = 16, 8 batches in flight (85 multiplies per batch)
ec2k_batch_tb: PASS
=== ec2k_step_tb ===
ec2k_step_tb: 200 steps in 2051 clk with 16 slots (10 multiplies per step)
ec2k_step_tb: PASS
=== ec2k_walker_tb ===
ec2k_walker_tb: 64 distinguished points from 797 steps in 19194 clk (32 walks, batches of 8, 4 in flight)
ec2k_walker_tb: PASS
=== ec2k_axil_tb ===
ec2k_axil_tb: 64 distinguished points, 797 steps through 2 engine(s) of 16 walks, 7672 status polls
ec2k_axil_tb: PASS

$ make rate
ec2k_batch_tb: 6144 steps in 34042 clk = 554/100 clk per step, W = 16, 8 batches in flight (85 multiplies per batch)
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
stage 5..8     27 products of 17 x 17: 11 in DSPs (4 clocks), 16 schoolbook (1, padded)
stage 9        level 3 post-combine
stage 10       level 2 post-combine
stage 11       level 1 post-combine                261-bit product
stage 12       r = to_onb(h)                       c-powers -> gamma
```

(With `MUL_DSP_LEAVES = 0` the leaf is one clock and the latency 10.)

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

`prep` is 991 XOR2 per operand, `to_onb` 3167. Binary-field arithmetic has
no carries, so the DSP's adder is useless to it — but its multiplier is
not, see "Leaves in DSPs" below.

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

### Leaves in DSPs

Once the engine count was bound by LUTs (128 engines are 83% of them with
the all-LUT multiplier) the device's 9 024 DSP48E2 blocks, all empty,
were the resource left. An integer product of two packed polynomials
counts the AND terms a GF(2) product would XOR:

```
A = sum a_u 2^(3u),  B = sum b_v 2^(3v)   =>   A B = sum_k c_k 2^(3k),
c_k = #{(u, v) : u + v = k, a_u = b_v = 1}
```

With the coefficients three bits apart and at most five on one side,
`c_k <= 5` never carries into the next field, and bit `3k` of the integer
product is the GF(2) coefficient for free. The DSP48E2 multiplies 27 × 18
signed, so one block covers eight coefficients against five (24 × 15
unsigned); `gf2_dsp_leaf` does a 17-bit leaf as a 2 × 3 grid of them
(`a[0..15]` × `b[0..14]`), A/B, M and P registered inside the block, plus
the corner the grid misses (`a[16]` against all of `b`, `a[0..15]` against
`b[15..16]` — one to three AND terms per output bit, one LUT6) a clock
earlier, and one XOR of the parities into place. Nothing in the source is
Xilinx-specific: an `unsigned * unsigned` GHDL simulates and Vivado puts
in a DSP (`use_dsp = "yes"`).

Leaves are numbered in instantiation order and the first `MUL_DSP_LEAVES`
(11) are DSP leaves; every leaf must retire on the same clock, so the
largest subtrees with no DSP leaf get their output delayed by the three
clocks (flip-flops, not SRLs — a LUT per bit would have cost what the
leaves save). Synthesised (`cl_probe`, 2 engines, 3.0 ns), **measured**:

| | LUTs / engine | FFs | DSPs | slack |
|---|---|---|---|---|
| 27 LUT leaves | 7 974 – 8 106 | 8 062 | 0 | +1.18 ns (probe) / +1.24 (engine) |
| **11 DSP leaves** | **7 170 – 7 318** | 9 546 | **66** | +1.18 / +1.24, no DSP path in the 80 worst |

A DSP leaf is six DSPs and ~70 LUTs against ~140: −788 LUTs per engine
(−9.7%). 128 engines take 8 448 DSPs (94%); 136 fit 10 leaves each
(`DSP_LEAVES=10` in `build_afi.sh`). The price is three clocks of
latency (13 from 10), which 8 batches in flight no longer hide (5.38
clocks per step at 32 × 8); 16 do — see "The batched step unit".

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

- every memory is a **hard memory block** with exactly one writer and one
  read address: the leaf table (`x, y, j, tag, valid`), written only by
  the fill, and `d = x + sigma^j(x)` stored twice at fill rather than
  formed at each use, in block RAM; the product tree twice (the tree
  levels read two entries per product), written only by retire, in
  **UltraRAM** — the tree is the largest array (`2W` words per batch) and
  the VU47P has 960 UltraRAM blocks the design otherwise leaves empty,
  while block RAM is what bounds the number of engines: 17 tiles per
  engine with the tree in block RAM, 13 with it in UltraRAM, plus four
  URAM288 whose depth (4096) is mostly unused — the ports are the resource,
  not the bits. A read is three clocks — address register, array, output
  register — and the tree write is registered too, so every address net
  is a register the placer can put beside the block (the UltraRAM columns
  are further from an engine's logic than its block RAMs). Nothing is read at retire:
  the final multiply of a walk, `lam · (x + x3)`, has `y`, the tag and
  `x3 = lam^2 + lam + d` in hand when it issues, and they ride a shift
  register beside the multiplier (SRLs, about 300 LUTs) and meet the
  product on the way out. (An earlier version kept a second leaf table
  and an `x3` table for the retire side to read through a look-ahead copy
  of the multiplier's tag: 4 RAMB36 + 1 RAMB18 of the engine's 20 + 2,
  and block RAM is what bounds the number of engines in a device.) No
  word is read, by a multiply that will use it, on the clock it is
  written — a simulation-only assertion checks it — so the read-during-
  write behaviour of the RAMs never matters;
- the port addresses are registered a clock before the read and operands
  are read a clock before they are formed, so no path runs level/index
  arithmetic → RAM or RAM → `sigma^j` → XOR in one clock, and the retire
  path is a memory write plus one XOR — there is no `sigma` and no
  popcount in it;
- the Hamming weight, 131 bits wide, is taken over two clocks on input
  (22 groups of six bits, one LUT6 per output bit, then the sum, of which
  only `hw/2 mod 8` is used) and four on output (groups, sums of four
  groups, their sum, compare — the 22-way sum in one clock was the
  engine's worst path);
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
| 8 | 8 | 8.02 | 6.62 | 5.63 |
| 16 | 4 | 9.86 | 8.38 | 5.31 |
| 16 | 8 | 6.98 | 5.67 | 5.31 |
| 16 | 16 | 6.46 | **5.34** | 5.31 |
| 32 | 8 | 7.20 | 5.20 | 5.16 |

Clocks per step, multiplier latency 10, three-clock memory reads (with
two-clock reads the same table read 6.83 / 7.57 / 5.54 / 5.32 / 5.24;
with LUTRAM's one-clock reads 5.95 / 7.16 / 5.49 / 5.32 / 5.21). With 8
batches in flight the long run sits 0.36 above the bound: about 500
clocks of the 34 800 are the fill of the first batch and the drain of the
last, the rest is the inversion's eight dependent single-multiply bursts,
each now thirteen clocks from issue to the next issue, which eight
batches do not cover and sixteen nearly do. Memory per batch of `W` is
`8W` field elements as stored (the leaf table's `x` and `y`, `d` twice,
the tree twice), and a block RAM is 72 × 512 whatever the design asks
for, so a 131-bit table of 128 or 256 entries costs the same two RAMB36
and the UltraRAMs are 4096 deep: any geometry holding 256 or 512 walks
costs 8 RAMB36 + 1 RAMB18 + 4 URAM288 per step unit. With the 10-clock
multiplier the image's default was 32 × 8: the same 256 walks and the
same memory as 16 × 16, and the walker testbench with 512 walks ran
**5.16 clocks per step against 5.31** — the bound is `5 + 5/W`. With the
13-clock multiplier (DSP leaves) eight batches no longer cover the
inversion chain and 32 × 8 reads 5.38; **the image's default is now
32 × 16** (`cl_ecc2k130_defines.vh`), all 512 of the engine's walks in
the step unit, at **5.22** (16 × 16: 5.31, 64 × 4: 6.06; all measured in
the walker testbench with 512 walks, `-gLOG_W=5 -gLOG_NB=4`). The
testbenches' default stays 16 × 8.

Degenerate inputs (`d = 0`, i.e. `sigma^j(x) = x`) are not special-cased,
matching the client: the chain returns `1/0 = 0`, the product tree zeroes
and every step in the batch produces garbage, as the client's batched
inversion does. On the challenge curve the case has probability about
`2^-131` per step.

## The walker

`ec2k_walker` is the sequencer `hdl/ecc/` leaves as "not here". Walks live
in a ready FIFO of `(id, steps, x, y, dp)`: a walk is either inside the
step unit or waiting in the FIFO, and a completed step goes straight back
into the FIFO with its count plus one, flagged if its `x` has weight at most
`DP_WEIGHT`; a flagged walk reaching the head is reported on the `dp_*`
port with its id and step count and goes idle until the host loads a new
start point into that id. The step count rides through the step unit in
the tag rather than in a per-walk counter memory, so the retire path has
no read-modify-write of a RAM and there is no state indexed by walk id at
all. The FIFO is block RAM — 512 words of 304 bits, four RAMB36 and a
RAMB18 — read through a four-entry prefetch buffer that hides the two-clock
read, so the step unit sees a walk every clock it can take one. That is
the same division of labour as the GPU client, whose kernel reports
`(seed, endpoint)` and whose host owns restarts, the corpus and collision
resolution (`ecc2k130/README.md`, "Distinguished points and restarts").
`NWALK` should comfortably exceed the `W · 2^LOG_NB` walks the step unit
holds, so a full batch is always forming; when it does not (start-up, a
host slow to reload) the flush keeps things moving at reduced efficiency,
never a stall. The default is 512 walks (`ID_W = 9`) over 256 in the step
unit, and `ec2k_walker_tb -gRATE_CLK=40000`, which keeps every id loaded,
measures **5.29 clocks per step** through the walker (5.54 with 8 batches
in flight and 256 walks).

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
  an id after that id has reported. A load while `RUN` is clear is dropped.
- **Report.** A distinguished step goes back into the walker's ready FIFO
  flagged `dp`; when it reaches the head it is moved into the engine's
  report register instead of the step unit and held there, valid, until
  the register block acknowledges it. Reports then travel to a
  `2^DP_FIFO_W`-deep queue the host reads through `DP_ID`, `DP_STEPS`,
  `DP_X*`, `DP_Y*` and pops with `DP_POP`. Nothing in this chain can lose
  a report — a full queue or a slow host only delays them, with the
  walker's FIFO as the backlog — so `DROPPED` and `STATUS.DP_OVERFLOW`
  read 0 and exist for register-map compatibility.
- **Count.** `STEPS` (64-bit, high word latched on the low read) is the
  total steps of every engine; `DPS` the reports queued. `GEOM` reads back
  `ID_W`, `LOG_W`, `LOG_NB`, `DP_WEIGHT` and `NENG`, so the host can refuse
  an image built for a different campaign; `CLOCK` the engine clock in
  kHz, so the operator can see which image is loaded; `MAGIC` is
  `0x2C130001`.

The engines hang off a **spine**, one register stage per engine, chained,
because with 48 or 64 engines over the VU47P's three SLRs it is the
interconnect and not the arithmetic that breaks timing: the first 48-engine
build, with the load registers fanned out to every engine and a 48:1 mux
of reports back, placed at −2.5 ns slack on a 4 ns clock while the engine
alone had +2 ns. On the spine every wire runs register to register between
neighbours. Down it, a load `(gid, x, y)` shifts outward a stage per clock
and the stage whose engine it names keeps a copy and offers it until taken;
up it, a report slot `(gid, steps, x, y)` shifts inward and an engine drops
its report into an empty slot passing by, with the load acknowledgement
and a running sum of step pulses riding beside it (so `STEPS` needs no
popcount over the die either). The chain never stalls and never drops, by
credit: the block issues one credit per free queue slot not already spoken
for, an engine with a report waiting takes the first credit that passes
and only then inserts, and a credit that reaches the far end unused comes
back up to be reissued — queue occupancy plus outstanding credits never
exceeds the queue depth. It costs about 870 flip-flops and 300 LUTs per
engine, a few percent, and `NENG` clocks of latency on events that happen
a few hundred times a second.

`CTRL.RUN` holds the engines in reset while clear — the reset sweeps down
the spine a stage per clock and the block ignores the chain for
`2·NENG + 8` clocks after `RUN` rises — and `CTRL.CLEAR` zeroes the
counters and the queue. The register map with byte offsets is at the top of
`ec2k_axil.vhd`. `ec2k_axil_tb` replays the walker testbench's 64 walks
through the registers alone, spread over two engines (one, three and four
have been run), with an 8-deep queue so the credits are exercised, and
checks the same per-id bookkeeping.

`ec2k_axil_cdc` is the bridge that lets the block run on a clock of its
own: an AXI-Lite slave on the host's clock and a master on the engines',
one write and one read in flight, each a four-phase level handshake whose
flags cross through `ASYNC_REG` synchronisers and whose data is written a
clock before its flag and copied only after the flag has been seen. The
register names spell out what crosses (`cdc_*_meta`, `cdc_*_src`,
`cdc_*_cap`) so the constraints can be two lines. It costs 250 LUTs and
about 50 ns per access, and `ec2k_axil_tb -gCDC=true` runs the same test
through it with the engines on a 7 ns clock against the host's 10 ns
(`make test` runs both; ratios from 3 to 21 ns have been run).

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

Vivado 2025.2, `xcvu47p-fsvh2892-2-e`, out of context, two `ec2k_walker`
behind the register block and the clock bridge (512 walks, W = 16, 16
batches in flight):

| | LUTs | of which LUTRAM | FFs | RAMB36 | RAMB18 | URAM | DSP |
|---|---|---|---|---|---|---|---|
| `ec2k_walker` (whole engine) | **7 170 – 7 320** | 180 (+614 SRL) | 9 550 | 12 | 2 | 4 | 66 |
| ├ walker body (FIFO, prefetch buffer, held report) | 366 | 176 | 337 | 4 | 1 | 0 | 0 |
| └ `ec2k_batch_pipe` | ~6 950 | 4 | ~9 200 | 8 | 1 | 4 | 66 |
| &nbsp;&nbsp; ├ step unit body (scheduler, operand and retire stages, the final multiply's shift register) | 2 630 | 4 | 2 994 | 8 | 1 | 4 | 0 |
| &nbsp;&nbsp; └ `gf131_mul` (three-level Karatsuba, 11 DSP leaves) | ~4 330 (5 117 with LUT leaves) | 0 | ~6 100 | 0 | 0 | 0 | 66 |
| `ec2k_axil` own (queue and its head registers, AXI, spine head; two stages) | 890 | 356 | 2 391 | 0 | 0 | 0 | 0 |
| `ec2k_axil_cdc` | 33 – 253 (the rest is merged into the block's read mux) | 0 | 201 | 0 | 0 | 0 | 0 |

(The two engines differ by 150 LUTs from `-keep_equivalent_registers`
falling differently; the body rows are from a hierarchical synthesis of
one engine with LUT leaves, where the multiplier was 5 117 LUTs: 4 428 in
the Karatsuba tree, of which 3 800 in the 27 leaf products, and 689 in
`prep`, `to_onb` and the tag pipe.) The register block's stage on the
spine is roughly 300 LUTs and 870 FFs per engine; the 180 LUTRAM left in
an engine are the walker's four-word prefetch buffer. Before the final
multiply's companions moved into a shift register the engine was 8 100 –
8 420 LUTs, 322 SRL, 7 420 FF and **20 RAMB36 + 2 RAMB18**: 150 LUTs and
260 SRLs bought four RAMB36 tiles. With the tree in block RAM it was
8 260 – 8 570 LUTs, 7 660 FF and 16 RAMB36 + 2 RAMB18 (the routed 48-,
64-, 80- and 96-engine images); the UltraRAM tree with its address and
write registers costs 220 FFs and buys four more tiles, 32-walk batches
and the enable-free stages take 250 LUTs back (7 970 – 8 110), and the
eleven DSP leaves another 790 for 1 480 FFs and 66 DSPs.

Worst slack at a 3.0 ns clock is **+1.18 ns** for the whole probe (the
register block's read mux, one copy on the die) and **+1.24 ns** inside
an engine, and none of the eighty worst paths touches a memory: the
output weight's second-stage sum (1.58 ns) and the FIFO's write enables
(1.17 ns) lead. The routed 80- and 96-engine images had every worst path
on a control net into a stage's data registers — the insert-mux select,
the input stage's valid, the operand stage's valid, 310 – 440 loads
each, 2.8 – 2.9 ns of route — so the operand, multiplier, retire and
leaf-write stages now move data every clock with no enable at all, and
every remaining control net over 250 loads carries a fanout limit so its
driver is replicated among the loads (`report_high_fanout_nets` on the
probe: nothing per engine above 180 now). Earlier worst paths and what removed
them: the operand stage's level/index arithmetic into a LUTRAM read
(+0.83 ns; an address stage); the walker's step-counter read-modify-write
(1.6 ns; the count now rides in the tag); the 22-group weight sum in one
clock (1.74 ns; sums of four, then of six); the queue's pointer
subtraction into the credit compare (1.70 ns; occupancy is a counter); the
bridge's address into the read mux and the queue's LUTRAM read into it
(1.67 ns; the address is registered in the block, the queue head is read
into registers a clock ahead, and a write no longer overlaps its own
response); the retire side's batch-level update, a 16:1 mux of the level
array into a decrement and back (1.56 ns; the level rides in the
multiplier tag too, so retire writes the array and never reads it). The
multiplier is now 59% of an engine and its sixteen LUT leaf products
(~2.2k LUTs) the single largest item; the other eleven leaves are the
engine's 66 DSPs. The one distributed RAM left is the register block's
64-deep report queue (356 LUTRAM, one copy on the die); the UltraRAMs
hold only the product trees.

Why block RAM and not the distributed RAM the first revisions used: the
first 48-engine implementation (`aws/README.md`) placed at −2.5 ns on a
4 ns clock with a quarter of a millisecond of total negative slack, and
the paths that failed were the same in every engine — the write address
and write enable of each LUTRAM array, a net of 1 000 to 1 800 loads that
the placer spreads over every SLICEM of a clock region and that physical
optimisation cannot split, because the loads are memory primitives, not
logic. 4 800 LUTRAM per engine was never going to place tightly 48 times
over. A block RAM has the same fanout inside one hard block.

The VU47P has **1 303 680 LUTs** (the "2.85M" in the marketing sheet is
logic cells), **2 016 RAMB36**, **960 URAM288** and **9 024 DSP48E2**, so
an engine is 0.56% of the LUTs, 0.64% of the block RAM (13 tiles: 12
RAMB36 and two RAMB18 sharing one), 0.42% of the UltraRAM and 0.73% of
the DSPs. With the tree in block RAM (17
tiles) **48 engines was a third of the device** (the routed 48-engine
image of the 21-tile revision used 31% of the LUTs and 50% of the RAM),
64 was 42% of the LUTs and 54% of the RAM, 80 was 53% and 67% (the routed
80-engine image: 686 854 LUTs, 1 360 tiles) and 96 63% and 81% — block RAM
ran out first. With the tree in UltraRAM 96 engines is 63% of the LUTs,
62% of the block RAM and 40% of the UltraRAM, 112 is 73% / 72% / 47%, 128
is 83% / 83% / 53%, and the LUTs bound the count. With eleven leaves in
DSPs 128 engines is 74% of the LUTs, 83% of the block RAM, 53% of the
UltraRAM and 94% of the DSPs; 136 engines with ten leaves each is 80% /
88% / 57% / 90%. At 5.31 clocks per step and 333 MHz that is 63 M
steps/s per engine and **3.0 G steps/s at 48 engines, 4.0 G at 64, 5.0 G
at 80, 6.0 G at 96 — all four measured on the device**; at 5.22 clocks
per step (32 × 16) 63.9 M per engine, 8.2 G at 128 and 8.7 G at 136; at
the shell's 250 MHz the measured four would be 2.3, 3.0, 3.8 and 4.5 G.

For scale, the measured client rate on an RTX PRO 6000 Blackwell is
6.9 G iterations/s (`ecc2k130/RTX-PRO6000.md`), reached by bitslicing 32
walks per word and spending 8859 instructions per multiply. The GPU has no
carryless multiply; the FPGA is nothing but carryless multiplies. That is
the structural reason a binary Koblitz curve suits an FPGA where the prime
field secp256k1 (`docs/ecc_fpga_cost_model.md`) does not: there the FPGA
advantage was "roughly 2x, fragile"; here a multiplier is 4.3k LUTs and 66 DSPs
that count AND terms, the inversion costs eight multiplies, and a step
costs 5.2.

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

- **A second multiplier per step unit.** The scheduler issues one product
  per clock; doubling that means two multipliers behind one ready queue
  and a two-port tree. The same throughput comes for free from
  instantiating two step units, which is the plan.
- **A campaign on F2 instances.** The 48-, 64-, 80- and 96-engine images
  run on an `f2.6xlarge` at the predicted rates (3.01, 4.01, 5.02 and
  6.02 G steps/s, `aws/README.md`, "What came back"); what has not happened yet is a
  fleet of workers feeding the shared corpus for hours, which is `f2.sh
  up` once the account has the worker instance profile.
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
