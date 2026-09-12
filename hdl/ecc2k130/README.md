# hdl/ecc2k130 — FPGA engine for the ECC2K-130 walk

The ECC2K-130 iteration function in portable VHDL-2008: a pipelined
GF(2^131) multiplier in the permuted type-II optimal normal basis, a step
unit that computes `R + sigma^j(R)` around it, and a sequencer that runs
many walks through the step unit and reports distinguished points. It is
the hardware counterpart of the GPU client in [`ecc2k130/`](../../ecc2k130/):
same field, same basis, same iteration function, same distinguished-point
rule, and the same oracle — every vector the testbenches check against
comes from the client's own field model in `ecc2k130/codegen/`.

Reference part is the AMD Virtex UltraScale+ VU47P (AWS `f2.6xlarge`), the
same one `hdl/ecc/` and `hdl/sha1/` target, but there are no vendor
primitives: plain `ieee.std_logic_1164` and `numeric_std`, simulates under
GHDL, reads into any synthesis flow. Nothing here has been synthesised or
run on an FPGA; every number below is labelled measured (in simulation),
derived, or estimated.

## Files

| File | What it is |
|---|---|
| `gf131_pkg.vhd` | The field: coordinates, `sigma^k`, `sigma^j` selection, weight, the two constant basis-change maps built at elaboration, and an independent direct-product oracle |
| `gf131_mul.vhd` | GF(2^131) multiplier, II = 1, latency 11, no DSPs |
| `ec2k_step_pipe.vhd` | One walk step around one shared multiplier: inversion (8 multiplies) plus addition (2), many steps in flight |
| `ec2k_walker.vhd` | The rho sequencer: N walks, host load port, distinguished-point output |
| `gf131_mul_tb.vhd`, `ec2k_step_tb.vhd`, `ec2k_walker_tb.vhd` | Self-checking testbenches |
| `gf131_tb_pkg.vhd` | Hex helpers shared by the testbenches |
| `ecc2k_ref.py` | Reference model of every hardware algorithm, proved against `ecc2k130/codegen/field.py`; generates the vectors |
| `vectors_ecc2k130.txt` | Generated: 256 products, 64 steps, 40 walks to a distinguished point |
| `aws/` | Out-of-context synthesis script for the VU47P |

## Run

```bash
make            # analyse, elaborate and run the three testbenches
make check      # prove the hardware algorithms against the client's field model
make vectors    # regenerate the vector file
make cost       # gate counts of the multiplier from its constant matrices
```

Needs GHDL and Python 3 (the generator imports `ecc2k130/codegen`). Current
output:

```
=== gf131_mul_tb ===
gf131_mul_tb: PASS -- 256 multiplications, II=1, latency 11 clk
=== ec2k_step_tb ===
ec2k_step_tb: 64 steps in 665 clk with 16 slots (10 multiplies per step)
ec2k_step_tb: PASS
=== ec2k_walker_tb ===
ec2k_walker_tb: 40 distinguished points from 576 steps in 7639 clk (16 walks, 16 slots)
ec2k_walker_tb: PASS
```

The multiplier vectors include zero, one (which is the all-ones vector in
this basis), single basis elements and their squares, and the top and
bottom coordinates alone; each product is checked against the file and
against `gf_mul_ref`, the direct normal-basis convolution computed in
simulation, which shares nothing with the multiplier's route. The step and
walk vectors are random points of order `l` on the challenge curve itself.

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
stages 2..9    h ^= a' * b'[17-bit digit r] << 17r   8 rows
stage 10       r = to_onb(h)                       c-powers -> gamma
```

The rows are operand scanning over GF(2)[c] with no carries: a row is
131 × 17 AND gates and an XOR tree five deep, two or three LUT6 levels.
That is what keeps the clock short without splitting anything.

Gate count per multiplier, read off the constant matrices (`make cost`), **derived**:

| | AND | XOR2 | Widest XOR |
|---|---|---|---|
| prep, both operands | | 1982 | 66 inputs |
| 131 × 131 product | 17161 | 16900 | |
| to_onb | | 3167 | 44 inputs |
| **total** | **17161** | **22049** | |

No DSP48 is involved anywhere — binary-field arithmetic has no carries, so
the DSP's adder is useless to it and the multiplier is pure LUT fabric.
Folding three AND terms and their XORs into one LUT6 and the remaining XORs
five per LUT gives roughly **8k LUTs and 5k flip-flops per multiplier,
estimated**. Compare `hdl/ecc/`: the secp256k1 multiplier is 256 DSPs plus
~10k LUTs for the same II = 1.

The `tag` rides alongside and the datapath never looks at it, so one
multiplier serves any number of independent contexts back-to-back.

## The step unit

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
multiplies and nothing else**.

Ten dependent multiplies cannot fill an 11-deep pipeline, so the unit keeps
up to `2^SLOT_W` independent steps in flight, one per slot. Scheduling is
dataflow: when a product retires, the slot's next operands are formed and
the slot is pushed on a ready queue; each clock the head of the queue
issues. New requests issue only when the queue is empty, so in-flight work
has priority and the multiplier never idles while anything is ready. Slot
and chain state ride in the multiplier tag; results route back with no
matching logic.

Measured in simulation (`ec2k_step_tb`, 64 independent steps back-to-back):

| Slots | Clocks | Clocks per step |
|---|---|---|
| 8 | 1051 | 16.4 |
| 16 | 665 | **10.4** |
| 32 | 665 | 10.4 |

With `MUL_LATENCY + 2 = 13` or more steps in flight the multiplier is
saturated and a step costs exactly its ten multiplies. The default of 16
slots is the smallest power of two that does it.

Degenerate inputs (`d = 0`, i.e. `sigma^j(x) = x`) are not special-cased,
matching the client: the chain returns `1/0 = 0` and the step produces
garbage, as the client does. On the challenge curve the case has
probability about `2^-131` per step.

## The walker

`ec2k_walker` is the sequencer `hdl/ecc/` leaves as "not here". Walks live
in a ready FIFO of `(id, x, y)`: a walk is either inside the step unit or
waiting in the FIFO, and a completed step goes straight back into the FIFO
unless its `x` has weight at most `DP_WEIGHT`, in which case it is reported
on the `dp_*` port with its id and step count and the walk goes idle until
the host loads a new start point into that id. That is the same division
of labour as the GPU client, whose kernel reports `(seed, endpoint)` and
whose host owns restarts, the corpus and collision resolution
(`ecc2k130/README.md`, "Distinguished points and restarts").

`ec2k_walker_tb` plays host: it loads 40 start points sixteen at a time,
checks every report's id, step count and point against the trajectory the
client's reference model walked from the same start (under a loose cutoff
of 56 so walks report within a few dozen steps), and reloads. Walks
interleave arbitrarily inside the engine; the per-id bookkeeping makes the
check independent of that. The clocks-per-step figure the walker testbench
prints is not a throughput measurement: the walk population dwindles as
records run out, and the step testbench is the one that measures the
saturated rate.

## Capacity and what it would mean — estimated

Per VU47P (2.85M LUTs), at ~12k LUTs per step unit including its
multiplier, slot memories and `sigma` multiplexers: **150–200 step units**,
LUT-bound, no DSPs or BRAM binding. At 10 clocks per step and 300–400 MHz
that is 30–40 M steps/s per unit and **5–8 G steps/s per FPGA**.

For scale, the measured client rate on an RTX PRO 6000 Blackwell is
6.9 G iterations/s (`ecc2k130/RTX-PRO6000.md`), reached by bitslicing 32
walks per word and spending 8859 instructions per multiply. The GPU has no
carryless multiply; the FPGA is nothing but carryless multiplies. That is
the structural reason a binary Koblitz curve suits an FPGA where the prime
field secp256k1 (`docs/ecc_fpga_cost_model.md`) does not: there the FPGA
advantage was "roughly 2x, fragile"; here the multiplier costs no DSPs and
the inversion costs eight multiplies.

None of this is measured. `aws/` synthesises the three modules on an AWS
build instance and brings back LUT, FF and Fmax; the estimates above are
what those numbers replace. The likely critical path is not the multiplier
but the step unit's retire path — a slot-memory read, `gf_sigma_j`, XORs,
write-back in one clock — and splitting it costs one clock of step latency
and nothing in throughput.

## What is not here

- **Montgomery batching.** The GPU client shares one inversion across 32
  walks, so its inversion cost per step is `3 + 8/32` multiplies instead of
  8. Doing the same here would take a step from 10 multiplies to about 5,
  close to a 2x throughput gain for more slot state and a longer chain. It is the
  single highest-value change and it is not implemented.
- **Karatsuba.** One level (66 + 65) turns 17161 AND gates into 13068, a
  24% saving on the product; the client's generator found one level to be
  the sweet spot on the GPU. The rows are a drop-in place to put it.
- **The AWS shell.** Like `hdl/ecc/`, this is a bare datapath with a clock,
  a load port and a report port. A loadable FPGA image needs the `aws-fpga`
  custom-logic wrapper, AXI interfaces and a host program; see
  `hdl/ecc/aws/README.md` for what that involves.
- **Restart on the FPGA.** Building a fresh start point costs 128 point
  additions with `sigma^i(P)`, once per report — roughly one in `2^25`
  steps at the real cutoff. The client keeps that off the hot path too.
- **The trace-zero check on load.** The walker trusts the host to load
  points of order `l`.
