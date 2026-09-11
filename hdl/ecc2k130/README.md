# hdl/ecc2k130 — the rho sequencer for the ECC2K-130 attack

The arithmetic for this attack is not written by hand: it is emitted from the
same straight-line IR the CUDA client is generated from, by
[`ecc2k130/codegen/genverilog.py`](../../ecc2k130/codegen/genverilog.py), and
checked against the same field model. What lives here is the part that is
control rather than arithmetic — the sequencer around that datapath.

The cost model both halves belong to is
[`ecc2k130/FPGA-CEILING.md`](../../ecc2k130/FPGA-CEILING.md), which asks whether
an FPGA beats one GPU on this walk. Short answer: on speed it is roughly a wash,
on cost per solved instance about 2x, on energy 5–10x.

## Files

| File | What it is |
|---|---|
| `ecc_rho_core.v` | One rho core: walks circulate, one step retires per clock |
| `tb_ecc_rho_core.v` | Self-checking testbench against reference trajectories |
| `Makefile` | Builds and runs the testbenches under Icarus Verilog |

The generated modules it instantiates — `ecc_pre131`, `ecc_post131`,
`ecc_inv131`, `ecc_mul131`, `ecc_delay` — come from
`../../ecc2k130/generated/rtl/`.

## Run

```bash
make -C ../../ecc2k130 verilog    # emit the datapath and the vectors
make                              # build and run the testbenches
```

Needs Icarus Verilog (`brew install icarus-verilog`, `apt install iverilog`).
Current output:

```
=== tb_ecc_rho_core ===
PASS: 124 walks x 3 steps retired, II=1, ring latency 124 clk, 517 distinguished
      points, all trajectories match curves.Curve
```

## The core

```
load ->|                                                        |
       +--> ecc_pre --> ecc_inv --> [state RAM] --> ecc_post ----+
              3 clk       96 clk        1 clk         24 clk
```

`R' = sigma^j(R) + R` with `j = 3 + ((HW(x_R)/2) mod 8)`, the iteration function
of [`ecc2k130/include/walk.h`](../../ecc2k130/include/walk.h).

There is **no state machine**. Every generated block has initiation interval 1
and no back pressure, so the ring the published FPGA designs describe is
literally a ring: a walk re-enters `ecc_pre` 124 clocks after it left, 124 walks
are resident, and one step retires per clock.

Three things carry the design:

- **The tag.** Each block carries an opaque tag beside its operands and never
  looks at it, so one pipelined block serves many independent walks and results
  route back to the right walk with no matching logic and no stalls. Same idea
  as [`../ecc/ec_add_pipe.vhd`](../ecc/ec_add_pipe.vhd).
- **A RAM, not delay lines.** `ecc_post` needs x, y, d and e when `d^-1` arrives
  96 clocks later. Registers would cost 50,304 flip-flops per core; a RAM
  indexed by the tag costs two block RAMs, and flip-flops are the second
  scarcest resource on the target part.
- **Pass 1 is free.** Splitting the iteration the way `walk.h` does puts the
  weight tree, the three conditional squarings and the addends in `ecc_pre`
  (2,887 gates, no multiplication) and only the two multiplications in
  `ecc_post` (33,641 gates). Squaring is a coordinate permutation in this basis,
  so it is wiring.

## What is not here

- **Montgomery batching.** This core inverts every denominator directly: 10
  multiplications per step instead of the 5.125 that batching reaches, so it
  costs about twice the LUTs per unit of throughput that the cost model
  projects. It is correct, and it is the reference the batching engine has to
  stay bit-exact against.
- **Host IO.** The report FIFO's PCIe or HBM side, and the start-point
  generator. Building a start point is 128 point additions, so the published
  designs have the host supply them; this core does the same.
- **Multi-core.** One inverter is meant to serve many cores, which is why it is
  a port here rather than an instance. Instantiating N cores around one inverter
  is wiring, and is not done.
- **Synthesis.** Nothing here has been through a synthesis tool, so there are no
  LUT counts, no timing closure and no power figures — only gate counts and
  simulated latencies.
