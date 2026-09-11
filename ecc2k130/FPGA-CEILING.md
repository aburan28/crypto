# Can one FPGA beat one GPU on the ECC2K-130 walk?

[THROUGHPUT-CEILING.md](THROUGHPUT-CEILING.md) shows that no arithmetic reaches
15–20 B iterations/s on one RTX PRO 6000, because the walk's bit operations all
issue through one 64-lane integer pipe. An FPGA has no such pipe: the same bit
operations are laid down in space and every one of them executes on every
clock. That is the reason to ask the question, and it is a real structural
advantage — but it buys less than it sounds like.

**The answer: roughly a wash per chip, about 2x cheaper per solved instance,
and about 5–10x better per watt.** A VU47P is estimated here at 5–12 B
iterations/s against the GPU's measured 6.9 B/s. One FPGA does not clear the
GPU's 17.6 B/s instruction ceiling, and the arithmetic — five multiplications
and a batched inversion per step — binds both platforms equally.

Numbers below are labelled **measured** (from circuits in this repository,
verified against the field model), **derived** (arithmetic on measured values
and published part or design parameters), or **estimated** (engineering
judgement). Nothing here has run on an FPGA. There is no FPGA toolchain in this
environment: no synthesis, no place-and-route, no timing closure, no power
measurement.

## 1. Why this is not the secp256k1 answer

[`docs/ecc_fpga_cost_model.md`](../docs/ecc_fpga_cost_model.md) costs a Pollard
rho on secp256k1 and finds the FPGA advantage is "roughly 2x, fragile". That
conclusion does not transfer, because it is about a **prime** field: the
primitive is a 256-bit modular multiply, the FPGA realisation is DSP-bound, and
a GPU's integer multiply-add pipes are built for exactly that.

ECC2K-130 is a **binary** field. `GF(2^131)` multiplication is XOR and AND — no
carries, no DSPs — which is the regime of
[`docs/sha1_fpga_cost_model.md`](../docs/sha1_fpga_cost_model.md), where LUTs
bind and the FPGA wins by more than 10x. The GPU has to *emulate* those bit
operations with general-purpose word instructions, and
[THROUGHPUT-CEILING.md](THROUGHPUT-CEILING.md) measures the price: 27.1 T
lane-ops/s across the whole card, about 1,540 lane-ops per update, hence 17.6
B/s as an absolute bound with perfect scheduling and free memory.

So ECC2K-130 sits on the SHA-1 side of the repository's own dichotomy. The
three cases together:

| Workload | Primitive | FPGA binding resource | FPGA advantage |
|---|---|---|---|
| SHA-1 collisions | 32-bit add/rotate/bool | LUTs | >10x |
| **ECC2K-130 rho** | **`GF(2^131)` multiply (XOR/AND)** | **LUTs** | **~1x speed, ~2x cost, ~5-10x energy** |
| secp256k1 rho | 256-bit modular multiply | DSPs | ~2x |

The middle row does not reach the top row because the binary-field *multiply*
is not free the way a SHA-1 round is: it is thousands of gates, and five of
them are needed per iteration.

## 2. The walk as gates (measured)

[benchmarks/fpga-lut/lutmap.py](benchmarks/fpga-lut/lutmap.py) maps the
straight-line programs the CUDA client is generated from onto 6-input lookup
tables, with a priority-cut mapper (the standard FlowMap construction: enumerate
cuts, keep the best ten per node, choose a cover). Every multiplier is verified
against `field.Onb` on 64 random inputs before it is counted. Raw output is
[benchmarks/fpga-lut/result.json](benchmarks/fpga-lut/result.json).

| Circuit | 2-input gates | Gate depth | LUT6, depth-opt | LUT6, area-opt |
|---|---:|---:|---:|---:|
| multiply, direct ONB | 42,706 | 10 | 10,970 @ 6 levels | 10,970 @ 6 |
| multiply, conversion cut 131 | 36,001 | 202 | 14,539 @ 54 | 9,515 @ 85 |
| multiply, conversion cut 66 | 27,940 | 139 | 11,337 @ 36 | 7,661 @ 54 |
| multiply, conversion cut 33 | 22,119 | 109 | 8,902 @ 28 | 6,344 @ 40 |
| multiply, conversion cut 17 | 18,424 | 95 | 7,344 @ 25 | 5,659 @ 33 |
| multiply, conversion cut 9 | 16,494 | 91 | 6,787 @ 24 | **5,563 @ 32** |
| hamming weight (8 bits) | 268 | 21 | 79 @ 5 | 91 @ 13 |

"Conversion" is the Bernstein-Lange multiplier the client already uses: ONB to
optimal polynomial basis, Karatsuba down to the stated leaf, convert back.
"Direct ONB" is the quadratic normal-basis alternative, built here for
comparison from `gamma_i gamma_j = gamma_{i+j} + gamma_{i-j}`.

**Cross-check.** The direct ONB multiplier can be counted by hand: each of 131
outputs is a XOR of 131 terms `a_i & (b_u ^ b_v)`, each term needs three inputs,
so one LUT6 covers two terms — 66 LUTs per output, then a fan-in-6 XOR tree
(66 → 11 → 2 → 1), about 10,350 LUTs at 4–5 levels. The mapper reports 10,970 at
6 levels, 6% above the hand count. That is the expected direction and size of
error for a cut-based mapper and is the reason to trust the other rows.

**The mapper's counts are an upper bound.** It does not model `LUT_6_2`
dual-output packing, which real designs exploit heavily: Bernstein et al. report
a hand-packed `GF(2^113)` multiplier at **3,071 LUTs** and `GF(2^127)` at
**3,620 LUTs**, where this mapper would report roughly 40% more. Packing quality
is therefore worth about 40% of everything downstream, and it is the single
biggest uncertainty in this document.

**Depth is not the FPGA's problem; area is.** The GPU wants the fewest
instructions and does not care about depth, since a straight-line program has no
depth cost — hence leaf 9–17. An FPGA pays for gates in area and for depth in
pipeline registers, and the measurement says the register cost is affordable:

| Multiplier | LUTs | Pipeline at 3 LUT levels/stage | Flip-flops |
|---|---:|---|---:|
| conversion cut 17 | 5,659 | 9 stages | 8,182 |
| conversion cut 9 | 5,563 | 8 stages | 6,542 |
| direct ONB | 10,970 | 2 stages | 3,143 |

The VU47P carries two flip-flops per LUT (2,607,360 against 1,303,680), so
trading 2x the registers for half the LUTs is the right trade. **The conversion
multiplier the client already generates is also the right FPGA multiplier** —
the shallow quadratic one costs twice the area to save pipeline registers that
are not scarce. At two levels per stage rather than three, though, the design
becomes flip-flop-limited (§4), so three is the operating point.

## 3. What one core costs (derived)

The published designs are the calibration. Bernstein, Engels, Lange,
Niederhagen, Paar, Schwabe and Zimmermann
([eprint 2016/382](https://eprint.iacr.org/2016/382)) build an unrolled core
that computes **one iteration per clock**, with multipliers pipelined at II = 1
and independent walks interleaved — the same idea as `hdl/ecc/ec_add_pipe.vhd`.
Inversion is not unrolled: a batch of walks shares one high-latency inverter
through Montgomery's trick and a dual-buffer memory, which is the hardware form
of what `ECC_BATCH` does in the CUDA client. Their published area, post-place
and route:

| Design | Field | Multiplier | Multipliers | Total LUTs | Per core |
|---|---|---:|---:|---:|---:|
| 3 cores, XC6SLX150 | `GF(2^113)` | 3,071 | 16 | 63,388 | 21.1 K |
| 3 cores, XC6SLX150 | `GF(2^127)` | 3,620 | 16 | 72,919 | 24.3 K |
| 7 cores, XC7K325T-2 | `GF(2^113)` | 3,071 | 37 (derived) | 145 K | 20.7 K |

Multipliers are 79% of the area; everything else — batching, the shared
inverter, the iteration logic, IO — is the other 21%. Their measured multiplier
growth from m=113 to m=127 is a factor 1.179 for a factor 1.124 in field size,
an exponent of 1.41, which extrapolates to **3,781 LUTs at m=131** (derived).

Per m=131 core, at 5.33 multipliers per core (16 for 3 cores):

* with literature-grade packing: 5.33 x 3,781 / 0.79 = **25.5 KLUT** (derived)
* with this repository's measured, unpacked circuit: 5.33 x 5,563 / 0.79 =
  **37.5 KLUT** (derived from measurement)

ECC2K-130 has two structural advantages over those 113- and 127-bit targets,
neither of which is credited above: **squaring is free** (Frobenius is a
coordinate permutation in the ONB — wires, not logic, where a polynomial basis
needs a reduction), and the walk needs no fruitless-cycle escape machinery,
since `j` comes from an orbit-invariant Hamming weight. Against that, m=131 is
the largest field of the three and the conversions in and out of the polynomial
basis are real logic. Treating these as cancelling is an **estimate**.

## 4. What fits on the AWS F2 part (derived)

The reference part is the AMD Virtex UltraScale+ HBM **VU47P** on
`f2.6xlarge`, the same part the other two cost models target.

> **Correction.** Both sibling cost models describe this part as having "2.85M
> LUTs". That is its *System Logic Cells* count (2,851,800). The VU47P has
> **1,303,680 CLB LUTs**, 2,607,360 flip-flops, 9,024 DSP slices, 70.9 Mb of
> block RAM, 270 Mb of UltraRAM and 16 GB of HBM. For UltraScale+ the
> logic-cell number is about 2.19x the LUT count, so capacity derived from
> 2.85M is more than twice too optimistic. Every figure below uses 1,303,680.

Budget:

* AWS F2 shell overhead: ~18% of the fabric (**estimated**, from the F1 shell's
  published footprint) leaves ~1.07 M LUTs.
* Practical utilisation 70–80%. The published designs ran at 63–96% of LUTs but
  hit stability and heat walls at the top: a 4-core Spartan-6 design "did not
  run stably, producing incorrect results", and Wenger and Wolfger had to drop
  an ML605 design from 275 MHz to 165 MHz for heat. Call it 748–855 K LUTs for
  cores.

| Per-core area | Cores on a VU47P |
|---|---:|
| 25.5 KLUT (literature packing) | 29–34 |
| 37.5 KLUT (this repository's unpacked circuit) | 20–23 |

**Clock: 250–350 MHz (estimated).** The published cores run at 180 MHz on a
28 nm Kintex-7; 16 nm UltraScale+ is roughly 1.4–2x faster for the same
pipelined logic depth, and three LUT levels per stage is a conventional target.
Power at this utilisation is the risk, not logic delay — every stage of this
design toggles on every clock, which is exactly the condition that defeated the
published 4-core design.

**Throughput: 5–12 B iterations/s per VU47P, central estimate ~8 B/s**, at one
iteration per core-clock (20 cores x 250 MHz to 34 cores x 350 MHz).

Three secondary resources, checked so that the LUT bound is known to be the
binding one:

| Resource | Demand at 27 cores | Available | Binding? |
|---|---:|---:|---|
| Flip-flops, 3 LUT levels/stage | ~1.18 M | 2.61 M | no |
| Flip-flops, 2 LUT levels/stage | ~1.78 M | 2.61 M | nearly — use 3 |
| Walk state + batch buffers | ~8 Mb | 70.9 Mb BRAM + 16 GB HBM | no |
| Distinguished-point output | ~7 KB/s | PCIe | no |

Walk state is 262 bits per walk and a report is one 32-byte record per 2^25.27
iterations, so neither memory nor IO comes close — the same conclusion the
CUDA client reaches for different reasons. **LUTs bind, then flip-flops.** That
is the SHA-1 signature, not the secp256k1 one.

## 5. Against the GPU

| | RTX PRO 6000 | VU47P (f2.6xlarge) |
|---|---:|---:|
| Iterations/s | **6.905 B measured** (6.767 B collecting) | 5–12 B estimated, ~8 central |
| Ceiling with any known arithmetic | 17.6 B absolute, 11–13 realistic | not established |
| On-demand | $4.14/GPU-h (g7e.48xlarge) | $1.980/h |
| Spot | $1.31/GPU-h | $0.709/h |

Expected work is 2^60.9 = 2.15e18 iterations, so at the central estimates:

| | GPU-hours or FPGA-hours | Spot | On-demand |
|---|---:|---:|---:|
| RTX PRO 6000 at 6.77 B/s | 88,300 | **$116 k** | $366 k |
| VU47P at 8 B/s | 74,700 | **$53 k** | $148 k |
| VU47P at 5 B/s (pessimistic) | 119,500 | $85 k | $237 k |
| VU47P at 12 B/s (optimistic) | 49,800 | $35 k | $99 k |

Cost is a property of the work, not the fleet size; the fleet buys wall-clock,
exactly as [aws/README.md](aws/README.md) sets out. The FPGA is **about 2x
cheaper per solved instance on spot** and about 2.5x on on-demand, and roughly
half of that advantage is F2 spot pricing rather than silicon.

**Energy (estimated).** This is where the binary-field FPGA case is strongest
and where this document is weakest, since nothing here was measured. The
published data points are per-FPGA measurements: a Spartan-3 at 111 M
iterations/s in 5 W (22 M/W) and a Spartan-6 at 300 M iterations/s in 8 W
(37.5 M/W). One RTX PRO 6000 at 6.9 B/s on a several-hundred-watt board is
roughly 12 M/W. A VU47P at 8 B/s in a plausible 75–150 W is 53–107 M/W, so
**5–10x the GPU per joule** — consistent with the historical pattern and
unverified here.

## 6. The datapath, generated (measured)

[codegen/genverilog.py](codegen/genverilog.py) is a third back end for the
straight-line IR, after C (`Prog.emit`) and CNF (`Prog.emitCnf`): it emits the
same DAG as pipelined Verilog. `make verilog` writes `generated/rtl/`, needing
python3 and no FPGA toolchain.

| Module | Gates | Stages | Pipeline registers | Latency | II |
|---|---:|---:|---:|---:|---:|
| `ecc_mul131.v` | 16,494 | 12 | 10,765 | 12 clk | 1 |
| `ecc_hamming131.v` | 268 | 3 | 40 | 3 clk | 1 |
| `ecc_sigma131.v` | **0** | — | — | combinational | — |
| `ecc_pre131.v` | 2,887 | 3 | 1,612 | 3 clk | 1 |
| `ecc_post131.v` | 33,641 | 24 | 29,765 | 24 clk | 1 |
| `ecc_inv131.v` | 8 × `ecc_mul131` | — | — | 96 clk | 1 |

`ecc_pre131` and `ecc_post131` are the walk's two passes, split exactly as
`include/walk.h` splits them. The asymmetry is the whole point: **pass 1 costs
no multiplication at all** — weight tree, three conditional squarings, and the
XORs that form `d = x + sigma^j(x)` and `e = y + sigma^j(y)` — while pass 2 is
two multiplications, `lam = e·d^-1` and `lam·(x + x3)`. 2,887 gates against
33,641.

`ecc_inv131` is Itoh-Tsujii, derived from the same addition chain
`FieldBs::inv` uses: 8 multiplications for m − 1 = 130, which is optimal, with
every intermediate squaring a permutation. It is emitted structurally, as 8
instances of the verified multiplier plus the wiring, so it inherits the
multiplier's verification. One inverter is ~44,500 LUT6, about 5% of a VU47P,
and Montgomery's trick means a device needs very few — which is why unrolling
it is affordable here and was not on the published Spartan-6 designs, where 8
multipliers would have been most of the chip.

The multiplier is the leaf-9 conversion multiplier, so its 16,494 gates are the
same 16,494 that §2 maps to 5,563 LUT6: the RTL and the area analysis describe
one circuit rather than two. Operands carry a tag through the pipeline and the
datapath never inspects it, which is what lets one multiplier serve many
independent walks at II = 1 with no stall logic — the same arrangement
[`hdl/ecc/ec_add_pipe.vhd`](../hdl/ecc/ec_add_pipe.vhd) uses for secp256k1.

`ecc_sigma131.v` contains **no gates**. `sigma^j` is a permutation of ONB
coordinates, so the module is wire assignments and three 2:1 mux layers
selecting `j = 3 + ((hw/2) mod 8)`. That is §2's "squaring is free" made
literal, and it is a genuine ECC2K-130 advantage: on a polynomial basis the
same operation is a reduction, which is why the published 113- and 127-bit
designs carry squarer logic and this one does not.

**Verification.** `make check-verilog` reads the emitted `.v` files back from
disk with an evaluator that knows nothing about the generator, simulates them
as synchronous circuits and compares against `field.Onb`:

```
ecc_mul131.v: PASS -- 8 multiplications, II=1, latency 12 clk, checked against field.Onb
ecc_hamming131.v: PASS -- 8 weights, II=1, latency 3 clk
ecc_sigma131.v: PASS -- sigma^(3+sel) for all 8 selects, 4 elements, 0 gates
ecc_pre131.v: PASS -- 8 steps, II=1, latency 3 clk, weight and sigma^j addends match field.Onb
ecc_post131.v: PASS -- 8 additions, II=1, latency 24 clk, points match curves.Curve
walk step: PASS -- 4 iterations of sigma^j(R)+R through ecc_pre and ecc_post, every
                   result on the curve and equal to curves.Curve
```

A real simulator now checks the same modules independently: `make -C
hdl/ecc2k130` runs the generated testbenches and the sequencer's under Icarus
Verilog.

```
PASS: 16 multiplications, II=1, latency 12 clk
PASS: 16 inversions, II=1, latency 96 clk
PASS: 124 walks x 3 steps retired, II=1, ring latency 124 clk, 517 distinguished
      points, all trajectories match curves.Curve
```

Both paths agreeing matters, because they fail differently. The Python checker
found nothing a simulator would not, but it runs with no toolchain at all. The
simulator immediately found what the Python checker structurally could not: the
valid pipeline had no initialiser, so it powered up as X and claimed valid
before the pipeline had filled. Data was never wrong — every product matched —
but a downstream consumer would have latched garbage on the first 12 clocks.

It recovers the latency and the initiation interval from the simulated
waveform rather than trusting the module header, so the two properties an
interleaved datapath lives on are checked rather than asserted. This caught a
real bug that no shape or count check would have: the Frobenius table was
indexed by destination rather than source, wiring `sigma^-j` — a permutation,
of the right size, with the right gate count, computing the wrong Frobenius
power. `tb_ecc_mul131.v` is the same check as a Verilog testbench for Icarus,
Verilator or Vivado; it has **not** been run, because no simulator exists here.

**One number moves.** The emitted pipeline cuts every 8 gate levels, about 3
LUT levels, which is §4's operating point — but it needs 10,765 registers per
multiplier against the 6,542 §2 estimated, because a depth-balanced cut by gate
level is wider than a cut by LUT level. At 5.33 multipliers per core and 27
cores that is 1.55 M flip-flops of 2.61 M, or 59%. The design stays LUT-bound,
so §4's conclusion holds, but the flip-flop margin is thinner than estimated
and a 2-levels-per-stage variant is now clearly out of reach.

## 7. The sequencer

[`hdl/ecc2k130/ecc_rho_core.v`](../hdl/ecc2k130/ecc_rho_core.v) is one complete
rho core. It is hand written, because it is control rather than arithmetic, and
it turned out to need **no state machine at all**: every generated block has
II = 1 and no back pressure, so the ring the published designs describe is
literally a ring.

```
load ->|                                                        |
       +--> ecc_pre --> ecc_inv --> [state RAM] --> ecc_post ----+
              3 clk       96 clk        1 clk         24 clk
```

A walk re-enters `ecc_pre` 124 clocks after it left, so 124 walks are resident
and each advances one step per 124 clocks — while the core retires **one step
per clock**, which is the number that matters. Two details are load-bearing:

* **The in-flight state is a RAM, not a delay line.** `ecc_post` needs x, y, d
  and e at the moment `d^-1` arrives, 96 clocks later. Registers would cost
  4 × 131 × 96 = 50,304 flip-flops per core, and §4 shows flip-flops are the
  second-scarcest resource; a RAM indexed by the walk's tag costs 2 block RAMs.
* **The tag is the whole mechanism.** Each block carries an opaque tag beside
  the operands and never inspects it, so results route back to the right walk
  with no matching logic — the arrangement
  [`hdl/ecc/ec_add_pipe.vhd`](../hdl/ecc/ec_add_pipe.vhd) already uses.

Distinguished points leave through a report port and the host supplies
replacement start points, as in the published designs: building a start point is
128 point additions and has no business inside the ring.

`tb_ecc_rho_core.v` loads one start point per slot from `curves.Curve`, lets the
ring free-run, and checks every retired step against that walk's reference
trajectory, that a step retires on every clock once full, and that nothing is
reported distinguished unless its weight really is within the cutoff. 372 steps,
all matching.

**One honest gap.** This core inverts every denominator directly, so it costs
2 + 8 = **10 multiplications per step**, not the 5.125 that §3's area numbers
assume. Montgomery's trick is what closes that: batching B denominators costs
3 + 8/B each, and §3's 5.33 multipliers per core depends on it. So the core is
correct but currently buys its correctness at about **half** the throughput per
LUT that §4 and §5 project. The batching engine is the next piece, and this
core is the reference it has to stay bit-exact against.

**What else is not built.** Host IO (the report FIFO's PCIe or HBM side), the
start-point generator, and multi-core instantiation with a shared inverter. And
**nothing has been synthesised**, so every area and clock number here is still
derived rather than measured.

## 8. What would settle it, in order

1. ~~**A Verilog backend for `ir.Prog`.**~~ Done — §6. The datapath is now
   generated from the same IR the CUDA client is built from and checked against
   the same field model, so the "derived" numbers in §3 describe a design that
   exists rather than one that might.
2. **Synthesis.** [`hdl/ecc/aws/`](../hdl/ecc/aws/) already has Vivado scripts
   that bring back utilisation and timing from an AWS build instance. Running
   them on one generated multiplier replaces the two biggest uncertainties —
   packing quality (§2) and clock (§4) — with measurements, and costs one
   build-instance hour rather than a quarter of engineering.
3. **Power, on a real card.** The published designs were limited by heat, not
   logic. Nothing short of a loaded board settles this, and it is what decides
   whether the core count in §4 is 20 or 34.

Until synthesis has been run, the honest summary is that an FPGA is
**competitive with, not decisively better than, one GPU** for this walk, and
that the case for building one rests on cost per instance and energy rather
than on speed. The datapath now exists and is verified; the sequencer around
it, and the synthesis run that would replace the two largest estimates here
with measurements, do not. The GPU client exists, works, and is audited; an
FPGA client is still a quarter of engineering for a number that may not beat
it.

## 9. Reproducing

```bash
cd ecc2k130/benchmarks/fpga-lut && python3 lutmap.py --json result.json
```

```bash
cd ecc2k130 && make verilog && make check-verilog
```

Both need `python3` only — no FPGA toolchain, no GPU, about 15 seconds each.
The mapper verifies every multiplier against `field.Onb` before counting it,
and `check-verilog` verifies the emitted RTL against the same model, so a
circuit that does not compute a product in `GF(2^131)` fails rather than being
reported.

Sources for the published figures: Bailey et al.,
[*Breaking ECC2K-130*](https://eprint.iacr.org/2009/541) (Spartan-3, 33.6 M
it/s per chip, 5 W); Fan et al., FPL 2010 (111 M it/s on the same part);
Wenger and Wolfger, [*Harder, better, faster, stronger*](https://eprint.iacr.org/2015/143)
(Kintex-7 XC7K325T, 203,800 LUTs, 5 cores at 180 MHz, 900 M it/s); Bernstein,
Engels, Lange, Niederhagen, Paar, Schwabe and Zimmermann,
[*Faster elliptic-curve discrete logarithms on FPGAs*](https://eprint.iacr.org/2016/382)
(multiplier and core areas, 7 cores at 180 MHz, 1.26 B it/s, 8 W per Spartan-6).
Device figures are from the AMD UltraScale+ product selection guide; prices were
observed on 2026-09-11.
