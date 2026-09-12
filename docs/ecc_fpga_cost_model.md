# ECDLP on FPGA vs GPU: architecture and cost model

Companion to `hdl/ecc/` (the VHDL engine) and `gpu/ecc/` (the CUDA kernels).
Both implement the same attack — a distinguished-point Pollard rho over
secp256k1 — against the same reference implementation, so the comparison is
between two realisations of one algorithm rather than two different attacks.

The sibling document `docs/sha1_fpga_cost_model.md` does this for SHA-1
collisions, where FPGAs win by more than 10x. **The answer here is
different**, and the reason is instructive.

There is a third case, between the two: `ecc2k130/FPGA-CEILING.md` costs a rho
on a *binary* curve, `GF(2^131)`. Being carry-free it is LUT-bound like SHA-1
rather than DSP-bound like this document, but a field multiplication is still
thousands of gates and five are needed per step, so the FPGA lands at roughly
GPU speed, ~2x the cost-efficiency and ~5-10x the energy efficiency. The
binding resource, not the curve, is what decides which of the three answers
applies.

Numbers are labelled: **measured** (from simulation or static analysis in
this repository), **derived** (arithmetic on measured values and published
part parameters), **estimated** (engineering judgement, stated as such).
Nothing here has been run on an FPGA or a GPU.

## 1. The algorithm both platforms run

Distinguished-point parallel rho with the negation map. Per walk step:

| Work | Modular multiplies |
|---|---|
| Affine point addition given the inverse | 3 |
| Montgomery batch-inversion share | 3 |
| Amortised inversion, batch of W | 270/W |
| **Total at W = 64** | **10.2** |

Cost is denominated in 256-bit modular multiplies because on both platforms
everything else is free by comparison. The batching is what makes affine
coordinates viable: without it a step costs a full inversion, 276 multiplies.

Solving a 128-bit-security instance needs about sqrt(pi*n/4) ≈ 2^127.4 steps
with the negation map. Nothing in this document makes secp256k1 breakable;
the point is the cost *ratio* between platforms, which is what tells you
where to spend money for feasible-size problems (research instances, 60-90
bit challenge curves, or attacks on structurally weakened parameters).

## 2. FPGA realisation

`hdl/ecc/` implements the datapath. Measured in GHDL simulation against
vectors from `gpu/ecc/ecref.py`:

| Property | Value | Source |
|---|---|---|
| Modular multiplier latency | 12 clocks | measured, `fp_mul_tb` |
| Modular multiplier initiation interval | 1 clock | measured, `fp_mul_tb` (256 back-to-back multiplies) |
| Point addition throughput | 2.98 clocks | measured, `ec_add_tb` (64 additions issued in 191 clocks) |

The point adder reaches one addition per 3 clocks because it interleaves
independent walks through a single multiplier, cycling between the three
phases of the addition. That is the whole architectural idea: an elliptic
curve addition has a dependency chain three multiplies deep, so a lone
addition would leave a 12-stage pipeline 92% idle. Feed it 16 unrelated
walks and the pipeline runs full.

### Resource estimate per multiplier — derived

The multiplier is eight operand-scanning rows, each a 32 x 256 multiply
(eight 32x32 products) plus a 296-bit add, followed by three reduction
stages.

| Resource | Per multiplier | Basis |
|---|---|---|
| DSP48E2 | 256 | 64 products of 32x32, 4 DSPs each (27x18 primitive, no Karatsuba) |
| LUT | ~10k | estimated: 8 row adders of 296 bits, 3 reduction folds, control |
| FF | ~12k | estimated: 12 pipeline stages x (296 + 3x256) bits |

Karatsuba on the 32x32 products would bring 4 DSPs down to 3, saving 25% of
the DSPs at some LUT cost — worth doing, not done here.

### Per-part capacity — derived

Reference part is the AMD Virtex UltraScale+ VU47P on AWS `f2.6xlarge`, the
same part the SHA-1 model uses: 1,303,680 LUTs (2.85M *logic cells*), 9024
DSP, 84 Mb BRAM.  Section 7c records the correction.

- DSP-bound: 9024 / 256 = **35 multiplier pipelines**
- LUT check: 35 x 10k = 350k LUTs of 1.30M — not binding (see section 7c:
  the 2.85M figure quoted here for the part is its *logic-cell* count,
  not its LUT count)
- Context memory: 16 slots x 6 x 256 bits per engine, ~1.5 Mb total — not binding

DSPs bind, comfortably. That is the signature of a multiply-dominated
workload and the opposite of SHA-1, where LUTs bind.

### Clock rate — estimated

The critical path as written is a 296-bit adder in each row stage: roughly
37 CARRY8 levels, so **300 MHz is the honest estimate** for this RTL, not
the 500-600 MHz the DSP48E2 itself supports. Splitting each row adder into
two pipeline stages (16 row stages instead of 8, latency 20 instead of 12)
should reach 500 MHz and costs nothing in throughput, since II stays 1 and
the engine already tolerates latency by interleaving more walks. That is the
single highest-value change to the HDL and it is not implemented.

**Derived throughput**, 35 engines at 300 MHz, 10.2 multiplies per step:

    35 x 300e6 / 10.2 = 1.03 G steps/s per VU47P

At the 500 MHz the deeper pipeline should reach: 1.7 G steps/s.

## 3. GPU realisation

`gpu/ecc/` implements the same walk. The cost model per SM is:

    steps/s = SMs x clock x IPC_int / (instructions_per_mul x 10.2)

**Measured** (static, `gpu/ecc/ptx_stats.sh`): one modular multiply is 210
PTX instructions with hand-written carry chains, 472 without. SASS will be
somewhat fewer after ptxas scheduling; that has not been measured, because
disassembling requires `nvdisasm`, which was not available.

**Measured** (`ptxas -arch=sm_100`): the walk kernel uses 224 registers per
thread by default, giving 12.5% occupancy at 128-thread blocks, or 128
registers and 25% occupancy when constrained. Occupancy this low means the
integer pipes will not be saturated, so an instruction-issue bound is
optimistic.

An SM retires at most 4 instructions per clock. Taking 160 SASS
instructions per multiply as an **estimate**, one multiply costs at least
40 SM-clocks, so:

    steps/s ≈ SMs x clock / (40 x 10.2) = SMs x clock / 408

Fill in the SM count and clock of the specific part — `./bench` prints both
— rather than trusting a number here. For a 132-SM part at 1.7 GHz this
gives roughly 0.55 G steps/s, in the same order of magnitude as one VU47P.
**This is an estimate built on an estimate and should be replaced by
`./bench rho` output as soon as a GPU is available.**

## 4. Where this lands

| | SHA-1 collision search | ECDLP rho |
|---|---|---|
| Primitive | 80 rounds of 32-bit add/rotate/bool | 256-bit modular multiply |
| FPGA binding resource | LUTs | DSPs |
| Work per unit of silicon | very high (bit-level ops map perfectly to LUTs) | moderate (multipliers are what DSPs are for, but GPUs have many too) |
| FPGA advantage | >10x, robust | roughly 2x, fragile — estimated |

The SHA-1 result comes from bit-level operations that a LUT fabric does
essentially for free and a GPU emulates with general-purpose instructions.
Modular multiplication is not like that: it is exactly what a GPU's integer
multiply-add pipes are built for, and a modern GPU has thousands of them.
The FPGA still wins on multiplies per joule and per dollar-hour, but by a
factor that engineering effort on either side can erase.

**Practical reading**: for ECDLP work, write the GPU kernel first. It is a
week of work rather than a quarter, it is portable across rented hardware,
and the analysis above says the FPGA will not repay the difference unless
the campaign is long-running and power-limited. That is the reverse of the
recommendation in the SHA-1 document, and the reversal is the finding.

## 5. What would change the answer

- **Karatsuba plus a deeper adder pipeline** on the FPGA: 25% fewer DSPs and
  a 1.7x clock, together ~2.2x. This is the realistic FPGA upside and it is
  ordinary engineering, not research.
- **Multi-target rho**: a sqrt(t) speedup for t simultaneous targets, on both
  platforms equally. Does not move the ratio, moves the absolute cost.
- **Larger jump tables**: on the GPU these compete for shared memory and, on
  datacenter Blackwell, 228 KB per SM makes r = 2^11 practical (see
  `gpu/ecc/OPTIMIZATION_BLACKWELL.md`). On the FPGA the table lives in BRAM
  and is nearly free — a genuine, if small, structural FPGA advantage,
  because it lowers the fruitless-cycle rate for nothing.
- **A curve that is not secp256k1**: the special reduction is worth roughly
  40% here (210 PTX instructions against Montgomery CIOS's 352). On a
  generic prime both platforms fall back to Montgomery
  multiplication, and the FPGA loses slightly more, since its reduction
  stages are hard-wired to 2^256 = 2^32 + 977.

## 7. Index calculus: the two new kernels

Sections 1-5 cost *rho*. `gpu/macaulay/` and `gpu/semaev/` are the two
kernels an index-calculus run would want instead, and they fall on
opposite sides of the same LUT-versus-DSP line that decides everything
above. Neither is worth RTL, for different reasons.

### 7a. The Macaulay reduction over `F_p` — the FPGA loses outright

`gaudry_cubic` reduces one `226 x 286` matrix over `F_p` per residual
(`RESEARCH_RESIDUAL_WALKS.md` 11.6). The primitive is a **32-bit
Montgomery multiply**, `p < 2^31`.

| | FPGA (VU47P) | GPU (132 SM, 1.755 GHz) |
|---|---|---|
| Primitive | 32x32 multiply | 32x32 multiply |
| Per unit | ~4 DSP48E2 for the two multiplies of a REDC — derived | ~5 integer instructions — estimated |
| Units available | 9024 / 4 = **2,256** — derived | 132 x 128 = **16,896 lanes** — part spec |
| Clock | 500 MHz (short carry chains, unlike the 296-bit adders of section 2) — estimated | 1.755 GHz |
| Throughput | 1.1e12 mulmod/s — derived | 5.9e12 mulmod/s — derived |

**The GPU wins by about 5x**, and the reason is the same one section 4
gives: a small modular multiply is *exactly* what a GPU's integer
multiply-add pipes are built for, and there are 16,896 of them against
9,024 DSPs. The 256-bit case of section 2 was close only because a
256-bit multiply costs 256 DSPs, which throws away the FPGA's count
advantage; at 32 bits there is no such tax and the raw unit counts
decide.

Nothing in `hdl/ecc/` is reusable here either: `fp_mul_secp256k1.vhd` is
a 256-bit engine with reduction stages hard-wired to `2^256 - 2^32 -
977`. An `F_p` Macaulay engine would be new RTL for a fight it loses.

### 7b. The binary decomposition sweep — the FPGA is competitive, and that is all

`semaev_decomp` / `gpu/semaev` sweep pairs over `GF(2^n)` with `n = 3l`,
so `n <= 54` at the `l = 18` this would be built to reach. Carry-free,
so **LUTs bind** — the SHA-1 regime, not the ECDLP one.

`ecc2k130/FPGA-CEILING.md` already costs this regime at `m = 131` and
its measured scaling transfers directly. Bernstein et al.'s post-place
multiplier areas give a growth exponent of **1.41** in field size
(3,071 LUTs at `m = 113`, 3,620 at `m = 127`), so

    3,071 x (54/113)^1.41 = 1,085 LUTs per GF(2^54) multiplier   (derived)

VU47P carries **1,303,680 LUTs**. At 60% utilisation and the 79%
multiplier share `FPGA-CEILING.md` measures:

    0.60 x 1,303,680 x 0.79 / 1,085 = 569 multipliers             (derived)

At the 250-350 MHz that document estimates for this fabric, and the
**190 field multiplications per pair** measured in
`RESEARCH_SEMAEV_DECOMPOSITION.md`:

    569 x 300e6 / 190 = 0.90 G pairs/s per VU47P                  (derived)

The GPU side is where I had this badly wrong, and the correction is the
finding:

> **NVIDIA has a native carry-less multiply.** PTX 9.3 introduced
> `clmad`, documented for `sm_80` and later, and
> `ecc2k130/NATIVE-CARRYLESS.md` *measures* a 22.4% end-to-end gain from
> it on the ECC2K-130 client (7.110 -> 8.704 B walk updates/s, RTX PRO
> 6000, CUDA 13.3).

The first version of `gpu/semaev` asserted the opposite, built `gf_mul`
as an `n`-iteration shift-reduce loop, and rested its whole cost
argument on the gap. With `clmad`, one `GF(2^54)` multiply is two
carry-less multiply-adds plus a reduction — call it **12 instructions,
estimated** — instead of the ~324 the loop needs:

    16,896 lanes x 1.755 GHz / (12 x 190) = 13 G pairs/s          (derived)

So the FPGA is at roughly **0.07x** the GPU on this kernel, not the
order of magnitude the no-`clmad` premise implied. Even discounting the
GPU heavily for occupancy and for the reduction being more than 12
instructions, it does not reach parity.

That lands in the same place `FPGA-CEILING.md` lands for ECC2K-130 —
**competitive with, not decisively better than, one GPU** — and for the
same reason, with the binary-field primitive now native on both sides.

### 7c. An accounting correction to section 2

Section 2 states the VU47P as "2.85M LUTs". That is its **logic-cell**
count; the part has **1,303,680 LUTs** and 2,607,360 flip-flops, which
is the figure `ecc2k130/FPGA-CEILING.md` uses. The LUT budget in
section 2 is therefore 2.19x smaller than written.

The conclusion there is unaffected — 35 engines x 10k LUTs = 350k is
27% of 1.30M rather than 12% of 2.85M, so DSPs still bind, comfortably —
but the margin is narrower than it reads, and a Karatsuba variant that
trades DSPs for LUTs has less room than section 5 implies.

### 7d. Recommendation

**Do not write RTL for either kernel.** For the Macaulay reduction the
FPGA loses on unit count and there is nothing to reuse. For the binary
sweep the FPGA is competitive rather than better, and the one argument
that made it look better than that was mine and was wrong.

The ordering that follows is the same one section 4 reaches for rho, now
for both new kernels: write the GPU kernel, use `clmad`, and spend the
quarter an FPGA campaign would cost on the thing neither platform can
buy. Under `AGENTS.md` both kernels are **engineering** rows — they move
wall-clock, and `RESEARCH_RESIDUAL_WALKS.md` 11.7 puts the crossover with
rho past `2^230` regardless of how fast either phase runs.

## 6. Reproducing

```bash
cd hdl/ecc && make          # GHDL: multiplier and point-adder testbenches
cd hdl/ecc/aws && ./run_aws_synthesis.sh --dry-run ...   # Vivado on AWS: replaces
                            # the DSP and clock estimates above with measurements
cd gpu/ecc  && make test    # CPU verification of the CUDA headers
cd gpu/ecc  && ./ptx_stats.sh   # static instruction / occupancy figures
cd gpu/ecc  && make bench && ./bench selftest && ./bench rho   # needs a GPU
```
