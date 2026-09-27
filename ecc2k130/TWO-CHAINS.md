# 30 B/s and the two-chain kernel

Question: can the table-walk kernel on one RTX PRO 6000 be taken from the
verified **20.078 B/s** of [ONE-BLOCK-GEOMETRY.md](ONE-BLOCK-GEOMETRY.md) to
**30 B complete scalar updates per second**?

Answer, derived before anything was built (§1) and unchanged by what was
measured (§5): **no, on this GPU.** 30 B/s is 15.2 SM-clocks per update, and
the five field products of one affine addition alone occupy the carry-less
unit for 18.5. The whole kernel's 33.1 `CLMAD`s per update put its floor at
**22.3 B/s** at 100% of the unit, and the 20.08 B/s build already sits at
0.90 of it. The remaining tenth is the only thing an iteration-function
change can still buy on this SKU, and it is what this note built for: a
kernel in which every warp carries two independent Montgomery chains, so the
inversion and the forward pass — the two phases the phase profile found
idling the unit — overlap inside the warp instead of waiting for other warps
that are in the same phase. Same arithmetic, same walk, same products per
update, the same 1,480,482 distinguished points byte for byte; static cost
per update unchanged to within 0.5% (§3). Measured on one RTX PRO 6000
against the reference rebuilt in the same session (§5): **13.90 B/s against
19.98, 0.696×**, every paired repetition between 0.692 and 0.699. Each warp
did issue `CLMAD`s 1.42× more often — the overlap the design was for is real
— but the kernel runs half as many warps, and 1.42 is not 2. **Engineering
that did not pay**; the reference stays. §6 then asks the question the
floor makes unavoidable — is the carry-less unit this slow on every GPU? —
and measures the answer: **no.** On an H100 and a B200 a `CLMAD` costs 2.0
and 3.8 logic slots instead of 38. On a B200 the same kernel runs at
**11.1 B/s** bound by the logic pipe with the unit 96% idle; the two
ALU→`CLMAD` trades that lost or measured neutral on the 6000 win there,
**+37.9% stacked (15.37 B/s**, verified, the part's best measured rate), at
par with the ALU slots they remove; and 30 B/s is not below a floor but a
1.95× cut of ALU work away — an engineering question this note prices and
does not build. Per dollar the 6000 stays 2.7× ahead.

## 1. Boundary

Unit: billions of complete scalar updates per second (`finished:` line), one
RTX PRO 6000 Blackwell Server Edition, 188 SMs, 2.42 GHz under load.

**Floor (derived):** the carry-less unit. ONE-BLOCK-GEOMETRY §1 measured
`CLMAD.LO`/`CLMAD.HI` at **1.62 lane-CLMADs per SM-clock** on this card,
independent of occupancy and operand width, and the PTX ISA (9.3, §9.7.1.5)
has only the `.u64` form, so there is no narrower instruction to trade to.
A 131×131 product is 6 `CLMAD` by two-way Karatsuba (three 64×64 products,
`lo` and `hi` each), and an affine addition with batched inversion is five
products per update (three for Montgomery's trick, one for λ, one for the
new y), plus 8 products and 20 ONB squarings per 16-slot inversion:

| what is on the unit | `CLMAD` / update | SM-clocks / update at 100% | B/s at 2.42 GHz |
|---|---:|---:|---:|
| five products of one affine addition, nothing else | 30.0 | 18.5 | 24.6 |
| the 20 B/s kernel: 4.81 products + inversion products + inversion squarings | 33.1 | 20.4 | **22.3** |
| the same with the inversion's squarings on the ALU (`PACKED_ALU_SQR=1`) | 31.85 | 19.7 | 23.1 |
| **30 B/s budget** | **24.6** | **15.2** | 30.0 |

The budget row is below the product-only row: 30 B/s needs fewer than 4.1
products per update on this unit, and no affine formula gives fewer than
five (ITERATION-FUNCTION.md §1 lists the alternatives and why each costs
more). One caveat on the rate itself, so the floor's uncertainty is on the
page: the tree has two measurements of the unit on this SKU. ONE-BLOCK-
GEOMETRY's 1.62 lane-CLMADs per SM-clock was taken on the card with the
kernel's own instruction mix and spacing and is the number the 22.3 rests
on. `benchmarks/clmad-price/probe.cu`, re-run here in the same session as
§6 ([probe-6000.txt](benchmarks/fast-clmad/probe/probe-6000.txt)), gives
**1.69** for a stream of `CLMAD.lo` alone and **1.99** for a stream of
`lo`+`hi` pairs — the same 1.687 / 2.00 ITERATION-FUNCTION §3.1 recorded on
the 4500. If the kernel's mix ran at the pair rate the floor would be 27.4
B/s and the 20 B/s build at 0.73 of it; at the on-card 1.62 it is 22.3 and
0.90. 30 B/s is above both: at the pair rate the five products alone are
15.0 SM-clocks against a budget of 15.2, which leaves the inversion,
squarings and every other cost 0.2 SM-clocks. The other execution resources are not the constraint — Nsight on the
17.4 B/s build put the ALU pipe at 36–39% and issue at 51% (THROUGHPUT-29B.md
§5), and a software product on those pipes costs about 800 slot-equivalents
against 6 `CLMAD`s at 38 slot-equivalents each (THROUGHPUT-30B.md), so
moving products off the unit costs more issue than it frees. Tensor cores
need one operand shared across eight or sixteen lanes and no product in the
walk has one. **30 B/s on one RTX PRO 6000 is therefore not an engineering
target; it is 1.35× the floor.** Two RTX PRO 6000s at the measured 20.08 are
40 B/s, which the fleet tooling already does.

**Reference:** `make gpu-rtx-pro6000-20b`, **20.078 B/s** median
(ONE-BLOCK-GEOMETRY §2), rebuilt in the same session as every row below.

**Where the 10% is** (ONE-BLOCK-GEOMETRY §6, `PHASE_PROFILE=1`, per warp per
16-slot step): forward pass 72,700 cycles against 58,400 of unit share (84%
fed), inversion 30,000 against 20,700 (69%), reverse pass 88,700 against
87,500 (99%). The reverse pass is unit-bound; the forward pass and the
inversion are not, and §7 of that note says why nothing was done about the
inversion: "Nothing in the warp can overlap it; only the other three warps
can, and they are in their own inversions at about the same time."

## 2. What was built

`ECC_PACKED_CHAINS=2` ([include/packedkernels.cuh](include/packedkernels.cuh)):
every thread runs **two independent Montgomery chains** of `ECC_BATCH/2`
slots, chain A on slots `[0, B/2)` and chain B on `[B/2, B)`, interleaved
slot by slot. In the forward pass chain A's twelve `CLMAD`s are issued, chain
A's next selection runs while they are in the unit, chain B's twelve are
issued, chain A's two reductions land, chain B's next selection runs, chain
B's reductions land: the same software pipeline `TABLE_PIPE_SELECT` builds
across slots, built across chains, with no deferred-W copy. In the reverse
pass the two chains' pairs (`inv·d` first, then `inv·W`) and then the two
y-products are issued together. The two inversions run **link by link**
through `inv131x2` ([include/packed131.h](include/packed131.h)): the same
Itoh–Tsujii chain as `inv131` (powers 2, 4, 8, 16, 32, 64, 65, 130) with
every product a `mul131x2` — one out-of-line body computing two independent
products, twelve `CLMAD`s with no mutual dependence — and every Frobenius map
a `sigma131x2`. A warp cannot issue a `CLMAD` more often than once per ~63
cycles and each link waited for the previous one; two chains give each link
two independent dependency chains and ptxas two conversion networks to
interleave.

The batch doubles so that the inversion share does not: 256 threads × 32
slots in two chains of 16 is the same 1,540,096 walks, the same 78.5 MB
persisting blob inside the 80 MiB L2 window, and the same 8 products per 16
slots as the reference's 512 × 16 — it trades 16 warps of one chain for 8
warps of two. `mul131x2`, `inv131x2` and `sigma131x2` are held bit-identical
to `mul131`, `inv131` and `sigma131` on the host (`make test-packed`,
including under the GPU preset's Frobenius networks), so the walk, its tags
and its distinguished points are the reference's; the job script below also
checks that directly.

## 3. The single table

Static columns from `kernel_cost.py` with the receipts' toolchain (nvcc
13.3.73 from the pinned pip wheels, `sm_120`; `--slot-updates 2` for the
two-chain loops, whose iterations advance one slot of each chain). Static
counts both arms of the slot-0 branches; the dynamic `CLMAD` count is the
static one less the predicated arms, 33.1 for every row with one inversion
per 16 slots. Measured columns are medians of five alternating repetitions
of `--bench --steps 1024 --launches 64` at the automatic worker count on
one RTX PRO 6000 Blackwell Server Edition (Modal, driver 580.95.05, nvcc
13.3.73, 2026-09-22; [summary.json](benchmarks/two-chains/summary.json),
raw logs beside it). The paired column is this variant over the reference
*within the same repetition*, so a drifting card cannot manufacture a
ratio; "verified" is 300 device reports re-walked by the host reference with
0 dropped, and "same DPs" is the sorted set of every distinguished point the
binary produced on the forced common 1,540,096 walks, compared by hash to
the reference's set of 1,480,482.

| variant | threads × blocks / SM | batch | regs | spills | instr / update (static) | ALU slots / update (static) | `CLMAD` / update (static; dynamic) | B/s, median of 5 | paired / ref (min – max) | / 22.3 floor | verified | same DPs | class |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|---|---|---|
| **reference** `gpu-rtx-pro6000-20b` | 512 × 1 | 16 | 116 | 0 | 2,114.6 | 1,921.4 | 40.25; 33.1 | **19.981** (ONE-BLOCK-GEOMETRY: 20.078) | 1 | 0.896 | 300/300 | yes | reference |
| reference + `PACKED_ALU_SQR=1` | 512 × 1 | 16 | 116 | 0 | 2,141.5 | 1,953.8 | 39.00; 31.85 | 19.818 | 0.992 (0.990 – 0.996) | 0.889 | 300/300 | yes | engineering, did not pay (as in ONE-BLOCK-GEOMETRY §2) |
| **two chains** `gpu-rtx-pro6000-chains2` | **256 × 1** | **32 = 2 × 16** | 182 | 0 | 2,123.3 | 1,910.7 | 40.25; 33.1 | **13.903** | **0.696 (0.692 – 0.699)** | 0.623 | 300/300 | yes | **engineering, did not pay** |
| two chains | 384 × 1 | 32 = 2 × 16 | 166 | 0 | 2,125.8 | 1,909.5 | 40.25; 33.1 | 15.627 | 0.784 (0.782 – 0.786) | 0.701 | 300/300 | yes | geometry scout, did not pay (118 MB blob, over the L2 window; clock 2362 – 2400 MHz) |
| two chains | 512 × 1 | 16 = 2 × 8 | 128 | 0 | 2,481.2 | 2,244.2 | 44.50; 37.35 | 15.690 | 0.785 (0.781 – 0.789) | 0.704 (0.80 of its own 19.7 floor) | 300/300 | yes | geometry scout, did not pay (one more inversion per 16 updates) |
| two chains + `PACKED_ALU_SQR=1` | 256 × 1 | 32 = 2 × 16 | 182 | 0 | | | 39.00; 31.85 | 13.422 | 0.672 (0.669 – 0.675) | 0.602 | 300/300 | yes | engineering, did not pay |
| **30 B/s** | | | | | | | ≤ 24.6 dynamic | 30 | 1.50 | **1.35** | | | **below the floor** |

Reading the static columns: the two-chain kernel is a rescheduling. Its
instructions per update are within 0.5% of the reference's, its `CLMAD`
count is identical, and its SASS is 1.7× longer because both passes carry
two slots per iteration. What it changes is not in this table; it is the
fraction of the 20.4 SM-clocks of unit time per update during which the unit
is actually busy, and the card measured that at 0.62 against the reference's
0.90 (§5).

Reading the measured columns: every row below the reference is slower than
it, every two-chain row by 20 – 30%, and every row is bit-for-bit the same
walk (six binaries, one distinguished-point hash `6cb064cd…`). The reference
itself re-measured 0.5% under its receipt, on a card that ran it at
2400 – 2422 MHz and 562 – 596 W; the two-chain kernel at 256 × 32 held
2422 MHz at 434 – 468 W throughout — a kernel drawing 130 W less than the
reference on the same card is a kernel leaving the datapath idle, which is
the measurement in one number.

Per `AGENTS.md` §3 every row is **engineering**: the floor is the cost of
the same generic algorithm, the ratio to it is bounded above by 1, and no
row can be an advance. A row that lowers `CLMAD`/update without changing the
walk (`ALU_SQR`) moves the floor itself by 3.8%, which is the most any row
here can move it, and measured −0.8%.

## 4. Falsification target, declared before the run

The two-chain build replaces `gpu-rtx-pro6000-20b` as the RTX PRO 6000 build
iff, in one session on one RTX PRO 6000 with both binaries rebuilt from the
same tree and run alternating:

- its median over five repetitions is **≥ 1.03×** the reference's median and
  **every** paired repetition is above 1.0;
- 300 of 300 device reports re-walk on the host reference with 0 dropped,
  for it and for the reference;
- its distinguished-point set at `dpWeight = 48`, forced to the reference's
  1,540,096 walks and `--run-id 7`, is **byte-identical** (as a sorted set of
  32-byte records) to the reference's — the same walks, the same points, not
  merely correct ones.

Below 1.03× it is engineering that did not pay and the reference stays. Any
row at or above **30 B/s** would falsify §1's floor and is not expected; any
row above 22.3 would falsify the measured `CLMAD` rate. Inadmissible: a
different worker count for the verification rows, a different `dpWeight`,
quoting the static table as a rate, or counting the 512 × 16 scout's extra
inversion as anything but the cost it is.

## 5. Measurement

The job is [benchmarks/two-chains/gpujob.sh](benchmarks/two-chains/gpujob.sh):
inside `nvidia/cuda:13.3.1-devel-ubuntu24.04` it builds the six binaries of
§3 plus `PHASE_PROFILE=1` forms of the reference and the two-chain kernel,
verifies each (300 re-walks; the distinguished-point sets written and
compared), benches them alternating for `REPS` rounds with SM clock, power
and temperature sampled after every run, and prints the two phase profiles.
Three launchers run that same script on one RTX PRO 6000 and bring
`results/` back with a `launch.json` receipt; the receipts here are from the
first:

```sh
cd ecc2k130
modal run modal_job.py --job benchmarks/two-chains/gpujob.sh --out /tmp/two-chains --env REPS=5
RUNPOD_API_KEY=… python3 runpod_job.py --job benchmarks/two-chains/gpujob.sh --out /tmp/two-chains --env REPS=5
python3 aws/bench_job.py --job benchmarks/two-chains/gpujob.sh --out /tmp/two-chains --env REPS=5
```

(The EC2 launcher could not be used in this session: `RunInstances` returns
`Blocked` in every region while dry runs pass, and AWS Health carries two
open risk events on the account, `AWS_RISK_CREDENTIALS_EXPOSURE_SUSPECTED`
and `AWS_RISK_ACCOUNT_CONSOLE_COMPROMISE`, with a support case open. It is
recorded because a launcher that was never run is a launcher that was never
tested end to end.)

### 5.1 Against the target

§4 asked for ≥ 1.03× with every paired repetition above 1.0. Measured:
**0.696×**, every paired repetition between 0.692 and 0.699
(13.895 / 13.903 / 13.902 / 13.905 / 13.912 against 20.073 / 19.994 /
19.920 / 19.981 / 19.904). The correctness rows are all green — 300/300
re-walked with 0 dropped for every binary, and one distinguished-point hash
across all six — so this is a clean negative: the same walk, 30% slower.
**The reference stays.** The two geometry scouts (384 × 32 at 0.784, 512 × 16
at 0.785) bracket it from above and say the loss is not one bad geometry.

### 5.2 Why two chains lost

The phase profiles ([profile-ref-prof.txt](benchmarks/two-chains/raw/profile-ref-prof.txt),
[profile-c2-256x32-prof.txt](benchmarks/two-chains/raw/profile-c2-256x32-prof.txt)),
cycles per warp per step, and from them the mean spacing between one warp's
consecutive `CLMAD`s in each phase:

| | reference, 16 warps/SM × 16 slots | two chains, 8 warps/SM × 32 slots | per-warp gain |
|---|---:|---:|---:|
| forward pass | 69,639 cycles (12 `CLMAD`/slot → 363 cycles each) | 110,121 (287) | 1.27× |
| inversion | 31,087 (68 `CLMAD` → 457) | 41,724 (136 → 307) | 1.49× |
| reverse pass | 83,560 (18/slot → 290) | 108,103 (188) | 1.54× |
| **per update** | **11,518 warp-cycles** | **8,123** | **1.42×** |
| unit fed (unit share at this occupancy / measured) | 0.91 | 0.64 | |
| throughput | 19.92 (this run) | 13.68 | 0.687 measured; 1.42 × 8/16 = 0.709 predicted from the profile |

The mechanism the kernel was built for is there: with two chains a warp
issues `CLMAD`s 1.27× to 1.54× more often in every phase, the inversion
included (§7 of ONE-BLOCK-GEOMETRY said nothing in a warp could overlap it;
`inv131x2` does, by half). But it does so with half the warps, and the
per-warp gain needed to break even was **2.0×**. At 2 warps per scheduler
the unit is fed 64% of the time against 91% at 4; the warp that would have
covered the other's stall is not there. The forward pass is where the gap
is widest — 1.27× — because two chains give it two `tableSelectSlot`s per
iteration, ~1,800 issue-cycles of lookups and ALU during which twelve
`CLMAD`s from one chain are in flight and the other chain's twelve are
waiting on the reductions of the first; the selection is on the ALU and MIO
pipes, and a second copy of it in the same warp does not overlap with the
first copy. It is the largest phase in the two-chain kernel (42%) where it
was the second-largest (38%) in the reference.

The register file is the reason 2 warps per scheduler is all there is: two
chains need 182 registers, 16 warps × 182 × 32 = 93 K of the SM's 64 K
registers. 384 × 32 (166 registers, 12 warps) and 512 × 16 (128 registers,
16 warps of two 8-slot chains, one more inversion per 16 updates) are the
two ways to get more warps back, and both land at 0.78: the first outgrows
the persisting-L2 window and gives back clock (2362 – 2400 MHz), the second
pays 12.8% more unit time for the extra inversion and, at 0.80 of its own
19.7 floor, is still fed worse than the one-chain reference at 0.90.

What this measures, in the terms of §1: on a unit that one warp can drive to
at most one `CLMAD` per ~63 cycles and that needs one every 79 cycles per
scheduler to be full, four warps per scheduler with 1.5 `CLMAD`s in flight
each beat two warps with 3 in flight each, because the stalls that idle the
unit are not the `CLMAD` latency but the ALU/MIO work between products, and
that work only overlaps across warps. The one-chain kernel at 512 × 1 is the
right shape for this SKU, and its remaining tenth is not reachable by
rescheduling within the warp.

### 5.3 Classification

| change | class | evidence |
|---|---|---|
| `PACKED_CHAINS=2` at 256 × 32 | engineering, did not pay | 0.696× paired, 5/5 repetitions; 300/300; identical DPs; unit fed 0.64 vs 0.91 |
| `PACKED_CHAINS=2` at 384 × 32, 512 × 16 | geometry scouts, did not pay | 0.784×, 0.785× |
| `PACKED_ALU_SQR=1` on either kernel | engineering, did not pay | 0.992×, 0.965× relative to each kernel's own row |
| the 22.3 B/s floor and the 30 B/s verdict | unchanged | every row ≤ 0.90 of the floor; the reference re-measured at 0.896 |

`gpu-rtx-pro6000-20b` remains the RTX PRO 6000 build. `PACKED_CHAINS=2`
stays in the tree behind its knob (off by default) as the measured answer to
"overlap the inversion inside the warp", so that it is not rebuilt.

## 6. The GPU where the question changes

Everything above is about the RTX PRO 6000, and on it the carry-less unit's
rate is the whole story. Three facts pointed at where that rate comes from:
Nsight Compute attributes `CLMAD` to the **FP64 pipe** (THROUGHPUT-29B.md
§5: "FP64/clmad 71%"); the RTX PRO 6000's die executes FP64 at 1/64 of FP32,
2 lanes per SM-clock, and the measured `CLMAD` rate is 0.8 of that; and
NVIDIA's CUDA 13.3 announcement reports GHASH on a B200 at a rate that needs
at least five times this card's unit. The hypothesis was that on a
full-rate-FP64 part the unit is 15–30× faster and this kernel is no longer
carry-less-bound. It was tested with one probe and one kernel run.

### 6.1 The probe: `CLMAD` on three parts

[benchmarks/clmad-price/probe.cu](benchmarks/clmad-price/probe.cu), unchanged,
through `modal run modal_job.py --job benchmarks/clmad-price/gpujob.sh --gpu …`
on three GPUs in one session (20,000 rounds × 3 passes × 3 runs, best per
stream; [probe/summary.json](benchmarks/fast-clmad/probe/summary.json) and
the raw output beside it). Lanes per SM-clock at the SM clock sampled after
each pass:

| GPU (SMs, clock) | `LOP3` | `CLMAD.lo` stream | `CLMAD` lo+hi product stream | **`LOP3` slots per `CLMAD.lo`** | `IMAD.WIDE` | `POPC` | `LDS.U8` random |
|---|---:|---:|---:|---:|---:|---:|---:|
| RTX PRO 6000 Blackwell SE (188, 2.35 GHz) | 63.9 | 1.69 | 1.99 | **37.8** | 31.3 | 16.1 | 9.2 |
| **B200** (148, 1.965 GHz, `sm_100`) | 63.7 | **16.9** | **29.1** | **3.78** | 22.9 | 16.0 | 9.2 |
| **H100 80GB HBM3** (132, 1.98 GHz, `sm_90`) | 59.8 | **29.9** | **30.1** | **2.00** | 22.6 | 16.0 | 9.2 |

The logic pipe, `POPC` and the shared-memory pipe are the same per SM-clock
on all three. The carry-less unit is not: **10× (lo stream) to 15× (product
stream) faster per SM-clock on the B200, 18× / 15× on the H100.** The
hypothesis holds — `CLMAD` is a full-rate instruction on the datacenter dies
and a 1/64-rate one on the workstation die — and the 6000's floor of §1 is a
property of that die, not of the instruction.

### 6.2 The kernel on a B200

[benchmarks/fast-clmad/gpujob.sh](benchmarks/fast-clmad/gpujob.sh): the
20 B/s knob set compiled for `sm_100` (`PRO6000_ARCH` override), plus the
ALU→`CLMAD` trades that the 6000 refused because its unit was full, each
verified (300 re-walks, 0 dropped) and held to one distinguished-point hash
across all six binaries (1,165,426 points on the forced 1,212,416 walks;
`626a6fce…`), then benched alternating for five repetitions. One B200
(driver 580.95.05, 1965 MHz throughout, 420 – 524 W of 1000; nvcc 13.3.73;
[summary.json](benchmarks/fast-clmad/summary.json)). Boundary for this part,
derived from §6.1 before reading the rates: 148 SMs × 1.965 GHz is 291 G
SM-clocks/s; the unit at 29.1 lane-`CLMAD`s per SM-clock does the kernel's
33.1 `CLMAD`s per update in **1.14 SM-clocks** — the carry-less floor here is
**255 B/s**, and irrelevant; the logic pipe at 63.7 lanes per SM-clock needs
**30.4 SM-clocks** for the reference's 1,939 static ALU slots. **30 B/s is
9.7 SM-clocks per update, i.e. ≤ 618 ALU slots at 100% of the logic pipe,
with room for 282 `CLMAD`s per update on the unit.**

| variant on the B200 | ALU slots / update (static, `sm_100`) | `CLMAD` / update (dynamic) | B/s, median of 5 | paired / ref (min – max) | SM-clocks / update | logic pipe at 100% would allow ≤ | unit busy | verified | same DPs | class |
|---|---:|---:|---:|---|---:|---:|---:|---|---|---|
| shipping σ-walk, 256 × 2 ([B200.md](B200.md), 2026-09) | 2,324 | 33.1 | 8.822 | — | 33.0 | | | 3/3 validate | different walk | before |
| **reference** = 20 B/s knob set, `sm_100` | 1,939 | 33.1 | **11.146** | 1 | 26.1 | 1,662 slots | 4.4% | 300/300 | yes | reference (+26% on the shipping walk) |
| **+ `PACKED_TOP_CLMAD=1`** (3-bit correction on the unit; −15% on the 6000) | 1,731 | 59.1 | **12.645** | **1.134 (1.134 – 1.135)** | 23.0 | 1,465 | 8.8% | 300/300 | yes | **engineering, +13.4%** |
| + `PACKED_ALU_SQUARE=0` (λ² on `CLMAD`), first run | — | — | ~~11.145~~ | ~~1.000~~ | | | | 300/300 | yes | **accounting: not a measurement.** The binary's identity line reads `packed alu square: 1`: the recursive `gpu-rtx-pro6000-20b` recipe hardcoded the knob, so the command-line 0 never reached the build and this row timed the reference against itself. Fixed in the Makefile (the recipe now forwards `PACKED_ALU_SQUARE`); re-measured below. |
| + `TOP_CLMAD` + `ALU_SQUARE=0`, first run | — | — | ~~12.644~~ | ~~1.134~~ | | | | 300/300 | yes | accounting, the same defect: identical to the `TOP_CLMAD` row |
| + `PACKED_ONB_INV=1` (−4.8% on the 6000) | 1,883 | 36.1 | 11.839 | 1.062 (1.062 – 1.062) | 24.6 | | 5.0% | 300/300 | yes | engineering, +6.2% |
| + `PACKED_CHAINS=2`, 256 × 32 | 1,911 | 33.1 | 9.204 | 0.826 (0.825 – 0.826) | 31.6 | | 3.6% | 300/300 | yes | did not pay here either |
| the `ALU_SQUARE=0` rows, re-measured, and the stacked best row | | | | | | | | | | **§6.4** |

Three things the table says. **The hypothesis was right about the unit and
wrong about what follows.** The unit is idle 91 – 96% of the time on the
B200, and the kernel is bound by the logic pipe instead — the reference's
26.1 SM-clocks per update against 30.4 of static ALU means the pipe is at
86 – 100% — so the B200 runs the same kernel at 0.56× the 6000's rate: its
SMs are 0.79 as many at 0.81 the clock, with the same logic pipe per
SM-clock, and per warp-cycle it is 12% *slower* than the 6000 (12,887 against
11,518 warp-cycles per update in the phase profiles) because the 6000's
unit, slow as it is, was still hiding some ALU latency that here lands on
the pipe. **The trades flip sign with the die**, which is the point of the
hypothesis: `TOP_CLMAD` cost 15% on the 6000 and buys 13.4% here; `ONB_INV`
cost 4.8% and buys 6.2%; the two-chain kernel loses on both, so its loss is
not a unit-rate effect but the register-file one §5.2 describes. (The two
`ALU_SQUARE=0` rows of this first run said "neutral" with a precision —
0.9998, and 1.1343 against 1.1342 — that should have been the tell: they
were the same binaries as their controls, a Makefile knob that did not
propagate, caught in review from the identity lines in the frozen logs.
§6.4 has the honest rows: **+8.8%** alone, **+25.5%** stacked on
`TOP_CLMAD`, **+37.9%** with `ONB_INV` on top.) **And 30 B/s on a B200 is
not below a floor**, unlike the 6000. It is 618 ALU slots per update against
the ≤ 1,205 the best row of §6.4 spends, i.e. a 1.95× cut of the logic-pipe
work, with 282 `CLMAD`s per
update of unit capacity to move it onto. The reduction (78 slots × 6.3 per
update, 490 in all) by `CLMAD` against the modulus is priced at −290 ALU /
+19 `CLMAD` in THROUGHPUT-20B §4; the basis conversions and the inversion's
Frobenius networks (≈ 300) are GF(2)-linear maps that a `CLMAD`-based
formulation might absorb; the selection's 175 slots are lookups that could
move to the idle unit as polynomial evaluations. None of that is built.
Whether the sum reaches 618 is the open engineering question on the B200 —
and it is an engineering question, not a floor.

What it is not is a campaign move. At Modal's list prices the 6000 does 6.6
B/s per dollar-hour at 19.98 B/s and the B200 2.5 at 15.37: **2.7× worse per
dollar**, and B200.md's break-even (29.1 B/s) is exactly the 30 B/s that is
not built. The receipts for the B200 rows also supersede B200.md's 8.82 B/s
as that part's best measured rate (+74% at 15.37, the same walk family);
B200.md keeps its number as the before mark.

### 6.4 The B200 rows re-measured, knob asserted

[benchmarks/fast-clmad/gpujob-clsq.sh](benchmarks/fast-clmad/gpujob-clsq.sh)
rebuilds the four affected binaries after the Makefile fix, refuses to time
any binary whose own identity line does not say the `PACKED_ALU_SQUARE` it
was asked for, verifies each (300/300, 0 dropped, one DP hash `626a6fce…`
across four) and benches them alternating, five repetitions, same B200 SKU
and driver as §6.2 ([clsq-rerun/summary.json](benchmarks/fast-clmad/clsq-rerun/summary.json)).
[gpujob-best.sh](benchmarks/fast-clmad/gpujob-best.sh) then stacks the ONB
inversion on the best of those ([best/summary.json](benchmarks/fast-clmad/best/summary.json)).

| variant on the B200 | ALU slots / update (static, `sm_100`) | `CLMAD` / update (static) | B/s, median of 5 | paired / ref (min – max) | SM-clocks / update | logic pipe allows ≤ | verified | same DPs | class |
|---|---:|---:|---:|---|---:|---:|---|---|---|
| reference (this session) | 1,939 | 40.25 | 11.058 | 1 | 26.3 | 1,675 | 300/300 | yes | reference |
| + `ALU_SQUARE=0` (λ² on `CLMAD`) | 1,831 (−5.6%) | 45.25 | **12.033** | **1.088 (1.088 – 1.089)** | 24.2 | 1,540 | 300/300 | yes | **engineering, +8.8%** (neutral on the 6000) |
| + `TOP_CLMAD=1` | 1,731 (−10.7%) | 66.25 | 12.595 | 1.139 (1.138 – 1.139) | 23.1 | 1,471 | 300/300 | yes | engineering, +13.9% (−15% on the 6000) |
| + `TOP_CLMAD=1` + `ALU_SQUARE=0` | 1,626 (−16.1%) | 71.25 | 13.876; 13.942 in the next session | 1.255 (1.253 – 1.256); 1.250 | 21.0 | 1,335 | 300/300 | yes | engineering, +25.5% |
| **+ `TOP_CLMAD=1` + `ALU_SQUARE=0` + `ONB_INV=1`** | **1,554 (−19.9%)** | 76.25 | **15.372** (ref that session 11.153) | **1.379 (1.378 – 1.379)** | **18.9** | 1,205 | 300/300 | yes | **engineering, +37.9%; the B200's best measured rate, 0.77 of the 6000's 19.98** |
| **30 B/s on a B200** | **≤ 618** | ≤ 282 | 30 | 2.69 | 9.7 | | | | **a 1.95× ALU cut away** |

The rate now moves with the static ALU column and with nothing else: −5.6%
slots buys +8.8%, −10.7% buys +13.9%, −16.1% buys +25.5%, −19.9% buys
+37.9%, while the `CLMAD` column rises 89% and the unit stays under 12%
busy. That is the signature of a kernel bound by one pipe with the other
empty, and it is the opposite of the 6000, where the same changes measured
−0.8%, −15% and −4.8% because that card's unit was the full one. The same
source, two dies, two different binding pipes, trades of opposite sign. The
gains run somewhat ahead of the static slot cut (−19.9% of slots for a
1.379× rate is what a −27.5% dynamic cut would give), which the profile
explains: `TOP_CLMAD` and `ONB_INV` both shorten the *dependent* part of a
product and of the inversion chain, and at four warps per scheduler the
kernel pays for latency as well as for issue.

### 6.5 What would have to be true for 30 B/s on one GPU

| part | binding pipe at the best row | best row | 30 B/s needs | status |
|---|---|---:|---|---|
| RTX PRO 6000 | carry-less unit (1.6 – 2.0 lane-`CLMAD`/SM-clock) | 19.98 | fewer than 4.1 products per update, or a faster unit | **below the floor**; two cards do 40 |
| B200 | logic pipe (63.7 lanes/SM-clock) | 15.37 here; **19.40** after the automatic sweep of [AUTOSWEEP.md](AUTOSWEEP.md) | ≤ 618 ALU slots per update, from ≤ 1,205 (≤ 955 after the sweep) | a 1.95× ALU cut onto an idle unit here, 1.55× after the sweep; every trade so far has paid at par or better, none of the remaining ones is built |
| H100 | logic pipe (59.8 lanes/SM-clock), by the probe; kernel not run | — | ≤ 521 ALU slots per update | as the B200, with 11% less pipe |

The per-GPU answer to the question at the top is therefore **no on every
part measured**, for two different reasons, and only one of them is a
floor. What the B200 rows add is the direction the remaining engineering
would take on a full-rate part: every logic-pipe instruction that a
`CLMAD` can replace is worth its slot count in rate, and the tree's
remaining priced candidates — the reduction against the modulus (−290 slots
for +19 `CLMAD`, THROUGHPUT-20B §4), the conversions and Frobenius networks
as `CLMAD`-based linear maps, the selection lookups as polynomial
evaluations — sum to something in the region of the 600 slots that 618
needs, and none has been built or measured. Per dollar the 6000 is still the
campaign's part: 6.6 B/s per dollar-hour against the B200's 2.5 at 15.37.
