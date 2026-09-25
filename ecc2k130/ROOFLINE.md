# Percent of peak: where the 20 B/s table walk sits on the RTX PRO 6000

Question: what fraction of the RTX PRO 6000's throughput does the Pollard rho
kernel (`make gpu-rtx-pro6000-20b`, the table walk `R ← R + ε·σᵏ(T_h)`,
**20.078 B/s**, [ONE-BLOCK-GEOMETRY.md](ONE-BLOCK-GEOMETRY.md)) actually use,
which unit is the wall, and how close can engineering take it to 30 B/s?

Answer:

- **The busiest unit is the ALU (logic) pipe, at 91% of its measured peak.**
  Each scalar update issues 1,320 ALU lane-instructions (LOP3, SHF, PRMT,
  IADD). At 64 lanes per SM-clock that is 20.6 of the 22.6 SM-clocks an
  update takes. The carry-less unit is 73% busy (33.1 `CLMAD`s at 2.0 lanes
  per SM-clock), issue slots 59%, the FMA pipe 13%, LSU 22% (§2).
- **Speed of light for this algorithm on this card is 27.4 B/s**: the
  carry-less unit 100% busy at 33.125 `CLMAD`s per update. The 20 B/s build
  is at **0.73** of it, and at **0.67** of 30 B/s.
- **The 22.3 B/s "carry-less floor" this tree has quoted since
  ONE-BLOCK-GEOMETRY is not a floor.** It assumes 1.62 lane-`CLMAD`s per
  SM-clock. A measured build ran faster than that rate allows: the
  three-limb-Karatsuba sweep arm issues 54.375 `CLMAD`s per update and took
  32.7 SM-clocks, so the unit sustained **at least 1.654** inside a real
  kernel (§3). This is an accounting correction; it moves no measured rate.
- **30 B/s is not an engineering target on one card** (§5). It needs at most
  30.3 `CLMAD`s per update with the carry-less unit never idle. That is
  exactly five products per update with a free inversion. It also needs the
  integer work cut by 27–33%.
- **Built here: two integer cuts, both bit-exact.** One is a shared-memory
  table for the polynomial square of λ (`PACKED_SQUARE_TABLE=1`). The other is
  a polynomial-basis Itoh–Tsujii inversion (`PACKED_INV_POLY=1`, or `=2` to
  keep one out-of-line copy of its product). Together they remove 86–87 ALU
  lane-instructions per update. **Predicted** +7.0%, to about 21.5 B/s; **not
  yet measured**: no GPU was reachable from the session that built them. [benchmarks/roofline/gpujob.sh](benchmarks/roofline/gpujob.sh)
  measures them in one run, together with the two pipe rates the model
  leaves open (§4).

Every number below is either a receipt already in the tree, cited, or a
count from the tools this note adds. Predictions are labelled as such.
Scope: these are `--bench` rates of the table walk. As built, that walk is
not the campaign default: its cycle rule lets fruitless cycles through
([WALK-CONSTANT.md](WALK-CONSTANT.md)), and nothing here changes that.

## 1. The tool: dynamic instruction counts per pipe

[`roofline.py`](roofline.py) builds the walk kernel with exactly the `-D`
flags the Makefile gives the client (read from `make -n`). It then
disassembles it with `nvdisasm -gi`, which gives every instruction's full
inlining chain, and builds the control-flow graph. Each conditional branch
gets a probability:

- a loop back-edge gets `1 − 1/trips`, with the trip count read from the
  source loop;
- a forward branch gets the probability of the source region it enters. The
  distinguished-point report and the overdue restart get about 0, the
  reverse pass's slot-0 arm gets 1/B, and so on (`REGION_RULES`).

The tool solves the visit equations exactly and prices each instruction on
the pipe that executes it, at a measured lane rate. Unlike
`kernel_cost.py`, it counts what a warp actually issues per update, not
both arms of every branch.

Two checks bind every count in this note:

- **CLMAD self-check.** The dynamic `CLMAD` count must equal what the
  arithmetic implies: 33.125 for the 20 B/s build, that is
  (5·16 − 3 + 8)/16 products × 6 plus 20/16 inversion squarings. All 46
  sweep builds in §3 pass. Before this note's fixes, two cases did not (§8).
- **SASS identity.** The standalone walk build is the client's kernel.
  Against `cuobjdump` of the full client, the 20 B/s build matches in all
  4,920 instructions and the build with both knobs in all 6,856. The only
  differences are in how the two disassemblers print branch targets.

```
./roofline.py --measured 20.078                              # the table in §2
./roofline.py --make-var PACKED_SQUARE_TABLE=1 --make-var PACKED_INV_POLY=1
./roofline.py --ops 30 --functions 30 --branches             # where the work is
./roofline_calibrate.py                                      # §3, every sweep receipt
```

Frozen outputs:

- [benchmarks/roofline/ref-20b.txt](benchmarks/roofline/ref-20b.txt),
  [both-20b.txt](benchmarks/roofline/both-20b.txt) and
  [both2-20b.txt](benchmarks/roofline/both2-20b.txt);
- [calibration.txt](benchmarks/roofline/calibration.txt) and
  [calibration.json](benchmarks/roofline/calibration.json).

The toolchain was nvcc, ptxas and nvdisasm 13.3.73 from NVIDIA's pip wheels,
`sm_120`.

## 2. The 20 B/s build against every pipe

Unit: work per complete scalar update, in lane-instructions. The clock is
188 SMs at 2.415 GHz, the median SM clock sampled after the timed runs of
this build in [two-chains](benchmarks/two-chains/summary.json). At that
clock, 20.078 B/s is **22.61 SM-clocks per update** per SM.

| unit | peak, lanes/SM-clk (receipt) | work/update | SM-clk/update | ceiling B/s | busy at 20.078 |
|---|---|---:|---:|---:|---:|
| ALU: LOP3, SHF, PRMT, IADD, … | 64 (LOP3 62.1–63.9, IADD3 63.0, SHF 62.5: hardware-limits, fast-clmad probes) | 1,319.6 | 20.62 | 22.0 | **91.2%** |
| carry-less (`CLMAD`, FP64 pipe) | 2.00 (1/64 of FP32; lo+hi probe stream 1.99) | 33.125 | 16.56 | **27.4** | 73.2% |
| issue (4 schedulers) | 128 | 1,698.9 | 13.27 | 34.2 | 58.7% |
| LSU (LDS, LDG, STG, …) | 16 (`ld.shared` 15.9) | 80.0 | 5.00 | 90.8 | 22.1% |
| FMA: IMAD; `.WIDE`/`.HI` count twice | 64 (IMAD 62.3) | 188.6 | 2.95 | 154 | 13.0% |
| XU: POPC, BREV, … | 16 (POPC 16.1) | 9.75 | 0.61 | 745 | 2.7% |

The kernel is within 9% of the ALU pipe's measured peak, and that pipe is the
tightest ceiling of the ones §3 leaves standing. The carry-less unit is the
largest ceiling that no integer cut can move: 27.4 B/s for this arithmetic
at 33.125 `CLMAD`s per update.

The integer work by source function: the innermost inlined frame, lane-instructions per
update, ALU and FMA together
([ref-20b.txt](benchmarks/roofline/ref-20b.txt)):

| function | 20 B/s build | + both knobs (`INV_POLY=1`) | what it is |
|---|---:|---:|---|
| `reducePolynomial131` | 440.6 | 402.1 | reduction of a 262-bit product mod the dense 19-term β-basis modulus, ~70 per reduction |
| `product131` | 372.1 | 367.9 | Karatsuba fold and the 3-bit top correction on the ALU |
| `spread32alu` | 109.0 | 52.0 | polynomial squaring of λ (bit spread) |
| `clmul64` | 94.8 | 98.4 | operand moves around the `CLMAD`s |
| `toPolynomial131` | 85.6 | 46.2 | ONB → β-basis transform |
| `fromPolynomialReduced131` + `fromPolynomialProduct131` | 103.6 | 76.4 | β-basis → ONB transforms |
| `sigmaInvNetwork131` + `sigmaWalkNetwork131` | 79.5 | 79.4 | Frobenius networks |
| table walk selection (`twPhase`, `twPivot`, `twAddend`, `twByte`, `twDenominator`, `tableSelectSlot`, `fusedSelect`) | 215.4 | 215.5 | next-point selection and addend fetch |
| `squarePolynomialTable131` | — | 125.0 | the table square, including its 65 `LDS` |

The reduction is the largest single item: 29% of the integer work.

## 3. Which peaks are real: calibration against 46 measured builds

A pipe peak is falsifiable. A build that ran faster than a model's floor
proves that model's peak too low. [`roofline_calibrate.py`](roofline_calibrate.py)
rebuilds all 46 builds of the two frozen automatic sweeps
([autosweep](benchmarks/autosweep/): one card per sweep, arms alternated
with their base). It counts every build and computes each one's efficiency
under six pipe models:

`η = model floor / measured SM-clocks per update`, which must be ≤ 1.

On the RTX PRO 6000 (21 arms; the mean and spread are over the 15 builds in
the 512×1, batch-16 geometry):

| model | max η | mean η | sd | verdict |
|---|---:|---:|---:|---|
| M1 `CLMAD` 1.62 only (ONE-BLOCK-GEOMETRY §1) | **1.027** | 0.894 | 0.041 | **refuted** by the three-limb Karatsuba arm |
| M2 `CLMAD` 2.00 only | 0.832 | 0.724 | 0.033 | consistent |
| M3 ALU 64 + `CLMAD` 1.62 | **1.027** | 0.908 | 0.037 | **refuted** (same arm) |
| **M4 ALU 64 + `CLMAD` 2.00** | 0.920 | **0.895** | **0.025** | consistent; flattest |
| M5 one integer datapath (ALU + FMA) 69 + `CLMAD` 2.00 | 0.977 | 0.944 | 0.036 | consistent |
| M6 ALU 128 (Nsight's denominator) + `CLMAD` 2.00 | 0.832 | 0.724 | 0.033 | consistent; `CLMAD` binds everywhere |

**The refutation.** The arm `PACKED_KARAT3=1 PACKED_TOP_CLMAD=0` issues 10
`CLMAD`s per product, so 54.375 per update. Its self-check passes. It
screened at **13.894 B/s**, 0.696× its base. Even at the part's maximum
2.43 GHz that is at most 32.88 SM-clocks per update. So the carry-less unit
executed at least **1.654 lane-`CLMAD`s per SM-clock** inside a real
kernel. The 1.62 rate would need 33.56. The margin is 2.1% at the most
favourable clock. The arm is a screen, not a verified run on this card.
Its arithmetic is bit-exact, and the same knob was verified on the B200
([AUTOSWEEP.md](AUTOSWEEP.md)). The probe streams agree: 1.69 for
`CLMAD.lo` alone, 1.99 for lo+hi pairs
([probe-6000.txt](benchmarks/fast-clmad/probe/probe-6000.txt)).

ONE-BLOCK-GEOMETRY does not record the probe behind 1.62. The walk issues
its `CLMAD`s as lo+hi pairs on the same operands, the pattern that ran at
1.99.

**Among the survivors,** M4 explains the knob-to-knob differences best. Its
η is flattest across same-geometry builds that move one pipe's work at a
time. The discriminating arm is `PACKED_ALU_SQR=1`: +1.9% ALU work and
−3.8% `CLMAD`s.

- M6 says that arm is `CLMAD`-bound and predicts +3.9%.
- M4 predicts −1.9%, and M5 predicts −2.1%.
- It measured −0.4%.

Nsight's "ALU 36–39%" on the 17.4 B/s build
([THROUGHPUT-29B.md](THROUGHPUT-29B.md) §5) reads as M6's denominator.
Every probe in the tree that isolates the pipe measures 62–64 lanes per
SM-clock, half that. §4's probe exists to settle which one is physical.

The B200 sweep (25 builds: base, 20 arms, 4 greedy combinations; `CLMAD`
at 29.1) is bound by integer work everywhere. Its best build, `gpu-b200-19b` at 19.40 B/s, is at η = 0.944
of one 69-lane integer datapath and 0.882 of a 63.7-lane ALU.

## 4. The open rates, and the job that settles them

Two pipe rates move the percentages above, and the tree has not pinned
either down.

1. **Are the ALU and FMA pipes separate on `sm_120` (M4), or one integer
   datapath (M5)?** Two receipts pull in different directions:
   - The hardware-limits probe measured LOP3 + IMAD interleaved at 69.3 lanes
     per SM-clock against 62–63 for each alone, which reads as shared.
   - The same probe put FFMA at 54.7, far under the part's 128 FP32 lanes, so
     it had a ceiling of its own.

   If the pipes are separate, shifts can move from the 91%-busy ALU pipe to
   the 13%-busy FMA pipe (IMAD.SHL) and the integer ceiling rises. If they
   are shared, the kernel is at 96% of it.
2. **What does the carry-less unit sustain in the walk's own pattern?** The
   pattern is independent lo/hi pairs, an RZ addend, and about 45 integer
   instructions per `CLMAD`. That rate is somewhere between the ≥ 1.654 of
   §3 and the 1.99 of the pair probe, and it sets the speed of light between
   22.7 and 27.4 B/s.

[`benchmarks/roofline/pipes.cu`](benchmarks/roofline/pipes.cu) times each
instruction class alone and in mixes: LOP3, SHF, PRMT, IADD3, IMAD, IMAD.HI,
IMAD.WIDE and FFMA; LOP3+IMAD, SHF+IMAD, LOP3+FFMA, and the walk's
8:4:2:1:1 integer proportions; and `CLMAD` as lo, hi, same-operand pairs,
dependent lo→hi, pairs beside 8-, 20- and 40-LOP3 chains, and a whole
128×128 Karatsuba. Each stream runs 16 independent chains per thread at
full occupancy.

[`sass_loops.py`](benchmarks/roofline/sass_loops.py) audits the SASS of
every timed loop. Three streams would have measured the wrong thing
without it:

- a Karatsuba middle product was hoisted out of the loop (32 `CLMAD`s per
  round instead of 96);
- FFMAs on thread-invariant data moved to the uniform datapath;
- a 32-bit `mad.hi` addend cost a MOV per multiply, because `sm_120`'s
  IMAD.HI takes a 64-bit addend.

All three are fixed. Every stream's loop is now the instructions its name
says.

[`benchmarks/roofline/gpujob.sh`](benchmarks/roofline/gpujob.sh), on one RTX
PRO 6000 under any of the tree's launchers, runs:

1. The probe, with the SASS audit of the binary that ran.
2. Seven builds, each priced by `roofline.py` with the container's own
   compiler: `ref`, `sqtab`, `invpoly`, `invpoly2`, `both` and `both2`, plus
   `fused2`. `fused2` is exploratory: both2 in the one-pass kernel
   (`TABLE_FUSED=1`).
3. For each build: 300 device reports re-walked by the host reference, and
   the distinguished-point set, which must be byte-identical across all
   seven builds at a forced common walk count.
4. `REPS` alternating rate rounds, with SM clock and power sampled.
5. Nsight Compute's per-pipe utilisation of `ref` and `both2`, where the host
   lets a container read the counters. ONE-BLOCK-GEOMETRY's host did not; a
   refusal is recorded, not fatal.

[`summarize.py`](benchmarks/roofline/summarize.py) freezes a run into
`summary.json`, with every arm's measured paired ratio beside its predicted
one.

```
modal run modal_job.py --job benchmarks/roofline/gpujob.sh --out /tmp/roofline-run1
# or, on EC2 or RunPod, which ship only the build tree plus --extra:
python3 aws/bench_job.py --job benchmarks/roofline/gpujob.sh --out /tmp/roofline-run1 \
    --extra benchmarks/roofline --extra roofline.py
python3 benchmarks/roofline/summarize.py /tmp/roofline-run1 benchmarks/roofline
```

## 5. What 30 B/s would take

At 2.415 GHz, 30 B/s is **15.13 SM-clocks per update**.

- **Carry-less unit.** At its architectural 2.0 lanes per SM-clock, 15.13
  SM-clocks is 30.3 `CLMAD`s.
  - One affine addition with Montgomery's trick is five products. At 6
    `CLMAD`s each, that is 30.0 before the inversion's share (8 products and 20
    squarings per batch).
  - At batch 16 the kernel issues 33.125; batch → ∞ tends to 30.
  - So 30 B/s needs the unit busy every cycle, and an inversion that costs
    nothing. If the walk's pattern runs no faster than the 1.654 of §3
    demonstrates, five products alone cap the card at 25.0 B/s.
  - No knob in the tree changes the six `CLMAD`s per product or the five
    products per update. [ITERATION-FUNCTION.md](ITERATION-FUNCTION.md) §1
    prices the alternative formulas, and each costs more.
- **Integer work.** At 64 lanes per SM-clock with the pipe never idle, 15.13
  SM-clocks is 968 ALU lane-instructions per update, against 1,320 now:
  −27%. At the 91% efficiency the kernel shows today, it is about 880:
  −33%.

Both conditions are required, and the first one cannot be met by
engineering. **30 B/s on one RTX PRO 6000 stays a two-card number.** Two
cards at 20.078 give 40 B/s, and the fleet tooling already runs them. The
reachable single-card range is set by the integer cut:

| integer work per update | ALU-bound rate (M4 at η 0.895–0.912) | carry-less ceiling |
|---:|---:|---:|
| 1,320 (20 B/s build) | 20.1 (measured) | 27.4 |
| 1,233 (+ both knobs) | 21.5 (**predicted**) | 27.4 |
| 1,060 | 24.5–25.0 (**extrapolated**; optimistic, since both pipes then bind together) | 27.4: from here the carry-less unit binds |

## 6. The single table

Unit: B complete scalar updates per second, one RTX PRO 6000.

- The **boundary** is the speed of light of §2: 27.4 B/s at 33.125
  `CLMAD`s per update. Its uncertainty is §4's rate: 22.7–27.4.
- The **reference** is `make gpu-rtx-pro6000-20b` at 20.078 B/s.

Every change in this note keeps the walk, the arithmetic's results, and the
operations per update. So none can be an advance, and each is engineering
or accounting.

| variant | B/s | / 27.4 SOL | / 20.078 | verified | class |
|---|---:|---:|---:|---|---|
| `gpu-rtx-pro6000-20b` (reference) | 20.078 measured | 0.733 | 1.000 | 300/300, ONE-BLOCK-GEOMETRY | reference |
| + `PACKED_SQUARE_TABLE=1` | 20.99 **predicted** (M4); 21.10 (M5) | 0.766 | 1.045 / 1.051 | 20,134 host cases bit-identical; GPU pending | engineering, unmeasured |
| + `PACKED_INV_POLY=1` | 20.53 **predicted** (M4); 20.40 (M5) | 0.749 | 1.022 / 1.016 | 2,133 host cases bit-identical; GPU pending | engineering, unmeasured |
| + `PACKED_INV_POLY=2` | 20.53 **predicted** (M4); 20.45 (M5) | 0.749 | 1.023 / 1.018 | as above | engineering, unmeasured |
| + both (table + `INV_POLY=1`) | 21.50 **predicted** (M4); 21.47 (M5) | 0.784 | 1.071 / 1.069 | as above | engineering, unmeasured |
| **+ both2 (table + `INV_POLY=2`)** | 21.48 **predicted** (M4); 21.51 (M5) | 0.784 | 1.070 / 1.071 | as above | engineering, unmeasured |
| fused2: both2 in the one-pass kernel (`TABLE_FUSED=1`), exploratory | 21.74 **predicted** (M4); 22.00 (M5) | 0.793 | 1.083 / 1.096 | as above | engineering, unmeasured; the fused kernel alone measured +0.2% against a +1.0% M4 prediction (AUTOSWEEP §4) |
| "22.3 B/s carry-less floor" (README, TWO-CHAINS §1) | — | — | — | refuted by a measured arm, §3 | **accounting** correction |
| 30 B/s objective | — | 1.095 | 1.494 | — | above the ceiling on one card |

How each prediction is formed: the reference's measured rate times the
ratio of the binding pipe's work, reference over arm, holding η at the
reference's value.

What the two knobs do:

- **`PACKED_SQUARE_TABLE=1`** squares λ as follows. The low 64 bits are
  spread on the ALU as before. The contribution of each 5-bit window of
  coefficients 66–130 is read from a 13 × 32-entry table of already-reduced
  squares, 8,320 bytes appended to the table walk's shared memory. The
  layout is structure-of-arrays, so a warp's 32 reads hit 32 banks. Shared
  memory per block goes from 48,732 to 57,052 bytes, which is still within
  the 64 KB carve-out at one block per SM. So the L1 keeps the 64 KB that
  ONE-BLOCK-GEOMETRY measured.
- **`PACKED_INV_POLY=1`** keeps the Itoh–Tsujii accumulator in both bases.
  The Frobenius powers stay in the ONB, where they are permutations, and the
  products stay in the β-basis. That saves the conversion that each ONB
  product performs on its way in.
- **`PACKED_INV_POLY=2`** is the same arithmetic with one out-of-line copy of
  the link's product and reduction, as the reference inverse keeps `mul131`
  out of line. The eight inlined links of `=1` grow the walk kernel from
  4,920 to 6,856 instructions and cost 10 registers; `=2` costs 272
  instructions and 1 register for the same dynamic work, ±0.2 lane-instructions
  per update.
- Registers per thread and kernel size in instructions:

  | build | registers | instructions |
  |---|---:|---:|
  | ref | 116 | 4,920 |
  | sqtab | 124 | 4,912 |
  | invpoly | 127 | 6,856 |
  | invpoly2 | 117 | 5,192 |
  | both | 127 | 6,856 |
  | both2 | 124 | 5,192 |

  All six builds run 512 threads per block with no spills.

## 7. Falsification targets, declared before the GPU run

For the knobs (`gpujob.sh`, REPS = 5). The primary arm is **both2**, chosen
now for its smaller code and register count at equal predicted work;
`both` is reported beside it and cannot substitute for it after the fact.

- **Success**: `both2` has a paired median ratio ≥ **1.040** against `ref`,
  every repetition > 1.0, 300/300 verified on all seven builds, and
  identical distinguished-point sets. `fused2` is exploratory and is
  reported, not scored. Then the two knobs join
  `gpu-rtx-pro6000-20b`, class engineering.
- **Partial**: 1.000 < ratio < 1.040. Then the knobs stay optional, and the
  gap to the 1.07 prediction goes on the page as the model's error on
  integer cuts.
- **Failure**: ratio ≤ 1.000, or any mismatch or DP difference. Then the
  knobs stay off. The model's "time follows ALU work" is refuted as a
  predictor of integer cuts. The likely suspects are LDS latency on the λ²
  path, the extra registers, and, for `=1`, instruction-cache pressure from
  a 39% larger kernel. The `invpoly` / `invpoly2` pair isolates that last
  one.

For the pipes (`pipes.cu`):

- **M4 over M5**: LOP3+IMAD 1:1 ≥ 1.5 × the faster single stream. Then
  moving SHFs onto the FMA pipe is the next integer lever.
- **M5**: that ratio < 1.2. Then the kernel is at 96% of its integer
  datapath, and only a cut in integer work moves it.
- **Speed of light**: `CLMAD` pairs beside a 40-LOP3 chain sustain ≥ 1.90
  lanes per SM-clock. Then 27.4 B/s stands. If they sustain ≤ 1.70, the SOL
  row becomes 22.7–23.3 B/s and the kernel is co-bound.

## 8. Tooling corrections made on the way

Class accounting. No earlier receipt changes.

- `roofline.py` initially dropped make variables passed to a preset target:
  the child make inherits them through `MAKEFLAGS` and the flag extraction
  did not. All four arms of §6 priced identically until it was fixed. The
  extracted flags are now checked equal to a real `make -n` for three
  overrides.
- With register spills (the 768-thread arm), ptxas closes the fruitless-cycle
  retry loop with an unconditional back-edge. The tool gave that loop's exit
  probability as 0, so all later flow vanished and the arm's CLMAD
  self-check read 0. The exit is now `1 − P_RETRY`, and all 46 builds pass.

## 9. Files

- [`roofline.py`](roofline.py), [`roofline_calibrate.py`](roofline_calibrate.py):
  the estimator and its calibration.
- [`benchmarks/roofline/`](benchmarks/roofline/): the probe (`pipes.cu`),
  its SASS audit (`sass_loops.py`), the GPU job, the summarizer, and the
  frozen calibration and reports.
- The knobs:
  - `PACKED_SQUARE_TABLE` and `PACKED_INV_POLY` in the Makefile;
  - `squarePolynomialTable131` / `fillSquareTable131` / `invPoly131` /
    `invPoly131Product` in `include/packed131.h`;
  - the table's place in `include/packedtablewalk.cuh`;
  - the call sites in `include/packedkernels.cuh`.
- Host tests in `src/testpacked.cpp`, run by `make test-packed-network`
  (CI runs it through `make certify-local`):
  - 20,134 table squares, bit-identical to the spread-then-reduce square and
    equal to the reference field's square;
  - 2,133 inverses, bit-identical to the conversion-sandwiched ONB inverse
    and equal to the reference field's inverse.
