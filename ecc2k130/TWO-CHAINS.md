# 30 B/s and the two-chain kernel

Question: can the table-walk kernel on one RTX PRO 6000 be taken from the
verified **20.078 B/s** of [ONE-BLOCK-GEOMETRY.md](ONE-BLOCK-GEOMETRY.md) to
**30 B complete scalar updates per second**?

Answer, derived before anything was built (§1): **no, on this GPU.** 30 B/s
is 15.2 SM-clocks per update, and the five field products of one affine
addition alone occupy the carry-less unit for 18.5. The whole kernel's 33.1
`CLMAD`s per update put its floor at **22.3 B/s** at 100% of the unit, and
the 20.08 B/s build already sits at 0.90 of it. The remaining tenth is the
only thing an iteration-function change can still buy on this SKU, and it is
what this note builds for: a kernel in which every warp carries two
independent Montgomery chains, so the inversion and the forward pass — the
two phases the phase profile found idling the unit — overlap inside the warp
instead of waiting for other warps that are in the same phase. Same
arithmetic, same walk, same products per update; static cost per update
unchanged to within 0.5% (§3). **Its rate is not yet measured** (§5 records
why and the exact command). The GPU that could change the answer is a
different one (§6).

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
more). The other execution resources are not the constraint — Nsight on the
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
of `--bench --steps 1024 --launches 64` at the automatic worker count, and
are **pending** (§5).

| variant | threads × blocks / SM | batch | SASS | regs | spills | instr / update (static) | ALU slots / update (static) | `CLMAD` / update (static; dynamic) | B/s | / 22.3 floor | / 20.078 | verified | class |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| **reference** `gpu-rtx-pro6000-20b` | 512 × 1 | 16 | 4,920 | 116 | 0 | 2,114.6 | 1,921.4 | 40.25; 33.1 | 20.078 (ONE-BLOCK-GEOMETRY) | 0.90 | 1.000 | 300/300 | reference |
| reference + `PACKED_ALU_SQR=1` | 512 × 1 | 16 | 5,352 | 116 | 0 | 2,141.5 | 1,953.8 | 39.00; 31.85 | *pending* (19.36 vs 19.43 at an earlier stage of the build, ONE-BLOCK-GEOMETRY §2) | | | | engineering |
| **two chains** `gpu-rtx-pro6000-chains2` | **256 × 1** | **32 = 2 × 16** | 8,320 | 182 | 0 | 2,123.3 | 1,910.7 | 40.25; 33.1 | *pending* | | | *pending* | engineering |
| two chains | 384 × 1 | 32 = 2 × 16 | 8,336 | 166 | 0 | 2,125.8 | 1,909.5 | 40.25; 33.1 | *pending* — 118 MB blob, over the L2 window | | | | geometry scout |
| two chains | 512 × 1 | 16 = 2 × 8 | 8,336 | 128 | 0 | 2,481.2 | 2,244.2 | 44.50; 37.35 | *pending* — one more inversion per 16 updates: floor 19.7 B/s | | | | geometry scout |
| two chains + `PACKED_ALU_SQR=1` | 256 × 1 | 32 = 2 × 16 | | | | | | 39.00; 31.85 | *pending* | | | | engineering |
| **30 B/s** | | | | | | | | ≤ 24.6 dynamic | 30 | **1.35** | 1.49 | | **below the floor** |

Reading the static columns: the two-chain kernel is a rescheduling. Its
instructions per update are within 0.5% of the reference's, its `CLMAD`
count is identical, and its SASS is 1.7× longer because both passes carry
two slots per iteration. What it changes is not in this table; it is the
fraction of the 20.4 SM-clocks of unit time per update during which the unit
is actually busy, and only the card measures that.

Per `AGENTS.md` §3 every row is **engineering**: the floor is the cost of
the same generic algorithm, the ratio to it is bounded above by 1, and no
row can be an advance. A row that lowers `CLMAD`/update without changing the
walk (`ALU_SQR`) moves the floor itself by 3.8%, which is the most any row
here can move it.

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

*Pending.* The job is [benchmarks/two-chains/gpujob.sh](benchmarks/two-chains/gpujob.sh):
inside `nvidia/cuda:13.3.1-devel-ubuntu24.04` it builds the six binaries of
§3 plus `PHASE_PROFILE=1` forms of the reference and the two-chain kernel,
verifies each (300 re-walks, the distinguished-point sets written and
compared), benches them alternating for `REPS` rounds with SM clock, power
and temperature sampled after every run, and prints the two phase profiles.
Three launchers run that same script on one RTX PRO 6000 and bring
`results/` back with a `launch.json` receipt:

```sh
cd ecc2k130
modal run modal_job.py --job benchmarks/two-chains/gpujob.sh --out /tmp/two-chains --env REPS=5
RUNPOD_API_KEY=… python3 runpod_job.py --job benchmarks/two-chains/gpujob.sh --out /tmp/two-chains --env REPS=5
python3 aws/bench_job.py --job benchmarks/two-chains/gpujob.sh --out /tmp/two-chains --env REPS=5
```

What stopped the run in the session that wrote this note, recorded so the
receipt can say where its numbers did not come from:

- **EC2**: `RunInstances` returns `Blocked` in every region ("This account is
  currently blocked and not recognized as a valid account") while dry runs
  pass; AWS Health carries two open risk events on the account,
  `AWS_RISK_CREDENTIALS_EXPOSURE_SUSPECTED` and
  `AWS_RISK_ACCOUNT_CONSOLE_COMPROMISE`, with a support case AWS opened.
  Resource creation stays blocked until the account owner answers that case.
  S3 still accepts the presigned transfers the launchers use.
- **Modal / RunPod**: no `MODAL_TOKEN_ID`/`MODAL_TOKEN_SECRET` and no
  `RUNPOD_API_KEY` were available to the agent environment. Both launchers
  were written and dry-checked (entry points, generated shell, SDK
  signatures) and wait on a token.
- **Datacenter parts**: the account's P-instance quota is 0 in every region,
  so §6 cannot be tested from it either.

When the job has run, this section gets the `finished:` rates, the paired
ratios, the identity check, the two phase profiles, and the receipts under
`benchmarks/two-chains/`; §3's pending cells are filled and §4 is applied as
written.

## 6. The GPU where the question changes

Everything above is about the RTX PRO 6000, and on it the carry-less unit's
rate is the whole story. Three facts point at where that rate comes from:

- Nsight Compute attributes `CLMAD` to the **FP64 pipe** (THROUGHPUT-29B.md
  §5: "FP64/clmad 71%").
- The RTX PRO 6000's die (GB202) executes FP64 at 1/64 of FP32: **2 lanes
  per SM-clock**. The measured `CLMAD` rate, 1.62 lanes per SM-clock, is
  0.81 of that.
- NVIDIA's CUDA 13.3 announcement reports GHASH — one GF(2¹²⁸) product per
  16-byte block, at least six `CLMAD`s — at 6.3 TB/s on a B200, and calls it
  memory-bound. That is ≥ 8 `CLMAD`s per SM-clock on 148 SMs, five times
  this card's unit, as a *lower* bound.

If `CLMAD` rides the FP64 datapath, then on a full-rate-FP64 part (A100,
H100, B200: 64 FP64 lanes per SM-clock) the unit is 20–30× faster than here
and this kernel is no longer carry-less-bound at all. It would be bound by
issue (1,700 instructions per update at 128 lanes per SM-clock is 13.3
SM-clocks) or by the ALU pipe, and every ALU→`CLMAD` trade that lost on this
card because the unit was full — `PACKED_TOP_CLMAD` (−313 ALU slots,
+19 `CLMAD`), the `CLMAD` squarings, a `CLMAD` reduction — becomes a
straightforward win there. Whether that adds up to 30 B/s on one B200 (148
SMs at ~1.9 GHz need ≤ 9.4 SM-clocks per update, i.e. about 1,200
instructions at 100% issue) is three unmeasured assumptions deep and is
stated here as a **hypothesis, not a result**. The B200 receipt in
[B200.md](B200.md) (8.82 B/s) is the *shipping* σ-walk at the 256 × 2
geometry, which was ALU-bound on both parts and so says nothing about the
B200's `CLMAD` rate.

The decisive experiment is one probe, not a kernel:
`benchmarks/clmad-price/probe.cu` on a B200 or H100 (Modal rents both;
`ECC_GPU=B200 modal run modal_job.py --job <a script that builds and runs the probe>`).
If it reports a `CLMAD` costing a few LOP3 slots instead of 38, the 30 B/s
question moves to that part and the levers are the ones this tree already
has behind knobs. If it reports ~38, the unit is a fixed-rate iterative
multiplier everywhere and the per-GPU answer is the one in §1 on every SKU.
