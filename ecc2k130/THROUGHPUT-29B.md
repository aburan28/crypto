# 29 B/s on one RTX PRO 6000, priced before the run

Question: can the table-walk kernel on this RTX PRO 6000 be pushed from the
measured **17.298 B/s** of the 20 B/s attempt
([THROUGHPUT-20B.md](THROUGHPUT-20B.md),
[benchmarks/throughput-20b-gpu/summary.json](benchmarks/throughput-20b-gpu/summary.json))
to **29 B complete scalar updates per second**?

This note states the boundary first (`AGENTS.md` §1). The levers below are
engineering of the Pollard walk, not an ECDLP exponent claim. Field products
stay 5.3125 per update. The campaign default stays the shipping walk until a
row clears its own target with a verified answer.

## 1. Boundary

Unit: billions of completed scalar updates per second (`finished:` line,
`codegen.benchreport.parseRate`). Same SKU, one GPU.

Pipe rates from ITERATION-FUNCTION.md §1 and the 20 B/s receipt's clocks
(188 SMs, 2.43 GHz max, logic pipe 64 lanes/SM-clock, `CLMAD.lo` 1.684
lanes/SM-clock):

| target | SM-clocks / update | ALU slots at 100% | clmad at 100% |
|---:|---:|---:|---:|
| 17.298 B/s (measured best, slot-unroll 2) | 26.4 | 1,690 | 44.5 |
| 20 B/s | 22.8 | 1,461 | 38.5 |
| 22 B/s (practical one-add line) | 20.7 | 1,327 | 34.9 |
| **29 B/s** | **15.8** | **1,008** | **26.5** |

**Floor (derived, not measured):** any rho walk that keeps the class structure
is one affine addition per step. At batch 16 that is 4.81 polynomial products
plus a batched inversion (5.3125 field products) and one squaring.
Karatsuba-on-`clmad` is 6 `clmad` per product. Squarings moved off `clmad`:

| component | clmad / update | ALU / update |
|---|---:|---:|
| 4.81 products + 0.5 inversion products | 28.9 | ≈ 1,090 |
| **29 B/s budget at 100%** | **26.5** | **1,008** |
| **29 B/s budget at 95%** | **25.2** | **958** |

Both pipes' 29 B/s budgets sit **below** the arithmetic floor. 29 / 26.5 ≈
1.09 on the carry-less unit even if every other cost is free; 1,090 / 1,008 ≈
1.08 on the ALU even if the walk is only the addition. This is the same shape
as ITERATION-FUNCTION.md's 28 B/s verdict.

**Reference:** table walk + byte pivot + pair-ILP + L2 persist + slot unroll 2
on this SKU, median **17.298 B/s**, 50,465,865,728 updates/sample, 300/300
reports replayed, 0 dropped.

**Class of every row below:** engineering. The ratio to the one-add floor
cannot fall below 1 without fewer products or a cheaper product.

## 2. What would have to move

17.298 → 29 is ×1.68. Instruction cuts that do not change the product count
cannot cross the floor. The leftover ALU in THROUGHPUT-20B.md that was not a
ceiling reduction was the 3-bit top-word correction (65 of `product131`'s 77
ALU, ~345 slots/update at 5.3125 products) and the inverse chain's eight
`mul131` calls, which convert to the polynomial basis and back even though
`inv131` is already in the ONB.

Those are the two arithmetic knobs this round adds, plus one scheduling knob
the 16-step Nsight DRAM share reopened:

| lever | priced effect | ceiling |
|---|---|---|
| `PACKED_TOP_HOIST=1` | rewrite the 3-bit correction with hoisted masks; bit-identical, 0 extra `clmad` | does not lower a pipe ceiling; SASS decides if it is fewer ALU than the k-loop |
| `PACKED_ONB_INV=1` | `inv131` uses the two-product ONB multiply instead of `mul131`'s convert-product-convert | +3 `clmad`/update (8 extra products / 16), ALU of the inverse chain should fall |
| `PACKED_SLOT_PIPELINE=1` | overlap the next compact-state load with the current product | scheduling; 0 static ALU |

`TABLE_ADDEND_GLOBAL=1` (selection in shared, addend in global, ~14 KB) is
wired so three blocks can fit the 100 KB SM, but Goal28 already measured
`MINBLOCKS=3` at 80 registers as slower on an issue-bound walk, so it is a
scout, not a 29 B/s path.

Inadmissible: a second GPU counted as one, a different SKU, dropping identity
or the 5.3125 product count, quoting a model as a rate, or scraping
boost-clock progress lines.

## 3. Acceptance, written before the run

Success is a verified median **> 29.0 B/s** on this RTX PRO 6000 Blackwell
Server Edition, identity matching the requested walk, every sample valid,
same 5.3125 field products per update, 300 device reports replayed with 0
dropped.

If every arm stays ≤ 29.0, the thread's answer is **no**, and the remaining
distance is the product and the reduction. A best row that beats 17.298
without clearing 29 is still engineering progress on the 20 B/s thread's
scoreboard, not a 29 B/s result.

## 4. Recipe

From `ecc2k130/`, GPU idle, CUDA 13.3.73 on `PATH`:

```sh
bash benchmarks/throughput-29b-gpu/run.sh
```

## 5. Measured, this SKU

Frozen: [benchmarks/throughput-29b-gpu/summary.json](benchmarks/throughput-29b-gpu/summary.json).
Unit: `finished:` B complete scalar updates/s. Each verified row replayed 300
reports with 0 dropped. 50,465,865,728 updates/sample. Class: engineering.

Round 1 priced the 3-bit hoist, ONB inverse, and slot pipeline against the
17.298 control (table + byte pivot + pair-ILP + L2 persist + unroll 2).
Nsight Compute on that control (sudo, `walk`, this card) is why round 2
moved the square: DRAM 11%, ALU pipe 36%, **FP64/`clmad` 71%**, 16 warps/SM.

| variant | median B/s | / 29 | / 22 B floor | / 17.298 | class |
|---|---:|---:|---:|---:|---|
| control (unroll 2) | 17.298 | 0.596 | 0.786 | 1.000 | reference |
| `TOP_HOIST` | 17.270 | 0.596 | 0.785 | 0.998 | engineering, did not pay |
| `ONB_INV` | 16.460 | 0.568 | 0.748 | 0.952 | engineering, did not pay |
| `SLOT_PIPELINE` | 16.463 | 0.568 | 0.748 | 0.952 | engineering, did not pay |
| hoist + ONB | 16.490 | 0.569 | 0.750 | 0.953 | engineering, did not pay |
| hoist + ONB + pipe | 16.498 | 0.569 | 0.750 | 0.954 | engineering, did not pay |
| **`ALU_SQUARE`** | **17.414** | **0.600** | **0.792** | **1.007** | **engineering** |

`ALU_SQUARE=1` is the only arm that beat the control: three alternating
samples 17.418 / 17.405 / 17.414 against a matched 17.298-class control
(17.296 / 17.308 / 17.298). 300/300 reports. Nsight after the change:
FP64 63%, ALU 39%, SM throughput 73%, issue slots 51%, IPC 2.07. The
squaring `clmad`s were on the busy pipe; the 3-bit hoist and ONB inverse
were not.

Inadmissible occupancy scouts (addend-global, 3 blocks at 80 registers,
128×5) all ran slower. `FROM_REDUCED` cut isolated `k_fromPoly` 61→48 and
measured a wash. Overlapping the 3-bit correction with `clmul128` still
cost 116 registers and 17.338 B/s.

Round 3 took the Nsight leftover (issue slots 51%, FP64 63%) as ILP and
occupancy. Matched 3-rep ALU_SQUARE control **17.403**. Every new arm
verified 300/300 and ran slower:

| variant | median B/s | notes |
|---|---:|---|
| `ALU_SQR` (ONB square on ALU) | 17.253 | one sample |
| `PAIR_CLMUL` (12 outstanding clmads) | 17.260 | 112 regs, one sample |
| `CLMUL_FLAT` (lo then hi) | 17.297 | one sample |
| `UNROLL_SLOTS=4` | 17.379 | 112 regs, 3-rep |
| `UNROLL_SLOTS=8` | 16.000 | 106 regs, 3-rep |
| combo of the three knobs | 17.214 | one sample |
| 192×3 + addend-global | 15.808 | 94 regs, 0 spill; SM clock 2340 MHz |

Nsight on the 17.41 kernel (sudo, `walk`): 3.44 active warps/scheduler of 12,
1.04 eligible, issue 52%. Dominant stall is `wait` (7.6%), then
`not_selected` (4.8%) and `math_pipe_throttle` (3.9%). Long scoreboard is
only 3.0%. The 192×3 scout did resident 3 blocks and still lost: extra
warps pulled the boost clock down and the addend missed LDS.

**Against the falsification target.** Success was a verified median > 29.0.
Measured best **17.414**. 29 B/s on one RTX PRO 6000 is not a result this
tree has. The one-add floor is still 22–25 B/s; this row is 0.792 of 22.
The campaign default stays the shipping walk. Knobs stay off by default.

The leftover is the 5.3125 products. Filling issue slots without cutting
`clmad` count does not move the ratio to the floor.
