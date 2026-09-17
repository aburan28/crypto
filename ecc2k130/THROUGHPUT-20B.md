# 20 B/s on one RTX PRO 6000, priced

Question: can the g7e client be pushed from the measured 16.56 B/s of the table
walk ([ITERATION-FUNCTION.md](ITERATION-FUNCTION.md) §6) to 20 B complete
scalar updates per second on one RTX PRO 6000?

Answer: **not with any lever this tree can price.** 20 B/s is above the
one-addition-per-step floor (≈ 22–23 B/s at 90% of either pipe, §1 of that
note), so unlike 28 B/s it is not forbidden — but it needs the ALU pipe to
shed **19–22% of its dynamic work** on top of the table walk while the
carry-less pipe sheds the per-update squaring, and every addressable item
this profile finds sums to about 3%. The reduction and the products are 60%
of the walk's ALU, and they are the floor.

The new evidence is a **per-routine static profile of the table walk**, made
with the exact compiler behind the tree's measured receipts (nvcc 13.3.73)
and no GPU. Everything below is from
[benchmarks/throughput-20b/](benchmarks/throughput-20b/); the toolchain that
produced it, which needs neither a card nor a system CUDA install, is in
[toolchain.txt](benchmarks/throughput-20b/toolchain.txt) there.

## 1. Reproduction first

`kernel_cost.py` run here on the shipping preset and the table walk reproduces
the committed receipts **to the decimal**:

| build | ALU slots / update | clmad / update | SASS | registers | receipt |
|---|---:|---:|---:|---:|---|
| shipping walk | 2,323.5 | 45.25 | 4,048 | 104 | [kernel-cost-shipping.txt](benchmarks/throughput-20b/kernel-cost-shipping.txt), matches [table-walk/raw/kernel-cost-legacy.txt](benchmarks/table-walk/raw/kernel-cost-legacy.txt) |
| table walk, `WALK_TABLE=1` | 1,830.7 | 45.25 | 3,544 | 102 | [kernel-cost-table.txt](benchmarks/throughput-20b/kernel-cost-table.txt), matches [table-walk/raw/kernel-cost-table-lut.txt](benchmarks/table-walk/raw/kernel-cost-table-lut.txt) |

So the static prices in this note are on the same footing as the ones the
6000 measurements were set against, not a second compiler's opinion.

## 2. The two budgets at 20 B/s

188 SMs, 2.40 GHz (the table walk ran the 6000 at 2385–2407 MHz), the logic
pipe at 64 lanes per SM-clock and `CLMAD` at 1/38.0 of it = 1.684 lanes per
SM-clock, both measured on the 6000 (§3.1, §4.5, §6.3 of ITERATION-FUNCTION):

| rate | SM-clocks / update | ALU slots at 100% | clmad at 100% |
|---:|---:|---:|---:|
| 16.56 B/s (measured, table walk) | 27.2 | 1,744 | 45.9 |
| 17.0 B/s (the §4.4 line) | 26.5 | 1,699 | 44.7 |
| **20 B/s** | **22.6** | **1,444** | **38.0** |

Against those, the table walk's static counts are 1,857 ALU slots
(`IMNMX` at its measured 1.81; 1,830.7 at `kernel_cost.py`'s 1.0) and 45.25
clmad. Static counts both arms of the slot-0 branches; taking the arms that do
not run (the slot-0 `mulPolynomial131` runs once per batch, the pair products
fifteen times) puts the dynamic count near **1,670 ALU and 37.7 clmad**, and
the measured rate bounds it from the other side at ≤ 1,744. At 16.56 B/s the
ALU is therefore at 96–100% and the carry-less unit at 82%, which is what
the neutral `ALU_SQUARE` measurement in §6.3 said: the ALU binds.

At 20 B/s the carry-less unit would be at 37.7 / 38.0 = **99%** with the
squaring where it is, so the squaring must move to the ALU (`ALU_SQUARE=1`,
−5 clmad, +100 ALU, measured neutral at 16.56 because the ALU binds). That
leaves clmad at 32.7 (86%) and the ALU needing to fall from ~1,770 to 1,444
at 100% or ~1,370 at a practical 95%: **−330 to −400 dynamic slots, 19–22%.**

## 3. Where the 1,857 go

`kernel_attribution.py` (new; §5) attributes every `__noinline__` callee as
one unit by its symbol, weighted by call sites, and the kernel body's own
instructions by source function. Table walk, per update, static
([attribution-table.txt](benchmarks/throughput-20b/attribution-table.txt)):

| routine | slots | of which | clmad | notes |
|---|---:|---|---:|---|
| `mulPolynomialPair131` ×2 | 622 | 311 per call = 2 × (77 product + 78 reduction) | 24 | prefix chain, reverse pass |
| `mulPolynomial131` ×2 | 322 | 161 per call; one is the slot-0 arm | 12 | λ·(x+x′), slot-0 λ |
| `mul131` ×0.5 | 180 | 360 per call: two `toPolynomial131`, product, `fromPolynomialProduct131` | 3 | the inversion's 8 products / 16 |
| squaring (`spread32p` + `reducePolynomial131` in the body) | 83 | 8 + 75 | 6.25 | λ² |
| Frobenius networks in the inversion, / 16 | 82 | three `sigmaInvNetwork131`, two `sigmaWalkNetwork131` | 0 | |
| **arithmetic** | **1,289** | | **45.25** | ITERATION-FUNCTION's floor row was ≈ 1,090 with the inversion at 268; this is the same work priced routine by routine |
| `fromPolynomialProduct131` (x → normal basis) | 115 | one per update, for `HW`, `h`, `k`, `p` | | required: the walk's invariants are normal-basis |
| pivot: 32 `IMNMX` (`math_functions.hpp`) + 33 `PRMT` in `twByte` + masks | ≈ 110 | 32 × 1.81 + 33 + ~20 | | **the one item this note moves**: 16 `IMNMX` + 17 `PRMT` with `TABLE_PIVOT_BYTES` |
| phase `twPhase` + its 17 `PRMT` | 72 | | | 17 `LDS.U8` besides |
| sign `twCoordinate` | 35 | | | 5 `LDS` besides |
| tags, cycle rule, history | 15 | | | |
| weight, DP test, report path, atomics | 66 | 22 + 44 | | the report path is static; it runs on ~2^−25 of updates |
| compact state load/store, `add131`, `toLimbs` | 55 | | | |
| kernel body: loop, addressing, guards | 89 | | | |
| **total** | **1,857** | | **45.25** | |

Reading it against §2: the arithmetic alone (1,289 static, ~1,110 dynamic) is
77% of the 1,444 the ALU affords at 20 B/s, before the conversion the walk
cannot avoid (115) and before any selection. The non-arithmetic remainder is
≈ 570 static, of which the DP report path (~60) does not run and the rest is
what makes the walk a class function. A 330–400 cut has to come out of the
arithmetic, and the arithmetic is one affine addition per step.

## 4. Every lever, priced

Static, per update, on top of the table walk. A lever that lowers one pipe's
ceiling below the current rate is listed because it keeps being proposed.

| lever | ALU | clmad | verdict |
|---|---:|---:|---|
| **`TABLE_PIVOT_BYTES=1`** (this note; §5) — pivot per byte of x, 17 lookups and 16 maxes for 33 and 32 | **−45** (−32 at `IMNMX` = 1) | 0 | 22 fewer `LDS`, same 102 registers, 0 spills; shared 48,732 B. If the ALU binds it is worth ~2.4%, i.e. **16.56 → ~16.96 B/s, within noise of the 17.0 line of §4.4**. Unmeasured on a card. |
| `ALU_SQUARE=1` | +100 | −5 | required for 20 B/s (clmad at 99% otherwise); neutral at 16.56 (§6.3) |
| `TOP_CLMAD=1` (TOP-CLMAD.md) | −313 | **+19.2** | 64.5 clmad → carry-less ceiling **11.8 B/s** at 100%. Never above 11.8 B/s. |
| one cross term on clmad (half of TOP_CLMAD) | −155 | +9.6 | 47.4 clmad → ceiling 16.1 B/s: below the current rate |
| reduction by `clmad` against the modulus's low terms | −290 | +19 | same wall as TOP_CLMAD |
| inversion chain kept in the polynomial basis (one `toPolynomial131` fewer per step) | ≈ 0 | 0 | `fromPolynomial131` of a 5-word input still costs ~100, so a step is ≥ the current 347; priced and dropped |
| Itoh–Tsujii in the polynomial basis (squarings by `spread`) | +630 | +41 | 130 squarings per inversion |
| batch 32 (halves the inversion share, −134) | −134 | −2 | measured **−34%** (BATCH-TUNING.md: 8.67 vs 13.2 B/s), state traffic |
| packed top words (part of `TABLE_PIVOT_BYTES`) | **+11** | 0 | `twCoordinate` 35 → 46 for the packed-`fromRow` address math; the price of the shared budget, already inside the −45 above |
| byte table for the phase | 0 | 0 | already bytes |
| `H = 4` (halves the table) | ~−10 | 0 | r-adding constant 1.125 vs 1.0625: +6% iterations for <1% rate |

Sum of everything that does not lower a ceiling below 16.56: **−45 static
ALU, 0 clmad**, against a need of −330 to −400 with the squaring moved.

## 5. What this change adds

**`TABLE_PIVOT_BYTES=1`** (Makefile knob; `-DECC_TABLE_PIVOT_BYTES=1`; off by
default) in [include/packedtablewalk.cuh](include/packedtablewalk.cuh):

- the pivot's max-L scan reads one byte of x per lookup (a 17 × 256-byte
  table) instead of one nibble (33 × 16): 17 `PRMT`/`LDS.U8`/`IMNMX` instead
  of 33;
- the byte table is 3,824 bytes larger than the nibble table and the walk's
  tables were at 48,508 of the 49,152 bytes that keep two blocks per SM, so
  the layout also packs the 3-bit top words: the table's `x | y << 3` tops as
  one byte per entry (two words per `k` at `H = 8` instead of one word per
  entry, −3,144 B) and `fromRow`'s top as a nibble (17 words for 131, −456 B).
  Net **48,732 bytes**;
- the selection primitives compile for the host (`TW_FN`), and
  **`src/testtablewalkhost.cpp`** (`make test-table-walk-host`, in `make test`)
  runs the device-side phase, pivot, sign, tag and addend from the packed
  buffer against the reference for 4,096 random subgroup points, both
  layouts: **0 mismatches on every row, cycle rule fired on 1,536**, the same
  counts the device probe reports. The default layout's SASS is unchanged by
  the refactor: 1,830.7 slots, 3,544 instructions, 102 registers before and
  after.

Static price, same compiler
([kernel-cost-table-pivot-bytes.txt](benchmarks/throughput-20b/kernel-cost-table-pivot-bytes.txt),
[attribution-table-pivot-bytes.txt](benchmarks/throughput-20b/attribution-table-pivot-bytes.txt)):

| | slots (`IMNMX` 1.81) | slots (`kernel_cost`) | clmad | SASS | `LDS` | `VIMNMX` | `PRMT` | regs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| table walk | 1,856.6 | 1,830.7 | 45.25 | 3,544 | 71 | 32 | 50 | 102 |
| + `TABLE_PIVOT_BYTES=1` | 1,811.7 | 1,798.7 | 45.25 | 3,488 | 49 | 16 | 34 | 102 |

**Nothing here is a device measurement.** What a card decides: whether the
byte lookups keep the MIO pipe's conflict-free behaviour the nibble table had
(16-byte rows span four banks; 256-byte rows do not, and `twByte` values are
random, so the 17 loads should spread — but the tree's rule is that a static
count is not a rate), and whether the 2.4% arrives. The command is the one
ITERATION-FUNCTION §6 used, with the knob added to the table build:

```bash
# on the g7e, both binaries alternating, six repetitions, per benchmarks/table-walk/gpujob.sh
make -B gpu ARCH="-gencode arch=compute_120,code=sm_120" BATCH=16 THREADS=256 MINBLOCKS=2 \
     PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3 \
     PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 \
     PACKED_DIRECT_REDUCE=1 PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256 \
     PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 \
     WALK_TABLE=1 TABLE_PIVOT_BYTES=1
make test-table-walk-cuda WALK_TABLE=1 TABLE_PIVOT_BYTES=1      # 0 mismatches required first
./ecc2k130 --curve 131 --packed --bench --steps 1024 --launches 32 --verify 0
```

If it lands at or above 17.0 B/s it crosses the line §4.4 set for the table
walk to become the campaign default; that decision, and the corpus fork it
implies, stay the campaign owner's (§4.2, §4.4).

**`kernel_attribution.py`** (new) is the profile tool of §3. It builds with
`-lineinfo`, disassembles with `nvdisasm -g`, and attributes callee bodies as
units and body instructions by source line. It deliberately does not attribute
lines *inside* `__noinline__` bodies: ptxas's markers there smear across the
header that inlined into them (a first version put a third of
`mulPolynomialPair131` under `sigmaWalkNetwork131` that way). Its
`kernel_cost.py`-equivalent total matches `kernel_cost.py` exactly (1,830.7).

## 6. What would change the answer

- **A cheaper reduction.** The polynomial basis is the OPB with modulus
  `0xd1d0d000d0000000d000000000000000d`, whose generated reduction is 78 ALU
  and runs 6.3 times per update: 490 slots, 26% of the walk. A pentanomial
  basis would reduce in ~40, but its conversions to and from the ONB — one
  per update for `x`, one per inversion product — lose the structure that
  makes the present ones ~100 each and become dense 131 × 131 maps at ~200
  ALU plus 66 `LDS` by nibble tables (byte tables do not fit). Net −240 in
  reductions against +200 or more in conversions: **not a route to −330.**
- **Fewer products per addition.** No affine formula does it (§1 of
  ITERATION-FUNCTION); the count is 4.81 and the weighted prefix already
  merged `e` into the chain.
- **The second GPU.** §5 of ITERATION-FUNCTION: per-GPU work above ~19 B/s
  buys nothing a second RTX PRO 6000 does not buy cheaper, and the shipping
  kernel already does 14+ B/s on each. 20 B/s of ECC2K-130 is a
  `g7e.2xlarge` and a third of another.
