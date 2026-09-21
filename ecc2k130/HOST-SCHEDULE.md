# Scheduling the generated straight-line code for the host register file

The generated `GF(2^131)` routines are emitted in the order the DAG was built.
That order was never chosen; it is whatever `codegen/ir.py` happened to append.
On a GPU it costs little, because a thread has 255 registers and `chooseLeaf`
already sized the leaf to fit them. On an x86 host the wide word is a `zmm`
register and there are 32 of them, so a routine with 254 simultaneously live
values spills, and a spilled value costs two memory instructions every time it
is touched.

This note asks whether reordering the same DAG -- same operations, same count,
same result -- recovers the spill traffic, and what that is worth end to end.

## 1. Boundary and unit

The unit is **instructions executed per `GF(2^131)` multiplication**, counted
from the disassembly of the built binary. Throughput in iterations/s is a
practicality note, never the metric (`AGENTS.md` §6).

*Floor.* A schedule may reorder the DAG but cannot change it: the operation
count after LOP3 fusion is invariant under scheduling. So

```
I_floor = opCount(leaf) + opCount(multPrep) + opCount(toOnb)
```

is a lower bound on instructions per multiply that no schedule can cross, and
every instruction above it is data movement the schedule is responsible for.
The reported ratio is `I / I_floor`, which is 1.0 for a routine that never
spills and grows with the spill traffic. For the shipping `m = 131` generator
`I_floor = 2835 + 331 + 903 = 4069`.

*Reference.* The shipping build: `SIZE_CHAIN = {131, 66}`, construction-order
emission, g++ 13.3 `-O3 -march=native`, AVX-512 (`Bits512`, 512 lanes).

## 2. Falsification target, declared in advance

The scheduler is a **success** if, on the `m = 131` AVX-512 host build:

* peak liveness of `mulLeaf` falls below **128** (half the register budget it
  was built for), and
* `I / I_floor` for the leaf falls to **1.5 or below** (from the reference
  measured in §3), and
* measured throughput does not regress: the paired median of the candidate is
  at or above the reference.

It is **abandoned** if the scheduler cannot bring the leaf's `I / I_floor`
below **1.8**, or if throughput regresses on the paired comparison despite the
instruction count falling.

Inadmissible: changing the operation count, the leaf size, the size chain, the
compiler, or the lane width between the two arms; reporting the leaf's
instruction ratio without the end-to-end throughput pair.

## 3. What the pass does

`ir.Prog.scheduleLive` scores two emission orders and keeps the cheaper, per
routine, by the number of C locals `emit()` allocates -- the same
first-fit-with-free-list the emitter runs, so the score is exactly what the
compiler is handed. `slotCount` is checked against `emit()`'s own count in
`codegen/testirschedule.py`.

* **Construction order**, what the generator emitted before this change. For a
  Karatsuba leaf this is the builder's subproduct-by-subproduct blocking,
  which is already a good schedule.
* **Sethi-Ullman depth-first order** (`suOrder`): roots one at a time, each
  subtree depth first, heavier operand first. Optimal on a tree, and the
  conversions are close to trees -- 131 independent output accumulations that
  hash-consing has partly merged.

A greedy "free the most registers" list schedule was tried first and is worse
than both (447 live in the m=131 leaf against construction order's 254): in a
schoolbook leaf every ready operation is an AND over array inputs, so nothing
frees anything and the tie-break chooses blind. It is not kept.

`--no-schedule` is the paired control and reproduces the previous headers byte
for byte, which is how the arms below were built.

## 4. Results

One `GF(2^131)` multiplication, g++ 13.3 `-O3 -march=native`, AVX-512
(`Bits512`, 512 lanes), instructions counted from `objdump` of the built
client. `I_floor` is the routine's post-fusion operation count, which no
schedule can cross.

| routine | `I_floor` | control `I` | scheduled `I` | control `I/I_floor` | scheduled `I/I_floor` | class |
|---|---|---|---|---|---|---|
| `mulLeaf` | 2835 | 5693 | 5693 | 2.008 | 2.008 | unchanged (construction order wins) |
| `toOnb` | 903 | 1918 | 1473 | 2.124 | **1.631** | engineering |
| `multPrep` | 331 | 745 | 806 | 2.251 | 2.435 | regression, see below |
| **per multiply** | **4069** | **8356** | **7972** | **2.054** | **1.959** | **engineering** |

Locals allocated, generator-reported: `toOnb` 233 -> 80, `multPrep` 116 -> 105,
`mulLeaf` 254 (unchanged), `hamming` 89 (unchanged).

**The target declared in §2 is missed.** It asked for the leaf below 128 live
and `I/I_floor <= 1.5`; the leaf is unchanged at 254 and 2.008, and the
per-multiply ratio only moves 2.054 -> 1.959. The abandon threshold (leaf
above 1.8) is also hit. The gain is real but it is 4.6% of instructions, an
**engineering** step by `AGENTS.md` §3, not an advance: the ratio to the floor
is essentially flat and the leaf, which is 70% of the multiply, did not move at
all.

**Slot count is an imperfect proxy.** `multPrep` cut locals 116 -> 105 and
emitted 15% *more* data movement. `toOnb` cut locals 233 -> 80 and emitted 44%
less. The pass optimises the proxy, not the spills, and it is right about the
sign only when the liveness change is large. No threshold is fitted here to
make `multPrep` pick the other arm: one sample is not enough to justify one,
and the net across the multiply is still -384 instructions.

**Throughput could not be resolved on the available host.** Paired, interleaved
runs on a shared 4-vCPU container, `--bench --steps 32`, one thread:

| arm | median M it/s (n=11) | range |
|---|---|---|
| control | 15.12 | 14.78 - 20.09 |
| scheduled | 15.69 | 15.31 - 19.44 |
| control, pinned, 192 launches (n=7) | 16.39 | 15.40 - 20.25 |
| scheduled, pinned, 192 launches (n=7) | 18.64 | 15.68 - 21.46 |

Paired sample deltas run from -22% to +31%. The medians favour the scheduled
arm by 3-16%, consistent in sign with the -4.6% instruction count, but the
container's noise is several times the effect and no throughput claim is made
from it. `AGENTS.md` §6 puts wall time in the practicality column anyway; the
instruction count is the result. A paired rerun on the 2.1 GHz Xeon of
README.md's table would settle the size.

## 5. The leaf is already well scheduled

Four root orders were scored on the m=131 leaf. None beats the order the
builder produces:

| order | locals |
|---|---|
| construction (Karatsuba blocking) | **254** |
| Sethi-Ullman, roots in natural order | 303 |
| Sethi-Ullman, roots reversed | 311 |
| Sethi-Ullman, roots by register need | 340 |
| Sethi-Ullman, roots by index | 424 |
| greedy "free the most registers" | 447 |

This is the DAG penalty: Karatsuba shares subproducts across many outputs, so
finishing one output first pins shared values live across the whole routine,
where the blocked order consumes each subproduct while it is hot.

The leaf's 2.008 therefore stands, and with it the 46% of the host multiply
that is spill traffic. Cutting it needs something other than a schedule --
fewer live values, not a better order for 254 of them.

## 6. Leaf size on x86, for the record

`Makefile`'s `REGS` note records this experiment on an M4 Pro with NEON words.
It replicates on AVX-512, paired medians of 7-9 interleaved runs on the same
container, same compiler and lane width:

| `REGS` | leaf | chain | M it/s | vs. shipping |
|---|---|---|---|---|
| 255 (shipping) | 66 | 131, 66 | **16.6** | -- |
| 128 | 33 | 131, 66, 33 | 15.4 | -6.8% |
| 32 | 9 | 131, 66, 33, 17, 9 | 10.7 | -36% |

A host-sized leaf remains measurably wrong, for the reason `chooseLeaf`'s
docstring gives: the C++ Karatsuba recursion costs more than the spills it
saves. The 32-register budget that x86 actually has is the worst of the three.

Also measured and not adopted: clang 18.1 against g++ 13.3 on the same source,
paired medians 14.39 against 15.36 M it/s, so clang is ~6% *behind*. The 3.02x
clang-over-gcc result in README.md's aarch64 table does not transfer; the host
multiply spills because 254 > 32, not because gcc mishandles it.

## 7. What this leaves

`VPCLMULQDQ` is present on the benchmark host and unused, and should stay that
way: the bitsliced leaf costs about 11 instructions per walk-multiply amortised
over 512 lanes, where a per-walk carryless 131-bit multiply needs roughly nine
`vpclmulqdq` plus reduction for one walk. The packed backend's native carryless
win does not carry to a host word.

`GFNI` is present and unused, and is the open lever. `vgf2p8affineqb` applies a
constant 8x8 `GF(2)` matrix per byte in one instruction, which is what
`multPrep` and `toOnb` are -- 1234 of the multiply's 4069 operations, the two
routines this note could only reschedule. That is a different pass, not a
schedule, and it is not attempted here.
