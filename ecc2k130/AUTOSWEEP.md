# Automatic sweep of the table walk on a B200

Question: [TWO-CHAINS.md](TWO-CHAINS.md) §6 found that on a B200 the
table-walk kernel is bound by the logic pipe with the carry-less unit
nearly idle, and that three knobs whose sign was negative on the RTX PRO
6000 are positive there. Every other knob in the tree, and the geometry,
was tuned on the 6000. What does the card itself pick when the whole knob
and geometry space is searched automatically?

Answer: a different point. On the B200 the sweep finds **19.40 B/s**, +26.3%
over the best hand-picked row, in three knobs the 6000 had rejected — one
fused pass per step, batch 32, and the three-limb Karatsuba product — and
the same sweep run on the RTX PRO 6000 finds nothing (§4): its one arm over
the line is a geometry that a warm card takes back in clock, and the B200's
winners are its largest losses. Both results are engineering; neither moves
a floor; 30 B/s is still nobody's.

## 1. Boundary

Unit: billions of complete scalar updates per second, `finished:` line, one
B200 (148 SMs, 1.965 GHz under load, `sm_100`).

**Floor:** the logic pipe. 63.7 `LOP3` lanes per SM-clock
([probe](benchmarks/fast-clmad/probe/summary.json)), so a kernel spending
`A` ALU slots per update cannot exceed `148 × 1.965 GHz × 63.7 / A`. The
carry-less unit at 29.1 lane-`CLMAD`s per SM-clock is not a floor at any
`CLMAD` count this kernel could reach (its 76 static `CLMAD`s per update are
2.6 SM-clocks against the 18.9 the kernel takes). **30 B/s needs `A ≤ 618`.**

**Reference:** the best row of TWO-CHAINS.md §6.4, the 20 B/s knob set for
`sm_100` plus `PACKED_TOP_CLMAD=1 PACKED_ALU_SQUARE=0 PACKED_ONB_INV=1`:
**15.372 B/s**, 1,554 static ALU slots, ≤ 1,205 dynamic, 18.9 SM-clocks per
update. It is the sweep's base point and is rebuilt and re-timed inside the
sweep, so every ratio below is against the same-session base.

**Target, declared before the run:** a combination that beats the base by
**≥ 2%** in the final alternating bench with every paired repetition above
1.0, 300/300 device reports re-walked, and the same distinguished-point set
as the base on the forced common walk count, is the new B200 build. A
single arm that beats it is reported as such; a sweep whose every arm loses
records that the 6000-tuned point is also the B200's, which would be a
finding about the kernel, not the card. Nothing here can be an advance: the
walk, the products per update and the floor are unchanged, and the ratio to
the floor is bounded above by 1.

## 2. Method

[benchmarks/autosweep/gpujob.sh](benchmarks/autosweep/gpujob.sh), through
`modal run --detach modal_job.py --job benchmarks/autosweep/gpujob.sh --gpu B200`:

1. **Star.** Build the base and twenty arms, each the base plus one change
   (`make gpu-preset KNOBS=…`, last assignment wins), and screen each with
   `--bench --steps 1024 --launches 32`, twice, all binaries alternating.
   The arms are every knob in the tree that the 6000 measured or rejected
   and that keeps the walk: geometry (256 × 2, 384 × 1, 640 × 1, 768 × 1;
   batch 32 at 256 and 512 threads), the product form (`KARAT3` in place of
   `TOP_CLMAD`, `CLMUL_FLAT`), scheduling (`TABLE_PIPE_SELECT`,
   `PACKED_CHAIN_FIRST`, `PACKED_INLINE_POLY` 0/1, `UNROLL_SLOTS=2`,
   `PACKED_PAIR_ILP` off), memory (`PACKED_L2_PERSIST` off, `TABLE_FUSED`
   with and without its load pipeline), and two small trades
   (`PACKED_FROM_REDUCED` off, `TABLE_PIVOT_BYTES` off, `PACKED_ALU_SQR`).
   Arms that change the walk (`TABLE_BRANCHES=16`) are excluded: their gain
   is in expected iterations, which a rate bench cannot see.
2. **Greedy.** Arms that beat the base by ≥ 0.5% screened, best first, are
   added one at a time; each combination is built and screened, and kept
   only if it beats the running best by 0.2%. Knobs that share the binding
   pipe do not add, and knobs that conflict fail to build and are dropped.
3. **Verify.** Base, best single arm, greedy result: 300 re-walks, 0 dropped
   required; distinguished-point sets on the forced 1,212,416 walks compared
   by hash.
4. **Final.** The same three, alternating, five repetitions of
   `--launches 64`, SM clock and power after each.

The arm list is the tree's knowledge; the ranking is the card's. Static
counts choose nothing here.

## 3. Measured on the B200

One B200 (Modal, driver 580.95.05, nvcc 13.3.73, 1965 MHz throughout,
2026-09-22; [summary.json](benchmarks/autosweep/b200/summary.json), raw logs
beside it). 21 builds screened, 6 combinations, 3 finalists; 50 minutes.

### 3.1 Star

Screened rate is the median of two alternating passes of `--launches 32`;
the two passes agree to 0.05% on every arm. Base: **15.372 B/s**.

| arm (base + …) | knobs | regs | spill | screened B/s | / base |
|---|---|---:|---:|---:|---:|
| **fused** — one pass per step | `TABLE_FUSED=1 TABLE_PIPE_SELECT=0` | 102 | 0 | **17.456** | **1.136** |
| fused, loads one slot ahead | `TABLE_FUSED=1 TABLE_FUSED_PIPE=1 TABLE_PIPE_SELECT=0` | 106 | 0 | 17.144 | 1.115 |
| batch 32, 8 warps | `BATCH=32 THREADS=256 MINBLOCKS=1` | 116 | 0 | 16.153 | 1.051 |
| **three-limb Karatsuba** in place of `TOP_CLMAD` | `PACKED_KARAT3=1 PACKED_TOP_CLMAD=0` | 113 | 0 | 16.062 | 1.045 |
| batch 32, 16 warps (124 MB blob, over the 83 MB L2 window) | `BATCH=32 THREADS=512 MINBLOCKS=1` | 116 | 0 | 15.909 | 1.035 |
| forward-pass software pipeline off | `TABLE_PIPE_SELECT=0` | 96 | 0 | 15.786 | 1.027 |
| chain product not first | `PACKED_CHAIN_FIRST=0` | 116 | 0 | 15.398 | 1.002 |
| pair-ILP off | `PACKED_PAIR_ILP=0` | 117 | 0 | 15.356 | 0.999 |
| slot unroll 2 | `UNROLL_SLOTS=2` | 118 | 0 | 15.341 | 0.998 |
| flat lo/hi issue | `PACKED_CLMUL_FLAT=1` | 114 | 0 | 15.307 | 0.996 |
| products out of line | `PACKED_INLINE_POLY=0` | 118 | 0 | 15.303 | 0.996 |
| 640 × 1 | `THREADS=640 MINBLOCKS=1` | 96 | 0 | 15.296 | 0.995 |
| 768 × 1 | `THREADS=768 MINBLOCKS=1` | 80 | 52 B | 15.237 | 0.991 |
| L2 persist off | `PACKED_L2_PERSIST=0` | 116 | 0 | 15.226 | 0.991 |
| 9-word inverse conversion | `PACKED_FROM_REDUCED=0` | 116 | 0 | 15.216 | 0.990 |
| 256 × 2 | `THREADS=256 MINBLOCKS=2` | 116 | 0 | 15.026 | 0.978 |
| 384 × 1 | `THREADS=384 MINBLOCKS=1` | 116 | 0 | 14.948 | 0.972 |
| single products inlined only | `PACKED_INLINE_POLY=1` | 114 | 0 | 14.850 | 0.966 |
| nibble pivot | `TABLE_PIVOT_BYTES=0` | 116 | 0 | 14.822 | 0.964 |
| ONB squarings on the ALU | `PACKED_ALU_SQR=1` | 116 | 0 | 14.621 | 0.951 |

Six arms clear the 0.5% line, and four of the six are knobs the 6000
measured as neutral or as losses: `TABLE_FUSED` was 19.22 against 19.43
there (ONE-BLOCK-GEOMETRY §5, "kept behind the knob because it is what a
larger batch would want"); batch 32 lost 34% there (BATCH-TUNING.md) and
12–15% in the same geometry (ONE-BLOCK-GEOMETRY §2); `TABLE_PIPE_SELECT`
was +0.4% there. The geometry the 6000 chose, 512 × 1 at batch 16, is
the B200's too among the one-chain block shapes — 256 × 2 is 2.2% slower
although this part has the shared memory for it, so the 6000's L1 story
is not what selects 512 × 1 here — but batch is where the two parts
differ: the B200's memory system (HBM3e, 8 TB/s) makes the state traffic
that sank batch 32 on the 6000 cheap.

### 3.2 Greedy

| step | combination | screened B/s | against | decision |
|---|---|---:|---:|---|
| combo1 | base + fused | 17.482 | 15.372 | keep |
| combo2 | combo1 + fused-pipe | 17.151 | 17.482 | drop |
| combo3 | combo1 + batch 32 at 256 threads | 17.753 | 17.482 | keep |
| combo4 | combo3 + `KARAT3` for `TOP_CLMAD` | 18.641 | 17.753 | keep |
| **combo5** | combo4 + batch 32 at 512 threads | **19.402** | 18.641 | **keep** |
| combo6 | combo5 + pipe-select off (already off in `fused`) | 19.397 | 19.402 | drop |

The result, `combo5`: the base plus `TABLE_FUSED=1 TABLE_PIPE_SELECT=0
PACKED_KARAT3=1 PACKED_TOP_CLMAD=0 BATCH=32 THREADS=512 MINBLOCKS=1`, 106
registers, no spills. Two things the greedy trace shows that the star
cannot: the fused pass and batch 32 *compound* (fused alone +13.6%, batch
32 alone +3.5% at 512 threads, together +24%, because the fused pass
halves the per-update state traffic that batch 32 doubles), and the
three-limb Karatsuba adds its full +4.9% on top — it removes the last of
the 3-bit correction's ALU by spending four more `CLMAD`s per product on a
unit that has them to spare.

### 3.3 Verify and final

| finalist | B/s, median of 5 | paired / base (min – max) | SM-clocks / update | logic pipe allows ≤ | verified | same DPs | class |
|---|---:|---|---:|---:|---|---|---|
| base (TWO-CHAINS.md §6.4 best row, rebuilt) | 15.365 | 1 | 18.9 | 1,206 | 300/300 | yes | reference |
| fused | 17.478 | 1.1375 (1.1372 – 1.1380) | 16.6 | 1,060 | 300/300 | yes | engineering |
| **combo5** = `make gpu-b200-19b` | **19.402** | **1.2627 (1.2623 – 1.2632)** | **15.0** | **955** | 300/300 | yes | **engineering; the B200 build** |
| RTX PRO 6000, `gpu-rtx-pro6000-20b`, same session series | 19.981 | | 20.4 (of a 22.3 B/s floor) | | 300/300 | — | reference, the other part |
| 30 B/s on a B200 | 30 | 1.95 | 9.7 | 618 | | | 1.55× ALU cut away |

Every finalist re-walked 300 of 300 device reports with 0 dropped, and the
three distinguished-point sets on the forced 1,212,416 walks hash
identically (`626a6fce…`, the same hash as every B200 binary in
TWO-CHAINS.md): the sweep changed the schedule and the batch, never the
walk. The final medians reproduce the screening medians to 0.1%
(19.402 against 19.402; 17.478 against 17.456), which is what the two
alternating screening passes were for.

**Against the target of §1:** +26.3% with every paired repetition in
[1.2623, 1.2632], 300/300, identical points. `combo5` is the B200 build,
as `make gpu-b200-19b`. Class **engineering**: same walk, same 5.3125
products per update, the logic-pipe floor unchanged, ratio to it bounded
by 1. What the sweep measured is that **the 6000-tuned point is not the
B200's**: 26% of the B200's rate was sitting in knobs the tree already had,
each rejected on the other die for reasons that do not hold on this one.

### 3.4 Where this leaves 30 B/s

| part | best build | B/s | binding pipe | 30 B/s needs | distance |
|---|---|---:|---|---|---|
| RTX PRO 6000 | `gpu-rtx-pro6000-20b` | 19.98 | carry-less unit, floor 22.3 | fewer than 4.1 products per update | **below the floor** |
| B200 | `gpu-b200-19b` | 19.40 | logic pipe, ≤ 955 dynamic ALU slots per update | ≤ 618 slots | a **1.55×** ALU cut, from 1.95× before the sweep |

The B200 has closed to 0.97× the 6000 at 2.06× the list price (3.1 B/s per
dollar-hour against 6.6), so it is still not the campaign's part, and 30
B/s is still nobody's. But the B200's remaining distance is engineering on
an idle unit — the reduction by `CLMAD` against the modulus (−290 slots for
+19 `CLMAD`, THROUGHPUT-20B §4) is the next priced lever and would take the
count to ≈ 665 — where the 6000's is a floor.

## 4. The same sweep on the RTX PRO 6000

The symmetric experiment, so that "the 6000-tuned point is the 6000's" is a
measurement and not an assumption: the same twenty arms around the 20 B/s
build itself (`BASE_KNOBS=`), same procedure, one RTX PRO 6000 (Modal,
2026-09-22; [summary.json](benchmarks/autosweep/rtx-pro-6000/summary.json)).
Base screened **19.961 B/s**.

| arm (20 B/s build + …) | screened / base | | arm | screened / base |
|---|---:|---|---|---:|
| 640 × 1 | **1.005** | | slot unroll 2 | 0.992 |
| fused | 1.002 | | fused, loads one slot ahead | 0.980 |
| pair-ILP off | 1.001 | | nibble pivot | 0.979 |
| flat lo/hi issue | 1.000 | | single products inlined only | 0.968 |
| forward-pass pipeline off | 0.999 | | products out of line | 0.960 |
| chain product not first | 0.998 | | L2 persist off | 0.939 |
| ONB squarings on the ALU | 0.996 | | 768 × 1 | 0.931 |
| 9-word inverse conversion | 0.995 | | 384 × 1 | 0.924 |
| | | | 256 × 2 | 0.908 |
| | | | batch 32, 8 warps | 0.861 |
| | | | batch 32, 16 warps | 0.851 |
| | | | three-limb Karatsuba for `TOP_CLMAD` | **0.696** |

One arm clears the 0.5% line, by 0.02%; the greedy therefore builds one
combination (the same arm) and stops. The four knobs that carried the B200
are the bottom of this table: batch 32 at −14 to −15%, and the three-limb
Karatsuba — four more `CLMAD`s per product on the die where a `CLMAD` costs
38 logic slots — at **−30%**, the largest loss of any arm on either card.
`TABLE_FUSED` is +0.2% here against +13.6% there. The sweep draws the two
dies' pipes as two orderings of the same list.

| finalist | B/s, reps 1 – 5 | median | paired / base per rep | SM clock, MHz | GPU °C | verified | same DPs | class |
|---|---|---:|---|---|---|---|---|---|
| base = `gpu-rtx-pro6000-20b` | 19.992 / 19.981 / 19.951 / 19.935 / 19.933 | **19.951** | 1 | 2422 throughout | 47 → 79 | 300/300 | yes | reference |
| 640 × 1 (96 registers) | 20.272 / 20.190 / 19.993 / 19.799 / 19.755 | 19.993 | 1.014, 1.010, 1.002, **0.993, 0.991** | 2422 → 2355 | 50 → 78 | 300/300 | yes | engineering, did not pay: fails "every repetition above 1.0" |

The 640-thread block is a per-clock gain of about 1.4% that the card takes
back in clock: 20 warps at 96 registers draw more per SM-clock than 16 at
116, and once the card is warm (it had run the sweep for an hour; 76 – 79 °C
in the last two repetitions) the 600 W limit costs the fuller kernel
60 – 70 MHz while the base holds 2422. Same shape as ITERATION-FUNCTION.md
§6.3 on the 4500. **The 20 B/s build stands on the RTX PRO 6000**, now
against an automatic search of the same space that moved the B200 by 26%.

## 5. Classification

| change | class | evidence |
|---|---|---|
| `TABLE_FUSED=1` on the B200 | engineering, +13.7% | paired 1.1375, 5/5; 300/300; identical DPs |
| `combo5` (fused + batch 32 + `KARAT3` + 512 threads) on the B200 | **engineering, +26.3%**, new B200 build | paired 1.2627, 5/5; 300/300; identical DPs |
| the 14 arms below the line | engineering that did not pay on this part | star table; each is a knob the 6000 selected for its own pipes |
| `kernel_cost.py` on the fused kernel | tooling gap, recorded | its loop finder sees one loop where the fused kernel has a step loop containing a slot loop; the fused rows' ALU column is the bound from the measured rate, not a static count |
| 640 × 1 on the RTX PRO 6000 | engineering, did not pay | +1.4% per clock, −0.9% once the card is warm; 2 of 5 paired repetitions below 1.0 |
| the 6000's other 19 arms | engineering that did not pay | star table of §4; the B200's four winners are the 6000's four largest losses |
