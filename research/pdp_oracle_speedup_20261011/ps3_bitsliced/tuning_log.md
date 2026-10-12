# PS-3 tuning log (interactive probes, kept as recorded; not the frozen measurement)

Protocol: [`../protocols/PS-3.md`](../protocols/PS-3.md) allows tuning before the frozen ABAB run and requires the rejected attempts to be kept.
Every figure below is **one** `pdp_speedup time --arm ARM --l 11 --targets 16 --passes 1 --min-secs 3` run pinned with `taskset -c 3` on the shared
container while another session's jobs were running (so **not** `isolated_bench` runs, no A/A, may be contended).  They steered the engineering; the frozen
numbers are in [`raw/times.jsonl`](raw/times.jsonl) and [`../ps3_bitsliced/README.md`](README.md).  The micro-benchmark rows are cycles at a nominal 2.8 GHz
from `pdp_speedup bs-micro --l 11` (33-plane product, dependent chain); the final-source rerun is [`raw/micro.jsonl`](raw/micro.jsonl).

| step | change | `bs64` | `bs256` | `bs512` | `group-scalar` (same session) | verdict |
|:--|:--|--:|--:|--:|--:|:--|
| 0 | first working const-generic build: lanes `u64`/AVX2/AVX-512, closures and iterator adaptors in `batch`/`stage_c`/`montgomery`, output-stationary product | 76.0 | 63.0 | 68.9 ns/pair | 31.9 | baseline of this task |
| 1 | callgrind on `bs256` showed `_mm256_xor_si256` and `Map::fold` executing as **out-of-line calls**: a closure inside a generic body is compiled as its own function without the caller's `#[target_feature]`, so intrinsics do not inline.  Removed every closure / iterator adaptor on lane values from the hot path (named `#[inline(always)]` functions, `Vec` pushes) | 76.1 | **33.3** | **41.2** | - | **kept** (`bs256` -47%, `bs512` -40%; `u64` has no intrinsics, unchanged) |
| 2 | column-by-column product, one register accumulator per output column (`colwise_mul`) | - | micro: 3,700-3,900 cyc (vs 950-1,170 output-stationary) | micro: 3,900-5,400 cyc (vs 2,450) | - | **rejected**: slower at every width; bounds checks and variable inner trip counts dominate (kept in the source as `colwise_mul`, micro mode `mul-columnwise`) |
| 3 | `T x T` register tile with an accumulator **array** (`tiled_mul`, T = 4 / 8) | 89.8 | 31.0 | 29.5 | 31.5 | superseded by 4: the disassembly of the standalone `mul_avx512` showed a load, a `vpternlogq` and a store per gate (the array is not promoted to registers before the loops unroll); standalone product 2,300-2,400 cyc against 1,350-1,550 when inlined into the micro loop |
| 4 | the same tiles with **named** accumulators written out gate by gate (`tile8_gates!`, `tile4_gates!`); `u64` lane back on the output-stationary product (tile 4 measured 1,040-1,240 cyc vs 920 for output-stationary) | 72.9 | **29.5** | **28.8** | 29.8 | **kept** (standalone AVX-512 product 1,380-1,420 cyc, AVX2 980-1,010, `u64` 920) |

Residual gaps identified and **not** closed (so the figures in the README are for this implementation, not a bound):

* **Product schedule.**  1,089 fused gates per product at `n = 33` take about 1,400 cycles at 512 lanes (1.3 cycles/gate): the zmm logic here retires at roughly one
  `vpternlogq` per cycle, not the two the plan's model assumes.  The AVX2 product takes 980-1,010 cycles for 2,178 uops (about 2.2 uops/cycle).
* **Non-product work.**  `bs-micro` mode `batch8` (one batch of 8 tasks) takes about 170k cycles at 512 lanes against 84 products x 1,420 = 119k: about 30% of the batch is
  outside the products (full-width `E` temporaries, zero-fills and copies of 2.5 KB arrays, the plane-wise adds that PS-2's program had dead-code-eliminated).  A lean
  `stage_c` that computes only planes `l..n` of `t` and `t~` and ORs them into the mask would remove most of it; not built.
* **Lane utilisation.**  A task is a (row, block) cell and a row only needs partners `i2 >= i1`, so wide registers waste the lower triangle: 91%, 78%, 65% of lane slots are valid pairs
  at 64, 256, 512 lanes for `l = 11` (`raw/utilisation.jsonl`).  Cost per lane slot is therefore 0.65 of the cost per valid pair at 512 lanes (arithmetic).
