# Where the table-v3 walk's 4x goes: false-hint probe latency under the block barrier

Question: the table walk measured 20.078 B/s under cycle rule v2
([ONE-BLOCK-GEOMETRY.md](ONE-BLOCK-GEOMETRY.md)) and 5.10 B/s under the exact
rule v3 with every scheduling repair applied
([CYCLE-ESCAPE-V3.md](CYCLE-ESCAPE-V3.md), `benchmarks/block-both2-confirm5`).
The v3 probe adds about 60 ALU lane-instructions per update on average, under
5% of the ALU budget in [ROOFLINE.md](ROOFLINE.md) §2. Why does it cost a
factor of four?

Answer: **the probe's cost is latency on the block's critical path, not
throughput, and 93% of it is spent on hints that are not cycles.** Every
number below is read from a receipt already in the tree; the one inference is
marked.

## 1. What the receipt says

[`benchmarks/cycle-probe-profile/counts.tsv`](benchmarks/cycle-probe-profile/counts.tsv),
long prefix, 3,154,116,608 updates on the B16 fast2 kernel:

| counter | value | share of hints |
|---|---:|---:|
| hints (`eccTagFruitless` fired) | 7,834,496 | 0.248% of updates |
| fast2 exact two-cycles | 2,876,286 | 36.7% |
| exact 4-cycles / 6-cycles | 2,282 / 14 | 0.03% |
| **open after eight steps (no cycle)** | **4,955,914** | **63.3%** |
| affine `next()` calls | 42,532,810 | 5.43 per hint |

The open-eight probes account for 4,955,914 × 8 = 39,647,312 of the 42,532,810
`next()` calls: **93.2% of all probe arithmetic resolves a false hint**, whose
outcome is "keep the raw step" — the step the walk was going to take anyway.

The hint predicate (`include/tablewalk.h:106-114`) fires when the new tag
negates any of the last four tags, or closes a τ-relation with the last three.
Only "negates the most recent tag" is a two-cycle; a match against tags 2–4 is
an `h`/`k`/`ε` coincidence with probability about `3/(2·8·131) ≈ 0.14%` per
step on the ECC2K-130 table, which is the false-hint rate the receipt shows.

## 2. Why 63% false hints costs 4x and not 5%

Per probe step (`include/cycleanchor_body.h:10-66`, `TwCycleOps` in
`include/packedtablewalk.cuh:362-375`), `next()` performs a full affine
addition: two `fromPolynomial131`, `add131`, **`inv131`** (eight `mul131`
links, each `mul131` two `toPolynomial131` plus the product plus
`fromPolynomialProduct131`, out of line via `ECC_BIG`), two more `mul131` and
one `sqr131`. An open-eight probe is eight of those in one serial dependency
chain, inside a `__noinline__` call that holds a 400-byte local frame
(`Point points[8]`, `tags[8]`).

The block schedule (`include/packedkernels.cuh:1189-1232`, the RTX table
preset's `TABLE_BLOCK_HINTS=1`) collects every hinted slot of the 512-thread
block into a shared queue, then resolves queue entry `q` on thread `q`. With
about 20 hints per block-step (512 threads × 16 slots × 0.248%), of which
about 13 are open-eight, the resolvers land in **one warp**, and the other
fifteen warps wait at the `__syncthreads` on line 1232. The probability that a
block-step has no open-eight hint is `e^{-13}`: the stall happens on every
step.

So the per-step cost is not the ALU work (≈ 0.0135 × 4,500 ≈ 60 lane-ops per
update, spread across the SM) but the **wall-clock length of the longest
probe in the block**: eight serial affine steps whose inversions are
latency-bound (ONE-BLOCK-GEOMETRY §7 measured ~30k cycles for the walk's own
batched inversion chain). Inference from the rates: a block-step is ≈ 191k
cycles at 20.08 B/s and ≈ 750k at 5.10 B/s, so ≈ 560k cycles per step are the
serialized probe, ≈ 70k per `next()`. The per-lane fallback without the
barrier (3.56 B/s, `benchmarks/batch-hints-fast2`) is slower still because
72% of warps then each carry an open-eight probe of their own, and the
serialized resolver kernel of `benchmarks/global-hints` (0.80×) shows the
same thing from the other side: moving the probe does not help while it stays
serial and on the critical path.

Everything else the repairs touched (split-forward reloads, two extra
barriers, 128 registers and the frame, de-pipelined selection) is a few
percent each and is hidden under the probe.

## 3. What would recover it

Ranked by expected gain; none is measured here.

1. **Deferred detection (rule v4).** On a hint, take the raw step and set a
   pending mark (hist bits 13–15 are free, `tablewalk.h:55-66`), keeping the
   start `x` in a cold per-slot field. For the next eight ordinary steps,
   which the batched kernel computes at full rate anyway, compare the walk's
   own new `x` against it (five-word compare on pending slots only). Only an
   exact return triggers anchor selection, from a per-slot ring of the last
   ≤ 8 `(weight, tag)` written while pending; the escape edge is taken at the
   next visit of the anchor. This removes all probe arithmetic for false
   hints and the barrier-serialized critical path with it; exact cycles need
   no affine arithmetic except `twCanonical` on weight ties. The seed to
   endpoint map should be unchanged (same exit vertex, same modified tag; the
   history after the exit differs, but history only delays discovery of a
   later cycle, by the v3 argument), while `iters` grows by at most one lap
   per anchor exit, so it is a **new rule identity**: `ref.h`, the v3 corpus
   header, `tablev3replay`, `test-cycle-escape` and `walkconstant` all move,
   and a finite-trail-guard edge case changes which trails are truncated.
   Estimate: most of the gap, 15–19 B/s with item 4.
2. **Exact v3, cheaper `next()`.** Keep probe points in the polynomial basis
   and use `invPoly131` / `mulPolynomial131` / `squarePolynomial131`
   (`include/packed131.h:1199-1227`, already built with `PACKED_INV_POLY=2`),
   converting `x` once per probe step for weight, phase and pivot; in
   `TwCycleOps::tag` the `toPolynomial131(p.y)` is redundant because `ε` is a
   bit of the ONB `y`. Bit-identical values, no map change. Estimate: 30–45%
   of the probe critical path, i.e. a modest rate gain while the structure of
   §2 remains.
3. **Exact v3 early termination.** Track the formal sum of the probe's table
   addends in `⊕_h Z[τ]/(τ²+τ+2)` and abort when the remaining steps cannot
   cancel it; equivalent to the exact return test unless a non-trivial
   `Z[τ]`-relation among table points exists (a solved DLP instance). Prunes
   only the last two or three steps of false hints (estimate −35% probe
   work), and must be mirrored in `ref.h`.
4. **Probe off the critical path, exact v3.** Stall a hinted slot for one
   step (forward pass fed `d = 1`, `e = 0`; reverse pass predicated no-write;
   `startIter` keeps `iters` exact) and resolve on a concurrent stream through
   the existing `selectGlobalHints` / `resolveGlobalHints`
   (`packedkernels.cuh:639-677`), polling a per-slot flag at the next step.
   Throughput cost ≈ 4.5% ALU, latency hidden. Medium risk: intra-launch
   memory ordering and a bounded-stall fallback.
5. **After 1 or 4:** revert to the 20 B/s forward structure (selection fused
   into the product loop under `TABLE_PIPE_SELECT`, one `x`/`y` load and one
   barrier per update), recovering the hidden 3–8%.

Compaction tweaks (sorting likely two-cycle hints into their own warps) do
not help: the open-eight lanes set the warp's length.

## 4. Boundary

Rates quoted are `--bench` rates of the table walk, which is not the campaign
default; the σ walk's 15.44 B/s ([RTX-PRO6000.md](RTX-PRO6000.md)) is
unaffected by any of this. Two rate figures exist for the exact-v3 block-queue
kernel, 5.10 B/s (1024-step launches, `block-both2-confirm5`) and 7.02 B/s
(95-step launches, `global-hints`); the difference is not resolved here and
does not change the §2 mechanism. Nothing in this note is a 26 B/s result or a
full-solve speedup.
