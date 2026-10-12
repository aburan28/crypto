# A large finite construction gain comes from packed direct multiplication

**Packed direct construction passes the predeclared 2x gate on both seed sets.**
On batch-64 changing-coefficient and degree-cycle workloads at n=6/8/10/12,
paired median ratios against the fastest of the five legacy arms are
**3.48–4.53x** in the primary run and **3.55–4.90x** in confirmation. The lowest
95% bootstrap lower bounds are **3.41x** and **3.47x**, respectively. Every one
of the eight declared comparisons passes in both runs.

This is an engineering gain in complete cold **matrix construction**, including
its own setup and output validation. It does not establish faster full solving,
index calculus, an asymptotic improvement or a cryptanalytic breakthrough. The
broader goal of a dramatic improvement to the complete algorithm remains open.

## What changed and what the controls establish

Both new arms build packed rows directly in fixed ambient coordinates, accumulate
actual support in an occupancy bitmap, and compact columns once. The packed
envelope also avoids scanning ineligible multipliers by selecting precomputed
degree-eligible views in numeric multiplier order.

The stronger packed-direct control retains only ambient coordinates and
multiplier lists. It evaluates the current products by XOR. It has no generator
envelopes, product-routing schedules or cached completed matrix. This arm is
faster than the packed envelope in all 16 attribution cells in both runs.
The experiment supports the combined change in ambient context reuse, packed
multiplication and compaction; it does not isolate the contribution of each
component individually, and does not credit the gain to symbolic reuse.

The packed envelope passes 45/48 incremental timing gates in each run, but
fails the universal 2x gate at n=12 on coefficient-changing inputs and fails
all 16 attribution gates against packed direct. Thus the frozen result's
`performance_gate: REJECTED` refers to **packed-envelope** promotion. The
separate `dramatic_gate.packed_direct: PASS` records the successful candidate.

## Cold timing table

Every value below is copied from the corresponding frozen `RESULT.md`. These
are milliseconds per cold batch of 64, pooled medians across four sizes, two
holdout seeds and fourteen balanced repetitions. Ratios in this table are
descriptive ratios of pooled medians against the unmodified direct constructor.
Acceptance instead uses paired per-size ratios against the **pointwise minimum
of all five legacy arms**; no favorable baseline is selected after measurement.
All rows have oracle equality PASS and the engineering classification.

| Variant | Primary coefficients (ms) | Primary degree cycle (ms) | Confirmation coefficients (ms) | Confirmation degree cycle (ms) | Primary direct / arm, coefficients | Confirmation direct / arm, coefficients |
|---|---:|---:|---:|---:|---:|---:|
| Direct construction | 1.448730 | 1.609520 | 1.571417 | 1.639750 | 1.000 | 1.000 |
| Verified layout reuse | 1.540104 | 1.661208 | 1.809438 | 1.637417 | 0.941 | 0.868 |
| Exact product schedule | 1.438396 | 1.623917 | 1.759542 | 1.614562 | 1.007 | 0.893 |
| Exact packed-matrix cache | 1.453021 | 1.601625 | 1.450812 | 1.613021 | 0.997 | 1.083 |
| Original support envelope | 1.836625 | 1.982729 | 1.847271 | 1.989396 | 0.789 | 0.851 |
| Packed support envelope | 0.586687 | 0.647167 | 0.775875 | 0.613354 | 2.469 | 2.025 |
| **Packed direct** | **0.319625** | **0.371125** | **0.332417** | **0.343000** | **4.533** | **4.727** |

The per-size cold paired ratios for packed direct, with 95% bootstrap intervals:

| Variables | Family | Primary ratio [interval] | Confirmation ratio [interval] |
|---:|---|---:|---:|
| 6 | coefficients | 4.444 [4.286, 4.573] | 4.538 [4.354, 4.611] |
| 6 | degree cycle | 4.038 [3.994, 4.078] | 4.283 [4.187, 4.340] |
| 8 | coefficients | 4.528 [4.400, 4.629] | 4.905 [4.700, 5.771] |
| 8 | degree cycle | 4.442 [4.341, 4.514] | 4.881 [4.558, 5.165] |
| 10 | coefficients | 4.250 [4.190, 4.332] | 3.879 [3.818, 3.931] |
| 10 | degree cycle | 4.455 [4.408, 4.536] | 4.323 [4.149, 4.399] |
| 12 | coefficients | 3.483 [3.407, 3.548] | 3.551 [3.469, 3.625] |
| 12 | degree cycle | 4.135 [4.060, 4.322] | 4.136 [3.994, 4.249] |

The exact matrix cache still wins on repeated identical inputs: primary pooled
median **0.165104 ms**, versus **0.423021 ms** for packed direct. Confirmation
is **0.167021 ms** versus **0.401292 ms**. Repeated-input and support-escape
guardrails remain in the raw data; neither is silently promoted into a win over
every possible workload.

## Cost, correctness and provenance

Packed direct retains a median **2,386 bytes** over these holdouts. The primary
original envelope retains **100,929 bytes** and the packed envelope **96,881
bytes**. These capacity-based numbers exclude allocator metadata. Whole-worker
peak RSS spans 2,293,760–4,915,200 bytes in the primary run and
2,342,912–4,964,352 bytes in confirmation; it includes all arms and references,
so it is not candidate-specific peak memory.

Each run has 256 cells, 25,088 batch-arm samples and 533,120 oracle-verified
matrix outputs. Together they contain **1,066,240 checked outputs**, without
failed or censored cells. Each envelope records 49,728 changed-input hits and
4,704 support-escape fallbacks per run. Packed direct accepts those escapes
because it needs only the unchanged ambient context, and independently matches
the oracle on them. Counts include measured repetitions.

Nine Rust tests include 8,192 exhaustive coefficient/degree/active-mask cases for
both packed arms, cross-word compaction, preserved row order, identity compaction,
current output caps, and an ambient basis wider than the allowed actual output.
Thirteen evidence tests replay both frozen runs, verify manifests and unchanged
worker source, check disjoint seed sets and unchanged gates, and reject corrupted
coefficient/context reuse claims. Each run's manifest covers 15 files.

`run_01` retains the first protocol and runner. `run_02` freezes the separate
confirmation protocol and the runner's added `--protocol` argument. Their worker
bytes are identical. No algorithm or acceptance threshold was tuned between
runs. Confirmation uses the same host and is not unaffiliated reproduction.
The bootstrap intervals describe repeated measurements on these fixed seeds.

## Limits and next unresolved question

The coordinate lookup has `2^n` entries and the ambient basis includes every
squarefree monomial of degree at most D. Here n is only 6–12 and D=3. These data
cannot justify extrapolating the lookup representation to large sparse systems.
A scalable coordinate lookup and a broader generic algebra workload need their
own frozen comparison. Nor do construction timings establish the fraction of
full solving time available to improve; that end-to-end effect remains unmeasured.

The prior exact-support and envelope studies remain immutable. The canonical
scoreboard now adds these construction rows without changing its full-method
verdict. No production solver integration, curve inputs or key-related work is
part of this experiment. Full-pipeline costs and rho ratios remain null.
