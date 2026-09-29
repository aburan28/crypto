# Is compact-orbit IC's low instruction throughput a cache-miss effect? A cachegrind check at n = 41

Follow-on to `RESEARCH_STRONG_RHO_LADDER_20260929.md` and
`RESEARCH_STRONG_RHO_SWEEP_PROTOCOL_20260929.md` (merged as PR #966). Written and
committed **before any cachegrind, native-timing or interference run of the IC or
rho arms on this host**; the host calibration below was run first because the
model's prices come from it.

## The question, and what it does and does not touch

The ladder note says, of the n = 53 cell: "IC's 1.36 GiB index runs at a lower
IPC, which is consistent with the memory-bound reading in PR #955 but was **not**
tested here (no cache simulation was run), so it stays a hypothesis." PR #966 and
`docs/ic/BOUNDARY_TARGETS.md` carry the same sentence. This note tests it, at
the smaller n = 41 cell where cachegrind is cheap.

**H_mem.** Compact-orbit IC retires instructions more slowly than the strong rho
(R3) at the frozen n = 41 cell *because* its data accesses miss the cache
hierarchy, and the misses account for most of the per-instruction time gap.

What the answer changes:

- It does **not** change the instruction-count verdict. IC/rho = 1.520 in retired
  instructions at this cell (sweep note); that ratio is the scored quantity and
  is hardware-neutral.
- It bears on the **wall-clock** gap (about 3.3× here against 1.52× in
  instructions) and on which IC engineering is worth trying: if H_mem holds, a
  cache-aware IC (prefetch, a smaller index) attacks the extra wall-clock factor
  but not the instruction factor; if it fails, the extra factor is core-bound
  (dependency chains, branch mispredicts) and memory work would not help.

## What I had seen before writing this (disclosure)

From the merged sweep at this cell (a different host, Xeon 2.1 GHz): retired
instructions IC 12,441,253,531 and rho R3 8,184,739,666; native wall IC 4.28 s
(index build 2.92 s, rank stage 0.25 s, targets 1.06 s) and rho 1.30 s;
peak RSS IC 338 MB, rho 41 MB; 2,664,777 root-table entries. So I already knew
IC ran at roughly 2.9 G instructions/s and R3 at roughly 6.3 G/s there.
I had **not** seen any cache-simulation output, any per-function miss count, any
interference behaviour, or any native timing on the host used here. The host
calibration below was the only thing measured before this note was committed, and
it touches neither arm.

One design change came from the calibration and is disclosed: I first intended
the two simulated last-level sizes to follow the host's nominal 1 MiB L2 and
33 MiB L3. The calibration showed that on this VM a random-access working set
above about 2-4 MiB already costs 90-130 ns (see the table), i.e. the nominal
33 MiB L3 gives essentially no random-access benefit. The simulated sizes were
therefore chosen to bracket the measured on-chip capacity (2 MiB and 8 MiB)
instead. This was decided before any arm was simulated.

## Frozen cell and instruments

- **Cell.** n = 41, a = 0, L = 1,024 targets, IC `K` = 255, corpus
  `n41-strong-sweep-L1024-v1`, batch seed 531310, rho `KIC_RHO_RUNG=3
  KIC_RHO_LANES=32 KIC_RHO_DP_BITS=4`, IC `construct:41:0:255 <scalars> 7 <out>`
  with the scalars file from the merged sweep cell
  (`cell_n41_L1024_K255/scalars.txt`, sha256 in `SHA256SUMS_inputs`). Identical
  to the cell in the sweep note; no re-tuning.
- **Binaries.** Rebuilt from the merged tree (sources unchanged from the sweep:
  `koblitz_rho_batch_ks_strong.rs` sha256 `3a4eb7b4…`, `koblitz_orbit_dlp_fast.rs`
  sha256 `8dea682a…`; binary hashes in `SHA256SUMS_inputs`, not byte-identical to
  the sweep's because the build path and host differ).
- **Host.** `host_manifest.txt`: 4 vCPU, Xeon @ 2.80 GHz nominal (effective clock
  below), `pclmulqdq`, 4 KiB pages (THP `madvise`, and nothing here calls
  `madvise`), all arms single-threaded.
- **Simulator.** `valgrind --tool=cachegrind --cache-sim=yes --branch-sim=yes`,
  `--I1=32768,8,64 --D1=32768,8,64`, `--LL=<size>,16,64` with `size` ∈ {2 MiB,
  8 MiB, 32 MiB}, both arms, whole process. The 32 MiB run (the nominal host L3)
  is exploratory and enters no decision.
- **Native timing.** Uninstrumented, pinned to core 0, five interleaved
  repetitions of (rho, IC); user + system CPU time of the child. Wall time is
  recorded alongside.
- **Interference test.** Three repetitions of each arm alone on core 0, and with
  three antagonists on cores 1-3, each running 8 independent random pointer chains
  over 256 MiB (`mem_calib antagonist`). It measures how much each arm slows when
  the shared last-level cache and memory system are contended.
- **Host calibration** (`mem_calib.c`, `calibrate.py`, three runs,
  `calibration_run{1,2,3}.json`, reduced by `summarize_calibration.py`): effective
  clock from a dependent add chain, random pointer-chase load-to-use latency by
  working-set size, memory-level parallelism from k independent chains over
  256 MiB, and the cost of a mispredicted branch.

### Calibration (measured; prices below come from it)

Median load-to-use latency in ns, min-max over the three runs:

| working set | 16 KiB | 512 KiB | 1 MiB | 2 MiB | 4 MiB | 8 MiB | 16 MiB | 64 MiB | 256 MiB | 1 GiB |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| ns | 1.24-1.34 | 6.3-6.6 | 15.8-18.4 | 24.3-27.1 | 92.6-128 | 139.6-146.2 | 153.4-156.2 | 167.8-187.5 | 230.3-253.9 | 378.8-400.6 |

Effective clock 3.21 GHz; 4-cycle L1 load-to-use (1.28 ns) as expected; peak
memory-level parallelism 11.05 (k = 10-16 chains); 7.49 ns per mispredicted
branch (about 24 cycles). The random-access curve leaves the on-chip tier between
2 and 4 MiB and is flat at 140-190 ns from 8 MiB to 64 MiB, then rises with page
walks. Working sets of 6 MiB and above pay TLB misses (4 KiB pages) that
`cachegrind` does not model; the prices are measured with the same page size the
arms use, so the page-walk cost rides in the price of a far miss.

Registered prices (nanoseconds over an L1 hit, `calibration_summary.json`):

| symbol | meaning | value |
|:--|:--|--:|
| `P_mid_low` | D1 miss served on chip, cheap end (min latency at 512 KiB − L1) | 4.98 |
| `P_mid_high` | D1 miss served on chip, expensive end (max latency at 2 MiB − L1) | 25.81 |
| `P_far_low` | LL miss, cheap end (min latency at 8 MiB − L1) | 138.36 |
| `P_far_high` | LL miss, expensive end (max latency at 256 MiB − L1) | 252.62 |
| `MLP` | memory-level parallelism, largest observed | 11.05 |
| `P_br` | ns per mispredicted branch (median) | 7.49 |

## Quantities and registered rules (fixed now)

Let `Ir_a`, `cpu_a` be an arm's retired instructions (cachegrind) and its median
native CPU time; `cpi_a = cpu_a / Ir_a` (ns per instruction). The throughput gap to
explain is `Δ = cpi_IC − cpi_rho`. For each arm, from cachegrind: `mid` = D1 misses
that hit the simulated LL, `far` = LL data misses (`DLmr + DLmw`).

- **Generous bound (F_high).** LL = 2 MiB; every miss paid at the expensive price
  with no overlap: `stall_hi = (mid·P_mid_high + far·P_far_high) / Ir`.
- **Stingy bound (F_low).** LL = 8 MiB; cheap prices and full overlap:
  `stall_lo = (mid·P_mid_low + far·P_far_low) / MLP / Ir`.
- `F_high = (stall_hi,IC − stall_hi,rho) / Δ` and `F_low = (stall_lo,IC − stall_lo,rho) / Δ`,
  the fraction of the per-instruction gap that modelled memory stalls can explain,
  at the two extremes of the model. (`cpi` and the stalls are in ns, so the ratio
  does not depend on the clock.) The simulator has no hardware prefetcher and no
  TLB, so `F_low` is not a strict lower bound on real stalls for prefetch-friendly
  streams; the interference test below is the native check on that blind spot.
- **Interference sensitivity.** `sens_a = median(CPU loaded) / median(CPU alone) − 1`.
  `S = sensitive` iff the alone-run spread (max/min − 1) is ≤ 0.05 for both arms,
  `sens_IC ≥ 0.10` and `sens_IC − sens_rho ≥ 0.08`; `S = insensitive` iff the spread
  test passes and `sens_IC < 0.05`; otherwise `ambiguous`.
- **Exploratory only:** branch-mispredict share `(Bcm_IC/Ir_IC − Bcm_rho/Ir_rho)·P_br / Δ`
  (cachegrind's predictor is a simple one and over-counts on modern cores); per-function
  miss attribution from `cg_annotate`.

### Gates (blocking unless noted)

- **G1** every native, cachegrind and interference run recovers/solves every one of
  the 1,024 targets for both arms, and the rho fixtures' planted scalars equal the
  scalars file the IC arm reads.
- **G2 (not blocking)** cachegrind `Ir` within 0.1 % of the sweep's callgrind `Ir`
  (IC 12,441,253,531; rho 8,184,739,666). If it fails the verdict is stated for
  the rebuilt binaries only and the sweep ratio is not cited for them.
- **G3** simulator sanity on controls (`mem_calib` under the same cachegrind
  configuration): for random pointer chains over 512 KiB (k = 1), 8 MiB (k = 1),
  256 MiB (k = 1) and 256 MiB (k = 10), the *native* extra time over the same
  loop at 16 KiB must lie in `[0.7·lower, 1.3·upper]`, where `lower`/`upper` are
  the stingy/generous model applied to the loop's own cachegrind counts (setup
  removed by differencing a long and a 1,000-step run). This is a deliberately
  loose check that the simulator's miss counting and level attribution work; it does
  not validate the prices. A sequential 256 MiB scan is also reported as the
  known prefetch blind spot, expected to run *below* the lower band, and enters no gate.
- **G4** LL data misses are non-increasing in LL size (2 → 8 → 32 MiB) for each
  arm, with 1 % slack.
- **G5** native CPU-time spread (max/min − 1) over the five repetitions ≤ 0.10 for
  each arm.

### Decision rule (registered)

If a blocking gate fails: **UNDETERMINED**. Else if `Δ ≤ 0.05·cpi_rho`: **no
throughput gap to explain**. Else:

- **H_mem SUPPORTED** iff `F_low ≥ 0.5`, *or* (`F_high ≥ 0.5` *and* `S = sensitive`).
- **H_mem REFUTED as the main explanation** iff `F_high < 0.25`, *or* (`S =
  insensitive` *and* `F_high < 0.5`).
- otherwise **UNDETERMINED**, with the interval `[F_low, F_high]` and `S` reported.

The two conditions cannot both hold (support needs `F_high ≥ 0.5`, refutation needs
`F_high < 0.5`).

**What each outcome licenses.** SUPPORTED means memory misses are a first-order
part of IC's throughput deficit at this cell on this host (by two independent
lines: a simulated miss budget and a native contention test); it does not show
that IC's wall time is memory-bound at other sizes, on other hardware, or at
n = 53, and says nothing about the instruction ratio. REFUTED means the extra
wall-clock factor is core-bound and memory-oriented IC engineering would not close
it at this cell. UNDETERMINED means the simulator and the contention test do not
jointly resolve it, and the intervals are all that is claimed.

**Prior (recorded now, not a finding).** I expect the index build (about 70 % of
IC's wall at this cell) to carry most of the far misses and the answer to land
in UNDETERMINED or SUPPORTED-via-contention rather than an outright `F_low ≥ 0.5`,
because the simulator has no overlap information and the stingy bound divides by
an MLP of 11.

**Inadmissible, stated in advance:** changing a price, a threshold or a simulated
size after seeing arm counts; adding, dropping or re-ordering repetitions after
seeing results; using wall time in place of CPU time; reporting `F_high` or
`F_low` alone; dropping a control that fails G3; treating a cachegrind
branch-mispredict count as a measurement of hardware mispredicts.

## Scope and honesty limits

- One cell (n = 41, L = 1,024, K = 255), one VM class, one implementation of each
  arm, IC as merged. Nothing here reaches n = 53 (where IC's index is 1.36 GiB) or
  larger, other L, aarch64, or m = 83 / ECC2K-130 transfer.
- Cachegrind simulates a two-level inclusive-style LRU cache with no prefetcher,
  no TLB, no out-of-order overlap and one thread; the prices and the MLP bound stand
  in for what it omits, and the contention test is the native cross-check. Neither
  is a hardware counter (`perf` is not available in this VM).
- The prices are measured on this VM with random pointer chases; a workload whose
  misses are partly sequential or partly overlapped in ways the chains do not
  reproduce is priced less accurately, which is exactly why the model is reported as
  an interval.
- Class (`AGENTS.md` §3): **accounting/measurement** on an explanation; it moves
  no boundary target and does not touch the instruction-count verdict.

---

<!-- Results are appended below by the commit that follows the runs. Nothing above this line is edited after seeing them. -->

## Amendment A1 (additive; written after the four registered stages ran, before any follow-up run)

**What I had seen when writing this.** All four registered stages had completed and
`analyze.py` had been run on them; their raw outputs and `analysis.txt` are committed
in the same commit as this amendment, so the order is checkable. The registered
outcome is:

- Gates: G1, G3, G4 pass; G2 passes (cachegrind `Ir` IC 12,442,73x,xxx and rho
  8,185,26x,xxx, within 0.013 % of the sweep's callgrind counts); **G5 fails** (native
  CPU-time spread over the five repetitions: IC 0.230, rho 0.122, against the
  registered ≤ 0.10; one IC repetition ran 18 % above the median on this shared VM).
- Model: `F_high = 1.91` (above 1, i.e. the no-overlap bound over-explains the whole
  gap and is vacuous as a bound), `F_low = 0.077`. Interference: IC slowed by −1.5 %
  and rho by +4.6 % with the antagonists, but the alone-run spread was 0.11 for both
  arms (> 0.05), so `S = ambiguous`.
- **Registered verdict: UNDETERMINED (a registered gate failed).** Nothing in this
  amendment changes it, and no result below can.

I had also looked at the raw simulated counts (IC 12.15 D1 misses per 1,000
instructions and 1.21 LL misses per 1,000 at 2 MiB; rho 6.63 and 0.147) and at the
`cg_annotate` per-function tables before writing this.

**Why an amendment.** Three weaknesses of the registered instrument are visible now,
and none of them is about the direction of the answer:

1. A five-repetition max/min noise gate on a shared VM is failed by one slow
   repetition; the medians may still be stable, but the registered stage cannot say.
2. The interference test has **no positive control**: "IC did not slow down" cannot be
   read unless the same antagonists demonstrably slow a memory-bound loop.
3. Neither the simulator nor the contention test isolates **address translation**:
   IC's 338 MB index on 4 KiB pages far exceeds the TLB's reach, and cachegrind has no
   TLB. Transparent huge pages (this VM: THP `madvise`; glibc 2.39 accepts
   `GLIBC_TUNABLES=glibc.malloc.hugetlb=1`) remove that term natively without
   touching the algorithm or its instruction count.

**Unregistered follow-ups (all exploratory; run once each unless a run dies on
infrastructure, in which case rerun that stage once and report both).**

- **E3 replicate.** Ten more interleaved repetitions of (rho, IC), core 0. Pooled
  median cpu-time over all fifteen repetitions gives `cpi_pooled`; spread is reported as
  max/min − 1 and as (Q3 − Q1)/median. `F_high` and `F_low` are recomputed with
  `cpi_pooled` and the same formulas and labelled post hoc.
- **E1 interference controls.** A random pointer chase over 256 MiB (positive control)
  and over 16 KiB (negative control), each alone on core 0 and with the same three
  antagonists, three interleaved repetitions. The registered `S` is *informative*
  only if the positive control slows by ≥ 0.10 and the negative control by ≤ 0.03; if
  not, `S` is uninformative regardless of IC's value.
- **E2 huge pages.** Seven interleaved repetitions of {rho, IC} × {4 KiB pages,
  `glibc.malloc.hugetlb=1`}, CPU time, core 0, sampling `AnonHugePages` from
  `/proc/<pid>/smaps_rollup`. *Valid* only if IC's sampled `AnonHugePages` reaches
  at least half of its peak RSS **and** a 256 MiB pointer-chase control runs at least
  1.3× faster with the tunable than without. With `d_a = (cpu_a,4K − cpu_a,THP) / Ir_a`
  (ns per instruction, medians), the **translation share** of the gap is
  `(d_IC − d_rho) / Δ`, with `Δ` from the 4 KiB medians of E3's pooled set.

**Reading rules for E2 (fixed now).** If E2 is valid and the translation share is
≥ 0.25, address-translation stalls that huge pages remove account for at least a
quarter of the per-instruction gap, natively measured — a lower bound on the memory
share, since huge pages leave the cache misses in place. If valid and the share is
< 0.05, translation stalls are not the mechanism (cache-miss latency remains
unresolved and stays an interval `[F_low, F_high]`). Anything else: unresolved.
If E2 is invalid, it is reported as not informative.

**Inadmissible, stated in advance:** using E1-E3 to change the registered verdict or
any threshold above; adding stages, repetitions or configurations to the follow-ups
after seeing their results; reporting only the direction that favours or disfavours
H_mem.

---

<!-- A1 results are appended below. Nothing above this line is edited after seeing them. -->
