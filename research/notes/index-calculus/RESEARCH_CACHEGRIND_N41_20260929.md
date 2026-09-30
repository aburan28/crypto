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

## A1 results, as run (added before Amendment A2; raw outputs `explore_*.json`, `work/`)

- **E3 replicate.** Ten more repetitions: CPU-time spread 0.256 (rho) and 0.262 (IC);
  this VM's run-to-run noise is about 25 %, so the registered five-repetition ≤ 0.10
  gate was unreachable here, not merely missed once.
- **E1 interference controls.** The positive control (pointer chase over 256 MiB) did
  **not** slow with the antagonists: median 244.2 ns/step alone, 220.1 ns loaded
  (−9.9 %); the negative control (16 KiB) 1.29 vs 1.28 ns. By the A1 rule (positive
  slowdown ≥ 0.10 required) the interference test is **uninformative on this VM**:
  IC's −1.5 % slowdown says nothing about its memory sensitivity.
- **E2 huge pages, CPU time (user + system).** Validity criteria met: IC had a
  sampled 272-276 MB in `AnonHugePages` (peak RSS 338 MB), and the 256 MiB chase ran
  1.33-1.53× faster with the tunable (221.8/226.4/219.8 → 156.4/170.3/143.4 ns).
  Yet **both arms got slower with huge pages**: median CPU 4.29 → 5.59 s for IC
  (ratio 0.77) and 1.19 → 1.38 s for rho (0.86). Read literally, the A1 rule ("valid
  and share < 0.05 ⇒ translation stalls are not the mechanism") would fire on a
  negative share. **I do not accept that reading**: the measured quantity is user +
  system CPU time, and first-touching 2 MiB pages costs kernel time (zeroing,
  compaction, and, in a VM, the host's faults) that A1 did not separate from user
  time. The rule was written without that confound in view. The literal outcome is
  recorded here, and the substantive answer is that E2 as run cannot separate the
  two effects.

## Amendment A2 (additive; after E1-E3, before the run it describes)

**What A2 changes.** Only the E2 measurement: repeat the huge-page test recording
`ru_utime` and `ru_stime` separately, and score on **user time**, the quantity in which
address-translation stalls appear (a page walk stalls the user thread; a fault handler
is kernel time). Everything else in A1 stands, including its validity criteria (the
chase control above already met the 1.3× criterion; `AnonHugePages` is sampled again).

**Rule (fixed now).** Seven interleaved repetitions of {rho, IC} × {4 KiB,
`glibc.malloc.hugetlb=1`}, core 0. Medians of user time: `d_a = (user_a,4K −
user_a,THP) / Ir_a` (ns per instruction; `Ir_a` from the registered cachegrind runs),
`Δ_user = user_IC,4K / Ir_IC − user_rho,4K / Ir_rho`, and the **user-time translation
share** `(d_IC − d_rho) / Δ_user`. Reading, as in A1: valid and share ≥ 0.25 ⇒
translation stalls removed by huge pages are at least a quarter of the per-instruction
user-time gap (a lower bound on the memory share); valid and share < 0.05 *and* the IC
system-time increase smaller than 5 % of its 4 KiB user time ⇒ translation is not the
mechanism (the last clause is what A1 lacked: if system time balloons, a small share is
uninformative); anything else, including a small share with a large system-time
increase, is **unresolved**. As before, no result here can change the registered
verdict.

---

<!-- A2 results are appended below. Nothing above this line is edited after seeing them. -->

## A2 results, as run (`explore_thp_a2.json`, `analysis_explore.txt`)

Valid: `AnonHugePages` reached 281 MB for IC (peak RSS 338 MB) and the 256 MiB chase
control ran 1.33-1.53× faster with the tunable (A1 run). Medians of seven repetitions:

| | user s, 4 KiB | user s, THP | sys s, 4 KiB | sys s, THP |
|:--|--:|--:|--:|--:|
| IC | 4.011 | 3.424 | 0.317 | 1.976 |
| rho R3 | 1.117 | 1.122 | 0.036 | 0.215 |

Huge pages cut IC's **user** time by 14.6 % (ratio 1.172) and left rho's unchanged
(0.995); both arms pay more **system** time in THP mode (page-fault path), which is
what made A1's user + system sums look like a slowdown. `d_IC = 0.0472` and
`d_rho = −0.0006` ns per instruction against `Δ_user = 0.1860`, so the **user-time
translation share is 0.257**, and the A2 rule fires its "≥ 0.25" branch.

Assessment of that label: it is a point estimate 0.007 above a threshold, from seven
repetitions whose user-time spread within a condition is 9-18 %, and no interval was
registered in A2. It is consistent with a translation share of about a quarter of the
user-time gap; it does not resolve it against 0.20 or 0.30. The other results of the
follow-ups: pooled over fifteen native repetitions the medians are IC 4.371 s and rho
1.237 s (IPC 0.89 and 2.06 at the effective clock), spread max/min 0.28 and 0.28 but
(Q3 − Q1)/median only 0.070 and 0.062, and the post hoc `F_high` / `F_low` recomputed
with the pooled medians are 1.911 / 0.077 (registered: 1.912 / 0.077). The registered
verdict stays **UNDETERMINED (G5)**.

---

## Registration R2 (a second registered attempt; written after A2, before its run)

**Why a second attempt.** Attempt 1 ended UNDETERMINED because (i) its noise gate
(max/min ≤ 0.10 over five repetitions) is unreachable on this VM (measured spread 0.28
over fifteen), (ii) its contention test has no power here (positive control did not
slow), and (iii) the one native mechanism test that worked, A2's huge-page user-time
share, landed at 0.257 with no interval. R2 repeats *only* the part that can be measured
natively and gives it a registered interval and gates. Carried over unchanged: the
frozen cell, the binaries (`SHA256SUMS_inputs`), the cachegrind counts and prices
(deterministic instruction and miss counts are not re-measured), and the simulator
interval `[F_low, F_high] = [0.077, 1.91]`. Dropped, and why: the antagonist
contention test (E1 shows it cannot detect memory contention on this VM).

**Design.** Fifteen blocks. Each block runs the four conditions {rho, IC} × {4 KiB,
`glibc.malloc.hugetlb=1`} once each, pinned to core 0, in an order drawn per block from a
fixed seed (20260930) to break any position effect. Recorded per run: user time, system
time, wall time, sampled `AnonHugePages`, correctness (G1 as before). Before the first
block and after the last, three alternating (4 KiB, THP) pairs of the 256 MiB
pointer-chase control.

**Estimator.** For each arm and page mode take the median user time over blocks; `d_a`,
`Δ_user` and the **translation share** `(d_IC − d_rho) / Δ_user` as in A2
(`Ir_a` from the registered cachegrind runs). Interval: percentile bootstrap over blocks
(10,000 resamples of the fifteen blocks with replacement, `random.Random(20260930)`),
2.5 % to 97.5 %.

**Gates (all must pass, else UNRESOLVED (gate)).**
- **R-G1** every run recovers every target for both arms.
- **R-V** validity: the median over THP-mode IC runs of the maximum sampled
  `AnonHugePages` is at least half of IC's peak RSS, and the median of the six chase
  control speedups is at least 1.3×.
- **R-N** noise: (Q3 − Q1)/median of the 4 KiB user times is ≤ 0.10 for each arm. This
  replaces G5's max/min form, which E3 showed is unreachable; the replacement is
  justified by E3 and applies to R2's new data only. Attempt 1's registered verdict is
  not reinterpreted under it.
- **R-D** drift: for each (arm, page mode), the median user time of blocks 1-7 and of
  blocks 9-15 differ by at most 10 %.

**Decision (registered).** With gates passed:
- **STRONG**: interval lower bound ≥ 0.25. Address-translation stalls that huge pages
  remove are, natively measured, at least a quarter of IC's per-instruction user-time
  deficit.
- **SUBSTANTIAL**: lower bound ≥ 0.10 (and not STRONG).
- **BELOW A QUARTER**: upper bound < 0.25 (may co-occur with SUBSTANTIAL: "between a tenth
  and a quarter").
- otherwise **UNRESOLVED**.

**What R2 can and cannot say about H_mem.** H_mem is SUPPORTED at this cell iff R2 is
STRONG (a natively measured lower bound on the memory share of at least a quarter) or
`F_low ≥ 0.5`. R2 **cannot refute** H_mem: the cache-latency term stays inside
`[0.077, 1.91]` and no native test in this note isolates it. Any outcome other than
STRONG leaves H_mem "partly supported by the measured translation share X, remainder
unresolved". A huge-page speed-up is also not a claim that IC would run faster with
THP: system time rises by more than the user time falls in this VM.

**Inadmissible:** changing seeds, block count, thresholds, estimator or gates after the
run; dropping blocks; reporting the point estimate without the interval.

---

<!-- R2 results are appended below. Nothing above this line is edited after seeing them. -->

## R2 results (2026-09-30; `r2.json`, `analysis_r2.txt`)

**Host.** R2 ran after the session was resumed on a **different VM instance** (kernel
`fc-v50`, uptime 2 minutes; the calibration, registered stages, A1 and A2 ran on
`fc-v37`); same CPU model string, effective clock 3.24 GHz. `host_manifest_r2.txt`. R2's
estimate uses only its own 4 KiB / THP contrasts and the deterministic `Ir`, so the change
does not enter it, but R2's user times are not comparable in absolute terms to A2's (IC
4 KiB user time 3.27 s against 4.01 s).

**Gates.** R-G1 pass (every run recovers every target); R-N pass ((Q3 − Q1)/median of
4 KiB user time 0.061 IC, 0.073 rho); R-D pass (block-half drift −5.5 %…−0.8 %);
**R-V fails**: the median of the six chase-control speedups is **1.24**, below the
registered 1.3 (pairs: 1.26, 0.88, 1.22, 1.18, 1.36, 1.32), although the other half of
R-V holds (IC's sampled `AnonHugePages` median 281 MB against a 338 MB peak RSS).
**Registered R2 verdict: UNRESOLVED (gate).**

**Numbers the gate withholds a label from** (medians over 15 blocks, user time):

| | 4 KiB s | THP s | ratio |
|:--|--:|--:|--:|
| IC | 3.266 | 2.828 | 1.155 |
| rho R3 | 0.983 | 0.951 | 1.035 |

Translation share **0.219, 95 % bootstrap interval [0.153, 0.322]** (`d_IC = 0.0352`,
`d_rho = 0.0040`, `Δ_user = 0.1424` ns per instruction). System time rises in THP
mode (IC 0.24 → 0.73 s) but wall time is flat (IC 3.79 vs 3.77 s).

## Conclusion of the check (what is and is not established)

1. **Registered outcome of attempt 1: UNDETERMINED (G5).** Registered outcome of R2:
   UNRESOLVED (R-V). Neither is relabelled here. The literal labels come from gates that
   fail narrowly (a noise gate the VM cannot meet; a control speedup of 1.24 against
   1.3), and the point of the substantive reading below is to say what the numbers do
   and do not support, not to rescue a label.
2. **H_mem is neither established nor refuted at n = 41.** The simulator interval for
   the fraction of IC's per-instruction time gap that cache-miss stalls can explain is
   `[0.077, 1.91]`, wide because cachegrind carries no overlap information, and the one
   native contention test has no power on this VM (E1). The native huge-page contrast
   is the only mechanism-isolating measurement: two runs on two VM instances give a
   translation share of 0.257 (A2, no interval) and 0.219 (R2, interval [0.153, 0.322],
   validity gate failed), i.e. address translation on 4 KiB pages plausibly accounts
   for about a fifth to a quarter of the user-time gap. That is a description of these
   two runs, not a finding: neither cleared its own bar. **Roughly three quarters of the
   per-instruction gap is unexplained by anything measured here.**
3. **Hardware-neutral facts (deterministic, cachegrind, Ir within 0.013 % of the
   sweep's callgrind).** IC has 12.15 D1 misses per 1,000 instructions against rho's
   6.63, and 1.21 against 0.147 last-level data misses per 1,000 at a 2 MiB LL
   (8.2×); at 32 MiB, 1.10 against 0.081 (13.5×). IC's last-level misses barely fall as
   the simulated cache grows (15.0 M at 2 MiB, 14.7 M at 8 MiB, 13.7 M at 32 MiB), which
   is what a random access into a several-hundred-MB table looks like. By function (2 MiB
   LL), `main` (the index build and probe loop, inlined) holds 55 % of IC's instructions,
   66 % of its far read misses and 99 % of its far write misses; `extract` holds 34 % of
   instructions and 34 % of far read misses. Instruction-cache misses (0.025 and 0.088
   per 1,000) and simulated branch mispredicts (2.91 and 2.39 per 1,000; 2 % of the gap
   at the calibrated cost) are minor.
4. **Why the verdict does not depend on the answer.** Wall time is `Ir × ns/instruction`.
   At this cell IC retires 1.52× rho's instructions; were IC's stalls removed entirely so
   that it ran at rho's ns per instruction, IC would still take 1.52× rho's time. No
   cache-level engineering of IC flips the verdict here; a candidate has to cut probes
   and index work, as PR #966 concluded. What the memory question decides is only how
   large the extra wall-clock factor (about 3.0-3.6× here against 1.52×) is and where
   it comes from.

## Limits

One cell (n = 41, L = 1,024, K = 255), one VM class observed on two instances, one
implementation of each arm, IC as merged. Nothing here reaches n = 53 (IC's 1.36 GiB
index), larger n, other L, aarch64 or m = 83 / ECC2K-130 transfer. Cachegrind has no
prefetcher or TLB; the prices are pointer-chase measurements on this VM; there is no
hardware-counter access. Registered gates failed on noise/controls; no gate or
threshold was edited, and the follow-ups (E1-E3, A2, R2) are explicitly outside the
attempt-1 verdict. Raw per-target streams from the follow-up stages were first committed
uncompressed (commit `ed9df0d13`, about 25 MB of history) and are gzip'd with a sha256
manifest from the next commit on.

## Superseded statements

The sentences "IC's 1.36 GiB index runs at a lower IPC, which is consistent with the
memory-bound reading in PR #955 but was **not** tested here" (ladder note and PR #966)
now have a test at n = 41: its outcome is the paragraph above, and the memory-bound
reading remains a partly supported hypothesis, not a result. It has not been tested at
n = 53.
