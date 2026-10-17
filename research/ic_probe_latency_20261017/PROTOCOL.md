# Protocol: index-calculus probe time against a pointer-chase null object

Frozen before any run below. The autoresearcher contract is
`EXP-ICPERF-9cd390` (hypothesis `H-ICPERF-b3d1a4`, proposal
`IDEA-20261015-b41f68`, goal `GOAL-ICPERF-e6b6a4`). This file, `run_step.sh`
and `sources.sha256` are its launcher. The protocol commit is the commit
that adds them, recorded in each run's `launch.tsv`.

## Question

The repository's folded pair-state index calculus at n = 83 spends most of
its online interval looking up a root table larger than any cache (crypto
`RESEARCH_IC_BOUNDARY_LEDGER.md` §23.5: lookups take 53–96% of the
interval). The proposal reads the IC online time as two factors:

- probes per target, which should be geometric;
- the host's random-access latency over the table's bytes.

A rho step is register arithmetic instead. If both readings hold, an IC/rho
online ratio carries the host's memory system as a factor, and a one-point
ratio has the spread of an exponential.

This study tests three of the proposal's four claims on one host:

- (1) geometric probe counts;
- (2) ns per probe against a pointer-chase null object;
- (4) 2 MiB pages on the root table.

Claim (3), the cross-host spread, needs a second host and is **not tested
here**.

This is a measurement and attribution study. It makes no IC arm faster or
slower and adds no gap-to-rho row.

## Workload

The workload is the archived K1798 dense arm of
`research/n83_dense_root_20261013`:

- **Curve.** Koblitz `a = b = 1` over F_2^83, subgroup order
  8569786107849059.
- **Base.** `construct:83:1:1798`: 1,798 signed-Frobenius columns.
- **Table.** Packed dense root table, 2^29 slots of 16 bytes (8 GiB).
- **Seeds.** Rank seed 29.
- **Targets.** The 64 public targets `targets.jsonl`, SHA-256 `e9ed7ddb…37cc`.
  The answers sidecar `answers.jsonl` (`6a7022f7…9181`) is read only by the
  claims step.
- **Environment.** `KIC_PACKED_ROOT_TABLE=1 KIC_PACKED_ROOT_TABLE_DENSE=1`,
  as archived. `KIC_ROOT_BLOOM=0` pins the archived run's setting: it ran
  without the Bloom prefilter, whose default has since changed to on.
  `KIC_RANK_THREADS=4` (this host has 4 vCPUs); `KIC_TARGET_THREADS=1`.
- **Online interval.** IC ns per probe is per target `target_ms × 10^6 /
  probes`. A run's value is the median over its 64 targets.

## Instruments

- **Null object** (`examples/ic_probe_latency.rs chase`):
  - **Layout.** One 64-byte node per cache line over an array of the
    table's byte size (8 GiB). Sattolo's shuffle links the nodes into a
    single cycle.
  - **What it times.** A dependent chase (latency) and an independent
    random gather (throughput with memory-level parallelism), 2^24
    accesses per repetition, 7 repetitions.
  - **Pages.** The array is advised `MADV_NOHUGEPAGE` (4 KiB) or
    `MADV_HUGEPAGE` (2 MiB) before it is first written.
  - **Recorded.** `AnonHugePages` is recorded.
- **Huge-page switch** (`examples/koblitz_orbit_dlp_fast.rs`):
  - **Switch.** `KIC_ROOT_TABLE_HUGEPAGE=1` advises the root slots
    `MADV_HUGEPAGE` before they are filled. Off by default, and the
    contents and probe order are unchanged.
  - **Summary fields.** The summary records `root_table_hugepage_requested`,
    `anon_huge_pages_bytes` and `vm_hwm_bytes` (peak RSS).
- **Compute control:**
  - **Walk.** `examples/n83_rho_optimization_eval.rs`, the archived
    `fold-L64-DP6-E0` policy (`rho_config.json`), seed 2026101301, on the
    archived primary target.
  - **Measured.** ns per walk step = `online_phase_ns.walk_and_collision /
    outcome.walk_steps`.
  - **Pages.** It allocates no advised memory, so the page mode cannot
    reach it.
- **Distribution** (`geomfit`):
  - **Statistics.** CV, median/mean, max/mean (against ln n + 1), and the
    Kolmogorov distance to Geometric(1/mean).
  - **Status.** The KS 5% reference 1.36/√n is descriptive, because the
    parameter is estimated.

## Host

This session's cloud container:

- Intel Xeon @ 2.10 GHz, family 6 model 207, 4 vCPUs, KVM;
- 15 GiB RAM;
- THP `madvise`;
- isolation L0, no PMU.

The archived K1798 run peaked at 8.7 GB RSS, so the IC fits with no other
large process running. No other build or test runs during the session.

## Trials, in order

| # | trial | step | what |
|--:|---|---|---|
| 0 | build | `build` | release build of the three examples at the protocol commit |
| 1 | stage0 | `stage0` | `geomfit` over the four archived 64-target panels (zero compute) |
| 2 | sweep | `sweep` | chase at 1 MiB, 16 MiB, 256 MiB, 1 GiB and 8 GiB × {4 KiB, 2 MiB} |
| 3 | rho-1 | `rho` | compute control, before |
| 4 | A1 | `ic-4k` | 8 GiB chase on 4 KiB pages, then K1798 on 4 KiB pages, 64 targets |
| 5 | B1 | `ic-2m` | the same on 2 MiB pages |
| 6 | A2 | `ic-4k` | repeat of A1 |
| 7 | B2 | `ic-2m` | repeat of B1 |
| 8 | rho-2 | `rho` | compute control, after |
| 9 | claims | `claims` | the claims below, computed in Rust from trials 1–8 |

Caps: IC 10,800 s per run, rho 600 s. A timed-out or
failed trial stays in the record with its exit status. No trial is re-run
to change an outcome. A re-run for infrastructure loss (the process died
before its first target) is allowed once and recorded.

## Claims

All are computed by `ic_probe_latency claims`, never by hand.

| id | claim | rule |
|---|---|---|
| G.scalars | every target of A1, B1, A2 and B2 recovers its published scalar | PASS / FAIL |
| G.probes | every target's probe count equals the archived K1798 dense panel's | PASS / FAIL |
| H1.\<base\> | probe counts are geometric, for each archived base and for A1 on this host | SUPPORT: CV ∈ [0.9, 1.1] and median/mean ∈ [0.6, 0.8]; FALSIFIED: CV < 0.6 or > 1.5; else NOT_SUPPORTED |
| L1 | median IC ns/probe (A1, A2) ÷ median 8 GiB dependent-chase ns/access (A1, A2) | SUPPORT: ∈ [1.0, 2.0]; FALSIFIED: < 0.5; else NOT_SUPPORTED |
| L1g | the same over the independent gather | descriptive |
| T1.applied | `AnonHugePages` ≥ 90% of the table bytes in B1, B2 and < 10% in A1, A2 | PASS / FAIL |
| T1.null | positive control: 8 GiB chase on 2 MiB ÷ on 4 KiB, B vs A runs | PASS if ≤ 0.9 |
| T1 | max(r1, r2), r_k = B_k ÷ A_k on IC ns/probe | gated on T1.applied and T1.null; SUPPORT: ≤ 0.67; FALSIFIED: both \|r − 1\| < 0.10; else NOT_SUPPORTED |
| R | rho ns/step, rho-2 vs rho-1 (max/min) | descriptive drift of the compute control |
| D | IC ns/probe ÷ rho ns/step, raw and ÷ chase latency, this host | descriptive |

The thresholds are the proposal's own predictions and falsification
conditions, copied before any run:

- H1 is its prediction (1).
- L1 is its prediction (2) with "falsified below 0.5".
- T1 is its prediction (4) with "< 10% falsifies the TLB mechanism".

The T1 gate is this protocol's addition. If the null object shows no
huge-page gain, or the kernel did not back the table with huge pages, T1
reads INSTRUMENT_CANNOT_SEE, never FALSIFIED.

How to read the result:

- **G.probes fails while G.scalars passes.** The on-host walk differs from
  the archived one, for example by host or ISA. The run is reported as a
  deviation. H1 and L1 still use the on-host probe counts.
- **G.scalars fails.** The run is invalid.
- **L1 is NOT_SUPPORTED above 2.** The probe carries a compute term beyond
  one latency, which this study does not separate (see "Not measured").
- **L1 is below 1 while L1g is near 1.** The lookups overlap, so the
  latency model is the wrong null.
- **H1 is NOT_SUPPORTED with CV in [0.6, 0.9).** That is the proposal's
  model with a lighter tail than geometric, a refinement.

## Not measured

The proposal's cache-resident control is a small table whose ns/probe would
be the compute term c0. It is not run, because it is not affordable here:

- the rank stage needs about K relations, each costing r/(A·#states)
  probes with #states ≈ 83K², so rank work grows as 1/K;
- a K = 256 smoke run was still in its rank stage after 12 minutes on this
  host;
- K = 32 would take hours.

c0 therefore stays unseparated, and L1's band [1.0, 2.0] admits a compute
term of up to one latency.

## Disclosure

While writing the launcher, the designer smoke-tested `geomfit` on two
archived panels (k1271_regular and k1798_dense) and saw their CV and
median/mean. A stand-in test of the `claims` step on archived files (archived
K1798 rows copied into four fake run directories) also printed the H1
verdicts of those two bases: NOT_SUPPORTED, with CV 0.81 and 0.92 and
median/mean 0.91 and 0.59. Stage 0 reads archived data, so its outcome was
fixed before this protocol. The thresholds above are the proposal's,
unchanged after seeing it. The
archived panels are public inputs the proposal itself cites. A 1 GiB chase
was also smoke-tested on both page sizes (195 and 151 ns/access). A K = 256
producer run was started on both page modes and stopped in its rank stage,
after `AnonHugePages` showed 254 MiB of its 256 MiB table on huge pages with
the switch on. No target was solved and no IC per-probe time was seen. No
smoke output is a run record.

## Commands

The autoresearcher trial plan calls `research/ic_probe_latency_20261017/run_step.sh STEP OUT`
from a checkout at the protocol commit. By hand:

```sh
R=<runs dir>; S=research/ic_probe_latency_20261017/run_step.sh
$S build $R/build; $S stage0 $R/stage0; $S sweep $R/sweep; $S rho $R/rho-1
$S ic-4k $R/a1; $S ic-2m $R/b1; $S ic-4k $R/a2; $S ic-2m $R/b2; $S rho $R/rho-2
$S claims $R/claims $R/stage0 $R/a1 $R/b1 $R/a2 $R/b2 $R/rho-1 $R/rho-2
```
