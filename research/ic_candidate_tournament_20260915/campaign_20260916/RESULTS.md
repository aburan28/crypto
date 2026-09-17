# Continued IC tournament results

Batch (16-target) parity verdict: **True**; selected implementation **combined_descent**.
Single-target verdict: **strictly below rho in both metrics on every cell**; selected implementation **tiny2** (round-0007), re-measured with the per-target descent certificate under the merged checker in round-0008 (retained, `beats_rho_strict`), and in round-0010 on a lean static executable under an evaluator that no longer forks itself: 0.5347 of rho's instructions and 0.8445 of its native time, `beats_rho_strict` on every cell (the strict-win record), after **tiny_batch1** reached parity in round-0006.

Completed three new tournaments with **4,704 profiled trials**, each paired with a fresh native run, followed by round 0006-batch16 (1,680 further profiled trials, incumbent retained, index-calculus admission enforced by a per-target descent certificate) and by the single-target rounds 0006–0010 (1,584, 1,488, 1,356, 1,356 and 1,440 profiled trials: parity, a strict win over rho, that win re-measured with the descent certificate, the same pipeline on a leaner shared job, and the strict win on every cell on a lean static executable under an evaluator that no longer forks itself). Every listed round passed its independent artifact/correctness audit. The 36-trial batch screen is separate development evidence.

## Profiled instruction cost relative to matched rho

| Round | Targets per cold job | Selected candidate | Instructions / rho | Paired 95% interval |
|---|---:|---|---:|---|
| [round-0002](../runs/round-0002/REPORT.md) | 1 | batch16 | 13.0349 | not calculated in original protocol |
| [round-0003b](../runs/round-0003b/REPORT.md) | 1 | combined_batch8 | 7.1494 | [6.563380843095217, 7.844932268218591] |
| [round-0004](../runs/round-0004/REPORT.md) | 1 | folded_lift_batch4 | 3.0429 | [2.688274518239304, 3.4833783647230367] |
| [round-0005-batch16](../runs/round-0005-batch16/REPORT.md) | 16 | combined_descent | 0.7195 | [0.6374471292425982, 0.8219836239443673] |
| [round-0006-batch16](../runs/round-0006-batch16/REPORT.md) | 16 | incumbent retained (combined_descent + descent certificate) | 0.7259 | [0.6410, 0.8372] |
| [round-0006](../runs/round-0006/REPORT.md) | 1 | tiny_batch1 | 0.8264 | [0.8032467836918035, 0.8488166193357406] |
| [round-0007](../runs/round-0007/REPORT.md) | 1 | tiny2 | 0.7376 | [0.7065693377992373, 0.7643545073767521] |
| [round-0008](../runs/round-0008/REPORT.md) | 1 | incumbent retained (tiny2 + descent certificate) | 0.7379 | [0.7144, 0.7572] |
| [round-0009](../runs/round-0009/REPORT.md) | 1 | incumbent retained (same pipeline; shared construction in single-word arithmetic for both arms) | 0.5511 | [0.5279, 0.5881] |
| [round-0010](../runs/round-0010/REPORT.md) | 1 | incumbent retained (same pipeline; static executable, C entry point, one CPUID; evaluator spawns without forking itself) | 0.5347 | [0.5106, 0.5645] |

## Native process time relative to matched rho

| Round | Targets per cold job | Selected candidate | Time / rho | Paired 95% interval |
|---|---:|---|---:|---|
| round-0002 | 1 | batch16 | diagnostic only | old polling-based timing excluded |
| round-0003b | 1 | combined_batch8 | 2.7579 | [2.408587066800431, 3.0866033786867497] |
| round-0004 | 1 | folded_lift_batch4 | 1.5213 | [1.4825577742485183, 1.5518235985269213] |
| round-0005-batch16 | 16 | combined_descent | 0.7796 | [0.6991620583548833, 0.8736653449888768] |
| round-0006-batch16 | 16 | incumbent retained (combined_descent + descent certificate) | 0.7853 | [0.7040, 0.8771] |

## Round 0006: index-calculus admission enforced, incumbent retained

Round 0006 ([pre-registration](ROUND6.md), [report](../runs/round-0006-batch16/REPORT.md))
ran on a different host and compiler (4 vCPU, rustc 1.94.1) from rounds
0002–0005, so only its within-round ratios are comparable. Two things changed
before measurement:

- **Only index-calculus algorithms compete with rho.** Every IC arm now
  reports, per target, the relation `[a]G + [b]Q = Σ P_i` its logarithm was
  derived from, and the frozen checker verifies it in the group and checks the
  scalar as its consequence under the verified column logs. All 1,680 receipts
  carry `certified_descents` equal to the target count. The incumbent is the
  round-0005 winner plus this certificate.
- **Seven profiled-bottleneck candidates**, each an exact source change with a
  preflight equivalence test: Itoh–Tsujii and binary-Euclid field inversion,
  one single-word curve per build and per verified batch, a lazily built
  `FieldStructure`, precomputed `λ^k`, single-word orbit tables, and their
  combination.

Development (four cells, three targets each, three repetitions; 324/324
verified), candidate/incumbent instruction ratios: combined 0.8026, euclid_inv
0.9142, lambda_table 0.9452, fast_orbits 0.9513, fast_curve_once 0.9603,
it_inv 0.9722, lazy_field 0.9732. Native-time ratios were much closer to one
(combined 0.93; the single mechanisms 0.96–1.01). Selection locked `combined`.

Confirmation and replay (60 fresh 16-target fixtures, five cells, three
repetitions; 1,080/1,080 verified): combined/incumbent **0.8018** instructions
(95% paired interval 0.7874–0.8146) and **0.9237** native time (0.9102–0.9365)
on confirmation; 0.8018 and 0.9136 on replay. Every cell was below 0.83 in
instructions and below 0.94 in native time. The pre-registered promotion gate
requires at most 0.80 on both metrics, so the challenger missed the
instruction threshold by 0.0018 and the native threshold by a wide margin:
**incumbent retained**. This is the outcome the pre-registration named as
plausible: 58% of the incumbent's instructions are shared with the rho arm
(target construction, general final verification, startup, reporting), and
native process time is dominated by fixed per-process cost the mechanisms do
not touch.

The retained incumbent reproduced the batch parity result on this host:
incumbent/rho **0.7259** instructions (0.6410–0.8372) and **0.7853** native
time (0.7040–0.8771) on confirmation, 0.7259 and 0.7942 on replay, every cell
at most 0.94 and 0.98 respectively; `rho_parity` true, `beats_rho` true.

Round 0006 is the first round without a promotion (one of the three the
protocol allows before an unsuccessful line stops). Evidence-based successor,
not yet run: the Euclid inversion dominated Itoh–Tsujii in every stage (0.914
versus 0.972), so a combination built on `euclid_inv` would be expected near
0.75 instructions against the incumbent; its native-time ratio would still be
expected near 0.89, so it could pass the instruction gate but not the native
gate as currently frozen. Whether the native-progress gate should stay on a
panel whose native time is dominated by fixed process cost is a protocol
decision to make before the next round, not during it. The mechanical
parameter-neighbour registry is in [round6-next-candidates.json](round6-next-candidates.json).

| round-0006 | 1 | tiny_batch1 | 0.9766 | [0.957057247571661, 0.9946203513598987] |
| round-0007 | 1 | tiny2 | 0.9461 | [0.9292926395602298, 0.9639483494225052] |
| round-0008 | 1 | incumbent retained (tiny2 + descent certificate) | 0.9355 | [0.9124, 0.9610] |
| round-0009 | 1 | incumbent retained (same pipeline; shared construction in single-word arithmetic for both arms) | 0.9406 | [0.9131, 0.9811] |
| round-0010 | 1 | incumbent retained (same pipeline; static executable, C entry point, one CPUID; evaluator spawns without forking itself) | 0.8445 | [0.8123, 0.8814] |

## Single-target parity (round-0006)

[Pre-registered](ROUND6-single-target.md) before measurement; seed 2026091606; 1,584 profiled trials with paired native runs, all verified; incumbent round-0004's `folded_lift_batch4`. The selected `tiny_batch1` runs the same bounded algorithm — the identical factor-base point set, three summands, a folded pair-sum table, a full-rank solve with every column logarithm certified in the group, a `[a]G + [b]Q` descent and the unchanged general-arithmetic final check — in single-word arithmetic: a Euclidean field inverse, a normal-basis orbit key, density-sized table rows, walked probes, no thread pool, and a directly written report ([patch](round6-tiny.patch), [winner record](WINNER-single-target.json)).

Confirmation and replay, one target per cold job, five cells:

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / incumbent (confirmation) | 0.2702 | [0.2320, 0.3073] | 0.6511 | [0.6236, 0.6801] |
| winner / rho (confirmation) | 0.8264 | [0.8032, 0.8488] | 0.9766 | [0.9571, 0.9946] |
| winner / rho (replay) | 0.8264 | [0.8032, 0.8488] | 0.9854 | [0.9641, 1.0062] |

Per-cell winner/rho instruction ratios 0.788–0.850 and native ratios 0.947–1.012; all within the 1.10 parity margin on both stages. The winner uses fewer instructions than rho with the whole interval below one. In native process time the point estimates are below one but the replay interval reaches 1.006, so the native result is parity, not a demonstrated speed advantage. The controls in the same round: `combined_descent` (the batch-panel winner's source, 0.949 of the incumbent), `fast_report` (serialisation only, 0.922) and `tiny_fulltable` (the table rule disabled, 0.2785 of the incumbent against 0.2711 for `tiny`; the rule's gain is on n13a0, where it builds two of seven rows).

Where the ratio comes from: on n13a0 rho's job is about 6.2 million instructions, of which curve construction, target hashing, general-arithmetic final verification and startup — shared by both arms — are about 4.3 million and the rho solve about 1.9 million, most of it the jump-table scalar multiplications and the general-arithmetic collision finish rather than the walk. The IC arm's own work has to fit in that solve plus ten percent of the job; on that cell it now takes about 0.7 million. The shipped rho's per-job cost on these cells is therefore mostly fixed setup; a rho specialised the same way has not been measured and would be cheaper. This is an implementation result in fixed-compiler Valgrind amd64 instructions and matched native time on five small Koblitz cells, class engineering; no arithmetic-complexity, family-wide or cryptographic-size claim follows.

## Beating rho on the single target (round-0007)

[Pre-registered](ROUND7-single-target.md) with `--objective rho`: promotion requires no regression against the round-0006 winner and, on confirmation and replay, candidate/rho upper paired 95% limits and every cell **below one** in both instructions and native process wall (`beats_rho_strict`). Seed 2026091607; 1,488 profiled trials with paired native runs, all verified; factor-base fingerprints equal to the incumbent's in every cell. The selected `tiny2` ([patch against the round-0006 winner source](round7-tiny2.patch), [winner record](WINNER-single-target.json)) keeps every obligation of round-0006 and changes how the work is done: a work-minimising rule builds one or two folded table rows instead of seven or eight (one row on n13–n19, two on n23), scalar multiplications run in López–Dahab projective coordinates or off a doubling table of `G` with all column certifications batched behind one inversion per bit, squaring goes through the carry-less multiplier, and the module carries its own reduction tables, nibble-table normal basis and a bit-exact scalar ChaCha12 sampler (tested draw for draw against `StdRng`), with no libm call or hashed container in the job.

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / incumbent (confirmation) | 0.8872 | [0.8493, 0.9242] | 0.9712 | [0.9531, 0.9862] |
| winner / rho (confirmation) | 0.7376 | [0.7066, 0.7644] | 0.9461 | [0.9293, 0.9639] |
| winner / rho (replay) | 0.7376 | [0.7066, 0.7644] | 0.9414 | [0.9215, 0.9611] |

Per-cell winner/rho: instructions 0.683–0.782, native 0.918–0.977 (n13a0 0.977 confirmation, 0.962 replay). The ablations in the same round: the row rule alone (`tiny2_rows`) 0.7667 of rho's instructions and 0.9469 of its native time on development; the arithmetic changes alone (`tiny2_arith`) 0.8162 and 0.9389; both together 0.7444 and 0.9269 (confirmation 0.7376 and 0.9461).

What the native figure means: in this virtual machine a complete cold job is about 7.5 ms of process wall, of which roughly 3.5 ms is process creation and loading and 0.5–2 ms the shared curve construction; a CPUID instruction traps at about 10 µs, so cold per-process costs (feature detection, first libm call, first-touch pages) weigh as much as the arithmetic. The winner's own work is about 0.2–0.5 ms per job against rho's 0.4–1.3 ms; after the shared part both arms pay, that is the 5–6% margin measured. The instruction ratio is the hardware-independent figure. Rho is the shipped implementation, whose per-job cost on these cells is mostly fixed setup and cold-process cost; a rho specialised the same way has not been measured and would be cheaper. Class engineering; no arithmetic-complexity, family-wide or cryptographic-size claim follows.

## Round 0008: the single-target win under the merged protocol

Rounds 0006 and 0007 were frozen with the evaluator of their day, before the index-calculus admission rule of round 0006-batch16 (a per-target descent relation, verified in the group) reached `main`; their frozen outputs carry no such relation and their audits stand under their own frozen checkers only. Round 0008 ([pre-registration](ROUND8-single-target.md), [report](../runs/round-0008/REPORT.md), `--objective rho`, seed 2026091608) re-measures the round-0007 winner with the certificate — `solve_target` returns the relation `[a]G + [b]Q = Σ P_i` it used and the worker writes it per solution ([patch](round8-tiny2_cert.patch)) — under the merged evaluator and checker, against a configuration control (`batch_trials: 4`, ratio 1.0019 to the incumbent, not promoted). All 1,356 profiled trials and paired native runs verified; every receipt carries `certified_descents = 1` and no degenerate descent; factor-base fingerprints equal the incumbent's in every cell. Decision: retained, and for the retained incumbent `beats_rho_strict` is true:

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / rho (confirmation) | 0.7379 | [0.7144, 0.7572] | 0.9355 | [0.9124, 0.9610] |
| winner / rho (replay) | 0.7379 | [0.7144, 0.7572] | 0.9414 | [0.9267, 0.9545] |

Per-cell winner/rho: instructions 0.695–0.764, native 0.902–0.967 (n13a0 0.967 confirmation, 0.949 replay). This is the single-target result the operation now stands on; the round-0007 figures agree with it within their intervals.

## Round 0009: the same pipeline on a leaner shared job

Every job pays, before either algorithm starts, for the curve construction and the public target, and at the end for the general-arithmetic final check; through round 0008 the construction ran in general big-integer arithmetic (3–11 million instructions, 0.5–1.8 ms per job for both arms). Round 0009 ([pre-registration](ROUND9-single-target.md), [report](../runs/round-0009/REPORT.md), `--objective rho`, seed 2026091609) puts that shared work on the tested single-word paths — the same abscissa sequence, lifts, tests and eigenvalue candidates, equal to the general arithmetic in a test on every cell and reproducing all 184 frozen fixtures of rounds 0007 and 0008 byte for byte ([patch](round9-fastcurve.patch)) — in the baseline executable that rho, the incumbent and the batch-4 control all run. The IC pipeline itself is unchanged from round 0008; the rho solve is the shipped implementation, untouched. All 1,356 trials verified with certificates; retained (control 1.005).

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / rho (confirmation) | 0.5511 | [0.5279, 0.5881] | 0.9406 | [0.9131, 0.9811] |
| winner / rho (replay) | 0.5511 | [0.5279, 0.5881] | 0.9323 | [0.9041, 0.9581] |

Median complete job: 4.33 million instructions and 6.1 ms for the winner against 8.1 million and 6.5 ms for rho (round 0008: 7.3 and 7.7 ms). Per-cell winner/rho instructions 0.528–0.625. The strict per-cell native gate was **not** met in this round: n13a0 came out at 1.017 on confirmation (0.960 on replay; the other cells 0.89–0.96 on both stages). On that cell the two arms differ by about 0.2 ms of their own work inside a 6 ms process, which is within the noise of twelve cases, so `beats_rho_strict` is false here while the aggregate native intervals are below one on both stages and `rho_parity` holds. Round 0008 remains the strict-win record; round 0009 shows what the instruction margin is once the shared job stops hiding it, and that cheaper shared work helps rho's wall time as much as the winner's.

## Round 0010: a lean job and a faithful clock — the strict win on every cell

Round 0009's one failing cell was decided by noise that neither arm produces: on `n13a0` the arms differ by 0.2 ms inside a 5.8 ms cold process whose spread (4.9–9.2 ms) came from the evaluator, which created every child with a `preexec_fn` and so made CPython `fork()` itself, charging a copy of its page tables to every native wall (re-creating `execute` with a 400 MiB parent: 23.7 ms per job against 2.9 ms without the fork). Round 0010 ([pre-registration](ROUND10-single-target.md), [report](../runs/round-0010/REPORT.md), `--objective rho`, seed 2026091610) changes the instrument and the job. The evaluator spawns without forking itself: the calling thread is pinned before the spawn so the child inherits the CPU, the memory and core caps are applied with `prlimit` while the child is still blocked on `stdin`, and the watchdog is armed before the timing window; the contract and every receipt record the spawn. The baseline executable, shared by every arm including rho ([patch](round10-lean.patch) on the round-0009 source), is statically linked through a sealed `.cargo/config.toml`, uses the C entry point instead of the standard runtime's start-up, and — the only IC-only change — issues one `CPUID` instead of the standard library's feature cache. It reproduces all 304 frozen round-0009 fixture outputs byte for byte; the shared `startup_and_input` phase fell from 385K to 84K instructions for both arms. Nothing about the algorithms changed. All 1,440 trials verified, every IC job's logarithm certified by its descent relation; retained (batch-4 control eliminated at selection; the ablation control `lean_stdprobe`, which restores the standard feature cache, measured 1.0002 of the incumbent's instructions and 0.999 [0.989, 1.010] of its native time on confirmation, 0.989 [0.974, 1.002] on replay: the probe's cost is not resolvable at this noise).

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / rho (confirmation) | 0.5347 | [0.5106, 0.5645] | 0.8445 | [0.8123, 0.8814] |
| winner / rho (replay) | 0.5347 | [0.5106, 0.5645] | 0.8665 | [0.8307, 0.9078] |

Per-cell winner/rho: instructions 0.503–0.594, native time 0.796–0.911 on confirmation and 0.812–0.941 on replay; `n13a0`, the cell that failed in round 0009, is 0.911 and 0.941. **`beats_rho_strict` holds**: every upper limit and every cell below one in both metrics on both stages, so round 0010 replaces round 0008 as the strict-win record with the same pipeline. Median complete job on confirmation: 4.03 million instructions and 2.41 ms for the winner against 7.68 million and 2.86 ms for rho (round 0009 under the old clock: 4.33 million and 6.1 ms against 8.1 million and 6.5 ms). Absolute times and instruction counts are not comparable with earlier rounds, whose native walls included the evaluator's fork and whose instruction counts included the dynamic loader; only within-round ratios are claimed, as before. The native ratio is still bounded below by the process creation both arms pay.

## Interpretation

Every ratio uses a fresh matched rho run in the same round. The 16-target panel charges all setup once to the complete job and solves every target; it is separate from the single-target result, and no ratio combines the two panels. Rho uses the existing per-target solver API on the same constructed curve. Additional cross-target rho optimizations, and a rho specialised like the round-0006 winner, have not been measured here.

The parity rule was declared before measurement: both candidate/rho upper paired 95% limits and every curve-cell ratio must be at most 1.10, in instructions and native time, on confirmation and replay. The full decisions contain the replay evidence.

These are implementation improvements in fixed-compiler Valgrind amd64 guest instructions and matched native wall time. Kernel/device and external-audit work are outside the instruction count. No arithmetic-complexity, broader-family, or cryptographic-size claim follows. The K-instruction rank floor is deliberately weak.

The successful mechanisms are exact arithmetic substitutions, folded pair-table construction, smaller relation batches, cheaper descent initialization/checking and, in round-0006, a single-word implementation of the whole pipeline. Every returned scalar still passes general worker verification and the independent Python checker.

## Reproduce and review

- [16-target winner source/configuration](WINNER.json) and [cumulative source patch](WINNER.patch).
- [Single-target winner record](WINNER-single-target.json) (round-0008: the round-0007 `tiny2` plus the descent certificate, [patch](round8-tiny2_cert.patch), [round8_candidates.py](round8_candidates.py)).
- Round-0010: the [lean-job patch](round10-lean.patch) applied to the round-0009 source for both arms (static link via `.cargo/config.toml`, C entry point, one CPUID), the [feature-cache ablation](round10-lean_stdprobe.patch), and [round10_candidates.py](round10_candidates.py).
- Round-0009: the [shared-construction patch](round9-fastcurve.patch) applied to the round-0008 source for both arms, and [round9_candidates.py](round9_candidates.py).
- Round-0007: the [patch against the round-0006 winner source](round7-tiny2.patch), the [row-rule](round7-tiny2_rows.patch) and [arithmetic](round7-tiny2_arith.patch) ablations, and [round7_candidates.py](round7_candidates.py) to rebuild the trees.
- Round-0006: the [parity patch against the round-0005 winner source](round6-tiny.patch), the [serialisation-only](round6-fast_report.patch) and [full-table](round6-tiny_fulltable.patch) controls, and [round6_candidates.py](round6_candidates.py).
- [Operation plan](PLAN.md), [single-target successor plan](ROUND4.md), [batch plan](ROUND5.md), [single-target parity pre-registration](ROUND6-single-target.md), [beat-rho pre-registration](ROUND7-single-target.md), [certificate re-measurement pre-registration](ROUND8-single-target.md), [leaner shared job pre-registration](ROUND9-single-target.md) and [lean job and faithful clock pre-registration](ROUND10-single-target.md).
- [Operating guide](../OPERATIONS.md) and [controller tests](controller-tests.json).
- [Arithmetic/table equivalence tests](preflight-next-tests.log) and [installed skill validation](skill-validation.json).
- [Retained failed build attempt](../runs/round-0003/prepare_failure.json); no measurements came from it.

Production library defaults were not changed. The frozen source, worker, configurations, raw profiles, receipts and replay evidence are retained under each linked round.
