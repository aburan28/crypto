# Continued IC tournament results

Batch (16-target) parity verdict: **True**; selected implementation **combined_descent**.
Single-target verdict: **strictly below rho in both metrics on every cell**; selected implementation **tiny2** (round-0007), re-measured with the per-target descent certificate under the merged checker in round-0008 (retained, `beats_rho_strict`), in round-0010 on a lean static executable under an evaluator that no longer forks itself (0.5347 of rho's instructions and 0.8445 of its native time, `beats_rho_strict` on every cell), in round-0011 with the general binary-field arithmetic in words for both arms (0.4796 of rho's instructions and 0.8653 of its native time, `beats_rho_strict` on every cell), in round-0012 as a musl static executable without relocations whose worker serves the heap from an arena, for both arms (0.4369 of rho's instructions and 0.7988 of its native time, `beats_rho_strict` on every cell), in round-0013 with one field inversion per scalar product for both arms, which halved rho's job because its setup is scalar products (0.7697 of rho's instructions and 0.9098 of its native time, `beats_rho_strict` on every cell: the strict-win record under the corrected baseline), in round-0014 with the worker built as one optimisation unit, which made both arms faster again and rho slightly more so (0.7818 and 0.9424), and in round-0015 on that same executable against fresh fixtures: 0.7372 and 0.9246, `beats_rho_strict` still on every cell, after **tiny_batch1** reached parity in round-0006. The last two figures are the same executable measured twice, which is the size of the fixture variation these ratios carry.

Completed three new tournaments with **4,704 profiled trials**, each paired with a fresh native run, followed by round 0006-batch16 (1,680 further profiled trials, incumbent retained, index-calculus admission enforced by a per-target descent certificate) and by the single-target rounds 0006–0011 (1,584, 1,488, 1,356, 1,356, 1,440 and 1,440 profiled trials: parity, a strict win over rho, that win re-measured with the descent certificate, the same pipeline on a leaner shared job, the strict win on every cell on a lean static executable under an evaluator that no longer forks itself, and the same again with the general binary-field arithmetic in words). Every listed round passed its independent artifact/correctness audit. The 36-trial batch screen is separate development evidence.

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
| [round-0011](../runs/round-0011/REPORT.md) | 1 | incumbent retained (same pipeline; the general binary-field arithmetic in words for both arms, one fast curve per job, tabulated irreducible, one-word eigenvalue root) | 0.4796 | [0.4498, 0.5217] |
| [round-0012](../runs/round-0012/REPORT.md) | 1 | incumbent retained (same pipeline; musl static executable without relocations and a bump-pointer arena allocator, both arms) | 0.4369 | [0.4069, 0.4815] |
| [round-0013](../runs/round-0013/REPORT.md) | 1 | incumbent retained (same pipeline; scalar products in López–Dahab coordinates for both arms, which halves rho's own job) | 0.7697 | [0.7148, 0.8370] |
| [round-0014](../runs/round-0014/REPORT.md) | 1 | incumbent retained (same pipeline; the worker built as one optimisation unit, both arms) | 0.7818 | [0.7200, 0.8601] |
| [round-0015](../runs/round-0015/REPORT.md) | 1 | incumbent retained (nothing shared changed; the round-0014 executable on fresh fixtures) | 0.7372 | [0.6628, 0.8326] |

## Native process time relative to matched rho

| Round | Targets per cold job | Selected candidate | Time / rho | Paired 95% interval |
|---|---:|---|---:|---|
| round-0002 | 1 | batch16 | diagnostic only | old polling-based timing excluded |
| round-0003b | 1 | combined_batch8 | 2.7579 | [2.408587066800431, 3.0866033786867497] |
| round-0004 | 1 | folded_lift_batch4 | 1.5213 | [1.4825577742485183, 1.5518235985269213] |
| round-0005-batch16 | 16 | combined_descent | 0.7796 | [0.6991620583548833, 0.8736653449888768] |
| round-0006-batch16 | 16 | incumbent retained (combined_descent + descent certificate) | 0.7853 | [0.7040, 0.8771] |
| round-0006 | 1 | tiny_batch1 | 0.9766 | [0.957057247571661, 0.9946203513598987] |
| round-0007 | 1 | tiny2 | 0.9461 | [0.9292926395602298, 0.9639483494225052] |
| round-0008 | 1 | incumbent retained (tiny2 + descent certificate) | 0.9355 | [0.9124, 0.9610] |
| round-0009 | 1 | incumbent retained (same pipeline; shared construction in single-word arithmetic for both arms) | 0.9406 | [0.9131, 0.9811] |
| round-0010 | 1 | incumbent retained (same pipeline; static executable, C entry point, one CPUID; evaluator spawns without forking itself) | 0.8445 | [0.8123, 0.8814] |
| round-0011 | 1 | incumbent retained (same pipeline; the general binary-field arithmetic in words for both arms, one fast curve per job, tabulated irreducible, one-word eigenvalue root) | 0.8653 | [0.8396, 0.8919] |
| round-0012 | 1 | incumbent retained (same pipeline; musl static executable without relocations and a bump-pointer arena allocator, both arms) | 0.7988 | [0.7668, 0.8319] |
| round-0013 | 1 | incumbent retained (same pipeline; scalar products in López–Dahab coordinates for both arms, which halves rho's own job) | 0.9098 | [0.8894, 0.9313] |
| round-0014 | 1 | incumbent retained (same pipeline; the worker built as one optimisation unit, both arms) | 0.9424 | [0.9126, 0.9715] |
| round-0015 | 1 | incumbent retained (nothing shared changed; the round-0014 executable on fresh fixtures) | 0.9246 | [0.9034, 0.9503] |

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

## Round 0011: the general arithmetic in words

After round 0010 the largest phase of the IC arm on every cell was the general-arithmetic final check both arms run on every recovered scalar (0.78M of its 1.87M instructions on `n13a0`, 2.83M of 6.82M on `n23a0`): the general field squared one bit at a time, reduced one bit position at a time and allocated per operation. Round 0011 ([pre-registration](ROUND11-single-target.md), [report](../runs/round-0011/REPORT.md), `--objective rho`, seed 2026091611) rewrites that arithmetic word-level — bit spreading for squares, carry-less products word by word, reduction a word of high bits at a time — with the previous implementation retained under `#[cfg(test)]` as the reference and a test comparing the two on 31 widths from 1 to 571 bits ([patch](round11-wordfield.patch) on the round-0010 source, shared by every arm including rho). The construction builds one single-word curve per job for the point count, generator search, eigenvalue check and target lift, tabulates the sparse irreducible polynomial (a test re-runs the search per degree and pins the table to it), takes the eigenvalue's square root in one word with the same Tonelli–Shanks steps as the big-integer routine (tested equal on every cell), and every field implementation shares one cached CPUID probe — the general `Gf2` constructor had been probing the standard library's feature cache for both arms, which is why round 0010's ablation of the IC arm's own probe measured nothing. The executable reproduces all 304 frozen round-0010 fixture outputs byte for byte; nothing about the algorithms changed. All 1,440 trials verified, every IC job's logarithm certified; retained. The IC-only challenger `wordfield_fastio` ([patch](round11-fastio.patch): the report's decimal numbers two digits at a time from a table, appended without a UTF-8 validation pass) measured 0.9651 [0.9614, 0.9692] of the incumbent's instructions on both final stages, inside the no-regression line, but its native upper limit was 1.0086 on confirmation (0.9982 point, cells 0.992–1.006; replay 0.9885 [0.9781, 0.9988]), so it was **not** promoted: a 3.5% instruction saving that the native clock cannot resolve at this noise. The batch-4 control was eliminated at selection (1.018 of the incumbent's instructions). The five pre-existing failures of the general pair-table tests are recorded in the pre-registration; they fail identically on the unmodified round-0010 snapshot.

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / rho (confirmation) | 0.4796 | [0.4498, 0.5217] | 0.8653 | [0.8396, 0.8919] |
| winner / rho (replay) | 0.4796 | [0.4498, 0.5217] | 0.8654 | [0.8454, 0.8883] |

Per-cell winner/rho: instructions 0.446–0.566, native time 0.823–0.912 on confirmation and 0.841–0.907 on replay. **`beats_rho_strict` holds** on every cell and both stages, so round 0011 is the strict-win record with the same pipeline. Median complete job on confirmation: 2.14 million instructions and 1.98 ms for the winner against 4.71 million and 2.30 ms for rho (round 0010, same clock and host: 4.03 million and 2.41 ms against 7.68 million and 2.86 ms). Rho's own solve also shrank where it uses the general arithmetic, which is why the native ratio moved less than the instruction ratio. The final check remains the largest IC phase on the largest cell; its remaining cost is the per-operation heap allocation of the general element type.

## Round 0012: the C runtime and the heap out of the way

After round 0011 two costs that were not the job's own were a large share of both arms. The executable was a static-pie glibc image, and on this virtual machine a Rust program that reads stdin and exits takes 1.17 ms from `posix_spawn` to reap against 0.26 ms for a static binary whose whole body is the `exit` system call: about 0.9 ms of every job, for both arms, was glibc relocating itself (`_dl_relocate_static_pie`, 59% of the `startup_and_input` phase), reading its tunables and probing the cache hierarchy with `CPUID`, which traps to the hypervisor; the same program linked statically against musl without relocations takes 0.36 ms. And the general field element allocates on every operation, so about 45% of the final check's instructions were `malloc`, `free` and `calloc`, for a process that asks for at most 321 KB in 5,700 calls (rho on `n23a0`; 115 KB in 1,353 calls for IC on `n13a0`), keeps at most 208 KB live and lives a few milliseconds. Round 0012 ([pre-registration](ROUND12-single-target.md), [report](../runs/round-0012/REPORT.md), `--objective rho`, seed 2026091612) selects the `x86_64-unknown-linux-musl` target with `crt-static` and `relocation-model=static` ([patch](round12-musl.patch), the cargo config only) and gives the worker a global allocator that bump-allocates from a 64 MiB zero region in the executable's `.bss`, frees only the most recent block, grows the most recent block in place and falls back to the system allocator if the region is exhausted ([patch](round12-arena.patch), three unit tests); the library source is byte for byte the round-0011 library, and rho runs in the same process with the same start-up and the same allocator. All four executables reproduce the 304 frozen round-0011 fixture outputs byte for byte. All 1,488 trials verified, every IC job's logarithm certified; retained. The IC-only challenger `arena_fastio` ([patch](round12-fastio.patch), round 0011's serialisation change rebased) measured 0.9541 [0.9487, 0.9602] of the incumbent's instructions on both final stages but a native upper limit of 1.0271 on confirmation (1.0078 point, cells 0.999–1.035; replay 1.0000 [0.9882, 1.0144]), so it was **not** promoted a second time: a 4.6% instruction saving in serialisation that the clock does not see. The two ablation controls behaved as pre-registered: `arena_glibc` (the arena worker on the round-0011 glibc link) measured 1.0589 [1.0345, 1.0864] of the incumbent's instructions and 1.7657 [1.6619, 1.8727] of its native time at selection, and `musl_sysalloc` (the musl executable on its own allocator) 1.4157 [1.3438, 1.4815] and 1.2278 [1.1551, 1.2882] at development; neither lever is enough alone, and the allocator is what makes the musl link affordable in instructions.

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / rho (confirmation) | 0.4369 | [0.4069, 0.4815] | 0.7988 | [0.7668, 0.8319] |
| winner / rho (replay) | 0.4369 | [0.4069, 0.4815] | 0.7914 | [0.7558, 0.8270] |

Per-cell winner/rho: instructions 0.404–0.530, native time 0.755–0.857 on confirmation and 0.741–0.857 on replay. **`beats_rho_strict` holds** on every cell and both stages, so round 0012 is the strict-win record with the same pipeline. Median complete job on confirmation: 1.67 million instructions and 1.155 ms for the winner against 3.94 million and 1.448 ms for rho (round 0011, same clock and host: 2.14 million and 1.98 ms against 4.71 million and 2.30 ms). The round moved costs both arms pay equally, so the native ratio moved more than either arm's own solve: `startup_and_input` fell from 0.085M to 0.011M instructions on every cell and the final check from 0.294M to 0.178M on `n13a0`. What remains of the IC job is now mostly its own: on `n13a0` the factor base and its tables (0.248M of 0.793M), the final check (0.178M), the curve construction (0.144M) and the log certificate (0.105M); on `n23a0` the factor base (0.821M of 2.98M), relation collection (0.741M), the final check (0.534M) and the construction (0.476M).

## Round 0013: one inversion per scalar product, and what it does to rho

A line-level profile of the round-0012 winner put most of the two shared phases in field inversions: every scalar product in the library — the general `scalar_mul` of the final check and the single-word `FastCurve::mul` of the generator search — was an affine double-and-add with one inversion per step. The same profile of the rho arm showed something the earlier rounds had not looked for: rho's solve begins every restart by drawing its jump table and its walk starts as `[a]G + [b]Q`, two of those affine products per jump and per walk, and on this panel that setup, not the walk, was most of what rho paid (6.39M of its 7.47M instructions on `n23a0`). Round 0013 ([pre-registration](ROUND13-single-target.md), [report](../runs/round-0013/REPORT.md), `--objective rho`, seed 2026091613) moves both products to López–Dahab projective coordinates with one inversion per product ([patch](round13-ld.patch) on the round-0012 source, shared by every arm including rho); every result is the same affine point, and the affine double-and-add is kept as the reference that tests compare against on every catalogue curve and on nine Koblitz cells. The executables reproduce all 304 frozen round-0012 fixture outputs byte for byte. All 1,440 trials verified, every IC job's logarithm certified; retained. **Rho gained more than the IC arm**, as pre-registered: the median complete job on confirmation is 1.28 million instructions and 1.156 ms for the winner against 1.72 million and 1.270 ms for rho (round 0012, same clock and host: 1.67 million and 1.155 ms against 3.94 million and 1.448 ms). The IC arm's own phases did not move; its shared phases fell (`final_verification` 0.178M → 0.100M and `curve_and_targets` 0.144M → 0.117M on `n13a0`), and rho's solve halved (`rho_solve` 6.39M → 2.89M on `n23a0`). The margins of rounds 0007–0012 therefore included rho's affine setup, and the ratios below are the ones to quote.

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / rho (confirmation) | 0.7697 | [0.7148, 0.8370] | 0.9098 | [0.8894, 0.9313] |
| winner / rho (replay) | 0.7697 | [0.7148, 0.8371] | 0.9206 | [0.8983, 0.9480] |

Per-cell winner/rho: instructions 0.684–0.860, native time 0.885–0.927 on confirmation and 0.899–0.965 on replay. **`beats_rho_strict` still holds** on every cell and both stages — the pre-registration expected `n23a0` to be lost natively from a two-fixture development estimate, and the sixty-fixture tournament says it is held at 0.927 and 0.965 — so round 0013 is the strict-win record under the corrected baseline. The IC-only challenger `ld_ic` (the orbit name of an abscissa found from the longest circular run of zeros, [patch](round13-canon.patch), plus the serialisation change of rounds 0011 and 0012, [patch](round13-fastio.patch)) measured 0.9221 [0.9106, 0.9341] of the incumbent's instructions on both final stages, inside the no-regression line, but its native upper limit was 1.0443 on confirmation (1.0199 point, cells 1.000–1.044; replay 0.9919 [0.9780, 1.0084]), so it was **not** promoted: an 8% instruction saving the native clock does not see. `ld_canon` alone measured 0.9758 [0.9505, 0.9939] at selection and was not the provisional challenger. What remains of the IC job is now its own: on `n23a0` the collection scans (0.852M of 2.61M), the factor base and its tables (0.821M) and the log certificate (0.270M); the shared phases are 0.437M.

## Round 0014: the worker as one optimisation unit

To find what the native clock spends, a copy of the round-0013 winner was instrumented with a monotonic clock and `getrusage` at every phase boundary. It found that most of the measured wall is not the job: on `n13a0` the phases sum to 324 µs of an 810 µs process wall, the rest being the spawn, the image page-in and the reap that the protocol charges to every trial and that both arms pay. It also found the worst-throughput phase to be a shared one — `curve_and_targets` at 1.10 ns per instruction against 0.17–0.29 for the job's own phases — because it holds the first `blake3` hash, whose runtime CPU-feature detection costs 30–37 µs on this machine, and about 20 of the job's 50 minor page faults. Two cures were tried in development and **dropped**, and the pre-registration records both: removing the runtime CPU probing (`blake3`'s `pure` feature and a compile-time `pclmulqdq`) took 3.4% of the instructions but cost 10 µs more on the clock, because portable `blake3` is dearer than the `CPUID` trap it removes; and pre-faulting with `madvise(MADV_POPULATE_READ|WRITE)` made the phases much faster (324 µs → 252 µs, the faults gone from them) while the two system calls cost about 380 µs, taking the process wall from 810 µs to 1119 µs.

What remained was the image itself. Round 0014 ([pre-registration](ROUND14-single-target.md), [report](../runs/round-0014/REPORT.md), `--objective rho`, seed 2026091614) builds the worker as **one optimisation unit** — fat LTO, a single codegen unit, dead sections dropped at link time ([patch](round14-lto.patch), shared by every arm including rho) — which takes the text from 2,749,393 to 2,067,481 bytes and the data from 167,592 to 111,080 without changing a line of the library or the worker, so every report byte is unchanged and all 304 frozen round-0013 fixture outputs reproduce byte for byte. The round also carries a tests-only fix ([patch](round14-arena-tests.patch)) to the arena tests, which took the process-global bump pointer for granted while the parallel test harness allocated from it too; the measured executable contains none of it. All 1,440 trials verified, every IC job's logarithm certified; retained.

**The round was pre-registered to make the headline ratio slightly worse, and it did.** Cross-crate inlining made both arms faster and took slightly more of rho's longer instruction stream than of the IC arm's, exactly as registered: the median complete job on confirmation is 1.27 million instructions and 1.143 ms for the winner against 1.64 million and 1.222 ms for rho (round 0013, same clock and host: 1.28 million and 1.156 ms against 1.72 million and 1.270 ms).

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / rho (confirmation) | 0.7818 | [0.7200, 0.8601] | 0.9424 | [0.9126, 0.9715] |
| winner / rho (replay) | 0.7818 | [0.7200, 0.8601] | 0.9284 | [0.9050, 0.9624] |

Per-cell winner/rho: instructions 0.712–0.913, native time 0.895–0.972 on confirmation and 0.907–0.989 on replay. **`beats_rho_strict` holds** on every cell and both stages, but at worse ratios than round 0013, so **the strict-win record stays with round 0013**: a round that makes both arms faster is not thereby a better result against rho, and the record tracks the ratio, not the calendar. The IC-only challenger `lto_ic` ([orbit naming](round14-canon.patch) from the longest circular run of zeros plus the [serialisation change](round14-fastio.patch)) measured 0.9226 [0.9108, 0.9347] of the incumbent's instructions on both final stages but a native upper limit of 1.0163 on confirmation (0.9947 point; replay 1.0022 [0.9825, 1.0256]), so it was **not** promoted for the fourth round running: an 8% instruction saving that four separate clocks have failed to resolve. `lto_canon` alone measured 0.9772 [0.9544, 0.9940] at selection with a native upper limit of 1.08 and was eliminated there.

## Round 0015: the scan, and what the native clock can see

Round 0015 ([pre-registration](ROUND15-single-target.md), [report](../runs/round-0015/REPORT.md), `--objective rho`, seed 2026091615) changed **nothing shared**: the baseline is the round-0014 winner source with an identical source-manifest hash, so the incumbent is the same executable measured again on fresh fixtures, and everything the challengers show is their own.

The round was designed from a line-level profile of that winner, which put the IC arm's own work on the largest cell in the block scan (1.10–1.76M instructions on `n23a0`) and the pair-table build (0.82M), and divided one scanned rest into the batched inversion of the block's denominators (25%), naming the rest's orbit (29%), and the field's reduction (15%). Two candidate changes were tried in development and **dropped, with their numbers recorded before the run**: reduction by folding a sparse modulus rather than by byte tables is **22–47% worse** (every irreducible here is a trinomial or pentanomial, but the table path needs only two or three lookups at these degrees, while folding a pentanomial costs four shift-mask-xor triples twice over), and indexing those tables through a slice of exactly 256 words so no bound is checked is a wash. The pair table's row count was also checked against measurement rather than assumption: building costs about 680 instructions per stored sum against about 380 per scanned rest, where the code's model assumes them equal, and the row count that model picks is still the one that minimises the expected total.

What survived went into two IC-only challengers, both byte-identical to the incumbent on all 304 frozen round-0014 fixtures: `scan` ([patch](round15-scan.patch)), the orbit name found from the longest circular run of zeros rather than by trying all `n − 1` rotations, together with the block scan keeping only each rest's **abscissa** and computing its ordinate for the one rest in a few hundred whose orbit is stored; and `scan_io` ([patch](round15-fastio.patch) on top), which adds the report serialisation change. All 1,440 trials verified, every IC job's logarithm certified; retained.

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / rho (confirmation) | 0.7372 | [0.6628, 0.8326] | 0.9246 | [0.9034, 0.9503] |
| winner / rho (replay) | 0.7372 | [0.6628, 0.8326] | 0.9157 | [0.8886, 0.9564] |

Per-cell winner/rho: instructions 0.646–0.901, native time 0.903–0.964 on confirmation and 0.885–0.988 on replay. **`beats_rho_strict` holds** on every cell and both stages. Median complete job on confirmation: 1.29 million instructions and 1.084 ms for the winner against 1.85 million and 1.186 ms for rho.

**The round's finding is about the gate, and it was predicted in advance.** `scan_io` reached **0.9059 [0.8802, 0.9272]** of the incumbent's instructions — the largest IC-only step this campaign has measured, 0.861 on `n23a0` — and was still **not promoted**, because its native upper limit was 1.0073 on confirmation and 1.0105 on replay. `scan` alone measured 0.9642 at selection with a native upper limit of 1.0119 and was eliminated there. The pre-registration set out the arithmetic before the run: the IC arm's own phases are about two thirds of its instructions, the job's own execution is about 60% of the measured process wall, so a 10% cut in those phases is about 4% of the wall against a paired interval half-width of 2.2%. Five rounds have now put an IC-only change of 3.5–9% of instructions against this gate and five have failed it. **An IC-only change worth less than about a tenth of the job cannot be resolved by a clock that spends a third of its measurement on process creation**, and future IC-only work should be judged on instructions with the native gate read as a no-regression check, or else pursued in steps large enough to clear the floor.

The same round also measures how much fixture draw moves these ratios: the round-0014 executable, unchanged, gave winner/rho 0.7818 instructions and 0.9424 native in round 0014 and 0.7372 and 0.9246 here. Differences of that size between recent rounds are therefore not rankings, which is why the strict-win record is left with round 0013 rather than moved.

## Round 0016: the panel was doing the work — `beats_rho_strict` fails on eight cells

Round 0016 ([pre-registration](ROUND16-single-target.md), [report](../runs/round-0016/REPORT.md), `--objective rho`, seed 2026091616) changed **no source at all**. Both arms are the byte-identical trees round 0015 sealed — the incumbent from `runs/round-0015/source`, `scan_io` from `runs/round-0015/source_candidates/scan_io/source`. What changed is the **panel**, for the first time since round 0006.

The round was designed from round-0015's own frozen receipts, and it cancelled the two rounds that were queued before spending anything on them. Re-running the evaluator's nested bootstrap over those receipts ([round16_resolution.py](round16_resolution.py), which reproduces the recorded interval to eight decimals) shows only the inner level responds to replication: `scan_io`'s native upper limit goes 1.00727 at twelve fixtures per cell to 1.00069 at a hundred — the standard profile, 8.3x the jobs — and converges to **0.99952**. Replication cannot promote it. Nor can removing measurement floor: subtracting a constant moves the point estimate the right way and the **upper limit the wrong way**, 1.0073 to 1.0113 at 404 us, because shrinking both denominators inflates the per-case log-ratio spread faster than it separates the arms. A floor fix would have flattered the headline and made the challenger gate harder.

That left the outer level, whose variance is `sigma_b^2 / C` exactly. A census of every `(degree, curve_a)` pair the worker admits ([round16_cell_census.py](round16_cell_census.py)) found **seventeen usable cells**, two never run: **n23a1**, whose subgroup of 4,196,903 is twice the largest cell in the panel, and **n31a0** at the highest usable degree. Degrees 21, 25 and 27 have no usable subgroup at either `a`, so 31 is the top of the family, not the next rung. `tournament.py` gained `--cells` and `--holdout-cells`, with defaults that reproduce the rounds 0006-0015 panel exactly; this round ran eight confirmation cells. All 2,142 trials verified.

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / rho (confirmation, 8 cells) | 0.8383 | [0.7479, 0.9420] | 0.9779 | [0.9339, **1.0288**] |
| winner / rho (replay, 8 cells) | 0.8383 | [0.7479, 0.9420] | 0.9482 | [0.9035, 0.9980] |
| winner / rho (confirmation, legacy 5 cells) | 0.7674 | [0.7013, 0.8559] | 0.9341 | [0.9118, 0.9653] |

**`beats_rho_strict` is FALSE**, for the first time since round 0013. The winner loses to rho outright at **n23a1** (instructions 1.0656, native 1.0586) and at **n29a1** (instructions 1.0309, native 1.0905), and natively at n31a0 (1.0187). `rho_parity` still holds, but n29a1's 1.0905 sits just under the 1.10 margin.

The third row is the control that settles what caused it. It is the **same frozen receipts** restricted to the rounds 0006-0015 panel ([round16_legacy_subset.py](round16_legacy_subset.py), validated by reproducing round-0015's recorded confirmation exactly), and on those five cells the strict win still passes on every cell. Same executables, same round, same seed: **the strict win was a property of the panel, not of the algorithm.** Every strict-win record in this document, round 0013 included, is scoped to that five-cell panel and must not be restated at family level.

The pre-registration's four predictions, checked in order. **One: `scan_io` refused again**, predicted "upper limit near 1.004" — observed 1.0034 on replay, while confirmation actually *passed* at 0.9947, the first time the native gate has been cleared in six rounds; the gate requires both. **Two: the non-algorithm remainder stays arm-dependent in every cell — NOT MET.** It held in seven of eight, reversing at n19a1 by 10.1 us against a ~690 us baseline; that is 1.5%, well inside the +-5% per-cell spread the A/A control shows for an executable compared against itself, so it is not evidence against the certificate-size explanation, but the prediction said "every cell" and it failed. The mechanism is supported by magnitude instead: the excess grows with cell size, +57.7, +77.6 and +84.4 us at n23a1, n29a1 and n31a0. **Three: the winner still beats rho panel-wide but by less than 0.9246 native** — confirmed and exceeded, 0.9779, far enough that the upper limit crossed one. **Four: the margin narrows with subgroup size** — direction confirmed, predictor wrong: n29a1 has a *small* subgroup (42,457) at high degree and is the worst cell of the eight.

That last miss is the round's most useful output. The factor base is sized `6 * degree` and so does not grow with the subgroup, while rho's cost grows as its square root. n23a1 has twice n23a0's subgroup at the same degree and the same 138-point factor base, which halves the chance a random point decomposes — and it is exactly where the winner flips from 0.898 to 1.066. Sizing the factor base by subgroup order rather than by degree is the open question round 0016 hands forward.

Median complete job on confirmation: 1.95 million instructions and 1.262 ms for the winner, 1.72 million and 1.242 ms for `scan_io`, against 2.00 million and 1.286 ms for rho.

## Interpretation

Every ratio uses a fresh matched rho run in the same round. The 16-target panel charges all setup once to the complete job and solves every target; it is separate from the single-target result, and no ratio combines the two panels. Rho uses the existing per-target solver API on the same constructed curve. Additional cross-target rho optimizations, and a rho specialised like the round-0006 winner, have not been measured here.

The parity rule was declared before measurement: both candidate/rho upper paired 95% limits and every curve-cell ratio must be at most 1.10, in instructions and native time, on confirmation and replay. The full decisions contain the replay evidence.

**Every strict-win statement above is scoped to the five-cell panel that produced it.** Rounds 0006 to 0015 all ran the cells `n13a0, n17a1, n19a0, n23a0` with `n19a1` added in confirmation and replay — a panel that was inherited, never chosen. Round 0016 ran the same executables on eight cells and `beats_rho_strict` came back false, with the winner losing to rho outright at `n23a1` and `n29a1`; the same receipts restricted to the legacy five still pass. Read the strict claims as "strictly below rho on those five cells", which is what was measured, and not as a statement about the Koblitz family or about every curve the worker admits.

These are implementation improvements in fixed-compiler Valgrind amd64 guest instructions and matched native wall time. Kernel/device and external-audit work are outside the instruction count. No arithmetic-complexity, broader-family, or cryptographic-size claim follows. The K-instruction rank floor is deliberately weak.

The successful mechanisms are exact arithmetic substitutions, folded pair-table construction, smaller relation batches, cheaper descent initialization/checking and, in round-0006, a single-word implementation of the whole pipeline. Every returned scalar still passes general worker verification and the independent Python checker.

## Reproduce and review

- [16-target winner source/configuration](WINNER.json) and [cumulative source patch](WINNER.patch).
- [Single-target winner record](WINNER-single-target.json) (round-0008: the round-0007 `tiny2` plus the descent certificate, [patch](round8-tiny2_cert.patch), [round8_candidates.py](round8_candidates.py)).
- Round-0011: the [word-level general arithmetic patch](round11-wordfield.patch) applied to the round-0010 source for both arms, the [serialisation challenger](round11-fastio.patch), and [round11_candidates.py](round11_candidates.py).
- Round-0012: the [musl static link](round12-musl.patch) and the [arena allocator](round12-arena.patch) applied to the round-0011 source for both arms, the [rebased serialisation challenger](round12-fastio.patch), and [round12_candidates.py](round12_candidates.py) (the two ablation controls are the same patches taken one at a time).
- Round-0013: the [López–Dahab scalar products](round13-ld.patch) applied to the round-0012 source for both arms, the IC-only [orbit naming](round13-canon.patch) and [serialisation](round13-fastio.patch) challengers, and [round13_candidates.py](round13_candidates.py).
- Round-0014: the [one-optimisation-unit build](round14-lto.patch) and the [arena test fix](round14-arena-tests.patch) applied to the round-0013 source for both arms, the IC-only [orbit naming](round14-canon.patch) and [serialisation](round14-fastio.patch) challengers, and [round14_candidates.py](round14_candidates.py).
- Round-0015: the IC-only [scan](round15-scan.patch) and [serialisation](round15-fastio.patch) challengers applied to the unchanged round-0014 source, and [round15_candidates.py](round15_candidates.py).
- Round-0016: no source change in either arm; the panel widened to eight cells. [round16_resolution.py](round16_resolution.py) (resolution budget and floor-removal tables from round-0015 receipts), [round16_cell_census.py](round16_cell_census.py) (the seventeen usable cells), [round16_floor_probe.py](round16_floor_probe.py) with [its patch](round16-floor-probe.patch) (wall decomposition), [round16_legacy_subset.py](round16_legacy_subset.py) (the legacy five-cell verdict from these receipts), and [round-0016-single-candidates.json](round-0016-single-candidates.json).
- Round-0010: the [lean-job patch](round10-lean.patch) applied to the round-0009 source for both arms (static link via `.cargo/config.toml`, C entry point, one CPUID), the [feature-cache ablation](round10-lean_stdprobe.patch), and [round10_candidates.py](round10_candidates.py).
- Round-0009: the [shared-construction patch](round9-fastcurve.patch) applied to the round-0008 source for both arms, and [round9_candidates.py](round9_candidates.py).
- Round-0007: the [patch against the round-0006 winner source](round7-tiny2.patch), the [row-rule](round7-tiny2_rows.patch) and [arithmetic](round7-tiny2_arith.patch) ablations, and [round7_candidates.py](round7_candidates.py) to rebuild the trees.
- Round-0006: the [parity patch against the round-0005 winner source](round6-tiny.patch), the [serialisation-only](round6-fast_report.patch) and [full-table](round6-tiny_fulltable.patch) controls, and [round6_candidates.py](round6_candidates.py).
- [Operation plan](PLAN.md), [single-target successor plan](ROUND4.md), [batch plan](ROUND5.md), [single-target parity pre-registration](ROUND6-single-target.md), [beat-rho pre-registration](ROUND7-single-target.md), [certificate re-measurement pre-registration](ROUND8-single-target.md), [leaner shared job pre-registration](ROUND9-single-target.md), [lean job and faithful clock pre-registration](ROUND10-single-target.md) and [general arithmetic in words pre-registration](ROUND11-single-target.md).
- [Operating guide](../OPERATIONS.md) and [controller tests](controller-tests.json).
- [Arithmetic/table equivalence tests](preflight-next-tests.log) and [installed skill validation](skill-validation.json).
- [Retained failed build attempt](../runs/round-0003/prepare_failure.json); no measurements came from it.

Production library defaults were not changed. The frozen source, worker, configurations, raw profiles, receipts and replay evidence are retained under each linked round.
