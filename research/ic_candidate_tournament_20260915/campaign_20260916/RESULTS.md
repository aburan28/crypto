# Continued IC tournament results

Batch (16-target) parity verdict: **True**; selected implementation **combined_descent**.
Single-target verdict: **strictly below rho in both metrics on every cell**; selected implementation **tiny2** (round-0007), re-measured with the per-target descent certificate under the merged checker in round-0008 (retained, `beats_rho_strict`), in round-0010 on a lean static executable under an evaluator that no longer forks itself (0.5347 of rho's instructions and 0.8445 of its native time, `beats_rho_strict` on every cell), in round-0011 with the general binary-field arithmetic in words for both arms (0.4796 of rho's instructions and 0.8653 of its native time, `beats_rho_strict` on every cell), in round-0012 as a musl static executable without relocations whose worker serves the heap from an arena, for both arms (0.4369 of rho's instructions and 0.7988 of its native time, `beats_rho_strict` on every cell), in round-0013 with one field inversion per scalar product for both arms, which halved rho's job because its setup is scalar products (0.7697 of rho's instructions and 0.9098 of its native time, `beats_rho_strict` on every cell: the strict-win record under the corrected baseline), in round-0014 with the worker built as one optimisation unit, which made both arms faster again and rho slightly more so (0.7818 and 0.9424), and in round-0015 on that same executable against fresh fixtures: 0.7372 and 0.9246, `beats_rho_strict` still on every cell, after **tiny_batch1** reached parity in round-0006. The last two figures are the same executable measured twice, which is the size of the fixture variation these ratios carry. Rounds 0016–0017 widened the panel to eight cells and round 0017 held `beats_rho_strict` there; round 0018b, the same executable under a fresh seed, did not, failing at `n23a1` alone. **Round 0019 settles it: `both` — round 0018's scan-block ceiling and representative column convention together — is PROMOTED with `beats_rho_strict` on all eight cells in both metrics on both final stages, at 0.6751 of rho's instructions and 0.8848 of its native time, under seed 2026092119 with forty confirmation fixtures at `n23a1` instead of twelve.** Round 0020 then **replicated that win on a second independent seed** (2026092120): `beats_rho_strict` on all eight cells again, at 0.6629 of rho's instructions and 0.8832 of its native time, with `n23a1` at 0.7948. Every eight-cell strict claim is scoped to its seeds as well as its panel, and the tournament classifies both rounds `engineering` rather than an advance: `S` fell, the ratio to the boundary at the largest cell did not.

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

## Round 0017: the certificate was the cost — `beats_rho_strict` on eight cells

Round 0017 ([pre-registration](ROUND17-single-target.md), [report](../runs/round-0017/REPORT.md), `--objective rho`, seed 2026091717, the round-0016 eight-cell panel) changed the **certificate**, not the algorithm. `orbits` is the round-0015 `scan_io` source with the factor base named by its orbit representatives ([patch](round17-orbits.patch)) and the checker regenerating the flat base with its own Frobenius and negation ([oracle patch](round17-oracle-orbits.patch), additive: both formats verify, both in one report is refused). `scan_io` rode unchanged as a control; `orbits_rows` added a re-weighted pair-table row rule ([patch](round17-rows.patch)). All 2,340 trials verified, every IC job's logarithm certified. **Promoted: `orbits`.**

| Comparison | Instructions | Paired 95% interval | Native time | Paired 95% interval |
|---|---:|---|---:|---|
| winner / rho (confirmation, 8 cells) | 0.6702 | [0.6050, 0.7443] | 0.8999 | [0.8755, 0.9266] |
| winner / rho (replay, 8 cells) | 0.6702 | [0.6050, 0.7443] | 0.8918 | [0.8705, 0.9166] |
| winner / rho (confirmation, legacy 5 cells) | 0.6424 | [0.5846, 0.7094] | 0.8882 | [0.8680, 0.9094] |
| winner / incumbent (confirmation) | 0.8319 | [0.8109, 0.8513] | 0.9417 | [0.9172, 0.9640] |
| winner / incumbent (replay) | 0.8319 | [0.8109, 0.8513] | 0.9287 | [0.9091, 0.9459] |

**`beats_rho_strict` holds on all eight cells**, on confirmation and replay, in both metrics — the first time on the panel that refuted round 0013's five-cell claim. Per-cell winner/rho: instructions 0.574 (n31a0) to 0.835 (n29a1); native 0.865 to 0.950 on confirmation and 0.858 to 0.941 on replay, the two cells that lost in round 0016 now at n23a1 0.950 / 0.933 and n29a1 0.947 / 0.941. Median complete job on confirmation: 1.54 million instructions and 1.139 ms for the winner against 2.12 million and 1.240 ms for rho.

The round first killed the hypothesis round 0016 handed forward — the base is quantised to whole orbits and the default already sits at the smallest size the sampler yields; sixteen orbits is worse at every cell ([sweep](round17_base_sweep.py)) — and two accounts of the native remainder: not page faults (rho takes more), and the stdout write target unresolvable in the noise. What it did establish is that the factor base as decimal strings was 9.7–15 KB of a 10–15.5 KB report and about 250,000 instructions of formatting at n29a1, 13% of the whole job there, for 464 points the checker can regenerate from eight. The format was admitted only after expanding 24 round-0016 reports byte for byte and verifying 16 in both formats to identical base hashes, relation counts and ranks.

Against its six pre-registered predictions, in order: **(1) confirmed** — on all 18 fixtures where both arms ran, `orbits` and `scan_io` produced identical relations, column logs, solutions and base hash, at 1,648 against 8,444 bytes. **(2) not met as written** — it predicted the native gain over `scan_io` would exceed the instruction gain at n29a1 and be nil at n23a1; n29a1 is a holdout cell, so the control never ran there, and at every development and selection cell the native gain was *smaller* than the instruction gain (n23a1: 0.967 instructions, 0.982 native; n31a0 no native gain at all). The output-size account of the arm-dependent remainder is weakened, not supported. **(3) confirmed** — n29a1 flipped, 0.947. **(4) FALSIFIED** — n23a1 was predicted not to flip and the strict gate to fail there; it flipped (0.950 / 0.933) and the gate passed everywhere. **(5) half and half** — `orbits_rows` stored two rows at n23a1 and was identical to `orbits` to four decimals at every other cell, exactly as the weighted rule's arithmetic said; but at n23a1 two rows cost **1.272×** the instructions of three, not 0.973×: the rule's build side was right and its scan side badly wrong, fewer rows needing far more probes than its coverage term allows for. It was not promoted. **(6) confirmed** — promotion over the incumbent on both stages.

Prediction 4's failure is the round's honest caveat. Development measured `orbits` at 1.005 of `scan_io`'s wall at n23a1 on four fixtures, and the cell flipped on twelve fresh ones under a new seed; the margin there (0.950 confirmation, upper limit for the panel 0.9266) is real on these fixtures and larger than the ±4% per-cell spread the A/A control shows, but cross-round fixture variation at a single cell is of the same order, and no reader should take 0.950 at n23a1 as a stable per-cell number. What is stable is the panel: eight cells, both metrics, both stages, upper limits 0.7443 and 0.9266.

## Round 0018: two levers that work, and a strict win that does not replicate

Round 0018 ([pre-registration](ROUND18-single-target.md), [report](../runs/round-0018b/REPORT.md),
`--objective rho`, seed 2026091818, the round-0016/0017 eight-cell panel) took
the n23a1 collection that round 0017 handed forward. It killed the hypothesis
it was given, measured two levers it found instead, ran twice because the
first run used the wrong baseline, and ended by refuting a claim of its own
campaign. **Retained: the incumbent.**

### What was killed before anything was built

Round 0017 found its `orbits_rows` arm badly wrong — two pair-table rows cost
1.272x three at n23a1 against a predicted 0.973x — which invites the
conclusion that the shipped row rule is wrong in the same direction. It is
not. Sweeping every row count at every cell
([round18_row_sweep.py](round18_row_sweep.py), with
[a probe patch](round18-rows-probe.patch) making the count settable and the
probe first checked byte-identical to the frozen worker) puts the shipped
choice at the measured optimum at **seven of eight cells, n23a1 included**;
the one miss is n23a0, where three rows beat two by 1.0%, inside the A/A
control's per-cell spread. The sweep reproduces round 0017's number from the
other side: two rows cost 1.307x three at n23a1 here against 1.272x there.
Cost of finding out: one probe build and about fifteen minutes.

### The two levers

`decompose` scans the base in blocks, each paying one batch inversion and
discarded on an early exit, so the rule sizes a block at about one expected
witness, `chunk = (1/hit + 1).clamp(8, 64)`. That ceiling binds at exactly two
cells — the rule asks 71 at n23a0 and 102 at n23a1 — and the obvious reading,
that the cap is too low, is wrong: at n23a1 the block the rule asks for costs
**1.040x** what 16 costs ([round18_block_sweep.py](round18_block_sweep.py)).
A single fixed block is a wash ([round18_block_grid.py](round18_block_grid.py):
the best constant on the development cells is 12, at 0.9910, and it makes
n13a0 and n29a1 worse), so **`block`** ([patch](round18-block.patch)) leaves the
rule's shape alone and lowers the ceiling to 16, chosen on the development
cells alone ([round18_ceiling_choice.py](round18_ceiling_choice.py)) and
recorded as marginal: 12, 16 and 24 sit within 0.1% of each other. What is not
marginal is that anything at or below 32 beats 64, and that a ceiling can only
lower a block, so no cell can regress.

**`column`** ([patch](round18-column.patch)) drops a cofactor multiplication
applied to a point that is already in the subgroup. `TinyIc::new` projects each
sampled point with `[h]` and closes its abscissa under Frobenius, so every
orbit representative is already in the order-`r` subgroup; forming the column
as `[h]rep` scaled every column by a constant and the row (`h*a`) and the
descent (`invmod(h*b)`, `sum - h*a`) each carried a matching constant.
Dropping all four leaves the same logarithms and saves `K` scalar
multiplications by the cofactor per job — nothing where `h` is 2 or 4, most of
the base phase where `h` is 12,646 (n29a1) or 1,492 (n31a0). The report
declares `column_convention: representative`, which the checker amendment
([oracle patch](round18-oracle-convention.patch), additive) reads: a base point
is then located by itself and its row reads `sum coeff*log == a`; a report
without the key is read exactly as before; an unknown value is refused. Under
the new convention every base point must carry a column, where the old one
excused a point whose `[h]` image is the identity — strictly stronger. The
amendment was admitted only after it was attacked at n29a1, where the two
conventions differ by a factor of 12,646: stripping the label, relabelling a
legacy report as `representative` and relabelling in the other direction are
each refused, as are a column logarithm, a relation scalar, a recovered
logarithm and an orbit representative's coordinates each moved by one; and 280
frozen round-0016 and round-0017 profiles sampled at random still verify
unchanged.

**`both`** is the two together. All 2,340 trials verified, every logarithm
certified, audit status VERIFIED over 2,340 receipts and 452 source files.

| both / incumbent | confirmation | replay |
|---|---|---|
| instructions | 0.9652 [0.9448, 0.9828] | 0.9652 [0.9448, 0.9828] |
| native wall | 0.9613 [0.9288, 0.9895] | 0.9739 [0.9475, 0.9973] |

Every cell inside 1.10; largest gains at n31a0 (0.9203) and n29a1 (0.9236),
smallest at n17a1 (0.9930). Against rho, `both` is below one at every cell in
both metrics **except n23a1 on instructions, at 1.0112**, and that single cell
is why nothing was promoted: the objective is `rho`, and a challenger must beat
matched rho strictly everywhere.

### The result that matters, which the round did not set out to find

The incumbent here is **byte-identical to round 0017's promoted winner**. On
this seed it measures 0.7115 [0.6206, 0.8188] against rho with n23a1 at
**1.0231**. Round 0017 measured the same executable at 0.6702 [0.6050, 0.7443]
with n23a1 at **0.950**. Instructions are deterministic, so the whole of that
difference is the fixture draw.

**Round 0017's eight-cell strict win does not replicate under a fresh seed**,
and n23a1 is where it breaks. Round 0017's own entry above had already said
that no reader should take 0.950 at n23a1 as a stable per-cell number; this is
that warning cashed out. `strict_win_record_eight_cells` still names round
0017, because those receipts are immutable and say what they say, but the
record is **seed-dependent** and is not a property of the executable or of the
panel. `rho_parity` does still hold here, in both metrics on both stages.

### Predictions, in the order they were checked

**2 confirmed.** `block` is a no-op where it provably cannot act: 1.00026,
1.00018, 1.00013 at n13a0, n17a1, n19a0 against a prediction of 1.000 within
0.001. **3 confirmed** — `both` over the incumbent at 0.9652, inside the
predicted 0.95–0.975, with the largest gains near 0.92 at n31a0 and n29a1 and
the smallest above 0.98 at n17a1. **4 confirmed on what could be tested**:
among the six development cells the ranking by `column`'s gain is exactly the
ranking by cofactor, all three `h=4` cells ahead of both `h=2` cells.
**6 confirmed, and its caveat earned** — both gates passed, but the replay
native upper limit cleared one by 0.27% (0.9973) and this round's own A/A
control returned a native panel interval of [0.9609, 1.034] against the
challenger's point estimate of 0.9613; the instruction side carries no such
caveat, reproducing to 1.0000034 [0.9999975, 1.0000117]. **7 falsified** —
`beats_rho_strict` does not hold on all eight cells.

**5, and part of 4, were not testable, and that is a fault in the
pre-registration rather than a result.** Single-lever arms run only in
development and selection, which exclude the holdout cells, so `block` at
n19a1 and `column` at n29a1 — the two cells those predictions name — were never
measured. Round 0017's prediction 2 hit the same wall for the same reason. A
prediction about a holdout cell has to be written about an arm that reaches
confirmation.

### The run that was thrown away, and kept

The first attempt ([runs/round-0018](../runs/round-0018)) ran all 2,340 trials
and verified every one, against the wrong incumbent: `prepare` was given
`runs/round-0017/source`, which is round 0017's *incumbent*, not its promoted
winner at `runs/round-0017/source_candidates/orbits/source` — the path this
campaign's own `next-proposal-single-target.json` names. Its incumbent worker
hashes to `7ca9953d…` against the winner's `e9f263b8…` and its source carries
no `factor_base_orbits` at all, so every arm-versus-incumbent number in it
conflates this round's levers with round 0017's certificate. Prediction 2 is
what caught it: the three cells where `block` cannot act read 0.876, 0.847 and
0.852 instead of 1.000.

The record is kept and superseded rather than deleted. The preflight that was
supposed to prevent exactly this ([round18_preflight.sh](round18_preflight.sh))
checked nine contract fields and all three candidate worker hashes, and never
checked the baseline those arms are measured against; it now verifies that the
incumbent hashes to the round-0017 winner and that the baseline source carries
the orbit certificate, and it refuses on either.

## Round 0019: the cell resolved, and the strict win promoted

Round 0018b left the campaign's central claim undecided rather than dead. Its
`both` arm — the scan-block clamp ceiling at 16 and the representative column
convention — beat the incumbent by 3.5% in instructions and failed the strict
rho gate at exactly one cell, `n23a1`, by 1.1%. Pooling rounds 0017 and 0018b,
24 independent fixtures for the same executable, put that cell's winner/rho
ratio at **0.8913 with a 95% band of [0.7802, 1.0182]**: an estimate straddling
the gate, on a cell whose per-case log spread is 0.329, where twelve fixtures
fail by sampling alone about **9%** of the time.

So round 0019 changed the instrument rather than the algorithm. It carries **no
new arm**: `tournament.py prepare` gained `--confirmation-cases cell=count`,
which is additive, raises only — a count below the profile floor is refused —
and leaves the estimator and both gates untouched. The allocation was computed
from frozen prior rounds alone by `round19_allocate.py` as the smallest count
whose one-sided failure probability is at most 1% in both metrics:
**forty fixtures at `n23a1`, twelve everywhere else, 124 cases against the flat
panel's 96, no cell measured less than round 0018b measured it.** More fixtures
move a cell's estimate toward its true value in whichever direction that lies,
so the allocation can resolve a cell but cannot buy it a pass.

**Result: `both` PROMOTED. `beats_rho_strict` true, both metrics, every cell,
both final stages.** Seed 2026092119, 2,646 trials, audit VERIFIED over 2,646
receipts and 452 source files.

| winner (`both`) / rho | confirmation | replay |
|---|---|---|
| instructions | 0.6751 [0.6298, 0.7265] | 0.6751 [0.6297, 0.7265] |
| native wall | 0.8848 [0.8624, 0.9098] | 0.8882 [0.8676, 0.9127] |

Worst cells: `n23a1` at **0.7498** in instructions, `n29a1` at **0.9427**
natively. Against the incumbent, `both` reads 0.9657 [0.9479, 0.9825] in
instructions and 0.9662 / 0.9787 natively across the two stages.

**The incumbent alone would not have passed.** Its `n29a1` native cell reads
1.0007 — the cell whose cofactor is 12,646, where round 0018's `column` lever
removes the redundant cofactor multiplication. The challenger is not merely
cheaper on average; it is what carries the panel's worst native cell under one.

### Predictions: four confirmed, two falsified, one split

**1 confirmed** — `both`/incumbent 0.9657 in instructions with both upper
limits below one on both stages. **3 confirmed** — `beats_rho_strict` holds.
**4 confirmed** — every cell's per-case spread landed within ±40% of the
pre-registered column, ratios 0.73 to 1.19; the variance model that sized the
allocation was right at all eight cells. **7 confirmed** — confirmation and
replay agree in instructions to 3.1e-5, as they must when two stages share
fixtures and Ir is deterministic.

**2 falsified.** `both`/rho at `n23a1` was predicted in [0.78, 1.00] with a
point value of 0.881, from rounds 0017 and 0018b pooled. It measured **0.7498**
— below the range, about three standard errors low. The miss is in the *mean*,
not the spread, which prediction 4 confirms held.

**5 falsified.** The per-case spread was predicted monotone in `r` across six
cells; `n19a1` (0.198) exceeds `n31a0` (0.184). One adjacent pair, by 0.014,
with the rest of the chain in order.

**6 split.** The same-degree contrast holds at degree 23 — winner/rho at
`n23a1` over `n23a0` is 1.1609, above the predicted 1.05 — and fails at degree
19, where 0.9870 is below the predicted 1.00. The crossover mechanism is
visible where the subgroup orders are large and is not resolvable at degree 19.

### What prediction 2's miss means, and the control it needed

`round19_seed_variance.py` combines the three rounds' per-cell log ratios by
inverse variance and tests their scatter against chi-square on two degrees of
freedom. Three of eight cells exceed p < 0.05 against 0.4 expected by chance:
`n13a0` (0.043), `n23a1` (0.044), and `n31a0` (0.0001, reading 0.5737, 0.6932,
0.7868 across the three seeds).

That needed a control before it could mean anything, because **`rho` is not a
fixed binary**: the tournament synthesises it from each round's baseline arm,
and round 0017's baseline was the pre-orbits incumbent (`7ca9953d…`) where
rounds 0018b and 0019 both use the orbits winner (`e9f263b8…`). A cross-round
ratio could have been reading a change of executable. Measured directly on
identical fixtures, three cases at each of eight cells, **the two binaries
return the same rho instruction count to within 2 parts in 10,000** (geometric
mean 0.99999, worst cell 1.0002). The orbits change never reached rho's code
path. The columns are comparable — and round 0018b's reading of its own
`n23a1` flip as the fixture draw now rests on a control rather than an
assumption.

**The test over-rejects here, and that matters more than the p-values.** rho's
cost is a collision time, so it is right-skewed with a long tail: a twelve-case
sample standard deviation underestimates the spread and a twelve-case mean is
not normal, and both push chi-square toward rejecting. With three seeds and
n=12 in two of them, a flagged cell is a question to put to a fourth seed, not
an established fact. Round 0019's n=40 column is the only well-resolved mean on
the table.

### How fragile is this win?

Measured on round 0019's own 124 cases, the per-cell failure probability of a
*fresh twelve-case seed* is at most **0.1%**, at `n23a1`; every other cell is
below 0.003% in both metrics. Round 0018b faced 9.1% at that cell. The
difference is not the allocation — it is that `both`'s margin at `n23a1` is now
25% rather than 1.1%, because round 0018b drew a hard fixture set there and
round 0019 did not. The cross-seed combined value for `both` at that cell is
about 0.78, so round 0018b's 1.0112 is the outlier of the three.

That is the honest scope: **the strict win is established for this seed and
panel with the binding cell resolved at forty fixtures, and the cross-seed
evidence puts its margin at roughly 22% rather than at the gate.** It is not
established that a per-cell ratio is a fixed property of the executable; three
cells say otherwise, weakly.

### What it is not

By the AGENTS.md §3 test this is **engineering**, and the tournament's own
decision record classifies it that way. `S` fell; the ratio to the boundary at
the largest cell did not. `round19_model.py` works the solver's own published
cost model: at its optimal factor-base size this pair-table collector is
`Θ(r^{2/3})` where rho is `Θ(r^{1/2})`, so the ratio grows as `r^{1/6}` and
must eventually exceed one. The same-degree contrast measures it rising 16.6%
per doubling of `r` at degree 23, and one further doubling past `n23a1` reads
1.040 — above rho.

The factor base cannot fix that. `round19_base_sweep.py` measured the cost over
the orbit count at every cell: it rises monotonically everywhere, and the base
the panel builds is the cheapest one reachable. That sweep also corrected the
model. Its first version varied the contract's `factor_base.points` from
`2·degree` to `12·degree` and measured a flat line, because
`build_subgroup_orbit_factor_base` samples orbits in batches of eight and stops
at the first rebuild that reaches the target. At degree 23, `points=1`,
`points=138` — the contract's value — and `points=368` all return eight orbits.
**The panel's factor base is seven or eight orbits at every cell and the
contract's `points` parameter has been inert since round 0002.**

## Round 0020: the replication, and what a fourth seed said about the flagged cells

Round 0019 promoted `both` with `beats_rho_strict` on all eight cells. That was
one seed, and this campaign's own failure mode is that a single-seed strict win
need not survive a fresh draw — round 0017 held the gate and round 0018b, the
same executable, did not. Round 0020 is the test round 0017 never got: identical
to round 0019 in panel, allocation, arms, objective and limits, changing only
the seed (2026092120).

The baseline stayed the **incumbent** rather than round 0019's promoted winner.
Promoting `both` into the baseline would make `both`/incumbent identically 1 and
destroy the quantity under test.

**Result: `beats_rho_strict` again, on a second independent seed. `both`
promoted.** 2,646 trials, audit VERIFIED over 2,646 receipts and 452 source
files.

| winner (`both`) / rho | confirmation | replay |
|---|---|---|
| instructions | 0.6629 [0.5916, 0.7426] | 0.6629 [0.5916, 0.7426] |
| native wall | 0.8832 [0.8568, 0.9129] | 0.8791 [0.8512, 0.9093] |

`n23a1` reads **0.7948** in instructions and 0.9322 natively, against 0.7498 and
0.8967 in round 0019 — and 1.0112 in round 0018b on twelve fixtures. Against the
incumbent, `both` reads 0.9638 [0.9440, 0.9822] in instructions and 0.9698
natively.

### Predictions: five confirmed, one reported rather than scored

**1 confirmed** — `beats_rho_strict` holds. **2 confirmed** — 0.6629 inside the
predicted [0.62, 0.74], 0.8832 inside [0.85, 0.92]. **3 confirmed** — `n23a1`
below 0.95 in instructions, at 0.7948. **4 confirmed** — `both`/incumbent 0.9638
inside [0.95, 0.98] with both upper limits below one on both stages. **6
confirmed** — the incumbent's worst native cell is `n29a1` as predicted.

**5 was pre-registered with no expected outcome**, which is what it means for
round 0019 to have recorded its flagged cells as open questions rather than
findings. Over four seeds:

| cell | 0017 | 0018b | 0019 | 0020 | combined | p (3 seeds) | p (4 seeds) |
|:--|--:|--:|--:|--:|--:|--:|--:|
| `n13a0` | 0.7146 | 0.7039 | 0.6932 | 0.7160 | 0.7058 | 0.0431 | **0.0204** |
| `n31a0` | 0.5737 | 0.6932 | 0.7868 | 0.6119 | 0.6612 | 0.0001 | **0.0001** |
| `n23a1` | 0.7765 | 1.0231 | 0.7647 | 0.8122 | 0.8022 | 0.0443 | 0.0953 |

**`n23a1` came off the list.** Adding a second forty-fixture column moved it from
0.044 to 0.095: its apparent seed-dependence was largely round 0018b's
twelve-fixture reading, which is the same diagnosis round 0019 made of the
original gate failure, now confirmed from the other side. `n13a0` tightened
rather than loosened, and `n31a0` is unchanged at 0.0001.

Two corrections to how round 0019 read this. The monotone rise at `n31a0` across
three seeds — 0.5737, 0.6932, 0.7868 — **is not a trend**; the fourth seed reads
0.6119. It is scatter. And the two cells that remain flagged sit at combined
0.7058 and 0.6612, **nowhere near the gate**, so whatever they are, they do not
threaten the strict win. The cell that could have threatened it is the one that
cleared.

The caveat on the test still stands and still matters: rho's cost is
right-skewed, twelve-case sample spreads underestimate, and chi-square
over-rejects here. Two flags out of eight cells at p < 0.05, on four seeds with
n=12 at six of them, is weak evidence of anything.

### Scope, unchanged

Two independent seeds now hold `beats_rho_strict` on the eight-cell panel with
the binding cell resolved at forty fixtures. That makes the result robust to the
fixture draw **on this panel, at these subgroup orders, for this executable**. It
says nothing about larger `r`, where the boundary of round 0019 §3 puts this
collector at `Θ(r^{2/3})` against rho's `Θ(r^{1/2})`. The classification is
**engineering** under AGENTS.md §3, in both rounds' own decision records.

## Round 0021: the crossover, measured — and the ceiling was eight bytes

Rounds 0019 and 0020 hold `beats_rho_strict` on the eight-cell panel under two
seeds, and both classify the gain **engineering** for one reason: the
`Θ(r^{2/3})` against `Θ(r^{1/2})` boundary was *derived*, its only measured
support a same-degree contrast between two panel cells.
[PROBE-degree-ceiling.md](PROBE-degree-ceiling.md) then found the cell that
would settle it — `n37a0`, `r = 230,603,167`, 55× the panel's largest — and
found that the promoted implementation could not run there. This round removes
that obstacle and takes the measurement.
[ROUND21-crossover.md](ROUND21-crossover.md) is the full write-up.

**No tournament ran, nothing was promoted, and `WINNER-single-target.json` is
unchanged.** This is a measurement round against the round-0020 winner source,
not a round of the tournament, so it has no `runs/round-0021/`, no stage
summaries and no evidence pack.

### The crossover

Sixty-four independent fixtures a cell, pooled over two seed streams
(20260922, 4242), every IC report checked by `oracle.py`
([round21_crossover.py](round21_crossover.py)):

| cell | r | IC/rho | 95% band | per-case sd(log) |
|:--|--:|--:|:--|--:|
| `n23a1` | 4,196,903 | **0.831** | [0.772, 0.894] | 0.300 |
| `n37a0` | 230,603,167 | **1.533** | [1.238, 1.899] | 0.872 |

**Both bands exclude one, in opposite directions: the collector crosses rho
between these two subgroup orders.** Measured on cells where both arms complete
and every answer is certified — not extrapolated. The ratio grows as
**r^0.153** across the two, where [round19_model.py](round19_model.py) derives
**r^{1/6} = r^0.167** for the balanced optimum. Two cells fix one rate and no
curvature, and the cells differ in degree as well as in `r`, so the agreement is
worth exactly what a two-point rate is worth — but it is the first measured
support the boundary has had beyond the panel. The `n23a1` value doubles as a
cross-check: 0.831 here against 0.7498 (round 0019) and 0.7948 (round 0020).

### The ceiling was eight bytes

`koblitz_tiny_ic` declared `MAX_DEGREE = 31`. Its field has always been one
`u64`; the ceiling lived in the pair table, which packed each stored sum's
coordinates into `u32`.
[round21-wide-pair-table.patch](round21-wide-pair-table.patch) widens those two
fields and lifts the three bounds that mirrored them — the module constant, the
worker's `(5..=31)` dispatch guard and `oracle.py`'s `Curve.__init__` — all to
61. Lifting one and not the others cost a build to discover, so
[round21_build.sh](round21_build.sh) now checks all three before it hands back a
binary.

The widening costs **0.07% panel-wide**, worst cell 0.15%, and **24 of 24
confirmation fixtures returned the same logarithms and the same factor-base
hash**:

| cell | n13a0 | n17a1 | n19a0 | n19a1 | n23a0 | n23a1 | n29a1 | n31a0 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| widened / promoted | 1.0011 | 1.0008 | 1.0015 | 1.0006 | 1.0003 | 1.0002 | 1.0004 | 1.0004 |

### A test that had been red since round 0018

Building the widened tree surfaced `complete_solve_verifies_in_general_arithmetic`
failing at `n13a0`, and it reproduces **identically on the unmodified round-0020
winner** — it is not this round's. The assertion carried the **cofactor**
convention (`h·(a + b·d)`) after round 0018's `column` patch moved the solver to
the **representative** convention. `oracle.py` was amended at the time; this test
was not, and has been red at every cell whose cofactor is not one ever since,
through the two rounds that promoted the arm that broke it.

It does not invalidate those rounds, for a specific reason: the failing assertion
is an internal cross-check between the tiny path and the general one, while the
assertion in the same test that the recovered logarithm equals the planted scalar
passed throughout, and every published result was verified by `oracle.py` — 2,646
receipts a round — against the planted secret. What it does expose is a real gap:
**the campaign verified through `oracle.py` and the tournament harness and never
ran `cargo test` on a candidate arm's own source.** Fixed here, in the patch.

### Withdrawn

Two claims of mine were wrong and are withdrawn. Before writing any code I
computed that `n37a0` would need 32 orbits and that 8 could not solve it inside
`max_trials`, and I let that shape the design: eight orbits solve it and are the
*best* configuration. `max_trials` bounds relation-producing trials; the quantity
in the cost model is the pair-table scan count. The same error produced the claim
that `n41a0` "needs ~39,922 orbits".

### What is not claimed

**`n41a0` yields no comparison.** At `r = 5.5·10¹¹` IC completes with a large
enough base (104 orbits, 3.8s) and rho does not complete at all, returning in
0.06s against the ~72,600 steps it would need. That is the frozen trial budget
cutting rho off, not rho being out-run, and a budget exhaustion is never evidence
about an algorithm. Recorded so the absence is on the record.

**The factor base was not tuned toward the answer.** The base is swept at
`n37a0` and every size reported — 8 orbits 1.533, 32 orbits 1.815, 40 orbits
1.884, 56 orbits 2.105 on the sixteen-fixture pass. Bigger is monotonically
worse, which is what round 0019's sweep found at all eight panel cells,
reproduced two decades of `r` higher; the headline uses the cheapest
configuration, which is the one the sampler builds by default.

**The sample size was arrived at the hard way.** Three fixtures read `n23a1` at
0.882, then 1.097 on a second draw; sixteen read `n37a0` at 1.782 [1.394, 2.278]
on one stream and **1.116 [0.578, 2.154]** on another — a band containing one,
which is no answer. That is round 0019's finding, that the per-case spread grows
with the cell because rho's collision search is a growing share of its cost,
holding two decades of `r` further out. Sixty-four fixtures bring both standard
errors under 0.11.

**This does not overturn rounds 0019 and 0020.** Their strict win stands exactly
as scoped: eight cells, subgroup orders 2·10³ to 4·10⁶, two seeds. This round
adds the measured cell beyond that scope, where the same collector loses. By
AGENTS.md §3 the campaign's gains remain **engineering**, with a measurement
rather than a derivation now behind the reason why.

## Round 0022: two checks round 0021 did not run, and both failed

[ROUND22-budget-and-curvature.md](ROUND22-budget-and-curvature.md) is the full
write-up. **Round 0021's crossing survives; both of its numbers do not.**

| cell | r | best base | IC/rho | 95% band |
|:--|--:|--:|--:|:--|
| `n23a1` | 4,196,903 | 8 orbits | **0.831** | [0.772, 0.894] |
| `n37a0` | 230,603,167 | 16 orbits | **1.396** | [1.244, 1.567] |
| `n43a1` | 4,644,189,029 | 24 orbits | **2.207** | [1.933, 2.519] |

64 fixtures a configuration over two seed streams, both arms required to
complete, every IC report checked by `oracle.py`
([round22_ladder.py](round22_ladder.py)). No tournament ran and nothing was
promoted.

### rho was charged while it was still running

`round21_crossover.py` checks that the IC arm completed and then measures both
arms, and `instructions()` reads callgrind's `Collected:` line rather than the
worker's JSON — so it cannot tell a solved rho from one that ran out of budget.
**Four of the 64 rho runs at `n37a0` were cut off and charged anyway**; 64 of 64
completed at `n23a1`.

My first correction said this inflated the ratio. That was reasoning rather than
measurement and it was wrong: an `incomplete` rho has exhausted its restarts, so
it did a great deal of work and produced nothing, where a completed run often
finds its collision early. Charging it **over**-charged the denominator. 1.533
was a *lower* bound; with every rho run required to finish, the same
default-base configuration reads 1.763 [1.547, 2.009].

### The factor base was optimal only on the panel

Round 0021 carried round 0019's "bigger is monotonically worse" from the
eight-cell panel to the new cells and reported the sampler's default single
batch. That holds at `n23a1`, where 8 orbits is best and the curve rises to
9.949 at 159 orbits. **It does not hold further out**
([round22_base_sweep.py](round22_base_sweep.py)):

| cell | 8 orbits | 16 | 24 | 32 | 40 | 48 |
|:--|--:|--:|--:|--:|--:|--:|
| `n23a1` | **0.831** | 1.043 | 1.337 | — | — | — |
| `n37a0` | 1.763 | **1.396** | 1.519 | 1.709 | 1.846 | — |
| `n43a1` | 5.794 | 3.006 | **2.207** | 2.261 | — | 2.449 |

At `n43a1` the default costs the candidate **2.6×**. Note the direction: the
default *overstated* the candidate's loss, which is the one direction a result
already going against the candidate will not be challenged on.

### The rate, at the third attempt

| step | r ratio | rate |
|:--|--:|:--|
| `n23a1` → `n37a0` | 55× | r^[0.130, 0.188] |
| `n37a0` → `n43a1` | 20× | r^[0.075, 0.187] |

Each step is a range over every base size in that cell's flat region — the
minimum plus every size whose band overlaps it — because the ladder must not
depend on choosing between sizes the data cannot separate. **The ranges overlap,
so no curvature is resolvable**, and the derived `r^{1/6} = r^0.167` lies inside
both.

The two earlier answers were artifacts. At the default base the rate appeared to
**accelerate** (r^0.188 then r^0.396) — a handicap widening with `r`. On a base
grid that skipped 16 orbits it appeared to **decelerate** on disjoint ranges
(r^[0.180, 0.199] then r^[0.060, 0.120]), with `n43a1`'s minimum on the edge of
the grid, which is how a minimum outside it announces itself. Filling the hole
moved `n37a0`'s best base from 32 orbits to 16 and its ratio from 1.709 to
1.396. Both spurious results were more interesting than the real one.

### What the trial cap is worth, and what the optima say

Raising `max_trials` from 4096 to 65,536 is a protocol change rather than an
extension — it is rho's iterations-per-restart — but
[round22_budget_effect.py](round22_budget_effect.py) bounds it at **14 ppm** on
either arm, both signs, with identical logarithms and factor bases everywhere;
`n23a1` reads 0.831 at both caps.

The model's rate survives and its other prediction does not.
`F = (c·#E·t/k)^{1/3}` predicts 30 and 83 orbits at the two new cells; measured
best bases are 16 and 24, growing about `r^0.157` rather than `r^{1/3}`. These
are not independent predictions — `r^{1/6}` is derived *from* that base — so the
agreement on the rate, at a base measurably not where the derivation puts it, is
a coincidence the model does not explain rather than a confirmation of it.

**The eight-cell panel is untouched.** Rounds 0019 and 0020 were measured at
4096 trials with a base that is optimal there, and both facts still hold. The
classification stays **engineering** under AGENTS.md §3.

## Interpretation

Every ratio uses a fresh matched rho run in the same round. The 16-target panel charges all setup once to the complete job and solves every target; it is separate from the single-target result, and no ratio combines the two panels. Rho uses the existing per-target solver API on the same constructed curve. Additional cross-target rho optimizations, and a rho specialised like the round-0006 winner, have not been measured here.

The parity rule was declared before measurement: both candidate/rho upper paired 95% limits and every curve-cell ratio must be at most 1.10, in instructions and native time, on confirmation and replay. The full decisions contain the replay evidence.

**Every strict-win statement above is scoped to the five-cell panel that produced it.** Rounds 0006 to 0015 all ran the cells `n13a0, n17a1, n19a0, n23a0` with `n19a1` added in confirmation and replay — a panel that was inherited, never chosen. Round 0016 ran the same executables on eight cells and `beats_rho_strict` came back false, with the winner losing to rho outright at `n23a1` and `n29a1`; the same receipts restricted to the legacy five still pass. Read the strict claims as "strictly below rho on those five cells", which is what was measured, and not as a statement about the Koblitz family or about every curve the worker admits. **Round 0018b then failed to replicate that restoration under a fresh seed**: the same executable, the
same panel and the same gates, seed 2026091818 instead of 2026091717, gives winner/rho 0.7115 [0.6206,
0.8188] on instructions with n23a1 at 1.0231 rather than 0.950, and `beats_rho_strict` false. Instructions
are deterministic, so that is the fixture draw alone. Read every eight-cell strict claim below as scoped to
its own seed as well as to its panel.

**Round 0017 then restored the strict claim on all eight cells** — winner/rho 0.6702 [0.6050, 0.7443] instructions and 0.8999 [0.8755, 0.9266] native on confirmation, every cell below one in both metrics on both stages — with the certificate named by orbit representatives. That claim is scoped to those eight cells and to this checker's second certificate format, both stated in the round's pre-registration; it is still not a statement about the family.

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
- Round-0017: `orbits` = the round-0015 `scan_io` source plus the [orbit-representative certificate](round17-orbits.patch), with the [checker amendment](round17-oracle-orbits.patch) applied to `oracle.py`; `orbits_rows` adds the [re-weighted row rule](round17-rows.patch); [round17_candidates.py](round17_candidates.py) rebuilds the trees, [round17_base_sweep.py](round17_base_sweep.py) is the base-size sweep, [round17_measure.py](round17_measure.py) the development measurement, and [round-0017-single-candidates.json](round-0017-single-candidates.json) the registry.
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
