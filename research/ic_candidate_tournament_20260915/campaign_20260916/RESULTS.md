# Continued IC tournament results

Batch (16-target) parity verdict: **True**; selected implementation **combined_descent**.
Single-target verdict: **strictly below rho in both metrics**; selected implementation **tiny2** (round-0007), after **tiny_batch1** reached parity in round-0006.

Completed three new tournaments with **4,704 profiled trials**, each paired with a fresh native run, followed by round 0006-batch16 (1,680 further profiled trials, incumbent retained, index-calculus admission enforced by a per-target descent certificate) and by the single-target rounds 0006 and 0007 (1,584 and 1,488 profiled trials, parity and then a strict win over rho). Every listed round passed its independent artifact/correctness audit. The 36-trial batch screen is separate development evidence.

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

## Interpretation

Every ratio uses a fresh matched rho run in the same round. The 16-target panel charges all setup once to the complete job and solves every target; it is separate from the single-target result, and no ratio combines the two panels. Rho uses the existing per-target solver API on the same constructed curve. Additional cross-target rho optimizations, and a rho specialised like the round-0006 winner, have not been measured here.

The parity rule was declared before measurement: both candidate/rho upper paired 95% limits and every curve-cell ratio must be at most 1.10, in instructions and native time, on confirmation and replay. The full decisions contain the replay evidence.

These are implementation improvements in fixed-compiler Valgrind amd64 guest instructions and matched native wall time. Kernel/device and external-audit work are outside the instruction count. No arithmetic-complexity, broader-family, or cryptographic-size claim follows. The K-instruction rank floor is deliberately weak.

The successful mechanisms are exact arithmetic substitutions, folded pair-table construction, smaller relation batches, cheaper descent initialization/checking and, in round-0006, a single-word implementation of the whole pipeline. Every returned scalar still passes general worker verification and the independent Python checker.

## Reproduce and review

- [16-target winner source/configuration](WINNER.json) and [cumulative source patch](WINNER.patch).
- [Single-target winner record](WINNER-single-target.json) (round-0007), its [patch against the round-0006 winner source](round7-tiny2.patch), the [row-rule](round7-tiny2_rows.patch) and [arithmetic](round7-tiny2_arith.patch) ablations, and [round7_candidates.py](round7_candidates.py) to rebuild the trees.
- Round-0006: the [parity patch against the round-0005 winner source](round6-tiny.patch), the [serialisation-only](round6-fast_report.patch) and [full-table](round6-tiny_fulltable.patch) controls, and [round6_candidates.py](round6_candidates.py).
- [Operation plan](PLAN.md), [single-target successor plan](ROUND4.md), [batch plan](ROUND5.md), [single-target parity pre-registration](ROUND6-single-target.md) and [beat-rho pre-registration](ROUND7-single-target.md).
- [Operating guide](../OPERATIONS.md) and [controller tests](controller-tests.json).
- [Arithmetic/table equivalence tests](preflight-next-tests.log) and [installed skill validation](skill-validation.json).
- [Retained failed build attempt](../runs/round-0003/prepare_failure.json); no measurements came from it.

Production library defaults were not changed. The frozen source, worker, configurations, raw profiles, receipts and replay evidence are retained under each linked round.
