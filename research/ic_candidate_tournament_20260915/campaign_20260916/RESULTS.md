# Continued IC tournament results

Batch parity verdict: **True**; selected implementation **combined_descent**.

Completed three new tournaments with **4,704 profiled trials**, each paired with a fresh native run, followed by round 0006 (1,680 further profiled trials with paired native runs, incumbent retained, index-calculus admission enforced by a per-target descent certificate). Every listed round passed its independent artifact/correctness audit. The 36-trial batch screen is separate development evidence.

## Profiled instruction cost relative to matched rho

| Round | Targets per cold job | Selected candidate | Instructions / rho | Paired 95% interval |
|---|---:|---|---:|---|
| [round-0002](../runs/round-0002/REPORT.md) | 1 | batch16 | 13.0349 | not calculated in original protocol |
| [round-0003b](../runs/round-0003b/REPORT.md) | 1 | combined_batch8 | 7.1494 | [6.563380843095217, 7.844932268218591] |
| [round-0004](../runs/round-0004/REPORT.md) | 1 | folded_lift_batch4 | 3.0429 | [2.688274518239304, 3.4833783647230367] |
| [round-0005-batch16](../runs/round-0005-batch16/REPORT.md) | 16 | combined_descent | 0.7195 | [0.6374471292425982, 0.8219836239443673] |
| [round-0006-batch16](../runs/round-0006-batch16/REPORT.md) | 16 | incumbent retained (combined_descent + descent certificate) | 0.7259 | [0.6410, 0.8372] |

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

## Interpretation

Every ratio uses a fresh matched rho run in the same round. The 16-target panel charges all setup once to the complete job and solves every target; it is separate from the single-target result. Rho uses the existing per-target solver API on the same constructed curve. Additional cross-target rho optimizations have not been measured here.

The parity rule was declared before measurement: both candidate/rho upper paired 95% limits and every curve-cell ratio must be at most 1.10, in instructions and native time, on confirmation and replay. The full decisions contain the replay evidence.

These are implementation improvements in fixed-compiler Valgrind amd64 guest instructions and matched native wall time. Kernel/device and external-audit work are outside the instruction count. No arithmetic-complexity, broader-family, or cryptographic-size claim follows. The K-instruction rank floor is deliberately weak.

The successful mechanisms are exact arithmetic substitutions, folded pair-table construction, smaller relation batches, and—where selected—cheaper descent initialization/checking. Every returned scalar still passes general worker verification and the independent Python checker.

## Reproduce and review

- [Complete winner source/configuration](WINNER.json) and [cumulative source patch](WINNER.patch).
- [Operation plan](PLAN.md), [single-target successor plan](ROUND4.md), and [batch plan](ROUND5.md).
- [Operating guide](../OPERATIONS.md) and [controller tests](controller-tests.json).
- [Arithmetic/table equivalence tests](preflight-next-tests.log) and [installed skill validation](skill-validation.json).
- [Retained failed build attempt](../runs/round-0003/prepare_failure.json); no measurements came from it.

Production library defaults were not changed. The frozen source, worker, configurations, raw profiles, receipts and replay evidence are retained under each linked round.
