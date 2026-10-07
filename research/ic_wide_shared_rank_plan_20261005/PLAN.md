# Plan: the index-calculus arm on the m = 83 gate curve in `ecbench` (wide shared-rank)

**Status: plan, registered 2026-10-05; no code, no measurement.**  Every
measurement below is pending.  Merging this plan does not close the
workstream; execution lands in linked follow-on PRs that update the
checklist at the end.

## 1. Why

AGENTS.md §8a requires m = 83 as the confidence gate for any index-calculus
improvement meant to transfer to ECC2K-130, and §12 requires new cross-method
measurements to go through `ecbench`.  Two of the three pieces exist:

- the wide `GF(2^m)` Koblitz group for `64 ≤ n ≤ 127` and the width-generic
  strong rho and claw (#1360, merged; `WideInstance`, `WideKoblitz`,
  `WideStrongRho`);
- the strong-rho reference on `icv1-f2m83-tm6151469093347-debefd74` (#1387,
  open: `rho.signed_frobenius_strong_escape`, part A calibrated on m = 53, 59,
  61, part B running).

The third does not: **no index-calculus method runs on a wide curve.**
`ecbench::methods::solve_wide` admits only `rho.signed_frobenius_strong` and
`claw.pair_table`; `dump_factor_base` and the factor-base plug-ins refuse a
wide instance ("no factor-base plug-in runs on a wide (n > 62) curve yet").
The shared-rank method (`ic.shared_rank`, the K8/K16 line of the n37 panels
and of `research/ecbench_n41_n53_shared_rank_20261005`) is built on
`koblitz_fast::FastPoint` (`x, y: u64`) through
`ic_framework::{shared_rank, plugins::CompactOrbitScanBase,
plugins::FrobeniusMitmOracle, stages}`.  Until that path is width-generic,
no IC figure can be read against the m = 83 reference, and the §8a gate
stays unestablished for every candidate on this line.

## 2. Hypothesis and prediction (to be falsified)

**H.** The counted cold index-calculus cost of the shared-rank method
(`compact-orbit-scan`, K8 and K16 columns) grows faster than rho's `√r`
from n = 37 to n = 83, so the counted IC/rho lower-bound quotient, 4.55×
(K8) and 5.03×–14.4× (K16/K42) at n = 37 (`research/ecbench_n37_*`),
is above 10× at n = 83.

**Prediction, from the n37–n61 phase exponents in
`docs/ic/leaderboard.json` and the n41/n53 counted panel when it lands:**
the reusable preparation (base, pair table, rank) dominates and its counted
cost fits `r^{a}` with `a > 1/2` over four or more sizes; the online
target-only cost stays below rho's.  *Falsified if* the counted cold
IC/rho quotient at n = 83 is below 2×, or if the fitted `a ≤ 1/2` over
n ∈ {37, 41, 53, 61, 67, 71, 73, 79, 83}.

Nothing here is a speed claim: counted units are lower bounds while native
work is unpriced, and wall time needs the L2 host that no session has yet
earned (`research/ecbench_n41_n53_shared_rank_20261005/L2_RUNBOOK.md`).

## 3. Frozen inputs

| item | value |
|:--|:--|
| curves | `icv1-f2m67-tm19346764963-82c84cca`, `icv1-f2m71-t48653080717-2cafaea3`, `icv1-f2m73-tm184271214331-9d25cc67`, `icv1-f2m79-tm420247971347-2a24b892`, `icv1-f2m83-tm6151469093347-debefd74` (registry, `K_0` over the registry moduli; m = 83 as §8a fixes it: `z^83 + z^45 + z^2 + z + 1`, order `4 · 2417851639230796216685689`) |
| arms | `ic.shared_rank` with `compact-orbit-scan:columns=8,raw_x_cap=1000000` and `columns=16`, `rank_max_trials=100000`, `target_max_attempts=512`, a frozen `rank_seed`; identical K16 control; reference `rho.signed_frobenius_strong_escape` with the `dp_bits` #1387 selects at m = 83 (and `rho.signed_frobenius_strong` where it fits, n ≤ 67) |
| targets | 8 public one-target workloads per curve (`target_kind: public`), fresh seeds, zero overlap with #1387's workloads unless pairing on them is chosen before the run and stated |
| rounds | counts: 1 measured round suffices (counts are deterministic and replayed); wall: not admissible on this Mac (L0) |
| caps | per-child timeout from `ecbench plan`; timeouts, failures and OOMs are kept and reported, never replaced |
| unit | GAE as `ecbench` charges it; IC totals are lower bounds while native work is unpriced |

## 4. Implementation plan (Rust; counts of existing ids must not change)

1. **Width-generic points in the IC framework.**  Introduce a point/field
   abstraction over `koblitz_fast::FastPoint` (u64) and the wide point of
   `koblitz_strong_rho::WideKoblitz` (u128 limbs), or a `Wide` variant of
   `InstanceCtx`, such that `CompactOrbitScanBase`, `FrobeniusMitmOracle`,
   `shared_rank::run_shared_rank_targets` and the rank/linear-algebra
   stages compile for both.  The narrow path must stay byte-identical in
   counts: CI replays every committed session (`ecbench verify --replay`),
   which is the guard.
2. **`ecbench` wiring.**  `solve_wide` gains `ic.shared_rank`;
   `dump_factor_base` accepts `compact-orbit-scan` on `Instance::Wide`.
   Decide, per README §9, whether the method id stays (`ic.shared_rank` on a
   wide instance was refused before, so no committed session's counts
   change) or a new id is registered; record the decision in the PR.
3. **Tests.**  Narrow-vs-wide agreement at n ≤ 61 on identical seeds
   (same base, same pair table, same rank rows, same scalar, same counts);
   a wide smoke at n = 67 that recovers and verifies a planted scalar.
4. **Measurement.**  One `ecbench.spec/v1` per curve under
   `research/ecbench_wide_shared_rank_<date>/sessions/`, run on this Mac for
   counts, `ecbench verify --replay-all`, committed with the receipt;
   the m = 83 cell paired with #1387's reference when it merges (cite its
   session id and receipt; until then the IC/rho cell at 83 is pending).
5. **Dashboard, same PR as the measurement (AGENTS.md §7a):** a point per
   curve in the counted-diagnostics series of `docs/ic/progress-timeline.json`,
   a panel on the scoreboard, the leaderboard `SOURCES` and the browser
   regenerated.

Rough size: the four framework files total ~4,600 lines with 58 narrow-type
sites in `shared_rank.rs` alone; this is days of work and a separate PR
from any measurement.

## 5. Success and stop conditions

- **Success of the implementation PR:** narrow/wide agreement tests pass;
  every committed session still replays identically; the n = 67 smoke
  recovers its scalar.
- **Success of the measurement PR:** every planned cell has a record
  (completed, failed or timed out), every completed cell replays, and the
  table reports the counted K8/rho, K16/rho and K16/K8 per curve with the
  pending m = 83 reference named as pending if #1387 has not merged.
- **Stop:** if the wide shared-rank preparation at n = 83 does not complete
  within 24 h of single-core work per cell, record the timeout and the
  phase it stopped in; do not shrink the cap after seeing it.

## 6. Cost accounting

Counts through `ecbench`'s counters (deterministic, host-independent).
Mac wall time recorded as L0 and labelled exploratory.  Build cost: a full
library rebuild per iteration (~7 min here); disk is the binding constraint
on this machine (see the n41/n53 panel's provenance).

## 7. Checklist (update in follow-on PRs)

- [ ] width-generic IC framework (implementation PR)
- [ ] `solve_wide` admits `ic.shared_rank`; `dump_factor_base` on wide
- [ ] narrow/wide agreement tests; n = 67 smoke
- [ ] counted sessions at n = 67, 71, 73, 79, 83 (K8, K16, control)
- [ ] m = 83 IC/rho cell paired with #1387's reference
- [ ] exponent fit over ≥ 4 sizes; dashboard and leaderboard updated
- [ ] L2 wall-time run (blocked on an authorized isolated host)
