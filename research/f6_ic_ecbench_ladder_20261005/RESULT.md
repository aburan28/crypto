# F6-IC on the E_0 ladder: closure capped at 256 points, no exponent signal, decision pending m = 31

**Status: the protocol's decision is not reached.** The m = 13, 19 and 23
sessions are complete and verified. The m = 31 size is unfinished: both
inherited F4 and F6-IC hit the 3,600-second cap on its first target, and
the session was then killed by a worker restart (see *What remains*). The
preregistered fallback size, m = 7, is degenerate: every arm, rho included,
exhausts on its 29-element subgroup. By `PROTOCOL.md` the round therefore
reports no decision. Everything below the decision line is a labelled
description of the three completed sizes, not a verdict.

Inputs: [`PROTOCOL.md`](PROTOCOL.md), [`AMENDMENT_1.md`](AMENDMENT_1.md)
(rounds are independent seeds), [`AMENDMENT_2.md`](AMENDMENT_2.md) (m = 31
stop rule) and the frozen `SPEC-n*.json`. The analysis is
[`ANALYSIS.json`](ANALYSIS.json), written by the native
`examples/f6_ladder_analyze.rs` from the sessions' `records.jsonl` and
`docs/ic/calibration.json`:

```sh
cargo build --release --example f6_ladder_analyze
target/release/examples/f6_ladder_analyze docs/ic/calibration.json \
  research/f6_ic_ecbench_ladder_20261005/ANALYSIS.json \
  research/f6_ic_ecbench_ladder_20261005/sessions/{n13,n19,n23,n31-killed-02}
```

## The table

The unit is `S = charged group-addition equivalents / √r`. IC `S` is a
**lower bound**: the Boolean solver's word XORs are counted but not
charged. F6-IC's geometric point additions *are* charged, at one addition
each. Each figure is a median over eight public one-target workloads. Each
workload's value is the mean of its three measured rounds (Amendment 1).
All 72 measured runs per arm and size verified `[k]G = Q`.

| m | `icv1` slug | base points | arm | `S` (lower bound) | `S` / floor | `S` / same-target rho (min–max) | solver word XORs | reductions | geometric additions | relations |
| ---: | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 13 | `icv1-f2m13-t181-515ee569` | 33 | strong rho | 14.71 | 59.9 | 1 | — | — | — | — |
| 13 | | | inherited F4 | 16.68 | 67.9 | **1.105** (0.720–1.246) | 6.01 × 10⁶ | 1,446 | 0 | 15.8 |
| 13 | | | F6-IC | 814.2 | 3,312 | **54.75** (37.81–84.22) | 4.63 × 10⁶ | 599 | 35,670 | 15.8 |
| 19 | `icv1-f2m19-t797-b6cf2467` | 139 | strong rho | 3.134 | 15.4 | 1 | — | — | — | — |
| 19 | | | inherited F4 | 4.951 | 24.4 | **1.582** (1.446–1.729) | 5.42 × 10⁸ | 21,451 | 0 | 63.7 |
| 19 | | | F6-IC | 6,398 | 31,470 | **2,091** (1,598–2,562) | 4.01 × 10⁸ | 8,784 | 2,312,972 | 63.7 |
| 23 | `icv1-f2m23-t5197-69e76b73` | 275 | strong rho | 1.046 | 5.66 | 1 | — | — | — | — |
| 23 | | | inherited F4 | 2.143 | 11.6 | **2.091** (1.806–2.160) | 1.06 × 10¹⁰ | 196,655 | 0 | 127.3 |
| 23 | | | F6-IC | 2.799 | 15.1 | **2.721** (2.302–2.811) | 1.06 × 10¹⁰ | 196,556 | 944 | 127.3 |
| 31 | `icv1-f2m31-tm90707-c95f16f5` | about 2,048 (predicted) | both IC arms | timeout at 3,600 s on target 1 | — | — | — | — | — | — |

The floors are `√(π/4m)`: 0.2458, 0.2033, 0.1848 and 0.1592. Strong rho's
`S` falls across the ladder (14.7, 3.13, 1.05) because its 32-lane setup,
whose scalar multiplications are charged at `1.5·log₂ r` additions each, is
amortised over a growing `√r`. That rho reference is the frozen one.

**Pinned-ratio sensitivity, m = 13 only** (the protocol's row; m = 31 has no
completed run). Pricing word XORs at `docs/ic/calibration.json`'s pinned
ratio for `icv1-f2m13-t181-515ee569` (0.002304 addition per XOR), the
median `S` is **326.8 for inherited F4 and 1,052.5 for F6-IC**. That ratio
was calibrated on another host, so this row is labelled as such and enters
no headline.

## Decision

`PROTOCOL.md`'s rule needs at least six complete targets at each of the
four sizes. m = 31 has none (one target observed, both arms capped).
Amendment 2's stop condition, three failed m = 31 targets, was not reached
before the session was killed. The fallback size m = 7 is undecidable: all
96 runs exhausted, because the walk and rho both cycle in a 29-element
group. **Decision: none.** `ANALYSIS.json`'s primary and first-round-only
readings both say "not decidable", so they agree.

## What the completed sizes show (descriptive, not the protocol decision)

1. **The closure stops at a hard-coded cap.** #1333's single-fixed-summand
   closure engages only when the base has at most 256 points
   (`koblitz_index_calculus.rs`, `F6GeometricGate`, the
   `fb.points.len() <= 256` test). The base has 33 points at m = 13 and
   139 at m = 19. There F6-IC refutes branches (350 and 3,853 refutations
   in the first record of each session) and cuts word XORs by
   log₂(W4/W6) = 0.377 and 0.433, about 1.30× and 1.35×. At m = 23 the
   base has 275 points, the closure is off, F6-IC records one refutation,
   and the two arms do the same Boolean work: median log₂ ratio 0.000. At
   m = 31 the base is about 2,048 points, so the closure would be off
   there too. The three-size OLS slope of log₂(W4/W6) on `m` is
   **−0.0342 per unit m** (bootstrap 95% −0.0349 to −0.0335). It is driven
   by the cap, not by a trend, and it is not the protocol's statistic,
   which needs four sizes. Nothing here points to an exponent gain. While
   the closure is on, the saving is a constant near 1.3×; past the cap it
   is nothing.
2. **In the counted unit the trade is a relabelling.** F6-IC lowers the
   solver's count (reductions 2.4×, word XORs about 1.3×) by doing
   geometric point additions. Charged at one addition each, those cost
   35,670 at m = 13 and 2.31 million at m = 19 against a whole F4
   relation phase of 720 and 1,769 GAE. Its lower-bound `S` rises 49× and
   1,292×. With word XORs priced at the pinned ratio (m = 13), F6-IC still
   costs 3.2× F4. By AGENTS.md §3 that is **relabelling**: the headline
   count fell and `S` rose.
3. **Wall time disagrees, and says why.** As a practicality note only
   (mixed L1/L2 runs on a 4-vCPU Linux x86-64 VM, partly overlapping other
   sessions), F6-IC's median solve wall was 1.93× faster than F4's at
   m = 13, 1.73× at m = 19 and 0.99× at m = 23. That matches #1333's
   1.3–1.6× on a Mac. The counted unit charges each of F6-IC's packed,
   batch-inverted additions as a full affine addition, while it leaves
   F4's word XORs uncharged. Under that convention F6-IC's real advantage
   reads as a loss. Settling it needs F6's geometric work priced at a
   measured, pinned ratio of its own. That is an accounting question, and
   no exponent changes with the answer.
4. **Against rho, IC loses at every size, and the gap widens.** Inherited
   F4's lower-bound `S` is 1.105×, 1.582× and 2.091× same-target strong rho
   at m = 13, 19 and 23. A lower bound above one proves IC is slower there,
   whatever the unpriced solver work costs. Per-target ranges overlap rho
   only at m = 13. The fitted slope of log₂ `S_lower` on `m` is −0.292 for
   F4 against −0.418 for rho, so on these sizes rho pulls away.

Two disclosures:

- One (target, round) pair at m = 19, workload `W6e35b77370b7` round 3, is
  a divergent query stream: F4 found 65 relations in 103 trials, F6-IC 66
  in 105. It stays in the medians and is listed in `ANALYSIS.json`.
- All runs were on one host class, x86-64 Linux in a cloud VM. Nothing
  here covers Arm64 or GPUs, and no claim transfers to m = 83 or
  ECC2K-130.

## Correctness and evidence

| session | ecbench id | records | status | audit receipt (SHA-256) |
| --- | --- | ---: | --- | --- |
| `sessions/n13` | `ECBS1h1c3104cd2c1d` | 96 | complete, 96 verified | `audits/n13-replay-all.json`: 72 of 72 replays identical, `58969db7…` |
| `sessions/n19` | `ECBS1h6c72efcb063f` | 96 | complete, 96 verified | `audits/n19-replay-all.json`: 72 of 72 identical, `021e0476…` |
| `sessions/n23` | `ECBS1hd8140f7200eb` | 96 | complete, 96 verified | `audits/n23-replay-4.json`: 4 of 4 identical, `ea3935c2…` |
| `sessions/n7` | `ECBS1h302934ff99ed` | 96 | complete, 0 verified (all exhausted) | `audits/n7-replay-12.json`: no verified run to replay, 0 problems, `2dd1f43c…` |
| `sessions/n23-interrupted-01` | `ECBS1h18a9874fbf72` | 8 | interrupted (relaunched detached: the 2-hour background limit) | — |
| `sessions/n31-interrupted-01` | `ECBS1h83ac48e3f4d3` | 2 | interrupted (same relaunch) | — |
| `sessions/n31-killed-02` | `ECBS1haeebe969c2a3` | 4 | killed by a worker restart, still marked `running` | — |

Every run was produced by one `ecbench` binary built from `f4994f9f`. A
rebuild was byte-identical. The interrupted m = 23 attempt's first seven records
match the complete session's first seven exactly: the same workloads,
seeds and `S`. Its eighth is the F4 run the SIGINT cut short, recorded as
`crashed`.

These sessions sit outside `research/ecbench_*`, so the `ecbench`
workflow's committed-evidence step does not replay them. A full m = 23
replay alone takes about 2.5 hours, beyond that job's budget. The
receipts above are this host's. A `vs_rho` claim would also need a
receipt from another host class, and no such claim is made.

## What remains

- **m = 31 on a persistent host.** This container cannot keep a 5–16-hour
  session alive: background jobs stop at two hours, the VM pauses when
  idle, and a worker restart killed the attempt. Running `SPEC-n31.json`
  on an independent runner (`ecbench-independent-runner`), under
  Amendment 2's stop rule, completes the protocol. From the cap mechanism
  above, the expectation is that m = 31 times out or, if it completes,
  shows log₂(W4/W6) ≈ 0. That expectation is extrapolation, not
  measurement.
- **Pricing F6's geometry.** A pinned ratio for packed, batched point
  additions, measured like the existing calibration entries, would turn
  point 3's wall-versus-count disagreement into one priced `S`.
- **The cap itself.** Lifting the 256-point limit would make a different
  candidate, with a new identity and its own protocol. Its cost per
  query, the fixed-summand enumeration over the base, grows with the
  base, so it is not obviously an exponent lever either.
