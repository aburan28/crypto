# Drastic prime conductor gap, round 3: paired targets

Measured on public synthetic prime-field curves, 2026-10-12 UTC. Classification: boundary control and counted-work accounting.

Round 3 adds the paired target law `planted_shared`: the same scalar `k` on both curves of a pair (their `r` is equal), so the per-index quotient removes the between-target variance of round 2's independent draws. The D = -3 pairs re-measure under it with 64 workloads; a fresh D = -4 (j = 1728, cofactor 2) pair measures the sqrt(2) rung with the new `rho.quartic`; a D = -19 pair at ell = 4067213 (32x the 131101 the programme started from) sweeps the gap under the theorem.

| family | ell | prediction | same-method paired (99%) | best-method paired (99%) |
|---|---|---|---|---|
| d3 | 1016839 | sqrt(3) = 1.7321 | **1.08144** (0.782885 to 1.47188) | **1.85905** (1.36982 to 2.50032) |
| d3 | 2033609 | sqrt(3) = 1.7321 | **0.931540** (0.696319 to 1.24097) | **1.91641** (1.39821 to 2.61342) |
| d4 | 1016839 | sqrt(2) = 1.4142 | **0.942538** (0.658469 to 1.40884) | **1.13648** (0.745565 to 1.95959) |
| d19 | 4067213 | 1 (the generic theorem at a 30x gap) | **1.01555** (0.646067 to 1.63847) | not run |

The D = -3 sqrt(3) confirms under the paired law: both best-method 99% intervals exclude one. The same-method controls sit at one everywhere (the pairing shows it directly: the crater and floor solve the same scalars for the same counted cost). The D = -19 cell holds the generic theorem at a 32x gap.

## The sqrt(2) rung under-realizes

The D = -4 pair's within-crater quotient reads 1.2058 [0.715 to 2.037] against the predicted sqrt(2) = 1.4142, and its best-method quotient 1.1365 [0.746 to 1.960] does not exclude one at 32 workloads: the quartic walk realizes about 1.2 of its fold where the cubic realizes its full sqrt(3). The walk's look-ahead is tuned for the negation map's 2-cycles; a four-class walk makes longer fruitless cycles likelier, and the counted rotations add one multiplication per canonicalisation. Walk tuning for A = 4 is the open engineering item; the mechanism is not in doubt (the D = -3 cubic and the rho-neg controls bracket it).

## Recorded execution and replay

- **d4-ell1016839-rho**, audit OK: true; verified records: 256; exact measured replays: 128. [audit](sessions/d4-ell1016839-rho.audit.json), [sealed session](sessions/d4-ell1016839-rho/session.json).
- **d4-ell1016839-auto**, audit OK: true; verified records: 64; exact measured replays: 32. [audit](sessions/d4-ell1016839-auto.audit.json), [sealed session](sessions/d4-ell1016839-auto/session.json).
- **d3-ell1016839-rho**, audit OK: true; verified records: 512; exact measured replays: 256. [audit](sessions/d3-ell1016839-rho.audit.json), [sealed session](sessions/d3-ell1016839-rho/session.json).
- **d3-ell1016839-cubic**, audit OK: true; verified records: 128; exact measured replays: 64. [audit](sessions/d3-ell1016839-cubic.audit.json), [sealed session](sessions/d3-ell1016839-cubic/session.json).
- **d3-ell2033609-rho**, audit OK: true; verified records: 512; exact measured replays: 256. [audit](sessions/d3-ell2033609-rho.audit.json), [sealed session](sessions/d3-ell2033609-rho/session.json).
- **d3-ell2033609-cubic**, audit OK: true; verified records: 128; exact measured replays: 64. [audit](sessions/d3-ell2033609-cubic.audit.json), [sealed session](sessions/d3-ell2033609-cubic/session.json).
- **d19-ell4067213-rho**, audit OK: true; verified records: 256; exact measured replays: 128. [audit](sessions/d19-ell4067213-rho.audit.json), [sealed session](sessions/d19-ell4067213-rho/session.json).

## Reproduce

`drastic_gap_round3_pairs --out pairs/pairs.json`; `drastic_gap_round3_specs --pairs pairs/pairs.json --round2 ../drastic_prime_conductor_gap_round2_20261012/pairs/pairs.json --out specs`; `run_sessions.sh`; analysis: `drastic_gap_round3_analyze --sessions research/drastic_prime_conductor_gap_round3_20261012/sessions --out research/drastic_prime_conductor_gap_round3_20261012/analysis`. Frozen seed 20261012; 20,000 resamples.

## Requirement status

| requirement | status | evidence or remaining gap |
|---|---|---|
| Higher-confidence intervals | measured | paired planted_shared law, 64 workloads; round-2 quotients re-measured |
| The sqrt(2) rung | measured, under-realized | rho.quartic realizes ~1.2 of sqrt(2); walk tuning open |
| The 32x gap sweep | measured | same-method = 1 at ell = 4067213 |
| Tier-1 Linux wall times | blocked | tools/isolated_bench.py needs a Linux host with PSI and taskset; none in this environment |
| External-host replay | blocked | same host class |
| Root-finding PDP across the gap | deferred | #1729's harness is not yet an ecbench arm |
| The cubic fold of the MITM table | frozen as round 4 | PROTOCOL.md E4 |

## Evidence and visual scope

[Protocol](PROTOCOL.md), [native tables](analysis.md), [machine-readable analysis](analysis.json), [pair search output](pairs/pairs.json), [report PDF](report.pdf). The scoreboard receives a bounded-count panel; no frontier or leaderboard promotion. The Conductor service was unreachable at session time; the check was attempted before editing.
